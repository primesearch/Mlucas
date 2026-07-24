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

int radix26_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-26 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-26 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n26, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
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
	double rt,it,r1,i1,r2,i2,r3,i3,r4,i4
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
	,aj1p00r,aj1p01r,aj1p02r,aj1p03r,aj1p04r,aj1p05r,aj1p06r,aj1p07r,aj1p08r,aj1p09r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r,aj1p18r,aj1p19r,aj1p20r,aj1p21r,aj1p22r,aj1p23r,aj1p24r,aj1p25r
	,aj1p00i,aj1p01i,aj1p02i,aj1p03i,aj1p04i,aj1p05i,aj1p06i,aj1p07i,aj1p08i,aj1p09i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i,aj1p18i,aj1p19i,aj1p20i,aj1p21i,aj1p22i,aj1p23i,aj1p24i,aj1p25i
	,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,temp,frac,scale;
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

/*...change n26 and n_div_wt to non-static to work around a gcc compiler bug. */
	n26   = n/26;
	n_div_nwt = n26 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n26)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/26 in radix26_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)26));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p01 = n26;
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

		bjmodnini=0;
		for(j=0; j < n26; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-26 final DIT pass is here.	*/

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
	cy22= 0;
	cy23= 0;
	cy24= 0;
	cy25= 0;

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
	bjmodn22= bjmodn21+bjmodnini-n; bjmodn22= bjmodn22+ ( (-(int)((uint32)bjmodn22>> 31)) & n);
	bjmodn23= bjmodn22+bjmodnini-n; bjmodn23= bjmodn23+ ( (-(int)((uint32)bjmodn23>> 31)) & n);
	bjmodn24= bjmodn23+bjmodnini-n; bjmodn24= bjmodn24+ ( (-(int)((uint32)bjmodn24>> 31)) & n);
	bjmodn25= bjmodn24+bjmodnini-n; bjmodn25= bjmodn25+ ( (-(int)((uint32)bjmodn25>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+26;
	co3=co2-26;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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
		!...gather the needed data (26 64-bit complex, i.e. 52 64-bit reals) and do 13 radix-2 transforms...
		*/
			aj1p00r = a[j1      ];		aj1p00i = a[j2    ];	/* x0a,b */
			rt      = a[j1  +p01];		it      = a[j2+p01];
			aj1p01r = aj1p00r-rt;		aj1p01i = aj1p00i-it;
			aj1p00r = aj1p00r+rt;		aj1p00i = aj1p00i+it;

			aj1p02r = a[j1  +p23];		aj1p02i = a[j2+p23];	/* x1a,b */
			rt      = a[j1  +p22];		it      = a[j2+p22];
			aj1p03r = aj1p02r-rt;		aj1p03i = aj1p02i-it;
			aj1p02r = aj1p02r+rt;		aj1p02i = aj1p02i+it;

			aj1p04r = a[j1  +p19];		aj1p04i = a[j2+p19];	/* x2a,b */
			rt      = a[j1  +p18];		it      = a[j2+p18];
			aj1p05r = aj1p04r-rt;		aj1p05i = aj1p04i-it;
			aj1p04r = aj1p04r+rt;		aj1p04i = aj1p04i+it;

			aj1p06r = a[j1  +p15];		aj1p06i = a[j2+p15];	/* x3a,b */
			rt      = a[j1  +p14];		it      = a[j2+p14];
			aj1p07r = aj1p06r-rt;		aj1p07i = aj1p06i-it;
			aj1p06r = aj1p06r+rt;		aj1p06i = aj1p06i+it;

			aj1p08r = a[j1  +p11];		aj1p08i = a[j2+p11];	/* x4a,b */
			rt      = a[j1  +p10];		it      = a[j2+p10];
			aj1p09r = aj1p08r-rt;		aj1p09i = aj1p08i-it;
			aj1p08r = aj1p08r+rt;		aj1p08i = aj1p08i+it;

			aj1p10r = a[j1  +p07];		aj1p10i = a[j2+p07];	/* x5a,b */
			rt      = a[j1  +p06];		it      = a[j2+p06];
			aj1p11r = aj1p10r-rt;		aj1p11i = aj1p10i-it;
			aj1p10r = aj1p10r+rt;		aj1p10i = aj1p10i+it;

			aj1p12r = a[j1  +p03];		aj1p12i = a[j2+p03];	/* x6a,b */
			rt      = a[j1  +p02];		it      = a[j2+p02];
			aj1p13r = aj1p12r-rt;		aj1p13i = aj1p12i-it;
			aj1p12r = aj1p12r+rt;		aj1p12i = aj1p12i+it;

			aj1p14r = a[j1  +p24];		aj1p14i = a[j2+p24];	/* x7a,b */
			rt      = a[j1  +p25];		it      = a[j2+p25];
			aj1p15r = aj1p14r-rt;		aj1p15i = aj1p14i-it;
			aj1p14r = aj1p14r+rt;		aj1p14i = aj1p14i+it;

			aj1p16r = a[j1  +p20];		aj1p16i = a[j2+p20];	/* x8a,b */
			rt      = a[j1  +p21];		it      = a[j2+p21];
			aj1p17r = aj1p16r-rt;		aj1p17i = aj1p16i-it;
			aj1p16r = aj1p16r+rt;		aj1p16i = aj1p16i+it;

			aj1p18r = a[j1  +p16];		aj1p18i = a[j2+p16];	/* x9a,b */
			rt      = a[j1  +p17];		it      = a[j2+p17];
			aj1p19r = aj1p18r-rt;		aj1p19i = aj1p18i-it;
			aj1p18r = aj1p18r+rt;		aj1p18i = aj1p18i+it;

			aj1p20r = a[j1  +p12];		aj1p20i = a[j2+p12];	/* x10a,b */
			rt      = a[j1  +p13];		it      = a[j2+p13];
			aj1p21r = aj1p20r-rt;		aj1p21i = aj1p20i-it;
			aj1p20r = aj1p20r+rt;		aj1p20i = aj1p20i+it;

			aj1p22r = a[j1  +p08];		aj1p22i = a[j2+p08];	/* x11a,b */
			rt      = a[j1  +p09];		it      = a[j2+p09];
			aj1p23r = aj1p22r-rt;		aj1p23i = aj1p22i-it;
			aj1p22r = aj1p22r+rt;		aj1p22i = aj1p22i+it;

			aj1p24r = a[j1  +p04];		aj1p24i = a[j2+p04];	/* x12a,b */
			rt      = a[j1  +p05];		it      = a[j2+p05];
			aj1p25r = aj1p24r-rt;		aj1p25i = aj1p24i-it;
			aj1p24r = aj1p24r+rt;		aj1p24i = aj1p24i+it;

	/*      ..and now do two radix-13 transforms.	*/

		/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
			t02=aj1p02r+aj1p24r;	t03=aj1p02i+aj1p24i;
			t14=aj1p02r-aj1p24r;	t15=aj1p02i-aj1p24i;

			t04=aj1p04r+aj1p22r;	t05=aj1p04i+aj1p22i;
			t16=aj1p04r-aj1p22r;	t17=aj1p04i-aj1p22i;

			t10=aj1p06r+aj1p20r;	t11=aj1p06i+aj1p20i;
			t22=aj1p06r-aj1p20r;	t23=aj1p06i-aj1p20i;

			t06=aj1p08r+aj1p18r;	t07=aj1p08i+aj1p18i;
			t18=aj1p08r-aj1p18r;	t19=aj1p08i-aj1p18i;

			t08=aj1p10r+aj1p16r;	t09=aj1p10i+aj1p16i;
			t20=aj1p16r-aj1p10r;	t21=aj1p16i-aj1p10i;		/* flip signs on as3 */

			t12=aj1p12r+aj1p14r;	t13=aj1p12i+aj1p14i;
			t24=aj1p12r-aj1p14r;	t25=aj1p12i-aj1p14i;

		/*** COSINE TERMS ***/

			t00 = t02+t06+t10;	t01 = t03+t07+t11;
			t02 = t02-t06   ;	t03 = t03-t07   ;
			t06 =     t06-t10;	t07 =     t07-t11;

			r1  = t04+t08+t12;	i1  = t05+t09+t13;
			t04 = t04    -t12;	t05 = t05    -t13;
			t08 =     t08-t12;	t09 =     t09-t13;

			t10 = t00+r1;		t11 = t01+i1;
			t00 = t00-r1;		t01 = t01-i1;
			t12 = t02-t08;		t13 = t03-t09;
			t02 = t02+t08;		t03 = t03+t09;
			t08 = t04+t06;		t09 = t05+t07;
			t06 = t04-t06;		t07 = t05-t07;

			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p00r;	i1 = t11*bp0 + aj1p00i;

			/* polypro(A,B) mod P1      */
			t00= t00*bp1;		t01= t01*bp1;

			aj1p00r += t10;		aj1p00i += t11;

			/* polypro(A,B) mod P2 = x^2-x+1: */
			t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
			t12 = t12*bp20;		t13 = t13*bp20;
			t08 = t08*bp21;		t09 = t09*bp21;
			t08 = t12 - t08;	t09 = t13 - t09;
			t04 = t04 - t12;	t05 = t05 - t13;

			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
			t02 = t02*bp30;		t03 = t03*bp30;
			t06 = t06*bp31;		t07 = t07*bp31;
			t06 = t02 - t06;	t07 = t03 - t07;
			t10 = t10 + t02;	t11 = t11 + t03;

			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t02= r1 +t00;		t03= i1 +t01;	/* a */
			t00= r1 -t00;		t01= i1 -t01;	/* b */
			r1 = t08+t06;		i1 = t09+t07;	/* c */
			t08= t08-t06;		t09= t09-t07;	/* d */
			t06= t04+t10;		t07= t05+t11;	/* e */
			t04= t04-t10;		t05= t05-t11;	/* f */
	#if(LO_ADD)
			t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
			t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t02-r1;		t11 = t03-i1;
			t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
			t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
			t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

			r1  = t00+t08;		i1  = t01+t09;
			t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
			t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
			t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
	#endif
		/*** SINE TERMS ***/

			r1 = t14-t18+t22;	i1 = t15-t19+t23;
			t14= t14-t22;		t15= t15-t23;
			t18= t18+t22;		t19= t19+t23;

			t00= t16-t20+t24;	t01= t17-t21+t25;
			t16= t16-t24;		t17= t17-t25;
			t20= t20+t24;		t21= t21+t25;

			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
			r1  = r1 *bq00;		i1  = i1 *bq00;
			t00 = t00*bq01;		t01 = t01*bq01;
			t24 = r1 - t00;		t25 = i1 - t01;
			t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;		i1  = t15+t19;
			t00 = t14*bq10;		t01 = t15*bq10;
			r2  = t18*bq12;		i2  = t19*bq12;
			r2  = t00 - r2;	i2  = t01 - i2;
			t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

			t14 = t14*bq11;		t15 = t15*bq11;
			t18 = t18*bq13;		t19 = t19*bq13;
			t18 = t14 - t18;	t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

			r1  = t16+t20;		i1  = t17+t21;
			r3  = t16*bq10;		i3  = t17*bq10;
			r4  = t20*bq12;		i4  = t21*bq12;
			r4  = r3 - r4;		i4  = i3 - i4;
			r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

			t16 = t16*bq11;		t17 = t17*bq11;
			t20 = t20*bq13;		t21 = t21*bq13;
			t20 = t16 - t20;	t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t00+t16;	t21 = t21+t01+t17;
			t16 = r2 -t16;		t17 = i2 -t17;
			t18 = t18+r4;		t19 = t19+i4;
			t14 = t14+r3;		t15 = t15+i3;

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

			aj1p12r=t12-t17;	aj1p12i=t13+t16;	/* C1 - S1 */
			aj1p24r=t10+t23;	aj1p24i=t11-t22;	/* C2 - S2 */
			aj1p10r=t02+t21;	aj1p10i=t03-t20;	/* C3 - S3 */
			aj1p22r=t04+t25;	aj1p22i=t05-t24;	/* C4 - S4 */
			aj1p08r=t08-t15;	aj1p08i=t09+t14;	/* C5 - S5 */
			aj1p20r=t06+t19;	aj1p20i=t07-t18;	/* C6 - S6 */
			aj1p06r=t06-t19;	aj1p06i=t07+t18;	/* C6 + S6 */
			aj1p18r=t08+t15;	aj1p18i=t09-t14;	/* C5 + S5 */
			aj1p04r=t04-t25;	aj1p04i=t05+t24;	/* C4 + S4 */
			aj1p16r=t02-t21;	aj1p16i=t03+t20;	/* C3 + S3 */
			aj1p02r=t10-t23;	aj1p02i=t11+t22;	/* C2 + S2 */
			aj1p14r=t12+t17;	aj1p14i=t13-t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t00 = t22+t18;		t01 = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
			t20 = t00+t18;		t21 = t01+t19;		/* S6 */
			t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

		/*...Inline multiply of sine parts by +-I ]nto finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/

			aj1p12r=t12-t15;	aj1p12i=t13+t14;	/* C1 - S1 */
			aj1p24r=t06+t23;	aj1p24i=t07-t22;	/* C2 - S2 */
			aj1p10r=t02+t17;	aj1p10i=t03-t16;	/* C3 - S3 */
			aj1p22r=t04+t25;	aj1p22i=t05-t24;	/* C4 - S4 */
			aj1p08r=t10-t19;	aj1p08i=t11+t18;	/* C5 - S5 */
			aj1p20r=t08+t21;	aj1p20i=t09-t20;	/* C6 - S6 */
			aj1p06r=t08-t21;	aj1p06i=t09+t20;	/* C6 + S6 */
			aj1p18r=t10+t19;	aj1p18i=t11-t18;	/* C5 + S5 */
			aj1p04r=t04-t25;	aj1p04i=t05+t24;	/* C4 + S4 */
			aj1p16r=t02-t17;	aj1p16i=t03+t16;	/* C3 + S3 */
			aj1p02r=t06-t23;	aj1p02i=t07+t22;	/* C2 + S2 */
			aj1p14r=t12+t15;	aj1p14i=t13-t14;	/* C1 + S1 */
	#endif

		/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
			t02=aj1p03r+aj1p25r;	t03=aj1p03i+aj1p25i;
			t14=aj1p03r-aj1p25r;	t15=aj1p03i-aj1p25i;

			t04=aj1p05r+aj1p23r;	t05=aj1p05i+aj1p23i;
			t16=aj1p05r-aj1p23r;	t17=aj1p05i-aj1p23i;

			t10=aj1p07r+aj1p21r;	t11=aj1p07i+aj1p21i;
			t22=aj1p07r-aj1p21r;	t23=aj1p07i-aj1p21i;

			t06=aj1p09r+aj1p19r;	t07=aj1p09i+aj1p19i;
			t18=aj1p09r-aj1p19r;	t19=aj1p09i-aj1p19i;

			t08=aj1p11r+aj1p17r;	t09=aj1p11i+aj1p17i;
			t20=aj1p17r-aj1p11r;	t21=aj1p17i-aj1p11i;		/* flip signs on as3 */

			t12=aj1p13r+aj1p15r;	t13=aj1p13i+aj1p15i;
			t24=aj1p13r-aj1p15r;	t25=aj1p13i-aj1p15i;

		/*** COSINE TERMS ***/

			t00 = t02+t06+t10;	t01 = t03+t07+t11;
			t02 = t02-t06   ;	t03 = t03-t07   ;
			t06 =     t06-t10;	t07 =     t07-t11;

			r1  = t04+t08+t12;	i1  = t05+t09+t13;
			t04 = t04    -t12;	t05 = t05    -t13;
			t08 =     t08-t12;	t09 =     t09-t13;

			t10 = t00+r1;		t11 = t01+i1;
			t00 = t00-r1;		t01 = t01-i1;
			t12 = t02-t08;		t13 = t03-t09;
			t02 = t02+t08;		t03 = t03+t09;
			t08 = t04+t06;		t09 = t05+t07;
			t06 = t04-t06;		t07 = t05-t07;

			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p01r;	i1 = t11*bp0 + aj1p01i;

			/* polypro(A,B) mod P1      */
			t00= t00*bp1;		t01= t01*bp1;

			aj1p01r += t10;		aj1p01i += t11;
			aj1p13r = aj1p01r;	aj1p13i = aj1p01i;

			/* polypro(A,B) mod P2 = x^2-x+1: */
			t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
			t12 = t12*bp20;		t13 = t13*bp20;
			t08 = t08*bp21;		t09 = t09*bp21;
			t08 = t12 - t08;	t09 = t13 - t09;
			t04 = t04 - t12;	t05 = t05 - t13;

			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
			t02 = t02*bp30;		t03 = t03*bp30;
			t06 = t06*bp31;		t07 = t07*bp31;
			t06 = t02 - t06;	t07 = t03 - t07;
			t10 = t10 + t02;	t11 = t11 + t03;

			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t02= r1 +t00;		t03= i1 +t01;	/* a */
			t00= r1 -t00;		t01= i1 -t01;	/* b */
			r1 = t08+t06;		i1 = t09+t07;	/* c */
			t08= t08-t06;		t09= t09-t07;	/* d */
			t06= t04+t10;		t07= t05+t11;	/* e */
			t04= t04-t10;		t05= t05-t11;	/* f */
	#if(LO_ADD)
			t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
			t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t02-r1;		t11 = t03-i1;
			t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
			t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
			t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

			r1  = t00+t08;		i1  = t01+t09;
			t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
			t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
			t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
	#endif
		/*** SINE TERMS ***/

			r1 = t14-t18+t22;	i1 = t15-t19+t23;
			t14= t14-t22;		t15= t15-t23;
			t18= t18+t22;		t19= t19+t23;

			t00= t16-t20+t24;	t01= t17-t21+t25;
			t16= t16-t24;		t17= t17-t25;
			t20= t20+t24;		t21= t21+t25;

			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
			r1  = r1 *bq00;		i1  = i1 *bq00;
			t00 = t00*bq01;		t01 = t01*bq01;
			t24 = r1 - t00;		t25 = i1 - t01;
			t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;		i1  = t15+t19;
			t00 = t14*bq10;		t01 = t15*bq10;
			r2  = t18*bq12;		i2  = t19*bq12;
			r2  = t00 - r2;	i2  = t01 - i2;
			t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

			t14 = t14*bq11;		t15 = t15*bq11;
			t18 = t18*bq13;		t19 = t19*bq13;
			t18 = t14 - t18;	t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

			r1  = t16+t20;		i1  = t17+t21;
			r3  = t16*bq10;		i3  = t17*bq10;
			r4  = t20*bq12;		i4  = t21*bq12;
			r4  = r3 - r4;		i4  = i3 - i4;
			r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

			t16 = t16*bq11;		t17 = t17*bq11;
			t20 = t20*bq13;		t21 = t21*bq13;
			t20 = t16 - t20;	t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t00+t16;	t21 = t21+t01+t17;
			t16 = r2 -t16;		t17 = i2 -t17;
			t18 = t18+r4;		t19 = t19+i4;
			t14 = t14+r3;		t15 = t15+i3;

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

			aj1p25r=t12-t17;	aj1p25i=t13+t16;	/* C1 - S1 */
			aj1p11r=t10+t23;	aj1p11i=t11-t22;	/* C2 - S2 */
			aj1p23r=t02+t21;	aj1p23i=t03-t20;	/* C3 - S3 */
			aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
			aj1p21r=t08-t15;	aj1p21i=t09+t14;	/* C5 - S5 */
			aj1p07r=t06+t19;	aj1p07i=t07-t18;	/* C6 - S6 */
			aj1p19r=t06-t19;	aj1p19i=t07+t18;	/* C6 + S6 */
			aj1p05r=t08+t15;	aj1p05i=t09-t14;	/* C5 + S5 */
			aj1p17r=t04-t25;	aj1p17i=t05+t24;	/* C4 + S4 */
			aj1p03r=t02-t21;	aj1p03i=t03+t20;	/* C3 + S3 */
			aj1p15r=t10-t23;	aj1p15i=t11+t22;	/* C2 + S2 */
			aj1p01r=t12+t17;	aj1p01i=t13-t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t00 = t22+t18;		t01 = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
			t20 = t00+t18;		t21 = t01+t19;		/* S6 */
			t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

		/*...Inline multiply of sine parts by +-I ]nto finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/

			aj1p25r=t12-t15;	aj1p25i=t13+t14;	/* C1 - S1 */
			aj1p11r=t06+t23;	aj1p11i=t07-t22;	/* C2 - S2 */
			aj1p23r=t02+t17;	aj1p23i=t03-t16;	/* C3 - S3 */
			aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
			aj1p21r=t10-t19;	aj1p21i=t11+t18;	/* C5 - S5 */
			aj1p07r=t08+t21;	aj1p07i=t09-t20;	/* C6 - S6 */
			aj1p19r=t08-t21;	aj1p19i=t09+t20;	/* C6 + S6 */
			aj1p05r=t10+t19;	aj1p05i=t11-t18;	/* C5 + S5 */
			aj1p17r=t04-t25;	aj1p17i=t05+t24;	/* C4 + S4 */
			aj1p03r=t02-t17;	aj1p03i=t03+t16;	/* C3 + S3 */
			aj1p15r=t06-t23;	aj1p15i=t07+t22;	/* C2 + S2 */
			aj1p01r=t12+t15;	aj1p01i=t13-t14;	/* C1 + S1 */
	#endif
/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 26 separate blocks of the A-array, we need 26 separate carries.	*/

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
			cmplx_carry_norm_errcheck(aj1p22r,aj1p22i,cy22,bjmodn22,22,prp_mult);
			cmplx_carry_norm_errcheck(aj1p23r,aj1p23i,cy23,bjmodn23,23,prp_mult);
			cmplx_carry_norm_errcheck(aj1p24r,aj1p24i,cy24,bjmodn24,24,prp_mult);
			cmplx_carry_norm_errcheck(aj1p25r,aj1p25i,cy25,bjmodn25,25,prp_mult);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-26 DIF pass is here:	*/
	#if PFETCH
	add0 = &a[j1];
	prefetch_p_doubles(add0);
	#endif
		/*...First radix-13 block uses aj1p[0:24:2] as inputs:	*/

			t02 =aj1p24r;		t03 =aj1p24i;
			r1  =aj1p02r;		i1  =aj1p02i;
			t14 =t02 -r1;		t15 =t03 -i1;
			t02 =t02 +r1;		t03 =t03 +i1;

			t04 =aj1p22r;		t05 =aj1p22i;
			r1  =aj1p04r;		i1  =aj1p04i;
			t16 =t04 -r1;		t17 =t05 -i1;
			t04 =t04 +r1;		t05 =t05 +i1;

			t10 =aj1p20r;		t11 =aj1p20i;
			r1  =aj1p06r;		i1  =aj1p06i;
			t22 =t10 -r1;		t23 =t11 -i1;
			t10 =t10 +r1;		t11 =t11 +i1;
	#if PFETCH
	addr = add0+p01;
	prefetch_p_doubles(addr);
	#endif
			t06 =aj1p18r;		t07 =aj1p18i;
			r1  =aj1p08r;		i1  =aj1p08i;
			t18 =t06 -r1;		t19 =t07 -i1;
			t06 =t06 +r1;		t07 =t07 +i1;

			t08 =aj1p16r;		t09 =aj1p16i;
			r1  =aj1p10r;		i1  =aj1p10i;
			t20 =r1 - t08;		t21 =i1 - t09;		/* flip signs on as3 */
			t08 =t08 +r1;		t09 =t09 +i1;

			t12 =aj1p14r;		t13 =aj1p14i;
			r1  =aj1p12r;		i1  =aj1p12i;
			t24 =t12 -r1;		t25 =t13 -i1;
			t12 =t12 +r1;		t13 =t13 +i1;
	#if PFETCH
	addr = add0+p02;
	prefetch_p_doubles(addr);
	#endif

		/*** COSINE TERMS ***/

			t00 = t02+t06+t10;	t01 = t03+t07+t11;
			t02 = t02-t06   ;	t03 = t03-t07   ;
			t06 =     t06-t10;	t07 =     t07-t11;

			r1  = t04+t08+t12;	i1  = t05+t09+t13;
			t04 = t04    -t12;	t05 = t05    -t13;
			t08 =     t08-t12;	t09 =     t09-t13;

			t10 = t00+r1;		t11 = t01+i1;
			t00 = t00-r1;		t01 = t01-i1;
			t12 = t02-t08;		t13 = t03-t09;
			t02 = t02+t08;		t03 = t03+t09;
			t08 = t04+t06;		t09 = t05+t07;
			t06 = t04-t06;		t07 = t05-t07;

			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p00r;	i1 = t11*bp0 + aj1p00i;

			/* polypro(A,B) mod P1      */
			t00= t00*bp1;		t01= t01*bp1;

			aj1p00r += t10;		aj1p00i += t11;
	#if PFETCH
	addr = add0+p03;
	prefetch_p_doubles(addr);
	#endif
			/* polypro(A,B) mod P2 = x^2-x+1: */
			t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
			t12 = t12*bp20;		t13 = t13*bp20;
			t08 = t08*bp21;		t09 = t09*bp21;
			t08 = t12 - t08;	t09 = t13 - t09;
			t04 = t04 - t12;	t05 = t05 - t13;

			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
			t02 = t02*bp30;		t03 = t03*bp30;
			t06 = t06*bp31;		t07 = t07*bp31;
			t06 = t02 - t06;	t07 = t03 - t07;
			t10 = t10 + t02;	t11 = t11 + t03;

			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t02= r1 +t00;		t03= i1 +t01;	/* a */
			t00= r1 -t00;		t01= i1 -t01;	/* b */
			r1 = t08+t06;		i1 = t09+t07;	/* c */
			t08= t08-t06;		t09= t09-t07;	/* d */
			t06= t04+t10;		t07= t05+t11;	/* e */
			t04= t04-t10;		t05= t05-t11;	/* f */
	#if PFETCH
	addr = add0+p04;
	prefetch_p_doubles(addr);
	#endif

	#if(LO_ADD)
			t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
			t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t02-r1;		t11 = t03-i1;
			t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
			t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
			t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

			r1  = t00+t08;		i1  = t01+t09;
			t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
			t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
			t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
	#endif
		/*** SINE TERMS ***/
	#if PFETCH
	addr = add0+p05;
	prefetch_p_doubles(addr);
	#endif
			r1 = t14-t18+t22;	i1 = t15-t19+t23;
			t14= t14-t22;		t15= t15-t23;
			t18= t18+t22;		t19= t19+t23;

			t00= t16-t20+t24;	t01= t17-t21+t25;
			t16= t16-t24;		t17= t17-t25;
			t20= t20+t24;		t21= t21+t25;

			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
			r1  = r1 *bq00;		i1  = i1 *bq00;
			t00 = t00*bq01;		t01 = t01*bq01;
			t24 = r1 - t00;		t25 = i1 - t01;
			t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */
	#if PFETCH
	addr = add0+p06;
	prefetch_p_doubles(addr);
	#endif
			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;		i1  = t15+t19;
			t00 = t14*bq10;		t01 = t15*bq10;
			r2  = t18*bq12;		i2  = t19*bq12;
			r2  = t00 - r2;	i2  = t01 - i2;
			t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

			t14 = t14*bq11;		t15 = t15*bq11;
			t18 = t18*bq13;		t19 = t19*bq13;
			t18 = t14 - t18;	t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;
	#if PFETCH
	addr = add0+p07;
	prefetch_p_doubles(addr);
	#endif
			r1  = t16+t20;		i1  = t17+t21;
			r3  = t16*bq10;		i3  = t17*bq10;
			r4  = t20*bq12;		i4  = t21*bq12;
			r4  = r3 - r4;		i4  = i3 - i4;
			r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

			t16 = t16*bq11;		t17 = t17*bq11;
			t20 = t20*bq13;		t21 = t21*bq13;
			t20 = t16 - t20;	t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t00+t16;	t21 = t21+t01+t17;
			t16 = r2 -t16;		t17 = i2 -t17;
			t18 = t18+r4;		t19 = t19+i4;
			t14 = t14+r3;		t15 = t15+i3;
	#if PFETCH
	addr = add0+p08;
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
	addr = add0+p09;
	prefetch_p_doubles(addr);
	#endif
			aj1p02r=t12-t17;	aj1p02i=t13+t16;	/* C1 - S1 */
			aj1p04r=t10+t23;	aj1p04i=t11-t22;	/* C2 - S2 */
			aj1p06r=t02+t21;	aj1p06i=t03-t20;	/* C3 - S3 */
			aj1p08r=t04+t25;	aj1p08i=t05-t24;	/* C4 - S4 */
			aj1p10r=t08-t15;	aj1p10i=t09+t14;	/* C5 - S5 */
			aj1p12r=t06+t19;	aj1p12i=t07-t18;	/* C6 - S6 */
			aj1p14r=t06-t19;	aj1p14i=t07+t18;	/* C6 + S6 */
			aj1p16r=t08+t15;	aj1p16i=t09-t14;	/* C5 + S5 */
			aj1p18r=t04-t25;	aj1p18i=t05+t24;	/* C4 + S4 */
			aj1p20r=t02-t21;	aj1p20i=t03+t20;	/* C3 + S3 */
			aj1p22r=t10-t23;	aj1p22i=t11+t22;	/* C2 + S2 */
			aj1p24r=t12+t17;	aj1p24i=t13-t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t00 = t22+t18;		t01 = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
			t20 = t00+t18;		t21 = t01+t19;		/* S6 */
			t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/
	#if PFETCH
	addr = add0+p09;
	prefetch_p_doubles(addr);
	#endif
			aj1p02r=t12-t15;	aj1p02i=t13+t14;	/* C1 - S1 */
			aj1p04r=t06+t23;	aj1p04i=t07-t22;	/* C2 - S2 */
			aj1p06r=t02+t17;	aj1p06i=t03-t16;	/* C3 - S3 */
			aj1p08r=t04+t25;	aj1p08i=t05-t24;	/* C4 - S4 */
			aj1p10r=t10-t19;	aj1p10i=t11+t18;	/* C5 - S5 */
			aj1p12r=t08+t21;	aj1p12i=t09-t20;	/* C6 - S6 */
			aj1p14r=t08-t21;	aj1p14i=t09+t20;	/* C6 + S6 */
			aj1p16r=t10+t19;	aj1p16i=t11-t18;	/* C5 + S5 */
			aj1p18r=t04-t25;	aj1p18i=t05+t24;	/* C4 + S4 */
			aj1p20r=t02-t17;	aj1p20i=t03+t16;	/* C3 + S3 */
			aj1p22r=t06-t23;	aj1p22i=t07+t22;	/* C2 + S2 */
			aj1p24r=t12+t15;	aj1p24i=t13-t14;	/* C1 + S1 */
	#endif

		/*...Second radix-13 block uses aj1p[1:25:2] as inputs:	*/
	#if PFETCH
	addr = add0+p10;
	prefetch_p_doubles(addr);
	#endif
			t02 =aj1p11r;		t03 =aj1p11i;
			r1  =aj1p15r;		i1  =aj1p15i;
			t14 =t02 -r1;		t15 =t03 -i1;
			t02 =t02 +r1;		t03 =t03 +i1;

			t04 =aj1p09r;		t05 =aj1p09i;
			r1  =aj1p17r;		i1  =aj1p17i;
			t16 =t04 -r1;		t17 =t05 -i1;
			t04 =t04 +r1;		t05 =t05 +i1;

			t10 =aj1p07r;		t11 =aj1p07i;
			r1  =aj1p19r;		i1  =aj1p19i;
			t22 =t10 -r1;		t23 =t11 -i1;
			t10 =t10 +r1;		t11 =t11 +i1;
	#if PFETCH
	addr = add0+p11;
	prefetch_p_doubles(addr);
	#endif
			t06 =aj1p05r;		t07 =aj1p05i;
			r1  =aj1p21r;		i1  =aj1p21i;
			t18 =t06 -r1;		t19 =t07 -i1;
			t06 =t06 +r1;		t07 =t07 +i1;

			t08 =aj1p03r;		t09 =aj1p03i;
			r1  =aj1p23r;		i1  =aj1p23i;
			t20 =r1 - t08;		t21 =i1 - t09;		/* flip signs on as3 */
			t08 =t08 +r1;		t09 =t09 +i1;

			t12 =aj1p01r;		t13 =aj1p01i;
			r1  =aj1p25r;		i1  =aj1p25i;
			t24 =t12 -r1;		t25 =t13 -i1;
			t12 =t12 +r1;		t13 =t13 +i1;
	#if PFETCH
	addr = add0+p12;
	prefetch_p_doubles(addr);
	#endif
			aj1p01r = aj1p13r;	aj1p01i = aj1p13i;

		/*** COSINE TERMS ***/

			t00 = t02+t06+t10;	t01 = t03+t07+t11;
			t02 = t02-t06   ;	t03 = t03-t07   ;
			t06 =     t06-t10;	t07 =     t07-t11;

			r1  = t04+t08+t12;	i1  = t05+t09+t13;
			t04 = t04    -t12;	t05 = t05    -t13;
			t08 =     t08-t12;	t09 =     t09-t13;

			t10 = t00+r1;		t11 = t01+i1;
			t00 = t00-r1;		t01 = t01-i1;
			t12 = t02-t08;		t13 = t03-t09;
			t02 = t02+t08;		t03 = t03+t09;
			t08 = t04+t06;		t09 = t05+t07;
			t06 = t04-t06;		t07 = t05-t07;

			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p01r;	i1 = t11*bp0 + aj1p01i;

			/* polypro(A,B) mod P1      */
			t00= t00*bp1;		t01= t01*bp1;

			aj1p01r += t10;		aj1p01i += t11;
	#if PFETCH
	addr = add0+p13;
	prefetch_p_doubles(addr);
	#endif
			/* polypro(A,B) mod P2 = x^2-x+1: */
			t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
			t12 = t12*bp20;		t13 = t13*bp20;
			t08 = t08*bp21;		t09 = t09*bp21;
			t08 = t12 - t08;	t09 = t13 - t09;
			t04 = t04 - t12;	t05 = t05 - t13;

			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
			t02 = t02*bp30;		t03 = t03*bp30;
			t06 = t06*bp31;		t07 = t07*bp31;
			t06 = t02 - t06;	t07 = t03 - t07;
			t10 = t10 + t02;	t11 = t11 + t03;

			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t02= r1 +t00;		t03= i1 +t01;	/* a */
			t00= r1 -t00;		t01= i1 -t01;	/* b */
			r1 = t08+t06;		i1 = t09+t07;	/* c */
			t08= t08-t06;		t09= t09-t07;	/* d */
			t06= t04+t10;		t07= t05+t11;	/* e */
			t04= t04-t10;		t05= t05-t11;	/* f */
	#if PFETCH
	addr = add0+p14;
	prefetch_p_doubles(addr);
	#endif

	#if(LO_ADD)
			t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
			t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t02-r1;		t11 = t03-i1;
			t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
			t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
			t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

			r1  = t00+t08;		i1  = t01+t09;
			t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
			t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
			t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
	#endif
		/*** SINE TERMS ***/
	#if PFETCH
	addr = add0+p15;
	prefetch_p_doubles(addr);
	#endif
			r1 = t14-t18+t22;	i1 = t15-t19+t23;
			t14= t14-t22;		t15= t15-t23;
			t18= t18+t22;		t19= t19+t23;

			t00= t16-t20+t24;	t01= t17-t21+t25;
			t16= t16-t24;		t17= t17-t25;
			t20= t20+t24;		t21= t21+t25;

			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
			r1  = r1 *bq00;		i1  = i1 *bq00;
			t00 = t00*bq01;		t01 = t01*bq01;
			t24 = r1 - t00;		t25 = i1 - t01;
			t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */
	#if PFETCH
	addr = add0+p16;
	prefetch_p_doubles(addr);
	#endif
			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;		i1  = t15+t19;
			t00 = t14*bq10;		t01 = t15*bq10;
			r2  = t18*bq12;		i2  = t19*bq12;
			r2  = t00 - r2;	i2  = t01 - i2;
			t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

			t14 = t14*bq11;		t15 = t15*bq11;
			t18 = t18*bq13;		t19 = t19*bq13;
			t18 = t14 - t18;	t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;
	#if PFETCH
	addr = add0+p17;
	prefetch_p_doubles(addr);
	#endif
			r1  = t16+t20;		i1  = t17+t21;
			r3  = t16*bq10;		i3  = t17*bq10;
			r4  = t20*bq12;		i4  = t21*bq12;
			r4  = r3 - r4;		i4  = i3 - i4;
			r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

			t16 = t16*bq11;		t17 = t17*bq11;
			t20 = t20*bq13;		t21 = t21*bq13;
			t20 = t16 - t20;	t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t00+t16;	t21 = t21+t01+t17;
			t16 = r2 -t16;		t17 = i2 -t17;
			t18 = t18+r4;		t19 = t19+i4;
			t14 = t14+r3;		t15 = t15+i3;
	#if PFETCH
	addr = add0+p18;
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
	addr = add0+p19;
	prefetch_p_doubles(addr);
	#endif
			aj1p03r=t12-t17;	aj1p03i=t13+t16;	/* C1 - S1 */
			aj1p05r=t10+t23;	aj1p05i=t11-t22;	/* C2 - S2 */
			aj1p07r=t02+t21;	aj1p07i=t03-t20;	/* C3 - S3 */
			aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
			aj1p11r=t08-t15;	aj1p11i=t09+t14;	/* C5 - S5 */
			aj1p13r=t06+t19;	aj1p13i=t07-t18;	/* C6 - S6 */
			aj1p15r=t06-t19;	aj1p15i=t07+t18;	/* C6 + S6 */
			aj1p17r=t08+t15;	aj1p17i=t09-t14;	/* C5 + S5 */
			aj1p19r=t04-t25;	aj1p19i=t05+t24;	/* C4 + S4 */
			aj1p21r=t02-t21;	aj1p21i=t03+t20;	/* C3 + S3 */
			aj1p23r=t10-t23;	aj1p23i=t11+t22;	/* C2 + S2 */
			aj1p25r=t12+t17;	aj1p25i=t13-t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t00 = t22+t18;		t01 = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
			t20 = t00+t18;		t21 = t01+t19;		/* S6 */
			t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/
	#if PFETCH
	addr = add0+p19;
	prefetch_p_doubles(addr);
	#endif
			aj1p03r=t12-t15;	aj1p03i=t13+t14;	/* C1 - S1 */
			aj1p05r=t06+t23;	aj1p05i=t07-t22;	/* C2 - S2 */
			aj1p07r=t02+t17;	aj1p07i=t03-t16;	/* C3 - S3 */
			aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
			aj1p11r=t10-t19;	aj1p11i=t11+t18;	/* C5 - S5 */
			aj1p13r=t08+t21;	aj1p13i=t09-t20;	/* C6 - S6 */
			aj1p15r=t08-t21;	aj1p15i=t09+t20;	/* C6 + S6 */
			aj1p17r=t10+t19;	aj1p17i=t11-t18;	/* C5 + S5 */
			aj1p19r=t04-t25;	aj1p19i=t05+t24;	/* C4 + S4 */
			aj1p21r=t02-t17;	aj1p21i=t03+t16;	/* C3 + S3 */
			aj1p23r=t06-t23;	aj1p23i=t07+t22;	/* C2 + S2 */
			aj1p25r=t12+t15;	aj1p25i=t13-t14;	/* C1 + S1 */
	#endif

	/*...and now do 13 radix-2 transforms:	*/

	#if PFETCH
	addr = add0+p20;
	prefetch_p_doubles(addr);
	#endif
			a[j1    ]=aj1p00r+aj1p01r;	a[j2    ]=aj1p00i+aj1p01i;
			a[j1+p01]=aj1p00r-aj1p01r;	a[j2+p01]=aj1p00i-aj1p01i;

			a[j1+p03]=aj1p02r+aj1p03r;	a[j2+p03]=aj1p02i+aj1p03i;
			a[j1+p02]=aj1p02r-aj1p03r;	a[j2+p02]=aj1p02i-aj1p03i;

	#if PFETCH
	addr = add0+p21;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p04]=aj1p04r+aj1p05r;	a[j2+p04]=aj1p04i+aj1p05i;
			a[j1+p05]=aj1p04r-aj1p05r;	a[j2+p05]=aj1p04i-aj1p05i;

			a[j1+p07]=aj1p06r+aj1p07r;	a[j2+p07]=aj1p06i+aj1p07i;
			a[j1+p06]=aj1p06r-aj1p07r;	a[j2+p06]=aj1p06i-aj1p07i;

	#if PFETCH
	addr = add0+p22;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p08]=aj1p08r+aj1p09r;	a[j2+p08]=aj1p08i+aj1p09i;
			a[j1+p09]=aj1p08r-aj1p09r;	a[j2+p09]=aj1p08i-aj1p09i;

			a[j1+p11]=aj1p10r+aj1p11r;	a[j2+p11]=aj1p10i+aj1p11i;
			a[j1+p10]=aj1p10r-aj1p11r;	a[j2+p10]=aj1p10i-aj1p11i;

	#if PFETCH
	addr = add0+p23;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p12]=aj1p12r+aj1p13r;	a[j2+p12]=aj1p12i+aj1p13i;
			a[j1+p13]=aj1p12r-aj1p13r;	a[j2+p13]=aj1p12i-aj1p13i;

			a[j1+p15]=aj1p14r+aj1p15r;	a[j2+p15]=aj1p14i+aj1p15i;
			a[j1+p14]=aj1p14r-aj1p15r;	a[j2+p14]=aj1p14i-aj1p15i;

	#if PFETCH
	addr = add0+p24;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p16]=aj1p16r+aj1p17r;	a[j2+p16]=aj1p16i+aj1p17i;
			a[j1+p17]=aj1p16r-aj1p17r;	a[j2+p17]=aj1p16i-aj1p17i;

			a[j1+p19]=aj1p18r+aj1p19r;	a[j2+p19]=aj1p18i+aj1p19i;
			a[j1+p18]=aj1p18r-aj1p19r;	a[j2+p18]=aj1p18i-aj1p19i;

	#if PFETCH
	addr = add0+p25;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p20]=aj1p20r+aj1p21r;	a[j2+p20]=aj1p20i+aj1p21i;
			a[j1+p21]=aj1p20r-aj1p21r;	a[j2+p21]=aj1p20i-aj1p21i;

			a[j1+p23]=aj1p22r+aj1p23r;	a[j2+p23]=aj1p22i+aj1p23i;
			a[j1+p22]=aj1p22r-aj1p23r;	a[j2+p22]=aj1p22i-aj1p23i;

			a[j1+p24]=aj1p24r+aj1p25r;	a[j2+p24]=aj1p24i+aj1p25i;
			a[j1+p25]=aj1p24r-aj1p25r;	a[j2+p25]=aj1p24i-aj1p25i;
		}

		jstart += nwt;
		jhi    += nwt;
		col += 26;
		co3 -= 26;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-26 forward DIF FFT of the first block of 26 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 26 outputs of (1);
!   (3) Reweight and perform a radix-26 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 26 elements and repeat (1-4).
*/
	t00 = cy25;
	cy25= cy24;
	cy24= cy23;
	cy23= cy22;
	cy22= cy21;
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
	cy00= t00;

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
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
	}
}

	if(fabs(cy00)+fabs(cy01)+fabs(cy02)+fabs(cy03)+fabs(cy04)+fabs(cy05)+fabs(cy06)+fabs(cy07)+fabs(cy08)+fabs(cy09)+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17)+fabs(cy18)+fabs(cy19)+fabs(cy20)+fabs(cy21)+fabs(cy22)+fabs(cy23)+fabs(cy24)+fabs(cy25) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix26_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix26_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-26 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix13_dif_pass for details on the radix-13 subtransforms.
*/
	int j,j1,j2;
	static int n26,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25, first_entry=TRUE;
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
	double r1,i1,r2,i2,r3,i3,r4,i4
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
	,aj1p00r,aj1p01r,aj1p02r,aj1p03r,aj1p04r,aj1p05r,aj1p06r,aj1p07r,aj1p08r,aj1p09r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r,aj1p18r,aj1p19r,aj1p20r,aj1p21r,aj1p22r,aj1p23r,aj1p24r,aj1p25r
	,aj1p00i,aj1p01i,aj1p02i,aj1p03i,aj1p04i,aj1p05i,aj1p06i,aj1p07i,aj1p08i,aj1p09i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i,aj1p18i,aj1p19i,aj1p20i,aj1p21i,aj1p22i,aj1p23i,aj1p24i,aj1p25i;

	if(!first_entry && (n/26) != n26)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n26=n/26;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n26;
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
	}

/*...The radix-26 pass is here.	*/

	for(j=0; j < n26; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*** 10/07/2004: NOTE: I actually wound up using the DIT-style radix-13 subtransforms here,
                 but due to the prime radix, some simple index permuting makes that OK.    ***/
	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
					=> 0,13,24,11,22, 9,20, 7,18, 5,16, 3,14, 1,12,25,10,23, 8,21, 6,19, 4,17, 2,15 modulo 26.
		I.e. start out with first 13-vector of indices {0,2,4,6,8,10,12,14,16,18,20,22,24},
		permute those according to {0, 2, 4, 6, 8,10,12,14,16,18,20,22,24}*25%26 => {0,24,22,20,18,16,14,12,10, 8, 6, 4, 2},
		then each is head of a length-2 list of indices with decrement 13.
	*/
	/*...First radix-13 block uses aj1p[0:24:2] as inputs:	*/

		aj1p00r = a[j1    ];	aj1p00i = a[j2    ];

		t02 =a[j1+p24];		t03 =a[j2+p24];
		r1  =a[j1+p02];		i1  =a[j2+p02];
		t14 =t02 -r1;		t15 =t03 -i1;
		t02 =t02 +r1;		t03 =t03 +i1;

		t04 =a[j1+p22];		t05 =a[j2+p22];
		r1  =a[j1+p04];		i1  =a[j2+p04];
		t16 =t04 -r1;		t17 =t05 -i1;
		t04 =t04 +r1;		t05 =t05 +i1;

		t10 =a[j1+p20];		t11 =a[j2+p20];
		r1  =a[j1+p06];		i1  =a[j2+p06];
		t22 =t10 -r1;		t23 =t11 -i1;
		t10 =t10 +r1;		t11 =t11 +i1;

		t06 =a[j1+p18];		t07 =a[j2+p18];
		r1  =a[j1+p08];		i1  =a[j2+p08];
		t18 =t06 -r1;		t19 =t07 -i1;
		t06 =t06 +r1;		t07 =t07 +i1;

		t08 =a[j1+p16];		t09 =a[j2+p16];
		r1  =a[j1+p10];		i1  =a[j2+p10];
		t20 =r1 - t08;		t21 =i1 - t09;		/* flip signs on as3 */
		t08 =t08 +r1;		t09 =t09 +i1;

		t12 =a[j1+p14];		t13 =a[j2+p14];
		r1  =a[j1+p12];		i1  =a[j2+p12];
		t24 =t12 -r1;		t25 =t13 -i1;
		t12 =t12 +r1;		t13 =t13 +i1;

	/*** COSINE TERMS ***/

		t00 = t02+t06+t10;	t01 = t03+t07+t11;
		t02 = t02-t06   ;	t03 = t03-t07   ;
		t06 =     t06-t10;	t07 =     t07-t11;

		r1  = t04+t08+t12;	i1  = t05+t09+t13;
		t04 = t04    -t12;	t05 = t05    -t13;
		t08 =     t08-t12;	t09 =     t09-t13;

		t10 = t00+r1;		t11 = t01+i1;
		t00 = t00-r1;		t01 = t01-i1;
		t12 = t02-t08;		t13 = t03-t09;
		t02 = t02+t08;		t03 = t03+t09;
		t08 = t04+t06;		t09 = t05+t07;
		t06 = t04-t06;		t07 = t05-t07;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p00r;	i1 = t11*bp0 + aj1p00i;

		/* polypro(A,B) mod P1      */
		t00= t00*bp1;		t01= t01*bp1;

		aj1p00r += t10;		aj1p00i += t11;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t08 = t08*bp21;		t09 = t09*bp21;
		t08 = t12 - t08;	t09 = t13 - t09;
		t04 = t04 - t12;	t05 = t05 - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
		t02 = t02*bp30;		t03 = t03*bp30;
		t06 = t06*bp31;		t07 = t07*bp31;
		t06 = t02 - t06;	t07 = t03 - t07;
		t10 = t10 + t02;	t11 = t11 + t03;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t02= r1 +t00;		t03= i1 +t01;	/* a */
		t00= r1 -t00;		t01= i1 -t01;	/* b */
		r1 = t08+t06;		i1 = t09+t07;	/* c */
		t08= t08-t06;		t09= t09-t07;	/* d */
		t06= t04+t10;		t07= t05+t11;	/* e */
		t04= t04-t10;		t05= t05-t11;	/* f */
#if(LO_ADD)
		t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
		t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t02-r1;		t11 = t03-i1;
		t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
		t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
		t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

		r1  = t00+t08;		i1  = t01+t09;
		t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
		t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
		t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
#endif
	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t00= t16-t20+t24;	t01= t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t00 = t00*bq01;		t01 = t01*bq01;
		t24 = r1 - t00;		t25 = i1 - t01;
		t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t00 = t14*bq10;		t01 = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t00 - r2;	i2  = t01 - i2;
		t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3 - r4;		i4  = i3 - i4;
		r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t00+t16;	t21 = t21+t01+t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3;		t15 = t15+i3;

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

		aj1p02r=t12-t17;	aj1p02i=t13+t16;	/* C1 - S1 */
		aj1p04r=t10+t23;	aj1p04i=t11-t22;	/* C2 - S2 */
		aj1p06r=t02+t21;	aj1p06i=t03-t20;	/* C3 - S3 */
		aj1p08r=t04+t25;	aj1p08i=t05-t24;	/* C4 - S4 */
		aj1p10r=t08-t15;	aj1p10i=t09+t14;	/* C5 - S5 */
		aj1p12r=t06+t19;	aj1p12i=t07-t18;	/* C6 - S6 */
		aj1p14r=t06-t19;	aj1p14i=t07+t18;	/* C6 + S6 */
		aj1p16r=t08+t15;	aj1p16i=t09-t14;	/* C5 + S5 */
		aj1p18r=t04-t25;	aj1p18i=t05+t24;	/* C4 + S4 */
		aj1p20r=t02-t21;	aj1p20i=t03+t20;	/* C3 + S3 */
		aj1p22r=t10-t23;	aj1p22i=t11+t22;	/* C2 + S2 */
		aj1p24r=t12+t17;	aj1p24i=t13-t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t00 = t22+t18;		t01 = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
		t20 = t00+t18;		t21 = t01+t19;		/* S6 */
		t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		aj1p02r=t12-t15;	aj1p02i=t13+t14;	/* C1 - S1 */
		aj1p04r=t06+t23;	aj1p04i=t07-t22;	/* C2 - S2 */
		aj1p06r=t02+t17;	aj1p06i=t03-t16;	/* C3 - S3 */
		aj1p08r=t04+t25;	aj1p08i=t05-t24;	/* C4 - S4 */
		aj1p10r=t10-t19;	aj1p10i=t11+t18;	/* C5 - S5 */
		aj1p12r=t08+t21;	aj1p12i=t09-t20;	/* C6 - S6 */
		aj1p14r=t08-t21;	aj1p14i=t09+t20;	/* C6 + S6 */
		aj1p16r=t10+t19;	aj1p16i=t11-t18;	/* C5 + S5 */
		aj1p18r=t04-t25;	aj1p18i=t05+t24;	/* C4 + S4 */
		aj1p20r=t02-t17;	aj1p20i=t03+t16;	/* C3 + S3 */
		aj1p22r=t06-t23;	aj1p22i=t07+t22;	/* C2 + S2 */
		aj1p24r=t12+t15;	aj1p24i=t13-t14;	/* C1 + S1 */
#endif

	/*...Second radix-13 block uses aj1p[1:25:2] as inputs:	*/

		aj1p01r = a[j1+p13];	aj1p01i = a[j2+p13];

		t02 =a[j1+p11];		t03 =a[j2+p11];
		r1  =a[j1+p15];		i1  =a[j2+p15];
		t14 =t02 -r1;		t15 =t03 -i1;
		t02 =t02 +r1;		t03 =t03 +i1;

		t04 =a[j1+p09];		t05 =a[j2+p09];
		r1  =a[j1+p17];		i1  =a[j2+p17];
		t16 =t04 -r1;		t17 =t05 -i1;
		t04 =t04 +r1;		t05 =t05 +i1;

		t10 =a[j1+p07];		t11 =a[j2+p07];
		r1  =a[j1+p19];		i1  =a[j2+p19];
		t22 =t10 -r1;		t23 =t11 -i1;
		t10 =t10 +r1;		t11 =t11 +i1;

		t06 =a[j1+p05];		t07 =a[j2+p05];
		r1  =a[j1+p21];		i1  =a[j2+p21];
		t18 =t06 -r1;		t19 =t07 -i1;
		t06 =t06 +r1;		t07 =t07 +i1;

		t08 =a[j1+p03];		t09 =a[j2+p03];
		r1  =a[j1+p23];		i1  =a[j2+p23];
		t20 =r1 - t08;		t21 =i1 - t09;		/* flip signs on as3 */
		t08 =t08 +r1;		t09 =t09 +i1;

		t12 =a[j1+p01];		t13 =a[j2+p01];
		r1  =a[j1+p25];		i1  =a[j2+p25];
		t24 =t12 -r1;		t25 =t13 -i1;
		t12 =t12 +r1;		t13 =t13 +i1;

	/*** COSINE TERMS ***/

		t00 = t02+t06+t10;	t01 = t03+t07+t11;
		t02 = t02-t06   ;	t03 = t03-t07   ;
		t06 =     t06-t10;	t07 =     t07-t11;

		r1  = t04+t08+t12;	i1  = t05+t09+t13;
		t04 = t04    -t12;	t05 = t05    -t13;
		t08 =     t08-t12;	t09 =     t09-t13;

		t10 = t00+r1;		t11 = t01+i1;
		t00 = t00-r1;		t01 = t01-i1;
		t12 = t02-t08;		t13 = t03-t09;
		t02 = t02+t08;		t03 = t03+t09;
		t08 = t04+t06;		t09 = t05+t07;
		t06 = t04-t06;		t07 = t05-t07;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p01r;	i1 = t11*bp0 + aj1p01i;

		/* polypro(A,B) mod P1      */
		t00= t00*bp1;		t01= t01*bp1;

		aj1p01r += t10;		aj1p01i += t11;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t08 = t08*bp21;		t09 = t09*bp21;
		t08 = t12 - t08;	t09 = t13 - t09;
		t04 = t04 - t12;	t05 = t05 - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
		t02 = t02*bp30;		t03 = t03*bp30;
		t06 = t06*bp31;		t07 = t07*bp31;
		t06 = t02 - t06;	t07 = t03 - t07;
		t10 = t10 + t02;	t11 = t11 + t03;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t02= r1 +t00;		t03= i1 +t01;	/* a */
		t00= r1 -t00;		t01= i1 -t01;	/* b */
		r1 = t08+t06;		i1 = t09+t07;	/* c */
		t08= t08-t06;		t09= t09-t07;	/* d */
		t06= t04+t10;		t07= t05+t11;	/* e */
		t04= t04-t10;		t05= t05-t11;	/* f */
#if(LO_ADD)
		t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
		t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t02-r1;		t11 = t03-i1;
		t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
		t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
		t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

		r1  = t00+t08;		i1  = t01+t09;
		t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
		t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
		t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
#endif
	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t00= t16-t20+t24;	t01= t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t00 = t00*bq01;		t01 = t01*bq01;
		t24 = r1 - t00;		t25 = i1 - t01;
		t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t00 = t14*bq10;		t01 = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t00 - r2;	i2  = t01 - i2;
		t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3 - r4;		i4  = i3 - i4;
		r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t00+t16;	t21 = t21+t01+t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3;		t15 = t15+i3;

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

		aj1p03r=t12-t17;	aj1p03i=t13+t16;	/* C1 - S1 */
		aj1p05r=t10+t23;	aj1p05i=t11-t22;	/* C2 - S2 */
		aj1p07r=t02+t21;	aj1p07i=t03-t20;	/* C3 - S3 */
		aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
		aj1p11r=t08-t15;	aj1p11i=t09+t14;	/* C5 - S5 */
		aj1p13r=t06+t19;	aj1p13i=t07-t18;	/* C6 - S6 */
		aj1p15r=t06-t19;	aj1p15i=t07+t18;	/* C6 + S6 */
		aj1p17r=t08+t15;	aj1p17i=t09-t14;	/* C5 + S5 */
		aj1p19r=t04-t25;	aj1p19i=t05+t24;	/* C4 + S4 */
		aj1p21r=t02-t21;	aj1p21i=t03+t20;	/* C3 + S3 */
		aj1p23r=t10-t23;	aj1p23i=t11+t22;	/* C2 + S2 */
		aj1p25r=t12+t17;	aj1p25i=t13-t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t00 = t22+t18;		t01 = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
		t20 = t00+t18;		t21 = t01+t19;		/* S6 */
		t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		aj1p03r=t12-t15;	aj1p03i=t13+t14;	/* C1 - S1 */
		aj1p05r=t06+t23;	aj1p05i=t07-t22;	/* C2 - S2 */
		aj1p07r=t02+t17;	aj1p07i=t03-t16;	/* C3 - S3 */
		aj1p09r=t04+t25;	aj1p09i=t05-t24;	/* C4 - S4 */
		aj1p11r=t10-t19;	aj1p11i=t11+t18;	/* C5 - S5 */
		aj1p13r=t08+t21;	aj1p13i=t09-t20;	/* C6 - S6 */
		aj1p15r=t08-t21;	aj1p15i=t09+t20;	/* C6 + S6 */
		aj1p17r=t10+t19;	aj1p17i=t11-t18;	/* C5 + S5 */
		aj1p19r=t04-t25;	aj1p19i=t05+t24;	/* C4 + S4 */
		aj1p21r=t02-t17;	aj1p21i=t03+t16;	/* C3 + S3 */
		aj1p23r=t06-t23;	aj1p23i=t07+t22;	/* C2 + S2 */
		aj1p25r=t12+t15;	aj1p25i=t13-t14;	/* C1 + S1 */
#endif

/*...and now do 13 radix-2 transforms:	*/

		a[j1    ]=aj1p00r+aj1p01r;	a[j2    ]=aj1p00i+aj1p01i;
		a[j1+p01]=aj1p00r-aj1p01r;	a[j2+p01]=aj1p00i-aj1p01i;

		a[j1+p03]=aj1p02r+aj1p03r;	a[j2+p03]=aj1p02i+aj1p03i;
		a[j1+p02]=aj1p02r-aj1p03r;	a[j2+p02]=aj1p02i-aj1p03i;

		a[j1+p04]=aj1p04r+aj1p05r;	a[j2+p04]=aj1p04i+aj1p05i;
		a[j1+p05]=aj1p04r-aj1p05r;	a[j2+p05]=aj1p04i-aj1p05i;

		a[j1+p07]=aj1p06r+aj1p07r;	a[j2+p07]=aj1p06i+aj1p07i;
		a[j1+p06]=aj1p06r-aj1p07r;	a[j2+p06]=aj1p06i-aj1p07i;

		a[j1+p08]=aj1p08r+aj1p09r;	a[j2+p08]=aj1p08i+aj1p09i;
		a[j1+p09]=aj1p08r-aj1p09r;	a[j2+p09]=aj1p08i-aj1p09i;

		a[j1+p11]=aj1p10r+aj1p11r;	a[j2+p11]=aj1p10i+aj1p11i;
		a[j1+p10]=aj1p10r-aj1p11r;	a[j2+p10]=aj1p10i-aj1p11i;

		a[j1+p12]=aj1p12r+aj1p13r;	a[j2+p12]=aj1p12i+aj1p13i;
		a[j1+p13]=aj1p12r-aj1p13r;	a[j2+p13]=aj1p12i-aj1p13i;

		a[j1+p15]=aj1p14r+aj1p15r;	a[j2+p15]=aj1p14i+aj1p15i;
		a[j1+p14]=aj1p14r-aj1p15r;	a[j2+p14]=aj1p14i-aj1p15i;

		a[j1+p16]=aj1p16r+aj1p17r;	a[j2+p16]=aj1p16i+aj1p17i;
		a[j1+p17]=aj1p16r-aj1p17r;	a[j2+p17]=aj1p16i-aj1p17i;

		a[j1+p19]=aj1p18r+aj1p19r;	a[j2+p19]=aj1p18i+aj1p19i;
		a[j1+p18]=aj1p18r-aj1p19r;	a[j2+p18]=aj1p18i-aj1p19i;

		a[j1+p20]=aj1p20r+aj1p21r;	a[j2+p20]=aj1p20i+aj1p21i;
		a[j1+p21]=aj1p20r-aj1p21r;	a[j2+p21]=aj1p20i-aj1p21i;

		a[j1+p23]=aj1p22r+aj1p23r;	a[j2+p23]=aj1p22i+aj1p23i;
		a[j1+p22]=aj1p22r-aj1p23r;	a[j2+p22]=aj1p22i-aj1p23i;

		a[j1+p24]=aj1p24r+aj1p25r;	a[j2+p24]=aj1p24i+aj1p25i;
		a[j1+p25]=aj1p24r-aj1p25r;	a[j2+p25]=aj1p24i-aj1p25i;

	}
}

/*
   0 (   0)  134.0000000000   118.0000000000   134.0000000000   118.0000000000
   1 (  13)  -26.0000000000   -10.0000000000   -26.0000000000   -10.0000000000
   2 (   1)   24.1310186456    15.4885860970    -6.4228791824     6.7397899092 ERR=   30.5538978279     8.7487961878
   3 (  14)   -6.4228791824     6.7397899092    24.1310186456    15.4885860970 ERR=   30.5538978279     8.7487961878
   4 (   2)  -26.0281658104   -17.2978161275   -26.0281658104   -17.2978161275
   5 (  15)   11.8690859868    -4.2511070234    11.8690859868    -4.2511070234
   6 (   3)  -13.8203590669   -16.8227432208    -6.0678817445   -18.4447434417 ERR=    7.7524773224     1.6220002209
   7 (  16)   -6.0678817445   -18.4447434417   -13.8203590669   -16.8227432208 ERR=    7.7524773224     1.6220002209
   8 (   4)   -2.8850529228    -1.4262732841    -2.8850529228    -1.4262732841
   9 (  17)  -38.5234312939    -4.6693525950   -38.5234312939    -4.6693525950
  10 (   5)  -13.3827842605     6.5589444051    -1.1540312470     8.7174502227 ERR=   12.2287530135     2.1585058177
  11 (  18)   -1.1540312470     8.7174502227   -13.3827842605     6.5589444051 ERR=   12.2287530135     2.1585058177
  12 (   6)   -8.9255627696     9.7676714182    -8.9255627696     9.7676714182
  13 (  19)   -6.2855442440     5.5760304149    -6.2855442440     5.5760304149
  14 (   7)   10.6801741389   -21.7505156543     4.8721379251   -16.6147035842 ERR=    5.8080362137     5.1358120701
  15 (  20)    4.8721379251   -16.6147035842    10.6801741389   -21.7505156543 ERR=    5.8080362137     5.1358120701
  16 (   8)   11.2605453848    -5.4557198696    11.2605453848    -5.4557198696
  17 (  21)    5.6647274160   -15.9879453528     5.6647274160   -15.9879453528
  18 (   9)    5.5156137849   -10.4641006451    18.1618857211    25.6068361388 ERR=   12.6462719361    36.0709367840
  19 (  22)   18.1618857211    25.6068361388     5.5156137849   -10.4641006451 ERR=   12.6462719361    36.0709367840
  20 (  10)   -1.8247176531   -11.5883567530    -1.8247176531   -11.5883567530
  21 (  23)  -10.1082280090   -16.0514483057   -10.1082280090   -16.0514483057
  22 (  11)  -10.5680463509    -4.8225641719     3.3919125554    10.9547240474 ERR=   13.9599589063    15.7772882193
  23 (  24)    3.3919125554    10.9547240474   -10.5680463509    -4.8225641719 ERR=   13.9599589063    15.7772882193
  24 (  12)   34.8473368801    -8.1871121941    34.8473368801    -8.1871121941
  25 (  25)  -14.3977538837     2.4244695699   -14.3977538837     2.4244695699
*/

/***************/

void radix26_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-26 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix13_dif_pass for details on the radix-13 subtransforms.
*/
	int j,j1,j2;
	static int n26,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25, first_entry=TRUE;
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
	double rt,it,r1,i1,r2,i2,r3,i3,r4,i4
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
	,aj1p00r,aj1p01r,aj1p02r,aj1p03r,aj1p04r,aj1p05r,aj1p06r,aj1p07r,aj1p08r,aj1p09r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r,aj1p18r,aj1p19r,aj1p20r,aj1p21r,aj1p22r,aj1p23r,aj1p24r,aj1p25r
	,aj1p00i,aj1p01i,aj1p02i,aj1p03i,aj1p04i,aj1p05i,aj1p06i,aj1p07i,aj1p08i,aj1p09i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i,aj1p18i,aj1p19i,aj1p20i,aj1p21i,aj1p22i,aj1p23i,aj1p24i,aj1p25i;

	if(!first_entry && (n/26) != n26)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n26=n/26;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n26;
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
	}

/*...The radix-26 pass is here.	*/

	for(j=0; j < n26; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (26 64-bit complex, i.e. 52 64-bit reals) and do 13 radix-2 transforms,	*/
	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
					=> 0,24,22,20,18,16,14,12,10, 8, 6, 4, 2,13,11, 9, 7, 5, 3, 1,25,23,21,19,17,15 modulo 26.
		I.e. start out with first pair of indices {0,13}, permute those according to
		{0,13}*25%26 = {0,13}, then each is head of a length-13 list of indices with decrement 2.

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25] contain
		x[0,13, 1,14, 2,15, 3,16, 4,17, 5,18, 6,19, 7,20, 8,21, 9,22,10,23,11,24,12,25], which get swapped to
		x[0,13,24,11,22, 9,20, 7,18, 5,16, 3,14, 1,12,25,10,23, 8,21, 6,19, 4,17, 2,15], which means the a-indices get swapped as
		a[0, 1,23,22,19,18,15,14,11,10, 7, 6, 3, 2,24,25,20,21,16,17,12,13, 8, 9, 4, 5].
	*/
		aj1p00r = a[j1      ];		aj1p00i = a[j2    ];	/* x0a,b */
		rt      = a[j1  +p01];		it      = a[j2+p01];
		aj1p01r = aj1p00r-rt;		aj1p01i = aj1p00i-it;
		aj1p00r = aj1p00r+rt;		aj1p00i = aj1p00i+it;

		aj1p02r = a[j1  +p23];		aj1p02i = a[j2+p23];	/* x1a,b */
		rt      = a[j1  +p22];		it      = a[j2+p22];
		aj1p03r = aj1p02r-rt;		aj1p03i = aj1p02i-it;
		aj1p02r = aj1p02r+rt;		aj1p02i = aj1p02i+it;

		aj1p04r = a[j1  +p19];		aj1p04i = a[j2+p19];	/* x2a,b */
		rt      = a[j1  +p18];		it      = a[j2+p18];
		aj1p05r = aj1p04r-rt;		aj1p05i = aj1p04i-it;
		aj1p04r = aj1p04r+rt;		aj1p04i = aj1p04i+it;

		aj1p06r = a[j1  +p15];		aj1p06i = a[j2+p15];	/* x3a,b */
		rt      = a[j1  +p14];		it      = a[j2+p14];
		aj1p07r = aj1p06r-rt;		aj1p07i = aj1p06i-it;
		aj1p06r = aj1p06r+rt;		aj1p06i = aj1p06i+it;

		aj1p08r = a[j1  +p11];		aj1p08i = a[j2+p11];	/* x4a,b */
		rt      = a[j1  +p10];		it      = a[j2+p10];
		aj1p09r = aj1p08r-rt;		aj1p09i = aj1p08i-it;
		aj1p08r = aj1p08r+rt;		aj1p08i = aj1p08i+it;

		aj1p10r = a[j1  +p07];		aj1p10i = a[j2+p07];	/* x5a,b */
		rt      = a[j1  +p06];		it      = a[j2+p06];
		aj1p11r = aj1p10r-rt;		aj1p11i = aj1p10i-it;
		aj1p10r = aj1p10r+rt;		aj1p10i = aj1p10i+it;

		aj1p12r = a[j1  +p03];		aj1p12i = a[j2+p03];	/* x6a,b */
		rt      = a[j1  +p02];		it      = a[j2+p02];
		aj1p13r = aj1p12r-rt;		aj1p13i = aj1p12i-it;
		aj1p12r = aj1p12r+rt;		aj1p12i = aj1p12i+it;

		aj1p14r = a[j1  +p24];		aj1p14i = a[j2+p24];	/* x7a,b */
		rt      = a[j1  +p25];		it      = a[j2+p25];
		aj1p15r = aj1p14r-rt;		aj1p15i = aj1p14i-it;
		aj1p14r = aj1p14r+rt;		aj1p14i = aj1p14i+it;

		aj1p16r = a[j1  +p20];		aj1p16i = a[j2+p20];	/* x8a,b */
		rt      = a[j1  +p21];		it      = a[j2+p21];
		aj1p17r = aj1p16r-rt;		aj1p17i = aj1p16i-it;
		aj1p16r = aj1p16r+rt;		aj1p16i = aj1p16i+it;

		aj1p18r = a[j1  +p16];		aj1p18i = a[j2+p16];	/* x9a,b */
		rt      = a[j1  +p17];		it      = a[j2+p17];
		aj1p19r = aj1p18r-rt;		aj1p19i = aj1p18i-it;
		aj1p18r = aj1p18r+rt;		aj1p18i = aj1p18i+it;

		aj1p20r = a[j1  +p12];		aj1p20i = a[j2+p12];	/* x10a,b */
		rt      = a[j1  +p13];		it      = a[j2+p13];
		aj1p21r = aj1p20r-rt;		aj1p21i = aj1p20i-it;
		aj1p20r = aj1p20r+rt;		aj1p20i = aj1p20i+it;

		aj1p22r = a[j1  +p08];		aj1p22i = a[j2+p08];	/* x11a,b */
		rt      = a[j1  +p09];		it      = a[j2+p09];
		aj1p23r = aj1p22r-rt;		aj1p23i = aj1p22i-it;
		aj1p22r = aj1p22r+rt;		aj1p22i = aj1p22i+it;

		aj1p24r = a[j1  +p04];		aj1p24i = a[j2+p04];	/* x12a,b */
		rt      = a[j1  +p05];		it      = a[j2+p05];
		aj1p25r = aj1p24r-rt;		aj1p25i = aj1p24i-it;
		aj1p24r = aj1p24r+rt;		aj1p24i = aj1p24i+it;

/*      ..and now do two radix-13 transforms.	*/

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		t02=aj1p02r+aj1p24r;	t03=aj1p02i+aj1p24i;
		t14=aj1p02r-aj1p24r;	t15=aj1p02i-aj1p24i;

		t04=aj1p04r+aj1p22r;	t05=aj1p04i+aj1p22i;
		t16=aj1p04r-aj1p22r;	t17=aj1p04i-aj1p22i;

		t10=aj1p06r+aj1p20r;	t11=aj1p06i+aj1p20i;
		t22=aj1p06r-aj1p20r;	t23=aj1p06i-aj1p20i;

		t06=aj1p08r+aj1p18r;	t07=aj1p08i+aj1p18i;
		t18=aj1p08r-aj1p18r;	t19=aj1p08i-aj1p18i;

		t08=aj1p10r+aj1p16r;	t09=aj1p10i+aj1p16i;
		t20=aj1p16r-aj1p10r;	t21=aj1p16i-aj1p10i;		/* flip signs on as3 */

		t12=aj1p12r+aj1p14r;	t13=aj1p12i+aj1p14i;
		t24=aj1p12r-aj1p14r;	t25=aj1p12i-aj1p14i;

	/*** COSINE TERMS ***/

		t00 = t02+t06+t10;	t01 = t03+t07+t11;
		t02 = t02-t06   ;	t03 = t03-t07   ;
		t06 =     t06-t10;	t07 =     t07-t11;

		r1  = t04+t08+t12;	i1  = t05+t09+t13;
		t04 = t04    -t12;	t05 = t05    -t13;
		t08 =     t08-t12;	t09 =     t09-t13;

		t10 = t00+r1;		t11 = t01+i1;
		t00 = t00-r1;		t01 = t01-i1;
		t12 = t02-t08;		t13 = t03-t09;
		t02 = t02+t08;		t03 = t03+t09;
		t08 = t04+t06;		t09 = t05+t07;
		t06 = t04-t06;		t07 = t05-t07;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p00r;	i1 = t11*bp0 + aj1p00i;

		/* polypro(A,B) mod P1      */
		t00= t00*bp1;		t01= t01*bp1;

		aj1p00r += t10;		aj1p00i += t11;
		a[j1  ] = aj1p00r;	a[j2] = aj1p00i;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t08 = t08*bp21;		t09 = t09*bp21;
		t08 = t12 - t08;	t09 = t13 - t09;
		t04 = t04 - t12;	t05 = t05 - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
		t02 = t02*bp30;		t03 = t03*bp30;
		t06 = t06*bp31;		t07 = t07*bp31;
		t06 = t02 - t06;	t07 = t03 - t07;
		t10 = t10 + t02;	t11 = t11 + t03;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t02= r1 +t00;		t03= i1 +t01;	/* a */
		t00= r1 -t00;		t01= i1 -t01;	/* b */
		r1 = t08+t06;		i1 = t09+t07;	/* c */
		t08= t08-t06;		t09= t09-t07;	/* d */
		t06= t04+t10;		t07= t05+t11;	/* e */
		t04= t04-t10;		t05= t05-t11;	/* f */
#if(LO_ADD)
		t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
		t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t02-r1;		t11 = t03-i1;
		t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
		t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
		t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

		r1  = t00+t08;		i1  = t01+t09;
		t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
		t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
		t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
#endif
	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t00= t16-t20+t24;	t01= t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t00 = t00*bq01;		t01 = t01*bq01;
		t24 = r1 - t00;		t25 = i1 - t01;
		t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t00 = t14*bq10;		t01 = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t00 - r2;	i2  = t01 - i2;
		t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3 - r4;		i4  = i3 - i4;
		r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t00+t16;	t21 = t21+t01+t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3;		t15 = t15+i3;

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

		a[j1+p12]=t12-t17;	a[j2+p12]=t13+t16;	/* C1 - S1 */
		a[j1+p24]=t10+t23;	a[j2+p24]=t11-t22;	/* C2 - S2 */
		a[j1+p10]=t02+t21;	a[j2+p10]=t03-t20;	/* C3 - S3 */
		a[j1+p22]=t04+t25;	a[j2+p22]=t05-t24;	/* C4 - S4 */
		a[j1+p08]=t08-t15;	a[j2+p08]=t09+t14;	/* C5 - S5 */
		a[j1+p20]=t06+t19;	a[j2+p20]=t07-t18;	/* C6 - S6 */
		a[j1+p06]=t06-t19;	a[j2+p06]=t07+t18;	/* C6 + S6 */
		a[j1+p18]=t08+t15;	a[j2+p18]=t09-t14;	/* C5 + S5 */
		a[j1+p04]=t04-t25;	a[j2+p04]=t05+t24;	/* C4 + S4 */
		a[j1+p16]=t02-t21;	a[j2+p16]=t03+t20;	/* C3 + S3 */
		a[j1+p02]=t10-t23;	a[j2+p02]=t11+t22;	/* C2 + S2 */
		a[j1+p14]=t12+t17;	a[j2+p14]=t13-t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t00 = t22+t18;		t01 = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
		t20 = t00+t18;		t21 = t01+t19;		/* S6 */
		t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

	/*...Inline multiply of sine parts by +-I ]nto finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p12]=t12-t15;	a[j2+p12]=t13+t14;	/* C1 - S1 */
		a[j1+p24]=t06+t23;	a[j2+p24]=t07-t22;	/* C2 - S2 */
		a[j1+p10]=t02+t17;	a[j2+p10]=t03-t16;	/* C3 - S3 */
		a[j1+p22]=t04+t25;	a[j2+p22]=t05-t24;	/* C4 - S4 */
		a[j1+p08]=t10-t19;	a[j2+p08]=t11+t18;	/* C5 - S5 */
		a[j1+p20]=t08+t21;	a[j2+p20]=t09-t20;	/* C6 - S6 */
		a[j1+p06]=t08-t21;	a[j2+p06]=t09+t20;	/* C6 + S6 */
		a[j1+p18]=t10+t19;	a[j2+p18]=t11-t18;	/* C5 + S5 */
		a[j1+p04]=t04-t25;	a[j2+p04]=t05+t24;	/* C4 + S4 */
		a[j1+p16]=t02-t17;	a[j2+p16]=t03+t16;	/* C3 + S3 */
		a[j1+p02]=t06-t23;	a[j2+p02]=t07+t22;	/* C2 + S2 */
		a[j1+p14]=t12+t15;	a[j2+p14]=t13-t14;	/* C1 + S1 */
#endif

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		t02=aj1p03r+aj1p25r;	t03=aj1p03i+aj1p25i;
		t14=aj1p03r-aj1p25r;	t15=aj1p03i-aj1p25i;

		t04=aj1p05r+aj1p23r;	t05=aj1p05i+aj1p23i;
		t16=aj1p05r-aj1p23r;	t17=aj1p05i-aj1p23i;

		t10=aj1p07r+aj1p21r;	t11=aj1p07i+aj1p21i;
		t22=aj1p07r-aj1p21r;	t23=aj1p07i-aj1p21i;

		t06=aj1p09r+aj1p19r;	t07=aj1p09i+aj1p19i;
		t18=aj1p09r-aj1p19r;	t19=aj1p09i-aj1p19i;

		t08=aj1p11r+aj1p17r;	t09=aj1p11i+aj1p17i;
		t20=aj1p17r-aj1p11r;	t21=aj1p17i-aj1p11i;		/* flip signs on as3 */

		t12=aj1p13r+aj1p15r;	t13=aj1p13i+aj1p15i;
		t24=aj1p13r-aj1p15r;	t25=aj1p13i-aj1p15i;

	/*** COSINE TERMS ***/

		t00 = t02+t06+t10;	t01 = t03+t07+t11;
		t02 = t02-t06   ;	t03 = t03-t07   ;
		t06 =     t06-t10;	t07 =     t07-t11;

		r1  = t04+t08+t12;	i1  = t05+t09+t13;
		t04 = t04    -t12;	t05 = t05    -t13;
		t08 =     t08-t12;	t09 =     t09-t13;

		t10 = t00+r1;		t11 = t01+i1;
		t00 = t00-r1;		t01 = t01-i1;
		t12 = t02-t08;		t13 = t03-t09;
		t02 = t02+t08;		t03 = t03+t09;
		t08 = t04+t06;		t09 = t05+t07;
		t06 = t04-t06;		t07 = t05-t07;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p01r;	i1 = t11*bp0 + aj1p01i;

		/* polypro(A,B) mod P1      */
		t00= t00*bp1;		t01= t01*bp1;

		aj1p01r += t10;		aj1p01i += t11;
		a[j1+p13] = aj1p01r;	a[j2+p13] = aj1p01i;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t04 =(t12+t08)*bp2B;	t05 =(t13+t09)*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t08 = t08*bp21;		t09 = t09*bp21;
		t08 = t12 - t08;	t09 = t13 - t09;
		t04 = t04 - t12;	t05 = t05 - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t02-t06)*bp3B;	t11 =(t03-t07)*bp3B;
		t02 = t02*bp30;		t03 = t03*bp30;
		t06 = t06*bp31;		t07 = t07*bp31;
		t06 = t02 - t06;	t07 = t03 - t07;
		t10 = t10 + t02;	t11 = t11 + t03;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t02= r1 +t00;		t03= i1 +t01;	/* a */
		t00= r1 -t00;		t01= i1 -t01;	/* b */
		r1 = t08+t06;		i1 = t09+t07;	/* c */
		t08= t08-t06;		t09= t09-t07;	/* d */
		t06= t04+t10;		t07= t05+t11;	/* e */
		t04= t04-t10;		t05= t05-t11;	/* f */
#if(LO_ADD)
		t02= t02- r1 + t04;	t03= t03-i1 + t05;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t02+ 3*r1;		t13= t03+ 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t04= t02- 3*t04;	t05= t03- 3*t05;	/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t00+ t08- t06;	t11= t01+ t09- t07;	/* C2 = x^5 = b + d - e		2 add        */
		t06= t10+ 3*t06;	t07= t11+ 3*t07;	/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t08= t10- 3*t08;	t09= t11- 3*t09;	/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t02-r1;		t11 = t03-i1;
		t12 = t02+r1 +r1 +t04;	t13 = t03+i1 +i1 +t05;	/* C1 */
		t02 = t10        +t04;	t03 = t11        +t05;	/* C3 */
		t04 = t10    -t04-t04;	t05 = t11    -t05-t05;	/* C4 */

		r1  = t00+t08;		i1  = t01+t09;
		t10 = t00-t08-t08-t06;	t11 = t01-t09-t09-t07;	/* C5 */
		t08 = r1     +t06+t06;	t09 = i1     +t07+t07;	/* C6 */
		t06 = r1         -t06;	t07 = i1         -t07;	/* C2 */
#endif
	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t00= t16-t20+t24;	t01= t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t00)*bq0B;	t23 =(i1 +t01)*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t00 = t00*bq01;		t01 = t01*bq01;
		t24 = r1 - t00;		t25 = i1 - t01;
		t22 = t22 - r1 - t00;	t23 = t23 - i1 - t01;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t00 = t14*bq10;		t01 = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t00 - r2;	i2  = t01 - i2;
		t00 = r1 *bq1A - t00;	t01 = i1 *bq1A - t01;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3 - r4;		i4  = i3 - i4;
		r3  = r1 *bq1A - r3;	i3  = i1 *bq1A - i3;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t00+t16;	t21 = t21+t01+t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3;		t15 = t15+i3;

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

		a[j1+p25]=t12-t17;	a[j2+p25]=t13+t16;	/* C1 - S1 */
		a[j1+p11]=t10+t23;	a[j2+p11]=t11-t22;	/* C2 - S2 */
		a[j1+p23]=t02+t21;	a[j2+p23]=t03-t20;	/* C3 - S3 */
		a[j1+p09]=t04+t25;	a[j2+p09]=t05-t24;	/* C4 - S4 */
		a[j1+p21]=t08-t15;	a[j2+p21]=t09+t14;	/* C5 - S5 */
		a[j1+p07]=t06+t19;	a[j2+p07]=t07-t18;	/* C6 - S6 */
		a[j1+p19]=t06-t19;	a[j2+p19]=t07+t18;	/* C6 + S6 */
		a[j1+p05]=t08+t15;	a[j2+p05]=t09-t14;	/* C5 + S5 */
		a[j1+p17]=t04-t25;	a[j2+p17]=t05+t24;	/* C4 + S4 */
		a[j1+p03]=t02-t21;	a[j2+p03]=t03+t20;	/* C3 + S3 */
		a[j1+p15]=t10-t23;	a[j2+p15]=t11+t22;	/* C2 + S2 */
		a[j1+p01]=t12+t17;	a[j2+p01]=t13-t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t00 = t22+t18;		t01 = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1;	t17 = t17+t17-i1;	/* S3 */
		t20 = t00+t18;		t21 = t01+t19;		/* S6 */
		t18 = t18+t18-t00;	t19 = t19+t19-t01;	/* S5 */

	/*...Inline multiply of sine parts by +-I ]nto finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p25]=t12-t15;	a[j2+p25]=t13+t14;	/* C1 - S1 */
		a[j1+p11]=t06+t23;	a[j2+p11]=t07-t22;	/* C2 - S2 */
		a[j1+p23]=t02+t17;	a[j2+p23]=t03-t16;	/* C3 - S3 */
		a[j1+p09]=t04+t25;	a[j2+p09]=t05-t24;	/* C4 - S4 */
		a[j1+p21]=t10-t19;	a[j2+p21]=t11+t18;	/* C5 - S5 */
		a[j1+p07]=t08+t21;	a[j2+p07]=t09-t20;	/* C6 - S6 */
		a[j1+p19]=t08-t21;	a[j2+p19]=t09+t20;	/* C6 + S6 */
		a[j1+p05]=t10+t19;	a[j2+p05]=t11-t18;	/* C5 + S5 */
		a[j1+p17]=t04-t25;	a[j2+p17]=t05+t24;	/* C4 + S4 */
		a[j1+p03]=t02-t17;	a[j2+p03]=t03+t16;	/* C3 + S3 */
		a[j1+p15]=t06-t23;	a[j2+p15]=t07+t22;	/* C2 + S2 */
		a[j1+p01]=t12+t15;	a[j2+p01]=t13-t14;	/* C1 + S1 */
#endif
	}
}

