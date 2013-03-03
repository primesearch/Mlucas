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
//#include "radix31.h"

/*
hanley wood: 1001 all time best selling home plans
hanley wood: 225 hillside home plans
chiras: the new ecological home
mcraven: building with stone
scutella / heberle: how to plan, contract and build your own home

Math Books donated:

bolloba's, modern graph theory
bolloba's, random graphs, 2nd ed
cercignani: the boltzmann equation and its applications
chorin / marsden, a mathematical introduction to fluid mechanics, 2nd ed
chorin / marsden, a mathematical introduction to fluid mechanics, 3rd ed
chung, spectral graph theory
courant / friedrichs: supersonic flow and shock waves
frisch, turbulence
gottlieb / orszag, numerical analysis of spectral methods
guckenheimer / holmes, nonlinear oscillations, dynamical systems, and bifurcations of vector fields
kevorkian / cole: perturbation methods in applied mechanics
lagerstrom, matched asymptotic expansions: ideas and techniques
reed / simon: vols 1-4
smoller, shock waves and reaction-diffusion equations
spekreijse, multigrid solution of the steady euler equations
swinney / gollub (eds), hydrodynamic instabilties and the transition to turbulence
weyl: symmetry
gardner: codes, ciphers and secret writing


Van Buskirk's original constants in left column, EWM-derived ones at right.
The derived ones arise from several common code modifications used to effect in-place
versions of certain code subsequences, such as:

[1] 2x2 const-matrix x vector in-place mul: Wish to compute the following using just 2 registers
which contain the variables x and y on input:

	/ a  b \   | x |   a*x + b*y
	|      | * |   | =
	\ c  d /   | y |   c*x + d*y ,

where a,b,c,d are precomputed constants. Rewrite the above results as

	a*[ x + (b/a)*y ]

	c*[ x + (d/c)*y ] = c*[ x + (d/c - b/a + b/a)*y ] = c*[ x + (b/a)*y + (d/c - b/a)*y]

then compute

	y *= (b/a)
	x += y
	y *= (a*d)/(b*c) - 1
	y += x
	x *= a
	y *= c

which has the same number of arithmetic operations (4 mul, 2 add) as the original sequence
and requires precomputation of 2 added constants:

	B := (b/a)
	D := (a*d)/(b*c) - 1 .

In KVB's radix-31 DFT implementation we encounter the above sequence for the following values
of a,b,c,d drawn from the left-column constants:

	dc3,4,8,3: call auxiliary consts B_c3483 and D_c3483
	dc5,6,-4,5:    "          "      B_c5645 and D_c5645
	dc9,-10,11,-9: "          "      B_c9019 and D_c9019

	dc15,12,18,-15:  "        "      B_c5285 and D_c5285
	dc16,-13,19,16:  "        "      B_c6396 and D_c6396
	dc17,14,20,-17:  "        "      B_c7407 and D_c7407

	dc24,21,27,-24:  "        "      C_c4174 and E_c4174 (Use C,E to avoid name collision with _c5285 and _c6396 terms above)
	dc25,-22,28,25:  "        "      C_c5285 and E_c5285
	dc26,23,29,-26:  "        "      C_c6396 and E_c6396

...and there is an analogous set of terms based on the sine (ds*) components, named [BB,C,D,E]_s****.

Another common mini-sequence is the Low-register-usage computation of

	a + b
	a - c
	a - b + c

If have 1 spare register x, can do with just 4 add/sub, with results in a,b,c:

	x = c
	c += a	// a + c
	a -= x	// a - c
	x -= c	// -a
	c -= b	// a - b + c
	b -= x	// a + b

If have 0 spare registers, best I have found so far needs 5 add/sub, 3 mul:

	a -= c	// a - c
	c *= 2	// 2*c
	c += a	// a + c
	c -= b	// a - b + c
	b *= 3	// 3*b
	b += c	// a + 2*b + c
	b += a	// 2*(a + b)
	b *= .5	// a + b

*/
#define MAT_2x2_MUL(__a,__B,__c,__D, __x,__y)\
{\
	__y *= __B;\
	__x += __y;\
	__y *= __D;\
	__y += __x;\
	__x *= __a;\
	__y *= __c;\
}

const double dc0    =-0.16666666666666666666666666666666700;	const double dc0_5th =-0.033333333333333333333333333333333333;
const double dc1    =-0.41848031916582736768981687098748646;	/* 046,1[0234],2[123] */
const double dc2    = 1.03708005397391977997160126838915860;	const double B_c3483 =-10.31938532050192055740019444;
const double dc3    = 0.06899615305609372142805240007781954;	const double D_c3483 =-1.039447157476621531648472621;
const double dc4    =-0.71199788901815727304919706300135852;	const double ndc4    = 0.71199788901815727304919706300135852;
const double dc5    = 0.49457655717669416346035981441598679;	const double B_c5645 =-1.957399490271361553945056740;
const double dc6    =-0.96808390091782605854354886831133920;	const double D_c5645 =-1.354874954742460424002602037;
const double dc7    = 0.32508216495576250692240420538780022;
const double dc8    = 0.16949439222093165653795560902818660;
const double dc9    =-0.29637372110299413755460095857226913;
const double dc10   = 0.04534684817389996223192362526889384;	const double B_c9019 = 0.1530056308809554651383615683;
const double dc11   =-0.34172056927689409978652458384116289;	const double D_c9019 =-6.668408724817368094978521568;
const double dc12   = 0.07310633583025898405268087482211124;
const double dc13   =-0.20167486904767495859952810665520680;	const double B_c5285 =-0.5605781319395146459386195222;
const double dc14   =-0.84201778467603236410459766557180865;	const double D_c5285 =-5.397694515310063808645584290;
const double dc15   =-0.13041239332209095847634591191071144;
const double dc16   =-0.45935884940861915763712570904324634;	const double B_c6396 =-0.4390355585993655658779099885;
const double dc17   = 0.48562853299381401878246608937483484;	const double D_c6396 =-8.169650232597136939166445099;
const double dc18   = 0.05290024156536825355103499921807512;
const double dc19   =-0.14593330617823442758015434230035912;	const double B_c7407 =-1.733872141912958941825760722;
const double dc20   =-0.60928979281795509222872426261479928;	const double D_c7407 =-1.459688060776070564000064511;
const double dc21   = 0.50714554918236722335933886804513386;
const double dc22   =-0.09614247372823604317748361539732969;	const double C_c4174 = 8.015643444687584573541352865;
const double dc23   =-0.57652834602834558094251739609128330;	const double E_c4174 =-1.021509017898617060353255097;
const double dc24   = 0.06326947458204131524249727162375559;
const double dc25   =-0.76932220806904896036979558584064568;	const double C_c5285 =-0.1249703605587412988133367679;
const double dc26   =-0.41013896184380627585461186107586186;	const double E_c5285 =-89.48778340022079724273433851;
const double dc27   = 0.36697396683700721180660373150156373;
const double dc28   =-0.06956934754225036504904405641510024;	const double C_c6396 = 1.405690265164092306971350691;
const double dc29   =-0.41717983028166295195804014713757581;	const double E_c6396 =-1.699387856677906977066853648;

const double ds0    = 0.92796072713833698701991188315309154;	const double ds0_5th = 0.18559214542766739740398237663061831;
const double ds1    = 0.75538558300353835543581057551249341;
const double ds2    =-0.71798870119142816504365819783899618;	const double B_s3483 =-0.1543847659293762455593100378;
const double ds3    =-1.14129841620797351615515640500369610;	const double D_s3483 = 77.72870185137802409062726308;
const double ds4    = 0.17619908884183581953599075099065958;	const double nds4    =-0.17619908884183581953599075099065958;
const double ds5    =-0.44789045813036076345265894471834662;	const double B_s5645 = 0.9451188506751776707508747134;
const double ds6    =-0.42330971501654535111149820716469997;	const double D_s5645 = 1.689563031430332839157319596;
const double ds7    =-0.54178961234959234550766744684833662;
const double ds8    = 0.09389915421923158205500850212998998;
const double ds9    = 0.10209749786491606368824206751661151;
const double ds10   = 0.36010442196019251577804178188101286;	const double B_s9019 =-3.527064124888175236105416334;
const double ds11   =-0.25800692409527645208979971436440136;	const double D_s9019 =-1.112194193764683738807670520;
const double ds12   =-0.68362751420088060523602271315212364;
const double ds13   = 1.42448046407360845685461718962915170;	const double B_s5285 = 1.568615537822747404332565815;
const double ds14   = 0.71705834223306922870900492103008074;	const double D_s5285 =-1.561648155256061432564325550;
const double ds15   =-0.43581585016667742322916851390840979;
const double ds16   = 0.44515550282226738565114480689187057;	const double B_s6396 =-3.199961485463981734674359059;
const double ds17   = 0.52796329829794124367930049003364155;	const double D_s6396 =-1.134960866988575539262620301;
const double ds18   =-0.49467751640467748805281786601391635;
const double ds19   = 1.03076374706570777841622291802646680;	const double B_s7407 = 1.358159448857025499908422489;
const double ds20   = 0.51886829082317972874576505841504815;	const double D_s7407 =-1.749196678153277082987222794;
const double ds21   = 0.82187006169806796277221650036735597;
const double ds22   =-0.79391727157204814808670251421336608;	const double C_s4174 =-4.143907256301630677016949990;
const double ds23   =-0.75874209018537673384980278323642215;	const double E_s4174 =-1.080478024630318207542706141;
const double ds24   =-0.19833215631171570317346941647854694;
const double ds25   =-0.43149551708674121082911820080227901;	const double C_s5285 =-1.839920092176650010161002721;
const double ds26   = 0.27843331019209505198390355841028348;	const double E_s5285 =-1.408224849824341570791003234;
const double ds27   = 0.59471076351191660168095301828858968;
const double ds28   =-0.57448393456065017248053353113326010;	const double C_s6396 =-2.725040655738747305153314734;
const double ds29   =-0.54903093419716620563980211536607990;	const double E_s6396 =-1.186102149729676553233822995;

const double c50    =-0.89442719099991587856366946749251043;
const double c51    =-0.27639320225002103035908263312687237;	const double c51d50m1 =-0.6909830056250525758977065828;	/* c51/c50 - 1 */
const double s52    = 0.61803398874989484820458683436563806;

/* This version implements a straight tangent-based radix-31 DFT (cf. the sample Forran-90 subroutine in fft31t.i90
by James Van Buskirk, ???) and only tries to minimize the number
of temporaries used in each of the major phases.
For a 16-register-optimized version which moreover mimics x86-style destructive (2-input, one overwritten with output)
register arithmetic, see the RADIX_31_DFT version of the same macro .

Output ordering is as for DIF, so for DIT caller must reverse order of last 12 B-output args in the macro call.

ASSUMES that the needed DC and DS-constants have been defined in the calling routine via inclusion of the radix31.h header file..
*/

/* C translation (macro-ized) of Van Buskirk`s Fortran-90 example code - This is much too temporary-heavy for practical use
(i.e. assembler implementation) but serves as a starting point for the optimizations we make on our way toward the latter: */
#define RADIX_31_JVB(\
__Ar00,__Ar01,__Ar02,__Ar03,__Ar04,__Ar05,__Ar06,__Ar07,__Ar08,__Ar09,__Ar10,__Ar11,__Ar12,__Ar13,__Ar14,__Ar15,__Ar16,__Ar17,__Ar18,__Ar19,__Ar20,__Ar21,__Ar22,__Ar23,__Ar24,__Ar25,__Ar26,__Ar27,__Ar28,__Ar29,__Ar30,\
__Ai00,__Ai01,__Ai02,__Ai03,__Ai04,__Ai05,__Ai06,__Ai07,__Ai08,__Ai09,__Ai10,__Ai11,__Ai12,__Ai13,__Ai14,__Ai15,__Ai16,__Ai17,__Ai18,__Ai19,__Ai20,__Ai21,__Ai22,__Ai23,__Ai24,__Ai25,__Ai26,__Ai27,__Ai28,__Ai29,__Ai30,\
__Br00,__Br01,__Br02,__Br03,__Br04,__Br05,__Br06,__Br07,__Br08,__Br09,__Br10,__Br11,__Br12,__Br13,__Br14,__Br15,__Br16,__Br17,__Br18,__Br19,__Br20,__Br21,__Br22,__Br23,__Br24,__Br25,__Br26,__Br27,__Br28,__Br29,__Br30,\
__Bi00,__Bi01,__Bi02,__Bi03,__Bi04,__Bi05,__Bi06,__Bi07,__Bi08,__Bi09,__Bi10,__Bi11,__Bi12,__Bi13,__Bi14,__Bi15,__Bi16,__Bi17,__Bi18,__Bi19,__Bi20,__Bi21,__Bi22,__Bi23,__Bi24,__Bi25,__Bi26,__Bi27,__Bi28,__Bi29,__Bi30\
)\
{\
	double xr0, xi0, xr1, xi1, xr2, xi2, xr3, xi3, xr4, xi4, xr5, xi5,\
	xr6, xi6, xr7, xi7, xr8, xi8, xr9, xi9, xra, xia, xrb, xib,\
	xrc, xic, xrd, xid, xre, xie, xrf, xif,\
	yr1, yi1, yr2, yi2, yr3, yi3, yr4, yi4, yr5, yi5, yr6, yi6,\
	yr7, yi7, yr8, yi8, yr9, yi9, yra, yia, yrb, yib, yrc, yic,\
	yrd, yid, yre, yie, yrf, yif,\
	xr16b, xi16b, xr49e, xi49e, xr27c, xi27c,\
	xr5af, xi5af, xr38d, xi38d, xr7_2, xi7_2,\
	xr8_3, xi8_3, xra_5, xia_5, xrb_6, xib_6,\
	xrd_8, xid_8, xre_9, xie_9, xr1_b, xi1_b,\
	xr2_c, xi2_c, xr4_e, xi4_e, xr5_f, xi5_f,\
	xr16b49e, xi16b49e, xr27c5af, xi27c5af,\
	xrsigma, xisigma, ursigma, uisigma,\
	ur1234, ui1234, ur12, ui12, ur34, ui34,\
	ur1, ui1, ur2, ui2,\
	ur3, ui3, ur4, ui4, ur5, ui5,\
	xra4_5e, xia4_5e, xrae_54, xiae_54,\
	xrd1_8b, xid1_8b, xrdb_81, xidb_81,\
	xra4d1_5e8b, xia4d1_5e8b, pr6, pi6, pr8a, pi8a,\
	pr8, pi8, pra, pia, prc, pic, pre, pie,\
	xrb5_6f, xib5_6f, xrbf_65, xibf_65,\
	xre2_9c, xie2_9c, xrec_92, xiec_92,\
	xrb5e2_6f9c, xib5e2_6f9c, pr7, pi7, pr9b, pi9b,\
	pr9, pi9, prb, pib, prd, pid, prf, pif,\
	ur6, ui6, ur7, ui7,\
	pr89, pi89, prcd, picd,\
	ur89, ui89, ur8, ui8, ur9, ui9,\
	urcd, uicd, urc, uic, urd, uid,\
	prab, piab, pref, pief,\
	urab, uiab, ura, uia, urb, uib,\
	uref, uief, ure, uie, urf, uif,\
	ur8a, ui8a, wr6, wi6, wr8eac, wi8eac,\
	wr8e, wi8e, wrac, wiac,\
	wr8_e, wi8_e, wra_c, wia_c,\
	wr8, wi8, wre, wie, wra, wia, wrc, wic,\
	ur9b, ui9b, wr7, wi7, wr9fbd, wi9fbd,\
	wr9f, wi9f, wrbd, wibd,\
	wr9_f, wi9_f, wrb_d, wib_d,\
	wr9, wi9, wrf, wif, wrb, wib, wrd, wid,\
	sr1, si1, sr2, si2, sr3, si3, sr4, si4, sr5, si5,\
	sr6, si6, sr7, si7, sr8, si8, sr9, si9, sra, sia,\
	srb, sib, src, sic, srd, sid, sre, sie, srf, sif,\
	yr16b, yi16b, yr49e, yi49e, yr27c, yi27c,\
	yr5af, yi5af, yr38d, yi38d, yr7_2, yi7_2,\
	yr8_3, yi8_3, yra_5, yia_5, yrb_6, yib_6,\
	yrd_8, yid_8, yre_9, yie_9, yr1_b, yi1_b,\
	yr2_c, yi2_c, yr4_e, yi4_e, yr5_f, yi5_f,\
	yr16b49e, yi16b49e, yr27c5af, yi27c5af,\
	yrsigma, yisigma, vrsigma, visigma,\
	vr1234, vi1234, vr12, vi12, vr34, vi34,\
	vr1, vi1, vr2, vi2,\
	vr3, vi3, vr4, vi4, vr5, vi5,\
	yra4_5e, yia4_5e, yrae_54, yiae_54,\
	yrd1_8b, yid1_8b, yrdb_81, yidb_81,\
	yra4d1_5e8b, yia4d1_5e8b, qr6, qi6, qr8a, qi8a,\
	qr8, qi8, qra, qia, qrc, qic, qre, qie,\
	yrb5_6f, yib5_6f, yrbf_65, yibf_65,\
	yre2_9c, yie2_9c, yrec_92, yiec_92,\
	yrb5e2_6f9c, yib5e2_6f9c, qr7, qi7, qr9b, qi9b,\
	qr9, qi9, qrb, qib, qrd, qid, qrf, qif,\
	vr6, vi6, vr7, vi7,\
	qr89, qi89, qrcd, qicd,\
	vr89, vi89, vr8, vi8, vr9, vi9,\
	vrcd, vicd, vrc, vic, vrd, vid,\
	qrab, qiab, qref, qief,\
	vrab, viab, vra, via, vrb, vib,\
	vref, vief, vre, vie, vrf, vif,\
	vr8a, vi8a, zr6, zi6, zr8eac, zi8eac,\
	zr8e, zi8e, zrac, ziac,\
	zr8_e, zi8_e, zra_c, zia_c,\
	zr8, zi8, zre, zie, zra, zia, zrc, zic,\
	vr9b, vi9b, zr7, zi7, zr9fbd, zi9fbd,\
	zr9f, zi9f, zrbd, zibd,\
	zr9_f, zi9_f, zrb_d, zib_d,\
	zr9, zi9, zrf, zif, zrb, zib, zrd, zid,\
	tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4, tr5, ti5,\
	tr6, ti6, tr7, ti7, tr8, ti8, tr9, ti9, tra, tia,\
	trb, tib, trc, tic, trd, tid, tre, tie, trf, tif;\
\
	xr0 = __Ar00;\
	xi0 = __Ai00;\
	xr1 = __Ar01 + __Ar30;\
	xi1 = __Ai01 + __Ai30;\
	yr1 = __Ar01 - __Ar30;\
	yi1 = __Ai01 - __Ai30;\
	xr2 = __Ar07 + __Ar24;\
	xi2 = __Ai07 + __Ai24;\
	yr2 = __Ar07 - __Ar24;\
	yi2 = __Ai07 - __Ai24;\
	xr3 = __Ar18 + __Ar13;\
	xi3 = __Ai18 + __Ai13;\
	yr3 = __Ar18 - __Ar13;\
	yi3 = __Ai18 - __Ai13;\
	xr4 = __Ar02 + __Ar29;\
	xi4 = __Ai02 + __Ai29;\
	yr4 = __Ar02 - __Ar29;\
	yi4 = __Ai02 - __Ai29;\
	xr5 = __Ar14 + __Ar17;\
	xi5 = __Ai14 + __Ai17;\
	yr5 = __Ar14 - __Ar17;\
	yi5 = __Ai14 - __Ai17;\
	xr6 = __Ar05 + __Ar26;\
	xi6 = __Ai05 + __Ai26;\
	yr6 = __Ar05 - __Ar26;\
	yi6 = __Ai05 - __Ai26;\
	xr7 = __Ar04 + __Ar27;\
	xi7 = __Ai04 + __Ai27;\
	yr7 = __Ar04 - __Ar27;\
	yi7 = __Ai04 - __Ai27;\
	xr8 = __Ar28 + __Ar03;\
	xi8 = __Ai28 + __Ai03;\
	yr8 = __Ar28 - __Ar03;\
	yi8 = __Ai28 - __Ai03;\
	xr9 = __Ar10 + __Ar21;\
	xi9 = __Ai10 + __Ai21;\
	yr9 = __Ar10 - __Ar21;\
	yi9 = __Ai10 - __Ai21;\
	xra = __Ar08 + __Ar23;\
	xia = __Ai08 + __Ai23;\
	yra = __Ar08 - __Ar23;\
	yia = __Ai08 - __Ai23;\
	xrb = __Ar25 + __Ar06;\
	xib = __Ai25 + __Ai06;\
	yrb = __Ar25 - __Ar06;\
	yib = __Ai25 - __Ai06;\
	xrc = __Ar20 + __Ar11;\
	xic = __Ai20 + __Ai11;\
	yrc = __Ar20 - __Ar11;\
	yic = __Ai20 - __Ai11;\
	xrd = __Ar16 + __Ar15;\
	xid = __Ai16 + __Ai15;\
	yrd = __Ar16 - __Ar15;\
	yid = __Ai16 - __Ai15;\
	xre = __Ar19 + __Ar12;\
	xie = __Ai19 + __Ai12;\
	yre = __Ar19 - __Ar12;\
	yie = __Ai19 - __Ai12;\
	xrf = __Ar09 + __Ar22;\
	xif = __Ai09 + __Ai22;\
	yrf = __Ar09 - __Ar22;\
	yif = __Ai09 - __Ai22;\
\
	xr16b = xr1 + xr6 + xrb;\
	xi16b = xi1 + xi6 + xib;\
	xr49e = xr4 + xr9 + xre;\
	xi49e = xi4 + xi9 + xie;\
	xr27c = xr2 + xr7 + xrc;\
	xi27c = xi2 + xi7 + xic;\
	xr5af = xr5 + xra + xrf;\
	xi5af = xi5 + xia + xif;\
	xr38d = xr3 + xr8 + xrd;\
	xi38d = xi3 + xi8 + xid;\
	xr7_2 = xr7 - xr2;\
	xi7_2 = xi7 - xi2;\
	xr8_3 = xr8 - xr3;\
	xi8_3 = xi8 - xi3;\
	xra_5 = xra - xr5;\
	xia_5 = xia - xi5;\
	xrb_6 = xrb - xr6;\
	xib_6 = xib - xi6;\
	xrd_8 = xrd - xr8;\
	xid_8 = xid - xi8;\
	xre_9 = xre - xr9;\
	xie_9 = xie - xi9;\
	xr1_b = xr1 - xrb;\
	xi1_b = xi1 - xib;\
	xr2_c = xr2 - xrc;\
	xi2_c = xi2 - xic;\
	xr4_e = xr4 - xre;\
	xi4_e = xi4 - xie;\
	xr5_f = xr5 - xrf;\
	xi5_f = xi5 - xif;\
	xr16b = xr16b - xr38d;\
	xi16b = xi16b - xi38d;\
	xr49e = xr49e - xr38d;\
	xi49e = xi49e - xi38d;\
	xr27c = xr27c - xr38d;\
	xi27c = xi27c - xi38d;\
	xr5af = xr5af - xr38d;\
	xi5af = xi5af - xi38d;\
	xr16b49e = xr16b + xr49e;\
	xi16b49e = xi16b + xi49e;\
	xr27c5af = xr27c + xr5af;\
	xi27c5af = xi27c + xi5af;\
	xrsigma = xr16b49e + xr27c5af;\
	xisigma = xi16b49e + xi27c5af;\
	ursigma = xrsigma + 5*xr38d;\
	uisigma = xisigma + 5*xi38d;\
/* DC componentt: */\
	__Br00 = xr0 + ursigma;\
	__Bi00 = xi0 + uisigma;\
	ur1234 = xr0 + dc0*xr38d + dc1*xrsigma;\
	ui1234 = xi0 + dc0*xi38d + dc1*xisigma;\
	ur12 = ur1234 + dc2*xr27c5af;\
	ui12 = ui1234 + dc2*xi27c5af;\
	ur34 = ur1234 + dc7*xr16b49e;\
	ui34 = ui1234 + dc7*xi16b49e;\
	ur1 = ur12 + dc3*xr49e + dc4*xr5af;\
	ui1 = ui12 + dc3*xi49e + dc4*xi5af;\
	ur2 = ur12 + dc5*xr16b + dc6*xr27c;\
	ui2 = ui12 + dc5*xi16b + dc6*xi27c;\
	ur3 = ur34 + dc8*xr49e + dc3*xr5af;\
	ui3 = ui34 + dc8*xi49e + dc3*xi5af;\
	ur4 = ur34 - dc4*xr16b + dc5*xr27c;\
	ui4 = ui34 - dc4*xi16b + dc5*xi27c;\
	ur5 = dc0*ursigma - (ur1 + ur2 + ur3 + ur4) + 5*xr0;\
	ui5 = dc0*uisigma - (ui1 + ui2 + ui3 + ui4) + 5*xi0;\
	xra4_5e = xra_5 + xr4_e;\
	xia4_5e = xia_5 + xi4_e;\
	xrae_54 = xra_5 - xr4_e;\
	xiae_54 = xia_5 - xi4_e;\
	xrd1_8b = xrd_8 + xr1_b;\
	xid1_8b = xid_8 + xi1_b;\
	xrdb_81 = xrd_8 - xr1_b;\
	xidb_81 = xid_8 - xi1_b;\
	xra4d1_5e8b = xra4_5e + xrd1_8b;\
	xia4d1_5e8b = xia4_5e + xid1_8b;\
	pr6 = xr7_2 + xra4d1_5e8b;\
	pi6 = xi7_2 + xia4d1_5e8b;\
	pr8a = c50*xr7_2 + c51*xra4d1_5e8b;\
	pi8a = c50*xi7_2 + c51*xia4d1_5e8b;\
	pr8 = pr8a + xrd1_8b;\
	pi8 = pi8a + xid1_8b;\
	pra = pr8a + xra4_5e;\
	pia = pi8a + xia4_5e;\
	prc = xrae_54 + s52*xrdb_81;\
	pic = xiae_54 + s52*xidb_81;\
	pre = s52*xrae_54 - xrdb_81;\
	pie = s52*xiae_54 - xidb_81;\
	xrb5_6f = xrb_6 + xr5_f;\
	xib5_6f = xib_6 + xi5_f;\
	xrbf_65 = xrb_6 - xr5_f;\
	xibf_65 = xib_6 - xi5_f;\
	xre2_9c = xre_9 + xr2_c;\
	xie2_9c = xie_9 + xi2_c;\
	xrec_92 = xre_9 - xr2_c;\
	xiec_92 = xie_9 - xi2_c;\
	xrb5e2_6f9c = xrb5_6f + xre2_9c;\
	xib5e2_6f9c = xib5_6f + xie2_9c;\
	pr7 = xr8_3 + xrb5e2_6f9c;\
	pi7 = xi8_3 + xib5e2_6f9c;\
	pr9b = c50*xr8_3 + c51*xrb5e2_6f9c;\
	pi9b = c50*xi8_3 + c51*xib5e2_6f9c;\
	pr9 = pr9b + xre2_9c;\
	pi9 = pi9b + xie2_9c;\
	prb = pr9b + xrb5_6f;\
	pib = pi9b + xib5_6f;\
	prd = xrbf_65 + s52*xrec_92;\
	pid = xibf_65 + s52*xiec_92;\
	prf = s52*xrbf_65 - xrec_92;\
	pif = s52*xibf_65 - xiec_92;\
	ur6 = dc9*pr6 + dc10*pr7;\
	ui6 = dc9*pi6 + dc10*pi7;\
	ur7 = dc11*pr6 + dc9*pr7;\
	ui7 = dc11*pi6 + dc9*pi7;\
	pr89 = pr8 + pr9;\
	pi89 = pi8 + pi9;\
	prcd = prc + prd;\
	picd = pic + pid;\
	ur89 = dc12*pr89 + dc15*prcd;\
	ui89 = dc12*pi89 + dc15*picd;\
	ur8 = ur89 + dc13*pr9 + dc16*prd;\
	ui8 = ui89 + dc13*pi9 + dc16*pid;\
	ur9 = ur89 + dc14*pr8 + dc17*prc;\
	ui9 = ui89 + dc14*pi8 + dc17*pic;\
	urcd = dc18*prcd - dc15*pr89;\
	uicd = dc18*picd - dc15*pi89;\
	urc = urcd - dc16*pr9 + dc19*prd;\
	uic = uicd - dc16*pi9 + dc19*pid;\
	urd = urcd - dc17*pr8 + dc20*prc;\
	uid = uicd - dc17*pi8 + dc20*pic;\
	prab = pra + prb;\
	piab = pia + pib;\
	pref = pre + prf;\
	pief = pie + pif;\
	urab = dc21*prab + dc24*pref;\
	uiab = dc21*piab + dc24*pief;\
	ura = urab + dc22*prb + dc25*prf;\
	uia = uiab + dc22*pib + dc25*pif;\
	urb = urab + dc23*pra + dc26*pre;\
	uib = uiab + dc23*pia + dc26*pie;\
	uref = dc27*pref - dc24*prab;\
	uief = dc27*pief - dc24*piab;\
	ure = uref - dc25*prb + dc28*prf;\
	uie = uief - dc25*pib + dc28*pif;\
	urf = uref - dc26*pra + dc29*pre;\
	uif = uief - dc26*pia + dc29*pie;\
	ur8a = ur8 + ura;\
	ui8a = ui8 + uia;\
	wr6 = ur6 + c50*ur8a;\
	wi6 = ui6 + c50*ui8a;\
	wr8eac = ur6 + c51*ur8a;\
	wi8eac = ui6 + c51*ui8a;\
	wr8e = wr8eac + ura;\
	wi8e = wi8eac + uia;\
	wrac = wr8eac + ur8;\
	wiac = wi8eac + ui8;\
	wr8_e = urc + s52*ure;\
	wi8_e = uic + s52*uie;\
	wra_c = s52*urc - ure;\
	wia_c = s52*uic - uie;\
	wr8 = wr8e + wr8_e;\
	wi8 = wi8e + wi8_e;\
	wre = wr8e - wr8_e;\
	wie = wi8e - wi8_e;\
	wra = wrac + wra_c;\
	wia = wiac + wia_c;\
	wrc = wrac - wra_c;\
	wic = wiac - wia_c;\
	ur9b = ur9 + urb;\
	ui9b = ui9 + uib;\
	wr7 = ur7 + c50*ur9b;\
	wi7 = ui7 + c50*ui9b;\
	wr9fbd = ur7 + c51*ur9b;\
	wi9fbd = ui7 + c51*ui9b;\
	wr9f = wr9fbd + urb;\
	wi9f = wi9fbd + uib;\
	wrbd = wr9fbd + ur9;\
	wibd = wi9fbd + ui9;\
	wr9_f = urd + s52*urf;\
	wi9_f = uid + s52*uif;\
	wrb_d = s52*urd - urf;\
	wib_d = s52*uid - uif;\
	wr9 = wr9f + wr9_f;\
	wi9 = wi9f + wi9_f;\
	wrf = wr9f - wr9_f;\
	wif = wi9f - wi9_f;\
	wrb = wrbd + wrb_d;\
	wib = wibd + wib_d;\
	wrd = wrbd - wrb_d;\
	wid = wibd - wib_d;\
	sr1 = ur1 + wrc;\
	si1 = ui1 + wic;\
	sr2 = ur3 - wr6 + wrd;\
	si2 = ui3 - wi6 + wid;\
	sr3 = ur5 - wr7;\
	si3 = ui5 - wi7;\
	sr4 = ur2 + wre;\
	si4 = ui2 + wie;\
	sr5 = ur4 - wr8 + wrf;\
	si5 = ui4 - wi8 + wif;\
	sr6 = ur1 - wr9;\
	si6 = ui1 - wi9;\
	sr7 = ur3 + wr6;\
	si7 = ui3 + wi6;\
	sr8 = ur5 + wr7 - wra;\
	si8 = ui5 + wi7 - wia;\
	sr9 = ur2 - wrb;\
	si9 = ui2 - wib;\
	sra = ur4 + wr8;\
	sia = ui4 + wi8;\
	srb = ur1 + wr9 - wrc;\
	sib = ui1 + wi9 - wic;\
	src = ur3 - wrd;\
	sic = ui3 - wid;\
	srd = ur5 + wra;\
	sid = ui5 + wia;\
	sre = ur2 + wrb - wre;\
	sie = ui2 + wib - wie;\
	srf = ur4 - wrf;\
	sif = ui4 - wif;\
\
	yr16b = yr1 + yr6 + yrb;\
	yi16b = yi1 + yi6 + yib;\
	yr49e = yr4 + yr9 + yre;\
	yi49e = yi4 + yi9 + yie;\
	yr27c = yr2 + yr7 + yrc;\
	yi27c = yi2 + yi7 + yic;\
	yr5af = yr5 + yra + yrf;\
	yi5af = yi5 + yia + yif;\
	yr38d = yr3 + yr8 + yrd;\
	yi38d = yi3 + yi8 + yid;\
	yr7_2 = yr7 - yr2;\
	yi7_2 = yi7 - yi2;\
	yr8_3 = yr8 - yr3;\
	yi8_3 = yi8 - yi3;\
	yra_5 = yra - yr5;\
	yia_5 = yia - yi5;\
	yrb_6 = yrb - yr6;\
	yib_6 = yib - yi6;\
	yrd_8 = yrd - yr8;\
	yid_8 = yid - yi8;\
	yre_9 = yre - yr9;\
	yie_9 = yie - yi9;\
	yr1_b = yr1 - yrb;\
	yi1_b = yi1 - yib;\
	yr2_c = yr2 - yrc;\
	yi2_c = yi2 - yic;\
	yr4_e = yr4 - yre;\
	yi4_e = yi4 - yie;\
	yr5_f = yr5 - yrf;\
	yi5_f = yi5 - yif;\
	yr16b = yr16b - yr38d;\
	yi16b = yi16b - yi38d;\
	yr49e = yr49e - yr38d;\
	yi49e = yi49e - yi38d;\
	yr27c = yr27c - yr38d;\
	yi27c = yi27c - yi38d;\
	yr5af = yr5af - yr38d;\
	yi5af = yi5af - yi38d;\
	yr16b49e = yr16b + yr49e;\
	yi16b49e = yi16b + yi49e;\
	yr27c5af = yr27c + yr5af;\
	yi27c5af = yi27c + yi5af;\
	yrsigma = yr16b49e + yr27c5af;\
	yisigma = yi16b49e + yi27c5af;\
	vrsigma = yrsigma + 5*yr38d;\
	visigma = yisigma + 5*yi38d;\
	vr1234 = ds0*yr38d + ds1*yrsigma;\
	vi1234 = ds0*yi38d + ds1*yisigma;\
	vr12 = vr1234 + ds2*yr27c5af;\
	vi12 = vi1234 + ds2*yi27c5af;\
	vr34 = vr1234 + ds7*yr16b49e;\
	vi34 = vi1234 + ds7*yi16b49e;\
	vr1 = vr12 + ds3*yr49e + ds4*yr5af;\
	vi1 = vi12 + ds3*yi49e + ds4*yi5af;\
	vr2 = vr12 + ds5*yr16b + ds6*yr27c;\
	vi2 = vi12 + ds5*yi16b + ds6*yi27c;\
	vr3 = vr34 + ds8*yr49e + ds3*yr5af;\
	vi3 = vi34 + ds8*yi49e + ds3*yi5af;\
	vr4 = vr34 - ds4*yr16b + ds5*yr27c;\
	vi4 = vi34 - ds4*yi16b + ds5*yi27c;\
	vr5 = ds0*vrsigma - (vr1 + vr2 + vr3 + vr4);\
	vi5 = ds0*visigma - (vi1 + vi2 + vi3 + vi4);\
	yra4_5e = yra_5 + yr4_e;\
	yia4_5e = yia_5 + yi4_e;\
	yrae_54 = yra_5 - yr4_e;\
	yiae_54 = yia_5 - yi4_e;\
	yrd1_8b = yrd_8 + yr1_b;\
	yid1_8b = yid_8 + yi1_b;\
	yrdb_81 = yrd_8 - yr1_b;\
	yidb_81 = yid_8 - yi1_b;\
	yra4d1_5e8b = yra4_5e + yrd1_8b;\
	yia4d1_5e8b = yia4_5e + yid1_8b;\
	qr6 = yr7_2 + yra4d1_5e8b;\
	qi6 = yi7_2 + yia4d1_5e8b;\
	qr8a = c50*yr7_2 + c51*yra4d1_5e8b;\
	qi8a = c50*yi7_2 + c51*yia4d1_5e8b;\
	qr8 = qr8a + yrd1_8b;\
	qi8 = qi8a + yid1_8b;\
	qra = qr8a + yra4_5e;\
	qia = qi8a + yia4_5e;\
	qrc = yrae_54 + s52*yrdb_81;\
	qic = yiae_54 + s52*yidb_81;\
	qre = s52*yrae_54 - yrdb_81;\
	qie = s52*yiae_54 - yidb_81;\
	yrb5_6f = yrb_6 + yr5_f;\
	yib5_6f = yib_6 + yi5_f;\
	yrbf_65 = yrb_6 - yr5_f;\
	yibf_65 = yib_6 - yi5_f;\
	yre2_9c = yre_9 + yr2_c;\
	yie2_9c = yie_9 + yi2_c;\
	yrec_92 = yre_9 - yr2_c;\
	yiec_92 = yie_9 - yi2_c;\
	yrb5e2_6f9c = yrb5_6f + yre2_9c;\
	yib5e2_6f9c = yib5_6f + yie2_9c;\
	qr7 = yr8_3 + yrb5e2_6f9c;\
	qi7 = yi8_3 + yib5e2_6f9c;\
	qr9b = c50*yr8_3 + c51*yrb5e2_6f9c;\
	qi9b = c50*yi8_3 + c51*yib5e2_6f9c;\
	qr9 = qr9b + yre2_9c;\
	qi9 = qi9b + yie2_9c;\
	qrb = qr9b + yrb5_6f;\
	qib = qi9b + yib5_6f;\
	qrd = yrbf_65 + s52*yrec_92;\
	qid = yibf_65 + s52*yiec_92;\
	qrf = s52*yrbf_65 - yrec_92;\
	qif = s52*yibf_65 - yiec_92;\
	vr6 = ds9*qr6 + ds10*qr7;\
	vi6 = ds9*qi6 + ds10*qi7;\
	vr7 = ds11*qr6 + ds9*qr7;\
	vi7 = ds11*qi6 + ds9*qi7;\
	qr89 = qr8 + qr9;\
	qi89 = qi8 + qi9;\
	qrcd = qrc + qrd;\
	qicd = qic + qid;\
	vr89 = ds12*qr89 + ds15*qrcd;\
	vi89 = ds12*qi89 + ds15*qicd;\
	vr8 = vr89 + ds13*qr9 + ds16*qrd;\
	vi8 = vi89 + ds13*qi9 + ds16*qid;\
	vr9 = vr89 + ds14*qr8 + ds17*qrc;\
	vi9 = vi89 + ds14*qi8 + ds17*qic;\
	vrcd = ds18*qrcd - ds15*qr89;\
	vicd = ds18*qicd - ds15*qi89;\
	vrc = vrcd - ds16*qr9 + ds19*qrd;\
	vic = vicd - ds16*qi9 + ds19*qid;\
	vrd = vrcd - ds17*qr8 + ds20*qrc;\
	vid = vicd - ds17*qi8 + ds20*qic;\
	qrab = qra + qrb;\
	qiab = qia + qib;\
	qref = qre + qrf;\
	qief = qie + qif;\
	vrab = ds21*qrab + ds24*qref;\
	viab = ds21*qiab + ds24*qief;\
	vra = vrab + ds22*qrb + ds25*qrf;\
	via = viab + ds22*qib + ds25*qif;\
	vrb = vrab + ds23*qra + ds26*qre;\
	vib = viab + ds23*qia + ds26*qie;\
	vref = ds27*qref - ds24*qrab;\
	vief = ds27*qief - ds24*qiab;\
	vre = vref - ds25*qrb + ds28*qrf;\
	vie = vief - ds25*qib + ds28*qif;\
	vrf = vref - ds26*qra + ds29*qre;\
	vif = vief - ds26*qia + ds29*qie;\
	vr8a = vr8 + vra;\
	vi8a = vi8 + via;\
	zr6 = vr6 + c50*vr8a;\
	zi6 = vi6 + c50*vi8a;\
	zr8eac = vr6 + c51*vr8a;\
	zi8eac = vi6 + c51*vi8a;\
	zr8e = zr8eac + vra;\
	zi8e = zi8eac + via;\
	zrac = zr8eac + vr8;\
	ziac = zi8eac + vi8;\
	zr8_e = vrc + s52*vre;\
	zi8_e = vic + s52*vie;\
	zra_c = s52*vrc - vre;\
	zia_c = s52*vic - vie;\
	zr8 = zr8e + zr8_e;\
	zi8 = zi8e + zi8_e;\
	zre = zr8e - zr8_e;\
	zie = zi8e - zi8_e;\
	zra = zrac + zra_c;\
	zia = ziac + zia_c;\
	zrc = zrac - zra_c;\
	zic = ziac - zia_c;\
	vr9b = vr9 + vrb;\
	vi9b = vi9 + vib;\
	zr7 = vr7 + c50*vr9b;\
	zi7 = vi7 + c50*vi9b;\
	zr9fbd = vr7 + c51*vr9b;\
	zi9fbd = vi7 + c51*vi9b;\
	zr9f = zr9fbd + vrb;\
	zi9f = zi9fbd + vib;\
	zrbd = zr9fbd + vr9;\
	zibd = zi9fbd + vi9;\
	zr9_f = vrd + s52*vrf;\
	zi9_f = vid + s52*vif;\
	zrb_d = s52*vrd - vrf;\
	zib_d = s52*vid - vif;\
	zr9 = zr9f + zr9_f;\
	zi9 = zi9f + zi9_f;\
	zrf = zr9f - zr9_f;\
	zif = zi9f - zi9_f;\
	zrb = zrbd + zrb_d;\
	zib = zibd + zib_d;\
	zrd = zrbd - zrb_d;\
	zid = zibd - zib_d;\
	tr1 = vr1 + zrc;\
	ti1 = vi1 + zic;\
	tr2 = vr3 - zr6 + zrd;\
	ti2 = vi3 - zi6 + zid;\
	tr3 = vr5 - zr7;\
	ti3 = vi5 - zi7;\
	tr4 = vr2 + zre;\
	ti4 = vi2 + zie;\
	tr5 = vr4 - zr8 + zrf;\
	ti5 = vi4 - zi8 + zif;\
	tr6 = vr1 - zr9;\
	ti6 = vi1 - zi9;\
	tr7 = vr3 + zr6;\
	ti7 = vi3 + zi6;\
	tr8 = vr5 + zr7 - zra;\
	ti8 = vi5 + zi7 - zia;\
	tr9 = vr2 - zrb;\
	ti9 = vi2 - zib;\
	tra = vr4 + zr8;\
	tia = vi4 + zi8;\
	trb = vr1 + zr9 - zrc;\
	tib = vi1 + zi9 - zic;\
	trc = vr3 - zrd;\
	tic = vi3 - zid;\
	trd = vr5 + zra;\
	tid = vi5 + zia;\
	tre = vr2 + zrb - zre;\
	tie = vi2 + zib - zie;\
	trf = vr4 - zrf;\
	tif = vi4 - zif;\
\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
	__Br01 = srf - tif;		__Bi01 = sif + trf;\
	__Br02 = src - tic;		__Bi02 = sic + trc;\
	__Br03 = sr8 + ti8;		__Bi03 = si8 - tr8;\
	__Br04 = sr9 - ti9;		__Bi04 = si9 + tr9;\
	__Br05 = sra - tia;		__Bi05 = sia + tra;\
	__Br06 = sr5 + ti5;		__Bi06 = si5 - tr5;\
	__Br07 = sre - tie;		__Bi07 = sie + tre;\
	__Br08 = sr6 - ti6;		__Bi08 = si6 + tr6;\
	__Br09 = sr1 - ti1;		__Bi09 = si1 + tr1;\
	__Br10 = sr7 - ti7;		__Bi10 = si7 + tr7;\
	__Br11 = sr4 + ti4;		__Bi11 = si4 - tr4;\
	__Br12 = sr2 + ti2;		__Bi12 = si2 - tr2;\
	__Br13 = srd + tid;		__Bi13 = sid - trd;\
	__Br14 = srb - tib;		__Bi14 = sib + trb;\
	__Br15 = sr3 + ti3;		__Bi15 = si3 - tr3;\
	__Br16 = sr3 - ti3;		__Bi16 = si3 + tr3;\
	__Br17 = srb + tib;		__Bi17 = sib - trb;\
	__Br18 = srd - tid;		__Bi18 = sid + trd;\
	__Br19 = sr2 - ti2;		__Bi19 = si2 + tr2;\
	__Br20 = sr4 - ti4;		__Bi20 = si4 + tr4;\
	__Br21 = sr7 + ti7;		__Bi21 = si7 - tr7;\
	__Br22 = sr1 + ti1;		__Bi22 = si1 - tr1;\
	__Br23 = sr6 + ti6;		__Bi23 = si6 - tr6;\
	__Br24 = sre + tie;		__Bi24 = sie - tre;\
	__Br25 = sr5 - ti5;		__Bi25 = si5 + tr5;\
	__Br26 = sra + tia;		__Bi26 = sia - tra;\
	__Br27 = sr9 + ti9;		__Bi27 = si9 - tr9;\
	__Br28 = sr8 - ti8;		__Bi28 = si8 + tr8;\
	__Br29 = src + tic;		__Bi29 = sic - trc;\
	__Br30 = srf + tif;		__Bi30 = sif - trf;\
/* Totals: 658 FADD, 234 FMUL. (Compare to 376 FADD, 88 FMUL for radix-32. (1.7x the ADDs per point, 2.6x the MUL per point...Ugh!) */\
}

/* Optimization Step 0: Separate Real and Imag-output instruction streams - This cuts the number of temps (aside from those we use
to store copies of the real parts of the inputs in order to support in-placeness) in half, since we re-use the same temps for bath halves: */
#define RADIX_31_DFT_0(\
__Ar00,__Ar01,__Ar02,__Ar03,__Ar04,__Ar05,__Ar06,__Ar07,__Ar08,__Ar09,__Ar10,__Ar11,__Ar12,__Ar13,__Ar14,__Ar15,__Ar16,__Ar17,__Ar18,__Ar19,__Ar20,__Ar21,__Ar22,__Ar23,__Ar24,__Ar25,__Ar26,__Ar27,__Ar28,__Ar29,__Ar30,\
__Ai00,__Ai01,__Ai02,__Ai03,__Ai04,__Ai05,__Ai06,__Ai07,__Ai08,__Ai09,__Ai10,__Ai11,__Ai12,__Ai13,__Ai14,__Ai15,__Ai16,__Ai17,__Ai18,__Ai19,__Ai20,__Ai21,__Ai22,__Ai23,__Ai24,__Ai25,__Ai26,__Ai27,__Ai28,__Ai29,__Ai30,\
__Br00,__Br01,__Br02,__Br03,__Br04,__Br05,__Br06,__Br07,__Br08,__Br09,__Br10,__Br11,__Br12,__Br13,__Br14,__Br15,__Br16,__Br17,__Br18,__Br19,__Br20,__Br21,__Br22,__Br23,__Br24,__Br25,__Br26,__Br27,__Br28,__Br29,__Br30,\
__Bi00,__Bi01,__Bi02,__Bi03,__Bi04,__Bi05,__Bi06,__Bi07,__Bi08,__Bi09,__Bi10,__Bi11,__Bi12,__Bi13,__Bi14,__Bi15,__Bi16,__Bi17,__Bi18,__Bi19,__Bi20,__Bi21,__Bi22,__Bi23,__Bi24,__Bi25,__Bi26,__Bi27,__Bi28,__Bi29,__Bi30\
)\
{\
	double __tr00,__tr01,__tr02,__tr03,__tr04,__tr05,__tr06,__tr07,__tr08,__tr09,__tr10,__tr11,__tr12,__tr13,__tr14,__tr15,__tr16,__tr17,__tr18,__tr19,__tr20,__tr21,__tr22,__tr23,__tr24,__tr25,__tr26,__tr27,__tr28,__tr29,__tr30;\
	double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf, y1, y2, y3, y4, y5, y6, y7, y8, y9, ya, yb, yc, yd, ye, yf, \
	x16b, x49e, x27c, x5af, x38d, x7_2, x8_3, xa_5, xb_6, xd_8, xe_9, x1_b, x2_c, x4_e, x5_f, x16b49e, x27c5af, xsigma, usigma, \
	u1234, u12, u34, u1, u2, u3, u4, u5, xa4_5e, xae_54, xd1_8b, xdb_81, xa4d1_5e8b, p6, p8a, p8, pa, pc, pe, xb5_6f, xbf_65, xe2_9c, xec_92, \
	xb5e2_6f9c, p7, p9b, p9, pb, pd, pf, u6, u7, p89, pcd, u89, u8, u9, ucd, uc, ud, pab, pef, uab, ua, ub, uef, ue, uf, u8a, \
	w6, w8eac, w8e, wac, w8_e, wa_c, w8, we, wa, wc, u9b, w7, w9fbd, w9f, wbd, w9_f, wb_d, w9, wf, wb, wd, \
	s1, s2, s3, s4, s5, s6, s7, s8, s9, sa, sb, sc, sd, se, sf, \
	y16b, y49e, y27c, y5af, y38d, y7_2, y8_3, ya_5, yb_6, yd_8, ye_9, y1_b, y2_c, y4_e, y5_f, y16b49e, y27c5af, ysigma, \
	vsigma, v1234, v12, v34, v1, v2, v3, v4, v5, ya4_5e, yae_54, yd1_8b, ydb_81, ya4d1_5e8b, \
	q6, q8a, q8, qa, qc, qe, yb5_6f, ybf_65, ye2_9c, yec_92, yb5e2_6f9c, q7, q9b, q9, qb, qd, qf, v6, v7, q89, qcd, \
	v89, v8, v9, vcd, vc, vd, qab, qef, vab, va, vb, vef, ve, vf, v8a, \
	z6, z8eac, z8e, zac, z8_e, za_c, z8, ze, za, zc, v9b, z7, z9fbd, z9f, zbd, z9_f, zb_d, z9, zf, zb, zd, \
	t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb, tc, td, te, tf;\
\
	__tr00 = __Ar00;\
	__tr01 = __Ar01;\
	__tr02 = __Ar02;\
	__tr03 = __Ar03;\
	__tr04 = __Ar04;\
	__tr05 = __Ar05;\
	__tr06 = __Ar06;\
	__tr07 = __Ar07;\
	__tr08 = __Ar08;\
	__tr09 = __Ar09;\
	__tr10 = __Ar10;\
	__tr11 = __Ar11;\
	__tr12 = __Ar12;\
	__tr13 = __Ar13;\
	__tr14 = __Ar14;\
	__tr15 = __Ar15;\
	__tr16 = __Ar16;\
	__tr17 = __Ar17;\
	__tr18 = __Ar18;\
	__tr19 = __Ar19;\
	__tr20 = __Ar20;\
	__tr21 = __Ar21;\
	__tr22 = __Ar22;\
	__tr23 = __Ar23;\
	__tr24 = __Ar24;\
	__tr25 = __Ar25;\
	__tr26 = __Ar26;\
	__tr27 = __Ar27;\
	__tr28 = __Ar28;\
	__tr29 = __Ar29;\
	__tr30 = __Ar30;\
\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms: */\
	x0 = __Ar00;\
	x1 = __Ar01 + __Ar30;\
	x2 = __Ar07 + __Ar24;\
	x3 = __Ar18 + __Ar13;\
	x4 = __Ar02 + __Ar29;\
	x5 = __Ar14 + __Ar17;\
	x6 = __Ar05 + __Ar26;\
	x7 = __Ar04 + __Ar27;\
	x8 = __Ar28 + __Ar03;\
	x9 = __Ar10 + __Ar21;\
	xa = __Ar08 + __Ar23;\
	xb = __Ar25 + __Ar06;\
	xc = __Ar20 + __Ar11;\
	xd = __Ar16 + __Ar15;\
	xe = __Ar19 + __Ar12;\
	xf = __Ar09 + __Ar22;\
\
	x16b = x1 + x6 + xb;\
	x49e = x4 + x9 + xe;\
	x27c = x2 + x7 + xc;\
	x5af = x5 + xa + xf;\
	x38d = x3 + x8 + xd;\
	x7_2 = x7 - x2;\
	x8_3 = x8 - x3;\
	xa_5 = xa - x5;\
	xb_6 = xb - x6;\
	xd_8 = xd - x8;\
	xe_9 = xe - x9;\
	x1_b = x1 - xb;\
	x2_c = x2 - xc;\
	x4_e = x4 - xe;\
	x5_f = x5 - xf;\
	x16b = x16b - x38d;\
	x49e = x49e - x38d;\
	x27c = x27c - x38d;\
	x5af = x5af - x38d;\
	x16b49e = x16b + x49e;\
	x27c5af = x27c + x5af;\
	xsigma = x16b49e + x27c5af;\
	usigma = xsigma + 5*x38d;\
/* DC component: */\
	__Br00 = x0 + usigma;\
	u1234 = x0 + dc0*x38d + dc1*xsigma;\
	u12 = u1234 + dc2*x27c5af;\
	u34 = u1234 + dc7*x16b49e;\
	u1 = u12 + dc3*x49e + dc4*x5af;\
	u2 = u12 + dc5*x16b + dc6*x27c;\
	u3 = u34 + dc8*x49e + dc3*x5af;\
	u4 = u34 - dc4*x16b + dc5*x27c;\
	u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0;\
	xa4_5e = xa_5 + x4_e;\
	xae_54 = xa_5 - x4_e;\
	xd1_8b = xd_8 + x1_b;\
	xdb_81 = xd_8 - x1_b;\
	xa4d1_5e8b = xa4_5e + xd1_8b;\
	p6 = x7_2 + xa4d1_5e8b;\
	p8a = c50*x7_2 + c51*xa4d1_5e8b;\
	p8 = p8a + xd1_8b;\
	pa = p8a + xa4_5e;\
	pc = xae_54 + s52*xdb_81;\
	pe = s52*xae_54 - xdb_81;\
	xb5_6f = xb_6 + x5_f;\
	xbf_65 = xb_6 - x5_f;\
	xe2_9c = xe_9 + x2_c;\
	xec_92 = xe_9 - x2_c;\
	xb5e2_6f9c = xb5_6f + xe2_9c;\
	p7 = x8_3 + xb5e2_6f9c;\
	p9b = c50*x8_3 + c51*xb5e2_6f9c;\
	p9 = p9b + xe2_9c;\
	pb = p9b + xb5_6f;\
	pd = xbf_65 + s52*xec_92;\
	pf = s52*xbf_65 - xec_92;\
	u6 = dc9*p6 + dc10*p7;\
	u7 = dc11*p6 + dc9*p7;\
	p89 = p8 + p9;\
	pcd = pc + pd;\
	u89 = dc12*p89 + dc15*pcd;\
	u8 = u89 + dc13*p9 + dc16*pd;\
	u9 = u89 + dc14*p8 + dc17*pc;\
	ucd = dc18*pcd - dc15*p89;\
	uc = ucd - dc16*p9 + dc19*pd;\
	ud = ucd - dc17*p8 + dc20*pc;\
	pab = pa + pb;\
	pef = pe + pf;\
	uab = dc21*pab + dc24*pef;\
	ua = uab + dc22*pb + dc25*pf;\
	ub = uab + dc23*pa + dc26*pe;\
	uef = dc27*pef - dc24*pab;\
	ue = uef - dc25*pb + dc28*pf;\
	uf = uef - dc26*pa + dc29*pe;\
	u8a = u8 + ua;\
	w6 = u6 + c50*u8a;\
	w8eac = u6 + c51*u8a;\
	w8e = w8eac + ua;\
	wac = w8eac + u8;\
	w8_e = uc + s52*ue;\
	wa_c = s52*uc - ue;\
	w8 = w8e + w8_e;\
	we = w8e - w8_e;\
	wa = wac + wa_c;\
	wc = wac - wa_c;\
	u9b = u9 + ub;\
	w7 = u7 + c50*u9b;\
	w9fbd = u7 + c51*u9b;\
	w9f = w9fbd + ub;\
	wbd = w9fbd + u9;\
	w9_f = ud + s52*uf;\
	wb_d = s52*ud - uf;\
	w9 = w9f + w9_f;\
	wf = w9f - w9_f;\
	wb = wbd + wb_d;\
	wd = wbd - wb_d;\
	s1 = u1 + wc;\
	s2 = u3 - w6 + wd;\
	s3 = u5 - w7;\
	s4 = u2 + we;\
	s5 = u4 - w8 + wf;\
	s6 = u1 - w9;\
	s7 = u3 + w6;\
	s8 = u5 + w7 - wa;\
	s9 = u2 - wb;\
	sa = u4 + w8;\
	sb = u1 + w9 - wc;\
	sc = u3 - wd;\
	sd = u5 + wa;\
	se = u2 + wb - we;\
	sf = u4 - wf;\
\
/* yi-terms: */\
	y1 = __Ai01 - __Ai30;\
	y2 = __Ai07 - __Ai24;\
	y3 = __Ai18 - __Ai13;\
	y4 = __Ai02 - __Ai29;\
	y5 = __Ai14 - __Ai17;\
	y6 = __Ai05 - __Ai26;\
	y7 = __Ai04 - __Ai27;\
	y8 = __Ai28 - __Ai03;\
	y9 = __Ai10 - __Ai21;\
	ya = __Ai08 - __Ai23;\
	yb = __Ai25 - __Ai06;\
	yc = __Ai20 - __Ai11;\
	yd = __Ai16 - __Ai15;\
	ye = __Ai19 - __Ai12;\
	yf = __Ai09 - __Ai22;\
\
	y16b = y1 + y6 + yb;\
	y49e = y4 + y9 + ye;\
	y27c = y2 + y7 + yc;\
	y5af = y5 + ya + yf;\
	y38d = y3 + y8 + yd;\
	y7_2 = y7 - y2;\
	y8_3 = y8 - y3;\
	ya_5 = ya - y5;\
	yb_6 = yb - y6;\
	yd_8 = yd - y8;\
	ye_9 = ye - y9;\
	y1_b = y1 - yb;\
	y2_c = y2 - yc;\
	y4_e = y4 - ye;\
	y5_f = y5 - yf;\
	y16b = y16b - y38d;\
	y49e = y49e - y38d;\
	y27c = y27c - y38d;\
	y5af = y5af - y38d;\
	y16b49e = y16b + y49e;\
	y27c5af = y27c + y5af;\
	ysigma = y16b49e + y27c5af;\
	vsigma = ysigma + 5*y38d;\
	v1234 = ds0*y38d + ds1*ysigma;\
	v12 = v1234 + ds2*y27c5af;\
	v34 = v1234 + ds7*y16b49e;\
	v1 = v12 + ds3*y49e + ds4*y5af;\
	v2 = v12 + ds5*y16b + ds6*y27c;\
	v3 = v34 + ds8*y49e + ds3*y5af;\
	v4 = v34 - ds4*y16b + ds5*y27c;\
	v5 = ds0*vsigma - (v1 + v2 + v3 + v4);\
	ya4_5e = ya_5 + y4_e;\
	yae_54 = ya_5 - y4_e;\
	yd1_8b = yd_8 + y1_b;\
	ydb_81 = yd_8 - y1_b;\
	ya4d1_5e8b = ya4_5e + yd1_8b;\
	q6 = y7_2 + ya4d1_5e8b;\
	q8a = c50*y7_2 + c51*ya4d1_5e8b;\
	q8 = q8a + yd1_8b;\
	qa = q8a + ya4_5e;\
	qc = yae_54 + s52*ydb_81;\
	qe = s52*yae_54 - ydb_81;\
	yb5_6f = yb_6 + y5_f;\
	ybf_65 = yb_6 - y5_f;\
	ye2_9c = ye_9 + y2_c;\
	yec_92 = ye_9 - y2_c;\
	yb5e2_6f9c = yb5_6f + ye2_9c;\
	q7 = y8_3 + yb5e2_6f9c;\
	q9b = c50*y8_3 + c51*yb5e2_6f9c;\
	q9 = q9b + ye2_9c;\
	qb = q9b + yb5_6f;\
	qd = ybf_65 + s52*yec_92;\
	qf = s52*ybf_65 - yec_92;\
	v6 = ds9*q6 + ds10*q7;\
	v7 = ds11*q6 + ds9*q7;\
	q89 = q8 + q9;\
	qcd = qc + qd;\
	v89 = ds12*q89 + ds15*qcd;\
	v8 = v89 + ds13*q9 + ds16*qd;\
	v9 = v89 + ds14*q8 + ds17*qc;\
	vcd = ds18*qcd - ds15*q89;\
	vc = vcd - ds16*q9 + ds19*qd;\
	vd = vcd - ds17*q8 + ds20*qc;\
	qab = qa + qb;\
	qef = qe + qf;\
	vab = ds21*qab + ds24*qef;\
	va = vab + ds22*qb + ds25*qf;\
	vb = vab + ds23*qa + ds26*qe;\
	vef = ds27*qef - ds24*qab;\
	ve = vef - ds25*qb + ds28*qf;\
	vf = vef - ds26*qa + ds29*qe;\
	v8a = v8 + va;\
	z6 = v6 + c50*v8a;\
	z8eac = v6 + c51*v8a;\
	z8e = z8eac + va;\
	zac = z8eac + v8;\
	z8_e = vc + s52*ve;\
	za_c = s52*vc - ve;\
	z8 = z8e + z8_e;\
	ze = z8e - z8_e;\
	za = zac + za_c;\
	zc = zac - za_c;\
	v9b = v9 + vb;\
	z7 = v7 + c50*v9b;\
	z9fbd = v7 + c51*v9b;\
	z9f = z9fbd + vb;\
	zbd = z9fbd + v9;\
	z9_f = vd + s52*vf;\
	zb_d = s52*vd - vf;\
	z9 = z9f + z9_f;\
	zf = z9f - z9_f;\
	zb = zbd + zb_d;\
	zd = zbd - zb_d;\
	t1 = v1 + zc;\
	t2 = v3 - z6 + zd;\
	t3 = v5 - z7;\
	t4 = v2 + ze;\
	t5 = v4 - z8 + zf;\
	t6 = v1 - z9;\
	t7 = v3 + z6;\
	t8 = v5 + z7 - za;\
	t9 = v2 - zb;\
	ta = v4 + z8;\
	tb = v1 + z9 - zc;\
	tc = v3 - zd;\
	td = v5 + za;\
	te = v2 + zb - ze;\
	tf = v4 - zf;\
\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
	__Br01 = sf - tf;\
	__Br02 = sc - tc;\
	__Br03 = s8 + t8;\
	__Br04 = s9 - t9;\
	__Br05 = sa - ta;\
	__Br06 = s5 + t5;\
	__Br07 = se - te;\
	__Br08 = s6 - t6;\
	__Br09 = s1 - t1;\
	__Br10 = s7 - t7;\
	__Br11 = s4 + t4;\
	__Br12 = s2 + t2;\
	__Br13 = sd + td;\
	__Br14 = sb - tb;\
	__Br15 = s3 + t3;\
	__Br16 = s3 - t3;\
	__Br17 = sb + tb;\
	__Br18 = sd - td;\
	__Br19 = s2 - t2;\
	__Br20 = s4 - t4;\
	__Br21 = s7 + t7;\
	__Br22 = s1 + t1;\
	__Br23 = s6 + t6;\
	__Br24 = se + te;\
	__Br25 = s5 - t5;\
	__Br26 = sa + ta;\
	__Br27 = s9 + t9;\
	__Br28 = s8 - t8;\
	__Br29 = sc + tc;\
	__Br30 = sf + tf;\
\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Replace __B[j] with __B[13-j] for j > 0: */\
	/***************/\
	/* xi-terms: */\
	x0 = __Ai00;\
	x1 = __Ai01 + __Ai30;\
	x2 = __Ai07 + __Ai24;\
	x3 = __Ai18 + __Ai13;\
	x4 = __Ai02 + __Ai29;\
	x5 = __Ai14 + __Ai17;\
	x6 = __Ai05 + __Ai26;\
	x7 = __Ai04 + __Ai27;\
	x8 = __Ai28 + __Ai03;\
	x9 = __Ai10 + __Ai21;\
	xa = __Ai08 + __Ai23;\
	xb = __Ai25 + __Ai06;\
	xc = __Ai20 + __Ai11;\
	xd = __Ai16 + __Ai15;\
	xe = __Ai19 + __Ai12;\
	xf = __Ai09 + __Ai22;\
\
	x16b = x1 + x6 + xb;\
	x49e = x4 + x9 + xe;\
	x27c = x2 + x7 + xc;\
	x5af = x5 + xa + xf;\
	x38d = x3 + x8 + xd;\
	x7_2 = x7 - x2;\
	x8_3 = x8 - x3;\
	xa_5 = xa - x5;\
	xb_6 = xb - x6;\
	xd_8 = xd - x8;\
	xe_9 = xe - x9;\
	x1_b = x1 - xb;\
	x2_c = x2 - xc;\
	x4_e = x4 - xe;\
	x5_f = x5 - xf;\
	x16b = x16b - x38d;\
	x49e = x49e - x38d;\
	x27c = x27c - x38d;\
	x5af = x5af - x38d;\
	x16b49e = x16b + x49e;\
	x27c5af = x27c + x5af;\
	xsigma = x16b49e + x27c5af;\
	usigma = xsigma + 5*x38d;\
/* DC component: */\
	__Bi00 = x0 + usigma;\
	u1234 = x0 + dc0*x38d + dc1*xsigma;\
	u12 = u1234 + dc2*x27c5af;\
	u34 = u1234 + dc7*x16b49e;\
	u1 = u12 + dc3*x49e + dc4*x5af;\
	u2 = u12 + dc5*x16b + dc6*x27c;\
	u3 = u34 + dc8*x49e + dc3*x5af;\
	u4 = u34 - dc4*x16b + dc5*x27c;\
	u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0;\
	xa4_5e = xa_5 + x4_e;\
	xae_54 = xa_5 - x4_e;\
	xd1_8b = xd_8 + x1_b;\
	xdb_81 = xd_8 - x1_b;\
	xa4d1_5e8b = xa4_5e + xd1_8b;\
	p6 = x7_2 + xa4d1_5e8b;\
	p8a = c50*x7_2 + c51*xa4d1_5e8b;\
	p8 = p8a + xd1_8b;\
	pa = p8a + xa4_5e;\
	pc = xae_54 + s52*xdb_81;\
	pe = s52*xae_54 - xdb_81;\
	xb5_6f = xb_6 + x5_f;\
	xbf_65 = xb_6 - x5_f;\
	xe2_9c = xe_9 + x2_c;\
	xec_92 = xe_9 - x2_c;\
	xb5e2_6f9c = xb5_6f + xe2_9c;\
	p7 = x8_3 + xb5e2_6f9c;\
	p9b = c50*x8_3 + c51*xb5e2_6f9c;\
	p9 = p9b + xe2_9c;\
	pb = p9b + xb5_6f;\
	pd = xbf_65 + s52*xec_92;\
	pf = s52*xbf_65 - xec_92;\
	u6 = dc9*p6 + dc10*p7;\
	u7 = dc11*p6 + dc9*p7;\
	p89 = p8 + p9;\
	pcd = pc + pd;\
	u89 = dc12*p89 + dc15*pcd;\
	u8 = u89 + dc13*p9 + dc16*pd;\
	u9 = u89 + dc14*p8 + dc17*pc;\
	ucd = dc18*pcd - dc15*p89;\
	uc = ucd - dc16*p9 + dc19*pd;\
	ud = ucd - dc17*p8 + dc20*pc;\
	pab = pa + pb;\
	pef = pe + pf;\
	uab = dc21*pab + dc24*pef;\
	ua = uab + dc22*pb + dc25*pf;\
	ub = uab + dc23*pa + dc26*pe;\
	uef = dc27*pef - dc24*pab;\
	ue = uef - dc25*pb + dc28*pf;\
	uf = uef - dc26*pa + dc29*pe;\
	u8a = u8 + ua;\
	w6 = u6 + c50*u8a;\
	w8eac = u6 + c51*u8a;\
	w8e = w8eac + ua;\
	wac = w8eac + u8;\
	w8_e = uc + s52*ue;\
	wa_c = s52*uc - ue;\
	w8 = w8e + w8_e;\
	we = w8e - w8_e;\
	wa = wac + wa_c;\
	wc = wac - wa_c;\
	u9b = u9 + ub;\
	w7 = u7 + c50*u9b;\
	w9fbd = u7 + c51*u9b;\
	w9f = w9fbd + ub;\
	wbd = w9fbd + u9;\
	w9_f = ud + s52*uf;\
	wb_d = s52*ud - uf;\
	w9 = w9f + w9_f;\
	wf = w9f - w9_f;\
	wb = wbd + wb_d;\
	wd = wbd - wb_d;\
	s1 = u1 + wc;\
	s2 = u3 - w6 + wd;\
	s3 = u5 - w7;\
	s4 = u2 + we;\
	s5 = u4 - w8 + wf;\
	s6 = u1 - w9;\
	s7 = u3 + w6;\
	s8 = u5 + w7 - wa;\
	s9 = u2 - wb;\
	sa = u4 + w8;\
	sb = u1 + w9 - wc;\
	sc = u3 - wd;\
	sd = u5 + wa;\
	se = u2 + wb - we;\
	sf = u4 - wf;\
\
/* yr-terms: */\
	y1 = __tr01 - __tr30;\
	y2 = __tr07 - __tr24;\
	y3 = __tr18 - __tr13;\
	y4 = __tr02 - __tr29;\
	y5 = __tr14 - __tr17;\
	y6 = __tr05 - __tr26;\
	y7 = __tr04 - __tr27;\
	y8 = __tr28 - __tr03;\
	y9 = __tr10 - __tr21;\
	ya = __tr08 - __tr23;\
	yb = __tr25 - __tr06;\
	yc = __tr20 - __tr11;\
	yd = __tr16 - __tr15;\
	ye = __tr19 - __tr12;\
	yf = __tr09 - __tr22;\
\
	y16b = y1 + y6 + yb;\
	y49e = y4 + y9 + ye;\
	y27c = y2 + y7 + yc;\
	y5af = y5 + ya + yf;\
	y38d = y3 + y8 + yd;\
	y7_2 = y7 - y2;\
	y8_3 = y8 - y3;\
	ya_5 = ya - y5;\
	yb_6 = yb - y6;\
	yd_8 = yd - y8;\
	ye_9 = ye - y9;\
	y1_b = y1 - yb;\
	y2_c = y2 - yc;\
	y4_e = y4 - ye;\
	y5_f = y5 - yf;\
	y16b = y16b - y38d;\
	y49e = y49e - y38d;\
	y27c = y27c - y38d;\
	y5af = y5af - y38d;\
	y16b49e = y16b + y49e;\
	y27c5af = y27c + y5af;\
	ysigma = y16b49e + y27c5af;\
	vsigma = ysigma + 5*y38d;\
	v1234 = ds0*y38d + ds1*ysigma;\
	v12 = v1234 + ds2*y27c5af;\
	v34 = v1234 + ds7*y16b49e;\
	v1 = v12 + ds3*y49e + ds4*y5af;\
	v2 = v12 + ds5*y16b + ds6*y27c;\
	v3 = v34 + ds8*y49e + ds3*y5af;\
	v4 = v34 - ds4*y16b + ds5*y27c;\
	v5 = ds0*vsigma - (v1 + v2 + v3 + v4);\
	ya4_5e = ya_5 + y4_e;\
	yae_54 = ya_5 - y4_e;\
	yd1_8b = yd_8 + y1_b;\
	ydb_81 = yd_8 - y1_b;\
	ya4d1_5e8b = ya4_5e + yd1_8b;\
	q6 = y7_2 + ya4d1_5e8b;\
	q8a = c50*y7_2 + c51*ya4d1_5e8b;\
	q8 = q8a + yd1_8b;\
	qa = q8a + ya4_5e;\
	qc = yae_54 + s52*ydb_81;\
	qe = s52*yae_54 - ydb_81;\
	yb5_6f = yb_6 + y5_f;\
	ybf_65 = yb_6 - y5_f;\
	ye2_9c = ye_9 + y2_c;\
	yec_92 = ye_9 - y2_c;\
	yb5e2_6f9c = yb5_6f + ye2_9c;\
	q7 = y8_3 + yb5e2_6f9c;\
	q9b = c50*y8_3 + c51*yb5e2_6f9c;\
	q9 = q9b + ye2_9c;\
	qb = q9b + yb5_6f;\
	qd = ybf_65 + s52*yec_92;\
	qf = s52*ybf_65 - yec_92;\
	v6 = ds9*q6 + ds10*q7;\
	v7 = ds11*q6 + ds9*q7;\
	q89 = q8 + q9;\
	qcd = qc + qd;\
	v89 = ds12*q89 + ds15*qcd;\
	v8 = v89 + ds13*q9 + ds16*qd;\
	v9 = v89 + ds14*q8 + ds17*qc;\
	vcd = ds18*qcd - ds15*q89;\
	vc = vcd - ds16*q9 + ds19*qd;\
	vd = vcd - ds17*q8 + ds20*qc;\
	qab = qa + qb;\
	qef = qe + qf;\
	vab = ds21*qab + ds24*qef;\
	va = vab + ds22*qb + ds25*qf;\
	vb = vab + ds23*qa + ds26*qe;\
	vef = ds27*qef - ds24*qab;\
	ve = vef - ds25*qb + ds28*qf;\
	vf = vef - ds26*qa + ds29*qe;\
	v8a = v8 + va;\
	z6 = v6 + c50*v8a;\
	z8eac = v6 + c51*v8a;\
	z8e = z8eac + va;\
	zac = z8eac + v8;\
	z8_e = vc + s52*ve;\
	za_c = s52*vc - ve;\
	z8 = z8e + z8_e;\
	ze = z8e - z8_e;\
	za = zac + za_c;\
	zc = zac - za_c;\
	v9b = v9 + vb;\
	z7 = v7 + c50*v9b;\
	z9fbd = v7 + c51*v9b;\
	z9f = z9fbd + vb;\
	zbd = z9fbd + v9;\
	z9_f = vd + s52*vf;\
	zb_d = s52*vd - vf;\
	z9 = z9f + z9_f;\
	zf = z9f - z9_f;\
	zb = zbd + zb_d;\
	zd = zbd - zb_d;\
	t1 = v1 + zc;\
	t2 = v3 - z6 + zd;\
	t3 = v5 - z7;\
	t4 = v2 + ze;\
	t5 = v4 - z8 + zf;\
	t6 = v1 - z9;\
	t7 = v3 + z6;\
	t8 = v5 + z7 - za;\
	t9 = v2 - zb;\
	ta = v4 + z8;\
	tb = v1 + z9 - zc;\
	tc = v3 - zd;\
	td = v5 + za;\
	te = v2 + zb - ze;\
	tf = v4 - zf;\
\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
	__Bi01 = sf + tf;\
	__Bi02 = sc + tc;\
	__Bi03 = s8 - t8;\
	__Bi04 = s9 + t9;\
	__Bi05 = sa + ta;\
	__Bi06 = s5 - t5;\
	__Bi07 = se + te;\
	__Bi08 = s6 + t6;\
	__Bi09 = s1 + t1;\
	__Bi10 = s7 + t7;\
	__Bi11 = s4 - t4;\
	__Bi12 = s2 - t2;\
	__Bi13 = sd - td;\
	__Bi14 = sb + tb;\
	__Bi15 = s3 - t3;\
	__Bi16 = s3 + t3;\
	__Bi17 = sb - tb;\
	__Bi18 = sd + td;\
	__Bi19 = s2 + t2;\
	__Bi20 = s4 + t4;\
	__Bi21 = s7 - t7;\
	__Bi22 = s1 - t1;\
	__Bi23 = s6 - t6;\
	__Bi24 = se - te;\
	__Bi25 = s5 + t5;\
	__Bi26 = sa - ta;\
	__Bi27 = s9 - t9;\
	__Bi28 = s8 + t8;\
	__Bi29 = sc - tc;\
	__Bi30 = sf - tf;\
\
/* Totals: 658 FADD, 234 FMUL. */\
}

/* Optimization Step 1: Analyze dependencies among temporaries to reduce needed number of same: */
#define RADIX_31_DFT_1(\
__Ar00,__Ar01,__Ar02,__Ar03,__Ar04,__Ar05,__Ar06,__Ar07,__Ar08,__Ar09,__Ar10,__Ar11,__Ar12,__Ar13,__Ar14,__Ar15,__Ar16,__Ar17,__Ar18,__Ar19,__Ar20,__Ar21,__Ar22,__Ar23,__Ar24,__Ar25,__Ar26,__Ar27,__Ar28,__Ar29,__Ar30,\
__Ai00,__Ai01,__Ai02,__Ai03,__Ai04,__Ai05,__Ai06,__Ai07,__Ai08,__Ai09,__Ai10,__Ai11,__Ai12,__Ai13,__Ai14,__Ai15,__Ai16,__Ai17,__Ai18,__Ai19,__Ai20,__Ai21,__Ai22,__Ai23,__Ai24,__Ai25,__Ai26,__Ai27,__Ai28,__Ai29,__Ai30,\
__Br00,__Br01,__Br02,__Br03,__Br04,__Br05,__Br06,__Br07,__Br08,__Br09,__Br10,__Br11,__Br12,__Br13,__Br14,__Br15,__Br16,__Br17,__Br18,__Br19,__Br20,__Br21,__Br22,__Br23,__Br24,__Br25,__Br26,__Br27,__Br28,__Br29,__Br30,\
__Bi00,__Bi01,__Bi02,__Bi03,__Bi04,__Bi05,__Bi06,__Bi07,__Bi08,__Bi09,__Bi10,__Bi11,__Bi12,__Bi13,__Bi14,__Bi15,__Bi16,__Bi17,__Bi18,__Bi19,__Bi20,__Bi21,__Bi22,__Bi23,__Bi24,__Bi25,__Bi26,__Bi27,__Bi28,__Bi29,__Bi30\
)\
{\
	double TMP01,TMP02;\
	double __tr01,__tr02,__tr03,__tr04,__tr05,__tr06,__tr07,__tr08,__tr09,__tr10,__tr11,__tr12,__tr13,__tr14,__tr15,\
	       __tr16,__tr17,__tr18,__tr19,__tr20,__tr21,__tr22,__tr23,__tr24,__tr25,__tr26,__tr27,__tr28,__tr29,__tr30;\
	double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;\
\
	/* To support in-place DFT, need to save copies of input real parts: */\
	__tr01 = __Ar01;\
	__tr02 = __Ar02;\
	__tr03 = __Ar03;\
	__tr04 = __Ar04;\
	__tr05 = __Ar05;\
	__tr06 = __Ar06;\
	__tr07 = __Ar07;\
	__tr08 = __Ar08;\
	__tr09 = __Ar09;\
	__tr10 = __Ar10;\
	__tr11 = __Ar11;\
	__tr12 = __Ar12;\
	__tr13 = __Ar13;\
	__tr14 = __Ar14;\
	__tr15 = __Ar15;\
	__tr16 = __Ar16;\
	__tr17 = __Ar17;\
	__tr18 = __Ar18;\
	__tr19 = __Ar19;\
	__tr20 = __Ar20;\
	__tr21 = __Ar21;\
	__tr22 = __Ar22;\
	__tr23 = __Ar23;\
	__tr24 = __Ar24;\
	__tr25 = __Ar25;\
	__tr26 = __Ar26;\
	__tr27 = __Ar27;\
	__tr28 = __Ar28;\
	__tr29 = __Ar29;\
	__tr30 = __Ar30;\
\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms: */\
	__Br00 = __Ar00;	/* x0; Free up by using output 0 in its place here */\
	x1 = __Ar01;\
	x2 = __Ar07;\
	x3 = __Ar18;\
	x4 = __Ar02;\
	x5 = __Ar14;\
	x6 = __Ar05;\
	x7 = __Ar04;\
	x8 = __Ar28;\
	x9 = __Ar10;\
	xa = __Ar08;\
	xb = __Ar25;\
	xc = __Ar20;\
	xd = __Ar16;\
	xe = __Ar19;\
	xf = __Ar09;		/* xr-terms: */\
	x1 += __Ar30;		/* x1 = __Ar01 + __Ar30 */\
	x2 += __Ar24;		/* x2 = __Ar07 + __Ar24 */\
	x3 += __Ar13;		/* x3 = __Ar18 + __Ar13 */\
	x4 += __Ar29;		/* x4 = __Ar02 + __Ar29 */\
	x5 += __Ar17;		/* x5 = __Ar14 + __Ar17 */\
	x6 += __Ar26;		/* x6 = __Ar05 + __Ar26 */\
	x7 += __Ar27;		/* x7 = __Ar04 + __Ar27 */\
	x8 += __Ar03;		/* x8 = __Ar28 + __Ar03 */\
	x9 += __Ar21;		/* x9 = __Ar10 + __Ar21 */\
	xa += __Ar23;		/* xa = __Ar08 + __Ar23 */\
	xb += __Ar06;		/* xb = __Ar25 + __Ar06 */\
	xc += __Ar11;		/* xc = __Ar20 + __Ar11 */\
	xd += __Ar15;		/* xd = __Ar16 + __Ar15 */\
	xe += __Ar12;		/* xe = __Ar19 + __Ar12 */\
	xf += __Ar22;		/* xf = __Ar09 + __Ar22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Br00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += __Br00 ;		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += __Br00;		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	__Br00 += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4);\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	__Br08= xe;	/* s6 = xe - x4 */\
	__Br14= x4;	/* sb = xe + x4 - xc */\
	__Br09= xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	__Br02= x5;	/* sc = x5 - xd */\
	__Br12= xd;	/* s2 = x5 + xd - x7 */\
	__Br10= x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	__Br15= x8;	/* s3 = x8 - x3 */\
	__Br03= x3;	/* s8 = x8 + x3 - x9 */\
	__Br13= x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	__Br04= xb;	/* s9 = xb - x1 */\
	__Br07= x1;	/* se = xb + x1 - xf */\
	__Br11= xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	__Br01= x2;	/* sf = x2 - xa */\
	__Br06= xa;	/* s5 = x2 + xa - x6 */\
	__Br05= x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1 = __Ai01;\
	x2 = __Ai07;\
	x3 = __Ai18;\
	x4 = __Ai02;\
	x5 = __Ai14;\
	x6 = __Ai05;\
	x7 = __Ai04;\
	x8 = __Ai28;\
	x9 = __Ai10;\
	xa = __Ai08;\
	xb = __Ai25;\
	xc = __Ai20;\
	xd = __Ai16;\
	xe = __Ai19;\
	xf = __Ai09;\
	x1 -= __Ai30;\
	x2 -= __Ai24;\
	x3 -= __Ai13;\
	x4 -= __Ai29;\
	x5 -= __Ai17;\
	x6 -= __Ai26;\
	x7 -= __Ai27;\
	x8 -= __Ai03;\
	x9 -= __Ai21;\
	xa -= __Ai23;\
	xb -= __Ai06;\
	xc -= __Ai11;\
	xd -= __Ai15;\
	xe -= __Ai12;\
	xf -= __Ai22;\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4);\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = __Br01;	x0 -= x2;	__Br01 = x0;	x2 *= 2;	x0 += x2;	__Br30 = x0;\
 	x0 = __Br02;	x0 -= x5;	__Br02 = x0;	x5 *= 2;	x0 += x5;	__Br29 = x0;\
 	x0 = __Br03;	x0 += x3;	__Br03 = x0;	x3 *= 2;	x0 -= x3;	__Br28 = x0;\
 	x0 = __Br04;	x0 -= xb;	__Br04 = x0;	xb *= 2;	x0 += xb;	__Br27 = x0;\
 	x0 = __Br05;	x0 -= x6;	__Br05 = x0;	x6 *= 2;	x0 += x6;	__Br26 = x0;\
 	x0 = __Br06;	x0 += xa;	__Br06 = x0;	xa *= 2;	x0 -= xa;	__Br25 = x0;\
 	x0 = __Br07;	x0 -= x1;	__Br07 = x0;	x1 *= 2;	x0 += x1;	__Br24 = x0;\
 	x0 = __Br08;	x0 -= xe;	__Br08 = x0;	xe *= 2;	x0 += xe;	__Br23 = x0;\
 	x0 = __Br09;	x0 -= xc;	__Br09 = x0;	xc *= 2;	x0 += xc;	__Br22 = x0;\
 	x0 = __Br10;	x0 -= x7;	__Br10 = x0;	x7 *= 2;	x0 += x7;	__Br21 = x0;\
 	x0 = __Br11;	x0 += xf;	__Br11 = x0;	xf *= 2;	x0 -= xf;	__Br20 = x0;\
 	x0 = __Br12;	x0 += xd;	__Br12 = x0;	xd *= 2;	x0 -= xd;	__Br19 = x0;\
 	x0 = __Br13;	x0 += x9;	__Br13 = x0;	x9 *= 2;	x0 -= x9;	__Br18 = x0;\
 	x0 = __Br14;	x0 -= x4;	__Br14 = x0;	x4 *= 2;	x0 += x4;	__Br17 = x0;\
 	x0 = __Br15;	x0 += x8;	__Br15 = x0;	x8 *= 2;	x0 -= x8;	__Br16 = x0;\
\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Swap +/- in final output sequence: */\
	/***************/\
	/* xi-terms: */\
	__Bi00 = __Ai00;	/* x0; Free up by using output 0 in its place here */\
	x1 = __Ai01;\
	x2 = __Ai07;\
	x3 = __Ai18;\
	x4 = __Ai02;\
	x5 = __Ai14;\
	x6 = __Ai05;\
	x7 = __Ai04;\
	x8 = __Ai28;\
	x9 = __Ai10;\
	xa = __Ai08;\
	xb = __Ai25;\
	xc = __Ai20;\
	xd = __Ai16;\
	xe = __Ai19;\
	xf = __Ai09;		/* xi-terms: */\
	x1 += __Ai30;		/* x1 = __Ai01 + __Ai30 */\
	x2 += __Ai24;		/* x2 = __Ai07 + __Ai24 */\
	x3 += __Ai13;		/* x3 = __Ai18 + __Ai13 */\
	x4 += __Ai29;		/* x4 = __Ai02 + __Ai29 */\
	x5 += __Ai17;		/* x5 = __Ai14 + __Ai17 */\
	x6 += __Ai26;		/* x6 = __Ai05 + __Ai26 */\
	x7 += __Ai27;		/* x7 = __Ai04 + __Ai27 */\
	x8 += __Ai03;		/* x8 = __Ai28 + __Ai03 */\
	x9 += __Ai21;		/* x9 = __Ai10 + __Ai21 */\
	xa += __Ai23;		/* xa = __Ai08 + __Ai23 */\
	xb += __Ai06;		/* xb = __Ai25 + __Ai06 */\
	xc += __Ai11;		/* xc = __Ai20 + __Ai11 */\
	xd += __Ai15;		/* xd = __Ai16 + __Ai15 */\
	xe += __Ai12;		/* xe = __Ai19 + __Ai12 */\
	xf += __Ai22;		/* xf = __Ai09 + __Ai22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Bi00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += __Bi00 ;		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += __Bi00;		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	__Bi00 += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4);\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	__Bi08= xe;	/* s6 = xe - x4 */\
	__Bi14= x4;	/* sb = xe + x4 - xc */\
	__Bi09= xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	__Bi02= x5;	/* sc = x5 - xd */\
	__Bi12= xd;	/* s2 = x5 + xd - x7 */\
	__Bi10= x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	__Bi15= x8;	/* s3 = x8 - x3 */\
	__Bi03= x3;	/* s8 = x8 + x3 - x9 */\
	__Bi13= x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	__Bi04= xb;	/* s9 = xb - x1 */\
	__Bi07= x1;	/* se = xb + x1 - xf */\
	__Bi11= xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	__Bi01= x2;	/* sf = x2 - xa */\
	__Bi06= xa;	/* s5 = x2 + xa - x6 */\
	__Bi05= x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1 = __tr01;\
	x2 = __tr07;\
	x3 = __tr18;\
	x4 = __tr02;\
	x5 = __tr14;\
	x6 = __tr05;\
	x7 = __tr04;\
	x8 = __tr28;\
	x9 = __tr10;\
	xa = __tr08;\
	xb = __tr25;\
	xc = __tr20;\
	xd = __tr16;\
	xe = __tr19;\
	xf = __tr09;\
	x1 -= __tr30;\
	x2 -= __tr24;\
	x3 -= __tr13;\
	x4 -= __tr29;\
	x5 -= __tr17;\
	x6 -= __tr26;\
	x7 -= __tr27;\
	x8 -= __tr03;\
	x9 -= __tr21;\
	xa -= __tr23;\
	xb -= __tr06;\
	xc -= __tr11;\
	xd -= __tr15;\
	xe -= __tr12;\
	xf -= __tr22;\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4);\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = __Bi01;	x0 += x2;	__Bi01 = x0;	x2 *= 2;	x0 -= x2;	__Bi30 = x0;\
 	x0 = __Bi02;	x0 += x5;	__Bi02 = x0;	x5 *= 2;	x0 -= x5;	__Bi29 = x0;\
 	x0 = __Bi03;	x0 -= x3;	__Bi03 = x0;	x3 *= 2;	x0 += x3;	__Bi28 = x0;\
 	x0 = __Bi04;	x0 += xb;	__Bi04 = x0;	xb *= 2;	x0 -= xb;	__Bi27 = x0;\
 	x0 = __Bi05;	x0 += x6;	__Bi05 = x0;	x6 *= 2;	x0 -= x6;	__Bi26 = x0;\
 	x0 = __Bi06;	x0 -= xa;	__Bi06 = x0;	xa *= 2;	x0 += xa;	__Bi25 = x0;\
 	x0 = __Bi07;	x0 += x1;	__Bi07 = x0;	x1 *= 2;	x0 -= x1;	__Bi24 = x0;\
 	x0 = __Bi08;	x0 += xe;	__Bi08 = x0;	xe *= 2;	x0 -= xe;	__Bi23 = x0;\
 	x0 = __Bi09;	x0 += xc;	__Bi09 = x0;	xc *= 2;	x0 -= xc;	__Bi22 = x0;\
 	x0 = __Bi10;	x0 += x7;	__Bi10 = x0;	x7 *= 2;	x0 -= x7;	__Bi21 = x0;\
 	x0 = __Bi11;	x0 -= xf;	__Bi11 = x0;	xf *= 2;	x0 += xf;	__Bi20 = x0;\
 	x0 = __Bi12;	x0 -= xd;	__Bi12 = x0;	xd *= 2;	x0 += xd;	__Bi19 = x0;\
 	x0 = __Bi13;	x0 -= x9;	__Bi13 = x0;	x9 *= 2;	x0 += x9;	__Bi18 = x0;\
 	x0 = __Bi14;	x0 += x4;	__Bi14 = x0;	x4 *= 2;	x0 -= x4;	__Bi17 = x0;\
 	x0 = __Bi15;	x0 -= x8;	__Bi15 = x0;	x8 *= 2;	x0 += x8;	__Bi16 = x0;\
/* Totals: 658 FADD, 234 FMUL. */\
}

/* Optimization Step 2: Expand the small inlined 2x2 mul-macros, try to interleave independent subsequences to reduce serialization,
and modify the parameters list to use pointer base/index-offsets form to reduce the numbers of macro parameters: */
#define RADIX_31_DFT(\
__A,__i00,__i01,__i02,__i03,__i04,__i05,__i06,__i07,__i08,__i09,__i10,__i11,__i12,__i13,__i14,__i15,__i16,__i17,__i18,__i19,__i20,__i21,__i22,__i23,__i24,__i25,__i26,__i27,__i28,__i29,__i30,/* Inputs: Base address plus 30 offsets */\
__B,__o00,__o01,__o02,__o03,__o04,__o05,__o06,__o07,__o08,__o09,__o10,__o11,__o12,__o13,__o14,__o15,__o16,__o17,__o18,__o19,__o20,__o21,__o22,__o23,__o24,__o25,__o26,__o27,__o28,__o29,__o30 /* Outputs: Base address plus 30 offsets */\
)\
{\
	double TMP01,TMP02;\
	double __tr01,__tr02,__tr03,__tr04,__tr05,__tr06,__tr07,__tr08,__tr09,__tr10,__tr11,__tr12,__tr13,__tr14,__tr15,\
	       __tr16,__tr17,__tr18,__tr19,__tr20,__tr21,__tr22,__tr23,__tr24,__tr25,__tr26,__tr27,__tr28,__tr29,__tr30;\
	double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;\
\
	/* To support in-place DFT, need to save copies of input real parts: */\
	__tr01 = *(__A+__i01);\
	__tr02 = *(__A+__i02);\
	__tr03 = *(__A+__i03);\
	__tr04 = *(__A+__i04);\
	__tr05 = *(__A+__i05);\
	__tr06 = *(__A+__i06);\
	__tr07 = *(__A+__i07);\
	__tr08 = *(__A+__i08);\
	__tr09 = *(__A+__i09);\
	__tr10 = *(__A+__i10);\
	__tr11 = *(__A+__i11);\
	__tr12 = *(__A+__i12);\
	__tr13 = *(__A+__i13);\
	__tr14 = *(__A+__i14);\
	__tr15 = *(__A+__i15);\
	__tr16 = *(__A+__i16);\
	__tr17 = *(__A+__i17);\
	__tr18 = *(__A+__i18);\
	__tr19 = *(__A+__i19);\
	__tr20 = *(__A+__i20);\
	__tr21 = *(__A+__i21);\
	__tr22 = *(__A+__i22);\
	__tr23 = *(__A+__i23);\
	__tr24 = *(__A+__i24);\
	__tr25 = *(__A+__i25);\
	__tr26 = *(__A+__i26);\
	__tr27 = *(__A+__i27);\
	__tr28 = *(__A+__i28);\
	__tr29 = *(__A+__i29);\
	__tr30 = *(__A+__i30);\
\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms: */\
	*(__B+__o00) = *(__A+__i00);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(__A+__i01);\
	x2  = *(__A+__i07);\
	x3  = *(__A+__i18);\
	x4  = *(__A+__i02);\
	x5  = *(__A+__i14);\
	x6  = *(__A+__i05);\
	x7  = *(__A+__i04);\
	x8  = *(__A+__i28);\
	x9  = *(__A+__i10);\
	xa  = *(__A+__i08);\
	xb  = *(__A+__i25);\
	xc  = *(__A+__i20);\
	xd  = *(__A+__i16);\
	xe  = *(__A+__i19);\
	xf  = *(__A+__i09);		/* xr-terms: */\
	x1 += *(__A+__i30);		/* x1 = __Ar01 + __Ar30 */\
	x2 += *(__A+__i24);		/* x2 = __Ar07 + __Ar24 */\
	x3 += *(__A+__i13);		/* x3 = __Ar18 + __Ar13 */\
	x4 += *(__A+__i29);		/* x4 = __Ar02 + __Ar29 */\
	x5 += *(__A+__i17);		/* x5 = __Ar14 + __Ar17 */\
	x6 += *(__A+__i26);		/* x6 = __Ar05 + __Ar26 */\
	x7 += *(__A+__i27);		/* x7 = __Ar04 + __Ar27 */\
	x8 += *(__A+__i03);		/* x8 = __Ar28 + __Ar03 */\
	x9 += *(__A+__i21);		/* x9 = __Ar10 + __Ar21 */\
	xa += *(__A+__i23);		/* xa = __Ar08 + __Ar23 */\
	xb += *(__A+__i06);		/* xb = __Ar25 + __Ar06 */\
	xc += *(__A+__i11);		/* xc = __Ar20 + __Ar11 */\
	xd += *(__A+__i15);		/* xd = __Ar16 + __Ar15 */\
	xe += *(__A+__i12);		/* xe = __Ar19 + __Ar12 */\
	xf += *(__A+__i22);		/* xf = __Ar09 + __Ar22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Br00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(__B+__o00);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(__B+__o00);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(__B+__o00) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(__B+__o08) = xe;	/* s6 = xe - x4 */\
	*(__B+__o14) = x4;	/* sb = xe + x4 - xc */\
	*(__B+__o09) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__B+__o02) = x5;	/* sc = x5 - xd */\
	*(__B+__o12) = xd;	/* s2 = x5 + xd - x7 */\
	*(__B+__o10) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__B+__o15) = x8;	/* s3 = x8 - x3 */\
	*(__B+__o03) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__B+__o13) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__B+__o04) = xb;	/* s9 = xb - x1 */\
	*(__B+__o07) = x1;	/* se = xb + x1 - xf */\
	*(__B+__o11) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(__B+__o01) = x2;	/* sf = x2 - xa */\
	*(__B+__o06) = xa;	/* s5 = x2 + xa - x6 */\
	*(__B+__o05) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(__A+__i01+RE_IM_STRIDE);\
	x2  = *(__A+__i07+RE_IM_STRIDE);\
	x3  = *(__A+__i18+RE_IM_STRIDE);\
	x4  = *(__A+__i02+RE_IM_STRIDE);\
	x5  = *(__A+__i14+RE_IM_STRIDE);\
	x6  = *(__A+__i05+RE_IM_STRIDE);\
	x7  = *(__A+__i04+RE_IM_STRIDE);\
	x8  = *(__A+__i28+RE_IM_STRIDE);\
	x9  = *(__A+__i10+RE_IM_STRIDE);\
	xa  = *(__A+__i08+RE_IM_STRIDE);\
	xb  = *(__A+__i25+RE_IM_STRIDE);\
	xc  = *(__A+__i20+RE_IM_STRIDE);\
	xd  = *(__A+__i16+RE_IM_STRIDE);\
	xe  = *(__A+__i19+RE_IM_STRIDE);\
	xf  = *(__A+__i09+RE_IM_STRIDE);\
	x1 -= *(__A+__i30+RE_IM_STRIDE);\
	x2 -= *(__A+__i24+RE_IM_STRIDE);\
	x3 -= *(__A+__i13+RE_IM_STRIDE);\
	x4 -= *(__A+__i29+RE_IM_STRIDE);\
	x5 -= *(__A+__i17+RE_IM_STRIDE);\
	x6 -= *(__A+__i26+RE_IM_STRIDE);\
	x7 -= *(__A+__i27+RE_IM_STRIDE);\
	x8 -= *(__A+__i03+RE_IM_STRIDE);\
	x9 -= *(__A+__i21+RE_IM_STRIDE);\
	xa -= *(__A+__i23+RE_IM_STRIDE);\
	xb -= *(__A+__i06+RE_IM_STRIDE);\
	xc -= *(__A+__i11+RE_IM_STRIDE);\
	xd -= *(__A+__i15+RE_IM_STRIDE);\
	xe -= *(__A+__i12+RE_IM_STRIDE);\
	xf -= *(__A+__i22+RE_IM_STRIDE);\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(__B+__o01);\
 	x0 -= x2;	\
 	*(__B+__o01) = x0;\
 	x2 *= 2;	\
 	x0 += x2;	\
 	*(__B+__o30) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__B+__o02);		x0 = *(__B+__o03);\
 	x2 -= x5;	 			x0 += x3;	\
 	*(__B+__o02) = x2;		*(__B+__o03) = x0;\
 	x5 *= 2;	 			x3 *= 2;	\
 	x2 += x5;	 			x0 -= x3;	\
 	*(__B+__o29) = x2;		*(__B+__o28) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers to reduce serialization (but stop at 4, more likely no help): */\
 	x0 = *(__B+__o04); 		x2 = *(__B+__o05); 		x3 = *(__B+__o06); 		x5 = *(__B+__o07);\
	x0 -= xb;				x2 -= x6;				x3 += xa;				x5 -= x1;	\
	*(__B+__o04) = x0;		*(__B+__o05) = x2;		*(__B+__o06) = x3;		*(__B+__o07) = x5;\
	xb *= 2;				x6 *= 2;				xa *= 2;				x1 *= 2;	\
	x0 += xb;				x2 += x6;				x3 -= xa;				x5 += x1;	\
	*(__B+__o27) = x0;		*(__B+__o26) = x2;		*(__B+__o25) = x3;		*(__B+__o24) = x5;\
\
 	x0 = *(__B+__o08); 		x2 = *(__B+__o09); 		x3 = *(__B+__o10); 		x5 = *(__B+__o11);\
	x0 -= xe;				x2 -= xc;				x3 -= x7;				x5 += xf;	\
	*(__B+__o08) = x0;		*(__B+__o09) = x2;		*(__B+__o10) = x3;		*(__B+__o11) = x5;\
	xe *= 2;				xc *= 2;				x7 *= 2;				xf *= 2;	\
	x0 += xe;				x2 += xc;				x3 += x7;				x5 -= xf;	\
	*(__B+__o23) = x0;		*(__B+__o22) = x2;		*(__B+__o21) = x3;		*(__B+__o20) = x5;\
\
 	x0 = *(__B+__o12); 		x2 = *(__B+__o13); 		x3 = *(__B+__o14); 		x5 = *(__B+__o15);\
	x0 += xd;				x2 += x9;				x3 -= x4;				x5 += x8;	\
	*(__B+__o12) = x0;		*(__B+__o13) = x2;		*(__B+__o14) = x3;		*(__B+__o15) = x5;\
	xd *= 2;				x9 *= 2;				x4 *= 2;				x8 *= 2;	\
	x0 -= xd;				x2 -= x9;				x3 += x4;				x5 -= x8;	\
	*(__B+__o19) = x0;		*(__B+__o18) = x2;		*(__B+__o17) = x3;		*(__B+__o16) = x5;\
\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Swap +/- in final output sequence: */\
	/***************/\
	/* xi-terms: */\
	*(__B+__o00+RE_IM_STRIDE) = *(__A+__i00+RE_IM_STRIDE);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(__A+__i01+RE_IM_STRIDE);\
	x2  = *(__A+__i07+RE_IM_STRIDE);\
	x3  = *(__A+__i18+RE_IM_STRIDE);\
	x4  = *(__A+__i02+RE_IM_STRIDE);\
	x5  = *(__A+__i14+RE_IM_STRIDE);\
	x6  = *(__A+__i05+RE_IM_STRIDE);\
	x7  = *(__A+__i04+RE_IM_STRIDE);\
	x8  = *(__A+__i28+RE_IM_STRIDE);\
	x9  = *(__A+__i10+RE_IM_STRIDE);\
	xa  = *(__A+__i08+RE_IM_STRIDE);\
	xb  = *(__A+__i25+RE_IM_STRIDE);\
	xc  = *(__A+__i20+RE_IM_STRIDE);\
	xd  = *(__A+__i16+RE_IM_STRIDE);\
	xe  = *(__A+__i19+RE_IM_STRIDE);\
	xf  = *(__A+__i09+RE_IM_STRIDE);		/* xi-terms: */\
	x1 += *(__A+__i30+RE_IM_STRIDE);		/* x1 = __Ai01 + __Ai30 */\
	x2 += *(__A+__i24+RE_IM_STRIDE);		/* x2 = __Ai07 + __Ai24 */\
	x3 += *(__A+__i13+RE_IM_STRIDE);		/* x3 = __Ai18 + __Ai13 */\
	x4 += *(__A+__i29+RE_IM_STRIDE);		/* x4 = __Ai02 + __Ai29 */\
	x5 += *(__A+__i17+RE_IM_STRIDE);		/* x5 = __Ai14 + __Ai17 */\
	x6 += *(__A+__i26+RE_IM_STRIDE);		/* x6 = __Ai05 + __Ai26 */\
	x7 += *(__A+__i27+RE_IM_STRIDE);		/* x7 = __Ai04 + __Ai27 */\
	x8 += *(__A+__i03+RE_IM_STRIDE);		/* x8 = __Ai28 + __Ai03 */\
	x9 += *(__A+__i21+RE_IM_STRIDE);		/* x9 = __Ai10 + __Ai21 */\
	xa += *(__A+__i23+RE_IM_STRIDE);		/* xa = __Ai08 + __Ai23 */\
	xb += *(__A+__i06+RE_IM_STRIDE);		/* xb = __Ai25 + __Ai06 */\
	xc += *(__A+__i11+RE_IM_STRIDE);		/* xc = __Ai20 + __Ai11 */\
	xd += *(__A+__i15+RE_IM_STRIDE);		/* xd = __Ai16 + __Ai15 */\
	xe += *(__A+__i12+RE_IM_STRIDE);		/* xe = __Ai19 + __Ai12 */\
	xf += *(__A+__i22+RE_IM_STRIDE);		/* xf = __Ai09 + __Ai22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Bi00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(__B+__o00+RE_IM_STRIDE);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(__B+__o00+RE_IM_STRIDE);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(__B+__o00+RE_IM_STRIDE) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(__B+__o08+RE_IM_STRIDE) = xe;	/* s6 = xe - x4 */\
	*(__B+__o14+RE_IM_STRIDE) = x4;	/* sb = xe + x4 - xc */\
	*(__B+__o09+RE_IM_STRIDE) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__B+__o02+RE_IM_STRIDE) = x5;	/* sc = x5 - xd */\
	*(__B+__o12+RE_IM_STRIDE) = xd;	/* s2 = x5 + xd - x7 */\
	*(__B+__o10+RE_IM_STRIDE) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__B+__o15+RE_IM_STRIDE) = x8;	/* s3 = x8 - x3 */\
	*(__B+__o03+RE_IM_STRIDE) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__B+__o13+RE_IM_STRIDE) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__B+__o04+RE_IM_STRIDE) = xb;	/* s9 = xb - x1 */\
	*(__B+__o07+RE_IM_STRIDE) = x1;	/* se = xb + x1 - xf */\
	*(__B+__o11+RE_IM_STRIDE) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(__B+__o01+RE_IM_STRIDE) = x2;	/* sf = x2 - xa */\
	*(__B+__o06+RE_IM_STRIDE) = xa;	/* s5 = x2 + xa - x6 */\
	*(__B+__o05+RE_IM_STRIDE) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = __tr01;\
	x2  = __tr07;\
	x3  = __tr18;\
	x4  = __tr02;\
	x5  = __tr14;\
	x6  = __tr05;\
	x7  = __tr04;\
	x8  = __tr28;\
	x9  = __tr10;\
	xa  = __tr08;\
	xb  = __tr25;\
	xc  = __tr20;\
	xd  = __tr16;\
	xe  = __tr19;\
	xf  = __tr09;\
	x1 -= __tr30;\
	x2 -= __tr24;\
	x3 -= __tr13;\
	x4 -= __tr29;\
	x5 -= __tr17;\
	x6 -= __tr26;\
	x7 -= __tr27;\
	x8 -= __tr03;\
	x9 -= __tr21;\
	xa -= __tr23;\
	xb -= __tr06;\
	xc -= __tr11;\
	xd -= __tr15;\
	xe -= __tr12;\
	xf -= __tr22;\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(__B+__o01+RE_IM_STRIDE);\
 	x0 += x2;	\
 	*(__B+__o01+RE_IM_STRIDE) = x0;\
 	x2 *= 2;	\
 	x0 -= x2;	\
 	*(__B+__o30+RE_IM_STRIDE) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__B+__o02+RE_IM_STRIDE);		x0 = *(__B+__o03+RE_IM_STRIDE);\
 	x2 += x5;						 	x0 -= x3;	\
 	*(__B+__o02+RE_IM_STRIDE) = x2;		*(__B+__o03+RE_IM_STRIDE) = x0;\
 	x5 *= 2;	 						x3 *= 2;	\
 	x2 -= x5;	 						x0 += x3;	\
 	*(__B+__o29+RE_IM_STRIDE) = x2;		*(__B+__o28+RE_IM_STRIDE) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers to reduce serialization, but stop at 4, more likely no help: */\
 	x0 = *(__B+__o04+RE_IM_STRIDE);		x2 = *(__B+__o05+RE_IM_STRIDE);		x3 = *(__B+__o06+RE_IM_STRIDE);		x5 = *(__B+__o07+RE_IM_STRIDE);\
	x0 += xb;							x2 += x6;							x3 -= xa;							x5 += x1;	\
	*(__B+__o04+RE_IM_STRIDE) = x0;		*(__B+__o05+RE_IM_STRIDE) = x2;		*(__B+__o06+RE_IM_STRIDE) = x3;		*(__B+__o07+RE_IM_STRIDE) = x5;\
	xb *= 2;							x6 *= 2;							xa *= 2;							x1 *= 2;	\
	x0 -= xb;							x2 -= x6;							x3 += xa;							x5 -= x1;	\
	*(__B+__o27+RE_IM_STRIDE) = x0;		*(__B+__o26+RE_IM_STRIDE) = x2;		*(__B+__o25+RE_IM_STRIDE) = x3;		*(__B+__o24+RE_IM_STRIDE) = x5;\
\
 	x0 = *(__B+__o08+RE_IM_STRIDE);		x2 = *(__B+__o09+RE_IM_STRIDE);		x3 = *(__B+__o10+RE_IM_STRIDE);		x5 = *(__B+__o11+RE_IM_STRIDE);\
	x0 += xe;							x2 += xc;							x3 += x7;							x5 -= xf;	\
	*(__B+__o08+RE_IM_STRIDE) = x0;		*(__B+__o09+RE_IM_STRIDE) = x2;		*(__B+__o10+RE_IM_STRIDE) = x3;		*(__B+__o11+RE_IM_STRIDE) = x5;\
	xe *= 2;							xc *= 2;							x7 *= 2;							xf *= 2;	\
	x0 -= xe;							x2 -= xc;							x3 -= x7;							x5 += xf;	\
	*(__B+__o23+RE_IM_STRIDE) = x0;		*(__B+__o22+RE_IM_STRIDE) = x2;		*(__B+__o21+RE_IM_STRIDE) = x3;		*(__B+__o20+RE_IM_STRIDE) = x5;\
\
 	x0 = *(__B+__o12+RE_IM_STRIDE);		x2 = *(__B+__o13+RE_IM_STRIDE);		x3 = *(__B+__o14+RE_IM_STRIDE);		x5 = *(__B+__o15+RE_IM_STRIDE);\
	x0 -= xd;							x2 -= x9;							x3 += x4;							x5 -= x8;	\
	*(__B+__o12+RE_IM_STRIDE) = x0;		*(__B+__o13+RE_IM_STRIDE) = x2;		*(__B+__o14+RE_IM_STRIDE) = x3;		*(__B+__o15+RE_IM_STRIDE) = x5;\
	xd *= 2;							x9 *= 2;							x4 *= 2;							x8 *= 2;	\
	x0 += xd;							x2 += x9;							x3 -= x4;							x5 += x8;	\
	*(__B+__o19+RE_IM_STRIDE) = x0;		*(__B+__o18+RE_IM_STRIDE) = x2;		*(__B+__o17+RE_IM_STRIDE) = x3;		*(__B+__o16+RE_IM_STRIDE) = x5;\
\
/* Totals: 658 FADD, 234 FMUL. */\
}

/****************************************************************************************************/
/***** Out-of-place versions of above (here we need separate specializations for DIF and DIT). ******/
/****************************************************************************************************/

/* Out-of-place DIF assumes outputs contiguous starting at output-base-address __out: */
#define RADIX_31_DIF(\
__A,__idx, /* Inputs: Base address plus 30 (index) offsets */\
__out /* Outputs: Base address plus 30 offsets */\
)\
{\
	double TMP01,TMP02;\
	double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;\
	double *Aim = __A + RE_IM_STRIDE;\
\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms: */\
*(__out)= *(__A+__idx[0x00]);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(__A+__idx[0x01]);\
	x2  = *(__A+__idx[0x07]);\
	x3  = *(__A+__idx[0x12]);\
	x4  = *(__A+__idx[0x02]);\
	x5  = *(__A+__idx[0x0E]);\
	x6  = *(__A+__idx[0x05]);\
	x7  = *(__A+__idx[0x04]);\
	x8  = *(__A+__idx[0x1C]);\
	x9  = *(__A+__idx[0x0A]);\
	xa  = *(__A+__idx[0x08]);\
	xb  = *(__A+__idx[0x19]);\
	xc  = *(__A+__idx[0x14]);\
	xd  = *(__A+__idx[0x10]);\
	xe  = *(__A+__idx[0x13]);\
	xf  = *(__A+__idx[0x09]);		/* xr-terms: */\
	x1 += *(__A+__idx[0x1E]);		/* x1 = __Ar01 + __Ar30 */\
	x2 += *(__A+__idx[0x18]);		/* x2 = __Ar07 + __Ar24 */\
	x3 += *(__A+__idx[0x0D]);		/* x3 = __Ar18 + __Ar13 */\
	x4 += *(__A+__idx[0x1D]);		/* x4 = __Ar02 + __Ar29 */\
	x5 += *(__A+__idx[0x11]);		/* x5 = __Ar14 + __Ar17 */\
	x6 += *(__A+__idx[0x1A]);		/* x6 = __Ar05 + __Ar26 */\
	x7 += *(__A+__idx[0x1B]);		/* x7 = __Ar04 + __Ar27 */\
	x8 += *(__A+__idx[0x03]);		/* x8 = __Ar28 + __Ar03 */\
	x9 += *(__A+__idx[0x15]);		/* x9 = __Ar10 + __Ar21 */\
	xa += *(__A+__idx[0x17]);		/* xa = __Ar08 + __Ar23 */\
	xb += *(__A+__idx[0x06]);		/* xb = __Ar25 + __Ar06 */\
	xc += *(__A+__idx[0x0B]);		/* xc = __Ar20 + __Ar11 */\
	xd += *(__A+__idx[0x0F]);		/* xd = __Ar16 + __Ar15 */\
	xe += *(__A+__idx[0x0C]);		/* xe = __Ar19 + __Ar12 */\
	xf += *(__A+__idx[0x16]);		/* xf = __Ar09 + __Ar22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Br00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(__out);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(__out);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(__out) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(__out+16) = xe;	/* s6 = xe - x4 */\
	*(__out+28) = x4;	/* sb = xe + x4 - xc */\
	*(__out+18) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__out+ 4) = x5;	/* sc = x5 - xd */\
	*(__out+24) = xd;	/* s2 = x5 + xd - x7 */\
	*(__out+20) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__out+30) = x8;	/* s3 = x8 - x3 */\
	*(__out+ 6) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__out+26) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__out+ 8) = xb;	/* s9 = xb - x1 */\
	*(__out+14) = x1;	/* se = xb + x1 - xf */\
	*(__out+22) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(__out+ 2) = x2;	/* sf = x2 - xa */\
	*(__out+12) = xa;	/* s5 = x2 + xa - x6 */\
	*(__out+10) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(Aim+__idx[0x01]);\
	x2  = *(Aim+__idx[0x07]);\
	x3  = *(Aim+__idx[0x12]);\
	x4  = *(Aim+__idx[0x02]);\
	x5  = *(Aim+__idx[0x0E]);\
	x6  = *(Aim+__idx[0x05]);\
	x7  = *(Aim+__idx[0x04]);\
	x8  = *(Aim+__idx[0x1C]);\
	x9  = *(Aim+__idx[0x0A]);\
	xa  = *(Aim+__idx[0x08]);\
	xb  = *(Aim+__idx[0x19]);\
	xc  = *(Aim+__idx[0x14]);\
	xd  = *(Aim+__idx[0x10]);\
	xe  = *(Aim+__idx[0x13]);\
	xf  = *(Aim+__idx[0x09]);\
	x1 -= *(Aim+__idx[0x1E]);\
	x2 -= *(Aim+__idx[0x18]);\
	x3 -= *(Aim+__idx[0x0D]);\
	x4 -= *(Aim+__idx[0x1D]);\
	x5 -= *(Aim+__idx[0x11]);\
	x6 -= *(Aim+__idx[0x1A]);\
	x7 -= *(Aim+__idx[0x1B]);\
	x8 -= *(Aim+__idx[0x03]);\
	x9 -= *(Aim+__idx[0x15]);\
	xa -= *(Aim+__idx[0x17]);\
	xb -= *(Aim+__idx[0x06]);\
	xc -= *(Aim+__idx[0x0B]);\
	xd -= *(Aim+__idx[0x0F]);\
	xe -= *(Aim+__idx[0x0C]);\
	xf -= *(Aim+__idx[0x16]);\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(__out+ 2);\
 	x0 -= x2;	\
 	*(__out+ 2) = x0;\
 	x2 *= 2;	\
 	x0 += x2;	\
 	*(__out+60) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__out+ 4);		x0 = *(__out+ 6);\
 	x2 -= x5;	 			x0 += x3;	\
 	*(__out+ 4) = x2;		*(__out+ 6) = x0;\
 	x5 *= 2;	 			x3 *= 2;	\
 	x2 += x5;	 			x0 -= x3;	\
 	*(__out+58) = x2;		*(__out+56) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers to reduce serialization (but stop at 4, more likely no help): */\
 	x0 = *(__out+ 8);		x2 = *(__out+10);		x3 = *(__out+12);		x5 = *(__out+14);\
	x0 -= xb;				x2 -= x6;				x3 += xa;				x5 -= x1;	\
	*(__out+ 8) = x0;		*(__out+10) = x2;		*(__out+12) = x3;		*(__out+14) = x5;\
	xb *= 2;				x6 *= 2;				xa *= 2;				x1 *= 2;	\
	x0 += xb;				x2 += x6;				x3 -= xa;				x5 += x1;	\
	*(__out+54) = x0;		*(__out+52) = x2;		*(__out+50) = x3;		*(__out+48) = x5;\
\
 	x0 = *(__out+16);		x2 = *(__out+18);		x3 = *(__out+20);		x5 = *(__out+22);\
	x0 -= xe;				x2 -= xc;				x3 -= x7;				x5 += xf;	\
	*(__out+16) = x0;		*(__out+18) = x2;		*(__out+20) = x3;		*(__out+22) = x5;\
	xe *= 2;				xc *= 2;				x7 *= 2;				xf *= 2;	\
	x0 += xe;				x2 += xc;				x3 += x7;				x5 -= xf;	\
	*(__out+46) = x0;		*(__out+44) = x2;		*(__out+42) = x3;		*(__out+40) = x5;\
\
 	x0 = *(__out+24);		x2 = *(__out+26);		x3 = *(__out+28);		x5 = *(__out+30);\
	x0 += xd;				x2 += x9;				x3 -= x4;				x5 += x8;	\
	*(__out+24) = x0;		*(__out+26) = x2;		*(__out+28) = x3;		*(__out+30) = x5;\
	xd *= 2;				x9 *= 2;				x4 *= 2;				x8 *= 2;	\
	x0 -= xd;				x2 -= x9;				x3 += x4;				x5 -= x8;	\
	*(__out+38) = x0;		*(__out+36) = x2;		*(__out+34) = x3;		*(__out+32) = x5;\
\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Swap +/- in final output sequence: */\
	/***************/\
	/* xi-terms: */\
*(__out+1) = *(Aim+__idx[0x00]);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(Aim+__idx[0x01]);\
	x2  = *(Aim+__idx[0x07]);\
	x3  = *(Aim+__idx[0x12]);\
	x4  = *(Aim+__idx[0x02]);\
	x5  = *(Aim+__idx[0x0E]);\
	x6  = *(Aim+__idx[0x05]);\
	x7  = *(Aim+__idx[0x04]);\
	x8  = *(Aim+__idx[0x1C]);\
	x9  = *(Aim+__idx[0x0A]);\
	xa  = *(Aim+__idx[0x08]);\
	xb  = *(Aim+__idx[0x19]);\
	xc  = *(Aim+__idx[0x14]);\
	xd  = *(Aim+__idx[0x10]);\
	xe  = *(Aim+__idx[0x13]);\
	xf  = *(Aim+__idx[0x09]);		/* xi-terms: */\
	x1 += *(Aim+__idx[0x1E]);		/* x1 = __Ai01 + __Ai30 */\
	x2 += *(Aim+__idx[0x18]);		/* x2 = __Ai07 + __Ai24 */\
	x3 += *(Aim+__idx[0x0D]);		/* x3 = __Ai18 + __Ai13 */\
	x4 += *(Aim+__idx[0x1D]);		/* x4 = __Ai02 + __Ai29 */\
	x5 += *(Aim+__idx[0x11]);		/* x5 = __Ai14 + __Ai17 */\
	x6 += *(Aim+__idx[0x1A]);		/* x6 = __Ai05 + __Ai26 */\
	x7 += *(Aim+__idx[0x1B]);		/* x7 = __Ai04 + __Ai27 */\
	x8 += *(Aim+__idx[0x03]);		/* x8 = __Ai28 + __Ai03 */\
	x9 += *(Aim+__idx[0x15]);		/* x9 = __Ai10 + __Ai21 */\
	xa += *(Aim+__idx[0x17]);		/* xa = __Ai08 + __Ai23 */\
	xb += *(Aim+__idx[0x06]);		/* xb = __Ai25 + __Ai06 */\
	xc += *(Aim+__idx[0x0B]);		/* xc = __Ai20 + __Ai11 */\
	xd += *(Aim+__idx[0x0F]);		/* xd = __Ai16 + __Ai15 */\
	xe += *(Aim+__idx[0x0C]);		/* xe = __Ai19 + __Ai12 */\
	xf += *(Aim+__idx[0x16]);		/* xf = __Ai09 + __Ai22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Bi00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(__out+1);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(__out+1);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(__out+1) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(__out+17) = xe;	/* s6 = xe - x4 */\
	*(__out+29) = x4;	/* sb = xe + x4 - xc */\
	*(__out+19) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__out+ 5) = x5;	/* sc = x5 - xd */\
	*(__out+25) = xd;	/* s2 = x5 + xd - x7 */\
	*(__out+21) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__out+31) = x8;	/* s3 = x8 - x3 */\
	*(__out+ 7) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__out+27) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__out+ 9) = xb;	/* s9 = xb - x1 */\
	*(__out+15) = x1;	/* se = xb + x1 - xf */\
	*(__out+23) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(__out+ 3) = x2;	/* sf = x2 - xa */\
	*(__out+13) = xa;	/* s5 = x2 + xa - x6 */\
	*(__out+11) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(__A+__idx[0x01]);\
	x2  = *(__A+__idx[0x07]);\
	x3  = *(__A+__idx[0x12]);\
	x4  = *(__A+__idx[0x02]);\
	x5  = *(__A+__idx[0x0E]);\
	x6  = *(__A+__idx[0x05]);\
	x7  = *(__A+__idx[0x04]);\
	x8  = *(__A+__idx[0x1C]);\
	x9  = *(__A+__idx[0x0A]);\
	xa  = *(__A+__idx[0x08]);\
	xb  = *(__A+__idx[0x19]);\
	xc  = *(__A+__idx[0x14]);\
	xd  = *(__A+__idx[0x10]);\
	xe  = *(__A+__idx[0x13]);\
	xf  = *(__A+__idx[0x09]);\
	x1 -= *(__A+__idx[0x1E]);\
	x2 -= *(__A+__idx[0x18]);\
	x3 -= *(__A+__idx[0x0D]);\
	x4 -= *(__A+__idx[0x1D]);\
	x5 -= *(__A+__idx[0x11]);\
	x6 -= *(__A+__idx[0x1A]);\
	x7 -= *(__A+__idx[0x1B]);\
	x8 -= *(__A+__idx[0x03]);\
	x9 -= *(__A+__idx[0x15]);\
	xa -= *(__A+__idx[0x17]);\
	xb -= *(__A+__idx[0x06]);\
	xc -= *(__A+__idx[0x0B]);\
	xd -= *(__A+__idx[0x0F]);\
	xe -= *(__A+__idx[0x0C]);\
	xf -= *(__A+__idx[0x16]);\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(__out+ 3);\
 	x0 += x2;	\
 	*(__out+ 3) = x0;\
 	x2 *= 2;	\
 	x0 -= x2;	\
 	*(__out+61) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__out+ 5);	x0 = *(__out+ 7);\
 	x2 += x5;			x0 -= x3;	\
 	*(__out+ 5) = x2;	*(__out+ 7) = x0;\
 	x5 *= 2;	 		x3 *= 2;	\
 	x2 -= x5;	 		x0 += x3;	\
 	*(__out+59) = x2;	*(__out+57) = x0;	/* x0,2,3,5 */\
/* Now take advantage of4 free registers to reduce serialization, but stop at 4, more likely no help: */\
 	x0 = *(__out+ 9);	x2 = *(__out+11);	x3 = *(__out+13);	x5 = *(__out+15);\
	x0 += xb;			x2 += x6;			x3 -= xa;			x5 += x1;	\
	*(__out+ 9) = x0;	*(__out+11) = x2;	*(__out+13) = x3;	*(__out+15) = x5;\
	xb *= 2;			x6 *= 2;			xa *= 2;			x1 *= 2;	\
	x0 -= xb;			x2 -= x6;			x3 += xa;			x5 -= x1;	\
	*(__out+55) = x0;	*(__out+53) = x2;	*(__out+51) = x3;	*(__out+49) = x5;\
\
 	x0 = *(__out+17);	x2 = *(__out+19);	x3 = *(__out+21);	x5 = *(__out+23);\
	x0 += xe;			x2 += xc;			x3 += x7;			x5 -= xf;	\
	*(__out+17) = x0;	*(__out+19) = x2;	*(__out+21) = x3;	*(__out+23) = x5;\
	xe *= 2;			xc *= 2;			x7 *= 2;			xf *= 2;	\
	x0 -= xe;			x2 -= xc;			x3 -= x7;			x5 += xf;	\
	*(__out+47) = x0;	*(__out+45) = x2;	*(__out+43) = x3;	*(__out+41) = x5;\
\
 	x0 = *(__out+25);	x2 = *(__out+27);	x3 = *(__out+29);	x5 = *(__out+31);\
	x0 -= xd;			x2 -= x9;			x3 += x4;			x5 -= x8;	\
	*(__out+25) = x0;	*(__out+27) = x2;	*(__out+29) = x3;	*(__out+31) = x5;\
	xd *= 2;			x9 *= 2;			x4 *= 2;			x8 *= 2;	\
	x0 += xd;			x2 += x9;			x3 -= x4;			x5 += x8;	\
	*(__out+39) = x0;	*(__out+37) = x2;	*(__out+35) = x3;	*(__out+33) = x5;\
\
/* Totals: 658 FADD, 234 FMUL. */\
}

/* Out-of-place DIT assumes inputs contiguous starting at input-base-address __in: */
#define RADIX_31_DIT(\
	__in, /* Inputs: Base address plus 30 offsets */\
	__B,__odx	/* Outputs: Base address plus 30 (index) offsets */\
)\
{\
	double TMP01,TMP02;\
	double x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc, xd, xe, xf;\
	double *Bim = __B + RE_IM_STRIDE;\
\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms: */\
	*(__B+__odx[0x0]) = *(__in);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(__in+ 2);\
	x2  = *(__in+14);\
	x3  = *(__in+36);\
	x4  = *(__in+ 4);\
	x5  = *(__in+28);\
	x6  = *(__in+10);\
	x7  = *(__in+ 8);\
	x8  = *(__in+56);\
	x9  = *(__in+20);\
	xa  = *(__in+16);\
	xb  = *(__in+50);\
	xc  = *(__in+40);\
	xd  = *(__in+32);\
	xe  = *(__in+38);\
	xf  = *(__in+18);		/* xr-terms: */\
	x1 += *(__in+60);		/* x1 = __Ar01 + __Ar30 */\
	x2 += *(__in+48);		/* x2 = __Ar07 + __Ar24 */\
	x3 += *(__in+26);		/* x3 = __Ar18 + __Ar13 */\
	x4 += *(__in+58);		/* x4 = __Ar02 + __Ar29 */\
	x5 += *(__in+34);		/* x5 = __Ar14 + __Ar17 */\
	x6 += *(__in+52);		/* x6 = __Ar05 + __Ar26 */\
	x7 += *(__in+54);		/* x7 = __Ar04 + __Ar27 */\
	x8 += *(__in+ 6);		/* x8 = __Ar28 + __Ar03 */\
	x9 += *(__in+42);		/* x9 = __Ar10 + __Ar21 */\
	xa += *(__in+46);		/* xa = __Ar08 + __Ar23 */\
	xb += *(__in+12);		/* xb = __Ar25 + __Ar06 */\
	xc += *(__in+22);		/* xc = __Ar20 + __Ar11 */\
	xd += *(__in+30);		/* xd = __Ar16 + __Ar15 */\
	xe += *(__in+24);		/* xe = __Ar19 + __Ar12 */\
	xf += *(__in+44);		/* xf = __Ar09 + __Ar22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Br00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(__B+__odx[0x0]);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(__B+__odx[0x0]);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(__B+__odx[0x0]) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(__B+__odx[0x08]) = xe;	/* s6 = xe - x4 */\
	*(__B+__odx[0x0E]) = x4;	/* sb = xe + x4 - xc */\
	*(__B+__odx[0x09]) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__B+__odx[0x02]) = x5;	/* sc = x5 - xd */\
	*(__B+__odx[0x0C]) = xd;	/* s2 = x5 + xd - x7 */\
	*(__B+__odx[0x0A]) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__B+__odx[0x0F]) = x8;	/* s3 = x8 - x3 */\
	*(__B+__odx[0x03]) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__B+__odx[0x0D]) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__B+__odx[0x04]) = xb;	/* s9 = xb - x1 */\
	*(__B+__odx[0x07]) = x1;	/* se = xb + x1 - xf */\
	*(__B+__odx[0x0B]) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(__B+__odx[0x01]) = x2;	/* sf = x2 - xa */\
	*(__B+__odx[0x06]) = xa;	/* s5 = x2 + xa - x6 */\
	*(__B+__odx[0x05]) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(__in+ 3);\
	x2  = *(__in+15);\
	x3  = *(__in+37);\
	x4  = *(__in+ 5);\
	x5  = *(__in+29);\
	x6  = *(__in+11);\
	x7  = *(__in+ 9);\
	x8  = *(__in+57);\
	x9  = *(__in+21);\
	xa  = *(__in+17);\
	xb  = *(__in+51);\
	xc  = *(__in+41);\
	xd  = *(__in+33);\
	xe  = *(__in+39);\
	xf  = *(__in+19);\
	x1 -= *(__in+61);\
	x2 -= *(__in+49);\
	x3 -= *(__in+27);\
	x4 -= *(__in+59);\
	x5 -= *(__in+35);\
	x6 -= *(__in+53);\
	x7 -= *(__in+55);\
	x8 -= *(__in+ 7);\
	x9 -= *(__in+43);\
	xa -= *(__in+47);\
	xb -= *(__in+13);\
	xc -= *(__in+23);\
	xd -= *(__in+31);\
	xe -= *(__in+25);\
	xf -= *(__in+45);\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(__B+__odx[0x01]);\
 	x0 -= x2;	\
 	*(__B+__odx[0x01]) = x0;\
 	x2 *= 2;	\
 	x0 += x2;	\
 	*(__B+__odx[0x1E]) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__B+__odx[0x02]);		x0 = *(__B+__odx[0x03]);\
 	x2 -= x5;	 					x0 += x3;	\
 	*(__B+__odx[0x02]) = x2;		*(__B+__odx[0x03]) = x0;\
 	x5 *= 2;	 					x3 *= 2;	\
 	x2 += x5;	 					x0 -= x3;	\
 	*(__B+__odx[0x1D]) = x2;		*(__B+__odx[0x1C]) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers to reduce serialization (but stop at 4, more likely no help): */\
 	x0 = *(__B+__odx[0x04]);		x2 = *(__B+__odx[0x05]);	x3 = *(__B+__odx[0x06]);	x5 = *(__B+__odx[0x07]);\
	x0 -= xb;						x2 -= x6;					x3 += xa;					x5 -= x1;	\
	*(__B+__odx[0x04]) = x0;		*(__B+__odx[0x05]) = x2;	*(__B+__odx[0x06]) = x3;	*(__B+__odx[0x07]) = x5;\
	xb *= 2;						x6 *= 2;					xa *= 2;					x1 *= 2;	\
	x0 += xb;						x2 += x6;					x3 -= xa;					x5 += x1;	\
	*(__B+__odx[0x1B]) = x0;		*(__B+__odx[0x1A]) = x2;	*(__B+__odx[0x19]) = x3;	*(__B+__odx[0x18]) = x5;\
\
 	x0 = *(__B+__odx[0x08]);		x2 = *(__B+__odx[0x09]);	x3 = *(__B+__odx[0x0A]);	x5 = *(__B+__odx[0x0B]);\
	x0 -= xe;						x2 -= xc;					x3 -= x7;					x5 += xf;	\
	*(__B+__odx[0x08]) = x0;		*(__B+__odx[0x09]) = x2;	*(__B+__odx[0x0A]) = x3;	*(__B+__odx[0x0B]) = x5;\
	xe *= 2;						xc *= 2;					x7 *= 2;					xf *= 2;	\
	x0 += xe;						x2 += xc;					x3 += x7;					x5 -= xf;	\
	*(__B+__odx[0x17]) = x0;		*(__B+__odx[0x16]) = x2;	*(__B+__odx[0x15]) = x3;	*(__B+__odx[0x14]) = x5;\
\
 	x0 = *(__B+__odx[0x0C]);		x2 = *(__B+__odx[0x0D]);	x3 = *(__B+__odx[0x0E]);	x5 = *(__B+__odx[0x0F]);\
	x0 += xd;						x2 += x9;					x3 -= x4;					x5 += x8;	\
	*(__B+__odx[0x0C]) = x0;		*(__B+__odx[0x0D]) = x2;	*(__B+__odx[0x0E]) = x3;	*(__B+__odx[0x0F]) = x5;\
	xd *= 2;						x9 *= 2;					x4 *= 2;					x8 *= 2;	\
	x0 -= xd;						x2 -= x9;					x3 += x4;					x5 -= x8;	\
	*(__B+__odx[0x13]) = x0;		*(__B+__odx[0x12]) = x2;	*(__B+__odx[0x11]) = x3;	*(__B+__odx[0x10]) = x5;\
\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Swap +/- in final output sequence: */\
	/***************/\
	/* xi-terms: */\
	*(Bim+__odx[0x0]) = *(__in+1);	/* x0; Free up by using output 0 in its place here */\
	x1  = *(__in+ 3);\
	x2  = *(__in+15);\
	x3  = *(__in+37);\
	x4  = *(__in+ 5);\
	x5  = *(__in+29);\
	x6  = *(__in+11);\
	x7  = *(__in+ 9);\
	x8  = *(__in+57);\
	x9  = *(__in+21);\
	xa  = *(__in+17);\
	xb  = *(__in+51);\
	xc  = *(__in+41);\
	xd  = *(__in+33);\
	xe  = *(__in+39);\
	xf  = *(__in+19);		/* xi-terms: */\
	x1 += *(__in+61);		/* x1 = __Ai01 + __Ai30 */\
	x2 += *(__in+49);		/* x2 = __Ai07 + __Ai24 */\
	x3 += *(__in+27);		/* x3 = __Ai18 + __Ai13 */\
	x4 += *(__in+59);		/* x4 = __Ai02 + __Ai29 */\
	x5 += *(__in+35);		/* x5 = __Ai14 + __Ai17 */\
	x6 += *(__in+53);		/* x6 = __Ai05 + __Ai26 */\
	x7 += *(__in+55);		/* x7 = __Ai04 + __Ai27 */\
	x8 += *(__in+ 7);		/* x8 = __Ai28 + __Ai03 */\
	x9 += *(__in+43);		/* x9 = __Ai10 + __Ai21 */\
	xa += *(__in+47);		/* xa = __Ai08 + __Ai23 */\
	xb += *(__in+13);		/* xb = __Ai25 + __Ai06 */\
	xc += *(__in+23);		/* xc = __Ai20 + __Ai11 */\
	xd += *(__in+31);		/* xd = __Ai16 + __Ai15 */\
	xe += *(__in+25);		/* xe = __Ai19 + __Ai12 */\
	xf += *(__in+45);		/* xf = __Ai09 + __Ai22 */\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;			/* [-]	xb_6 = xb - x6 */\
	xb += x1;			/* x16b = x1 + x6 + xb */\
	x1 -= x0;			/* x1_b = x1 - xb */\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;			/* [-]	xe_9 = xe - x9 */\
	xe += x4;			/* x49e = x4 + x9 + xe */\
	x4 -= x0;			/* x4_e = x4 - xe */\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;			/* x7_2 = x7 - x2 */\
	x2 += xc;			/* x27c = x2 + x7 + xc */\
	xc -= x0;			/* [-]	x2_c = x2 - xc */\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;			/* xa_5 = xa - x5 */\
	x5 += xf;			/* x5af = x5 + xa + xf */\
	xf -= x0;			/* [-]	x5_f = x5 - xf */\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;			/* [-]	x8_3 = x8 - x3 */\
	x8 += xd;			/* x38d = x3 + x8 + xd */\
	xd -= x0;			/* xd_8 = xd - x8 */\
\
	xb -= x8;			/* x16b = x16b - x38d */\
	xe -= x8;			/* x49e = x49e - x38d */\
	x2 -= x8;			/* x27c = x27c - x38d */\
	x5 -= x8;			/* x5af = x5af - x38d */\
	xb += xe;			/* x16b49e = x16b + x49e */\
	x2 += x5;			/* x27c5af = x27c + x5af */\
	x0 = xb;			/* x16b49e */\
	x0 += x2;			/* xsigma = x16b49e + x27c5af */\
	x8 *= 5;			/* 5*x38d */\
	/* Bi00 needed for final output; use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= dc0_5th;		/* dc0*x38d */\
	x0 *= dc1;			/* dc1*xsigma */\
	x8 += *(Bim+__odx[0x0]);		/* x0 + dc0*x38d */\
	x8 += x0;			/* u1234 = x0 + dc0*x38d + dc1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= dc2;			/* dc2*x27c5af */\
	x8 *= dc7;			/* dc7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + dc2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + dc7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(dc3,B_c3483,dc8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= dc3;\
	x5 *= dc8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(dc5,B_c5645,-dc4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= dc5;\
	x2 *= ndc4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= dc0_5th;		/* dc0*usigma/5 */\
	x0 += *(Bim+__odx[0x0]);		/* dc0*usigma/5 + x0 */\
	x0 *= 5;			/* dc0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = dc0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
/* DC component: */\
	x0 = TMP01;\
	*(Bim+__odx[0x0]) += x0;	/* x0 */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(dc9,B_c9019,dc11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= dc9;\
	x3 *= dc11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(dc15,B_c5285,dc18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= dc15;\
	x0 *= dc18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(dc16,B_c6396,dc19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= dc16;\
	x9 *= dc19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(dc17,B_c7407,dc20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= dc17;\
	x1 *= dc20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(dc24,C_c4174,dc27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= dc24;\
	x0 *= dc27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(dc25,C_c5285,dc28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= dc25;\
	x6 *= dc28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(dc26,C_c6396,dc29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= dc26;\
	x4 *= dc29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* xe + xc */\
	x4 += x0;	/* xe - xc + x4 */\
\
	*(Bim+__odx[0x08]) = xe;	/* s6 = xe - x4 */\
	*(Bim+__odx[0x0E]) = x4;	/* sb = xe + x4 - xc */\
	*(Bim+__odx[0x09]) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(Bim+__odx[0x02]) = x5;	/* sc = x5 - xd */\
	*(Bim+__odx[0x0C]) = xd;	/* s2 = x5 + xd - x7 */\
	*(Bim+__odx[0x0A]) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(Bim+__odx[0x0F]) = x8;	/* s3 = x8 - x3 */\
	*(Bim+__odx[0x03]) = x3;	/* s8 = x8 + x3 - x9 */\
	*(Bim+__odx[0x0D]) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(Bim+__odx[0x04]) = xb;	/* s9 = xb - x1 */\
	*(Bim+__odx[0x07]) = x1;	/* se = xb + x1 - xf */\
	*(Bim+__odx[0x0B]) = xf;	/* s4 = xb + xf */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* x2 + x6 */\
	xa += x0;	/* x2 - x6 + xa */\
\
	*(Bim+__odx[0x01]) = x2;	/* sf = x2 - xa */\
	*(Bim+__odx[0x06]) = xa;	/* s5 = x2 + xa - x6 */\
	*(Bim+__odx[0x05]) = x6;	/* sa = x2 + x6 */\
	/****************************************************************************************************************/\
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all dc-consts --> ds, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(__in+ 2);\
	x2  = *(__in+14);\
	x3  = *(__in+36);\
	x4  = *(__in+ 4);\
	x5  = *(__in+28);\
	x6  = *(__in+10);\
	x7  = *(__in+ 8);\
	x8  = *(__in+56);\
	x9  = *(__in+20);\
	xa  = *(__in+16);\
	xb  = *(__in+50);\
	xc  = *(__in+40);\
	xd  = *(__in+32);\
	xe  = *(__in+38);\
	xf  = *(__in+18);\
	x1 -= *(__in+60);\
	x2 -= *(__in+48);\
	x3 -= *(__in+26);\
	x4 -= *(__in+58);\
	x5 -= *(__in+34);\
	x6 -= *(__in+52);\
	x7 -= *(__in+54);\
	x8 -= *(__in+ 6);\
	x9 -= *(__in+42);\
	xa -= *(__in+46);\
	xb -= *(__in+12);\
	xc -= *(__in+22);\
	xd -= *(__in+30);\
	xe -= *(__in+24);\
	xf -= *(__in+44);\
\
	x0 = xb;\
	xb += x6;\
	x6 -= x0;\
	xb += x1;\
	x1 -= x0;\
	x0 = xe;\
	xe += x9;\
	x9 -= x0;\
	xe += x4;\
	x4 -= x0;\
	x0 = x2;\
	x2 += x7;\
	x7 -= x0;\
	x2 += xc;\
	xc -= x0;\
	x0 = x5;\
	x5 += xa;\
	xa -= x0;\
	x5 += xf;\
	xf -= x0;\
	x0 = x8;\
	x8 += x3;\
	x3 -= x0;\
	x8 += xd;\
	xd -= x0;		/* x0 */\
\
	xb -= x8;\
	xe -= x8;\
	x2 -= x8;\
	x5 -= x8;\
	xb += xe;\
	x2 += x5;\
	x0 = xb;\
	x0 += x2;\
	x8 *= 5;\
	/* Use TMP0* for spills */\
	TMP01 = x0 + x8;	/* usigma = xsigma + 5*x38d */\
	x8 *= ds0_5th;		/* ds0*x38d */\
	x0 *= ds1;			/* ds1*xsigma */\
	x8 += x0;			/* u1234 = ds0*x38d + ds1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= ds2;			/* ds2*x27c5af */\
	x8 *= ds7;			/* ds7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + ds2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + ds7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(ds3,B_s3483,ds8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= ds3;\
	x5 *= ds8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(ds5,B_s5645,-ds4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= ds5;\
	x2 *= nds4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= ds0;			/* ds0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = ds0*usigma - (u1 + u2 + u3 + u4) */\
/* No DC component in this half: */\
	x0 = x1;\
	x1 += xd;\
	xd -= x0;\
	x0 = x4;\
	x4 += xa;\
	xa -= x0;	/* x0 */\
	x4 += x1;\
	x7 += x4;	/* ~x7 = x7 + x4 */\
	x0 = x4;\
	x0 *= c51d50m1;\
	x0 += x7;\
	x0 *= c50;	/* c50*x7 + c51*x4 */\
	x4 -= x1;\
	x1 += x0;\
	x4 += x0;	/* x0 */\
	x0 = xd;\
	xd *= s52;\
	xd += xa;	/* xa + s52*xd */\
	xa *= s52;\
	xa -= x0;	/* s52*xa - xd *//* x0 */\
	x0 = x6;\
	x6 += xf;\
	xf -= x0;\
	x0 = x9;\
	x9 += xc;\
	xc -= x0;\
	x9 += x6;\
	x3 += x9;	/* ~x3 = x3 + x9 */\
	x0 = x9;\
	x0 *= c51d50m1;\
	x0 += x3;\
	x0 *= c50;	/* c50*x3 + c51*x9 */\
	x9 -= x6;\
	x6 += x0;\
	x9 += x0;	/* x0 */\
	x0 = xc;\
	xc *= s52;\
	xc += xf;	/* xf + s52*xc */\
	xf *= s52;\
	xf -= x0;	/* s52*xf - xc *//* x0 */\
	/* MAT_2x2_MUL(ds9,B_s9019,ds11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= ds9;\
	x3 *= ds11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(ds15,B_s5285,ds18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= ds15;\
	x0 *= ds18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(ds16,B_s6396,ds19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= ds16;\
	x9 *= ds19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(ds17,B_s7407,ds20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= ds17;\
	x1 *= ds20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the dc-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(ds24,C_s4174,ds27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= ds24;\
	x0 *= ds27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(ds25,C_s5285,ds28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= ds25;\
	x6 *= ds28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(ds26,C_s6396,ds29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= ds26;\
	x4 *= ds29;\
\
	xa += x8;	/* x8 */\
	x4 += x0;\
	x0 = xc;\
	x0 += xf;\
	x0 *= c50;\
	x7 += x0;	/* w6 */\
	x0 *= c51d50m1;\
	x0 += x7;\
	xf += x0;\
	xc += x0;\
	x0 = x6;\
	x6 *= s52;\
	x6 += x9;	/* x9 + s52*x6 */\
	x9 *= s52;\
	x9 -= x0;	/* s52*x9 - x6 *//* x0 */\
	x0 = x6;\
	x6 += xf;	/* w8 */\
	xf -= x0;	/* we */\
	x0 = x9;\
	x9 += xc;	/* wa */\
	xc -= x0;	/* wc */\
	x0 = xd;\
	x0 += xa;\
	x0 *= c50;\
	x3 += x0;	/* w7 */\
	x0 *= c51d50m1;\
	x0 += x3;\
	xa += x0;\
	xd += x0;\
	x0 = x4;\
	x4 *= s52;\
	x4 += x1;	/* x1 + s52*x4 */\
	x1 *= s52;\
	x1 -= x0;	/* s52*x1 - x4 *//* x0 */\
	x0 = x4;\
	x4 += xa;	/* w9 */\
	xa -= x0;	/* wf */\
	x0 = x1;\
	x1 += xd;	/* wb */\
	xd -= x0;	/* wd */\
	x8 = TMP02;\
/* Now form LHS-terms of the x +- y combos which will ultimately be output. */\
/* abc -> ec4: */\
	x0  = xe;\
	xe -= x4;	/* s6 = xe - x4 */\
	x4 -= xc;\
	xc += x0;	/* s1 = xe + xc */\
	x4 += x0;	/* sb = xe - xc + x4 */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* sc = x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* s7 = x5 + x7 */\
	xd += x0;	/* s2 = x5 - x7 + xd */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* s3 = x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* sd = x8 + x9 */\
	x3 += x0;	/* s8 = x8 - x9 + x3 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* s9 = xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* s4 = xb + xf */\
	x1 += x0;	/* se = xb - xf + x1 */\
/* abc -> 26a: */\
	x0  = x2;\
	x2 -= xa;	/* sf = x2 - xa */\
	xa -= x6;\
	x6 += x0;	/* sa = x2 + x6 */\
	xa += x0;	/* s5 = x2 - x6 + xa */\
/* In ASM, do xr and yi-terms and combine to get real parts of outputs,\
   then do xi and yr-terms and combine to get imaginary parts: */\
 	x0 = *(Bim+__odx[0x01]);\
 	x0 += x2;	\
 	*(Bim+__odx[0x01]) = x0;\
 	x2 *= 2;	\
 	x0 -= x2;	\
 	*(Bim+__odx[0x1E]) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(Bim+__odx[0x02]);	x0 = *(Bim+__odx[0x03]);\
 	x2 += x5;					x0 -= x3;	\
 	*(Bim+__odx[0x02]) = x2;	*(Bim+__odx[0x03]) = x0;\
 	x5 *= 2;	 				x3 *= 2;	\
 	x2 -= x5;	 				x0 += x3;	\
 	*(Bim+__odx[0x1D]) = x2;	*(Bim+__odx[0x1C]) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers toreduce serialization, but stop at 4, more likely no help: */\
 	x0 = *(Bim+__odx[0x04]);	x2 = *(Bim+__odx[0x05]);	x3 = *(Bim+__odx[0x06]);	x5 = *(Bim+__odx[0x07]);\
	x0 += xb;					x2 += x6;					x3 -= xa;					x5 += x1;	\
	*(Bim+__odx[0x04]) = x0;	*(Bim+__odx[0x05]) = x2;	*(Bim+__odx[0x06]) = x3;	*(Bim+__odx[0x07]) = x5;\
	xb *= 2;					x6 *= 2;					xa *= 2;					x1 *= 2;	\
	x0 -= xb;					x2 -= x6;					x3 += xa;					x5 -= x1;	\
	*(Bim+__odx[0x1B]) = x0;	*(Bim+__odx[0x1A]) = x2;	*(Bim+__odx[0x19]) = x3;	*(Bim+__odx[0x18]) = x5;\
\
 	x0 = *(Bim+__odx[0x08]);	x2 = *(Bim+__odx[0x09]);	x3 = *(Bim+__odx[0x0A]);	x5 = *(Bim+__odx[0x0B]);\
	x0 += xe;					x2 += xc;					x3 += x7;					x5 -= xf;	\
	*(Bim+__odx[0x08]) = x0;	*(Bim+__odx[0x09]) = x2;	*(Bim+__odx[0x0A]) = x3;	*(Bim+__odx[0x0B]) = x5;\
	xe *= 2;					xc *= 2;					x7 *= 2;					xf *= 2;	\
	x0 -= xe;					x2 -= xc;					x3 -= x7;					x5 += xf;	\
	*(Bim+__odx[0x17]) = x0;	*(Bim+__odx[0x16]) = x2;	*(Bim+__odx[0x15]) = x3;	*(Bim+__odx[0x14]) = x5;\
\
 	x0 = *(Bim+__odx[0x0C]);	x2 = *(Bim+__odx[0x0D]);	x3 = *(Bim+__odx[0x0E]);	x5 = *(Bim+__odx[0x0F]);\
	x0 -= xd;					x2 -= x9;					x3 += x4;					x5 -= x8;	\
	*(Bim+__odx[0x0C]) = x0;	*(Bim+__odx[0x0D]) = x2;	*(Bim+__odx[0x0E]) = x3;	*(Bim+__odx[0x0F]) = x5;\
	xd *= 2;					x9 *= 2;					x4 *= 2;					x8 *= 2;	\
	x0 += xd;					x2 += x9;					x3 -= x4;					x5 += x8;	\
	*(Bim+__odx[0x13]) = x0;	*(Bim+__odx[0x12]) = x2;	*(Bim+__odx[0x11]) = x3;	*(Bim+__odx[0x10]) = x5;\
\
/* Totals: 658 FADD, 234 FMUL. */\
}

/***************/
#if 0
int radix31_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
}
#endif
/***************/

void radix31_dif_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-31 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,arr_offsets[31];
	static int n31,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30, first_entry=TRUE;
	double tmp[62];

/*...initialize things upon first entry	*/

	if(!first_entry && (n/31) != n31)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}
	if(first_entry)
	{
		first_entry=FALSE;
		n31=n/31;

/*   constant index offsets for array load/stores are here.	*/

		p01= n31;
		p02= p01+p01;
		p03= p02+p01;
		p04= p03+p01;
		p05= p04+p01;
		p06= p05+p01;
		p07= p06+p01;
		p08= p07+p01;
		p09= p08+p01;
		p10= p09+p01;
		p11= p10+p01;
		p12= p11+p01;
		p13= p12+p01;
		p14= p13+p01;
		p15= p14+p01;
		p16= p15+p01;
		p17= p16+p01;
		p18= p17+p01;
		p19= p18+p01;
		p20= p19+p01;
		p21= p20+p01;
		p22= p21+p01;
		p23= p22+p01;
		p24= p23+p01;
		p25= p24+p01;
		p26= p25+p01;
		p27= p26+p01;
		p28= p27+p01;
		p29= p28+p01;
		p30= p29+p01;

		p01= p01+ ( (p01>> DAT_BITS) << PAD_BITS );
		p02= p02+ ( (p02>> DAT_BITS) << PAD_BITS );
		p03= p03+ ( (p03>> DAT_BITS) << PAD_BITS );
		p04= p04+ ( (p04>> DAT_BITS) << PAD_BITS );
		p05= p05+ ( (p05>> DAT_BITS) << PAD_BITS );
		p06= p06+ ( (p06>> DAT_BITS) << PAD_BITS );
		p07= p07+ ( (p07>> DAT_BITS) << PAD_BITS );
		p08= p08+ ( (p08>> DAT_BITS) << PAD_BITS );
		p09= p09+ ( (p09>> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
		p18= p18+ ( (p18>> DAT_BITS) << PAD_BITS );
		p19= p19+ ( (p19>> DAT_BITS) << PAD_BITS );
		p20= p20+ ( (p20>> DAT_BITS) << PAD_BITS );
		p21= p21+ ( (p21>> DAT_BITS) << PAD_BITS );
		p22= p22+ ( (p22>> DAT_BITS) << PAD_BITS );
		p23= p23+ ( (p23>> DAT_BITS) << PAD_BITS );
		p24= p24+ ( (p24>> DAT_BITS) << PAD_BITS );
		p25= p25+ ( (p25>> DAT_BITS) << PAD_BITS );
		p26= p26+ ( (p26>> DAT_BITS) << PAD_BITS );
		p27= p27+ ( (p27>> DAT_BITS) << PAD_BITS );
		p28= p28+ ( (p28>> DAT_BITS) << PAD_BITS );
		p29= p29+ ( (p29>> DAT_BITS) << PAD_BITS );
		p30= p30+ ( (p30>> DAT_BITS) << PAD_BITS );

		arr_offsets[0x00] = 0;
		arr_offsets[0x01] = p01;
		arr_offsets[0x02] = p02;
		arr_offsets[0x03] = p03;
		arr_offsets[0x04] = p04;
		arr_offsets[0x05] = p05;
		arr_offsets[0x06] = p06;
		arr_offsets[0x07] = p07;
		arr_offsets[0x08] = p08;
		arr_offsets[0x09] = p09;
		arr_offsets[0x0A] = p10;
		arr_offsets[0x0B] = p11;
		arr_offsets[0x0C] = p12;
		arr_offsets[0x0D] = p13;
		arr_offsets[0x0E] = p14;
		arr_offsets[0x0F] = p15;
		arr_offsets[0x10] = p16;
		arr_offsets[0x11] = p17;
		arr_offsets[0x12] = p18;
		arr_offsets[0x13] = p19;
		arr_offsets[0x14] = p20;
		arr_offsets[0x15] = p21;
		arr_offsets[0x16] = p22;
		arr_offsets[0x17] = p23;
		arr_offsets[0x18] = p24;
		arr_offsets[0x19] = p25;
		arr_offsets[0x1A] = p26;
		arr_offsets[0x1B] = p27;
		arr_offsets[0x1C] = p28;
		arr_offsets[0x1D] = p29;
		arr_offsets[0x1E] = p30;
	}

/*...The radix-31 pass is here.	*/

	for(j=0; j < n31; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	#if 0
		RADIX_31_DFT(a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
	/* outputs: */	,a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30);
	#else	/* Out-of-place version, needed for combining radix-31 with powers-of-2: */
		RADIX_31_DIF(\
			a+j1,arr_offsets,\
			tmp\
		);
		a[j1    ] = tmp[ 0];	a[j2    ] = tmp[ 1];
		a[j1+p01] = tmp[ 2];	a[j2+p01] = tmp[ 3];
		a[j1+p02] = tmp[ 4];	a[j2+p02] = tmp[ 5];
		a[j1+p03] = tmp[ 6];	a[j2+p03] = tmp[ 7];
		a[j1+p04] = tmp[ 8];	a[j2+p04] = tmp[ 9];
		a[j1+p05] = tmp[10];	a[j2+p05] = tmp[11];
		a[j1+p06] = tmp[12];	a[j2+p06] = tmp[13];
		a[j1+p07] = tmp[14];	a[j2+p07] = tmp[15];
		a[j1+p08] = tmp[16];	a[j2+p08] = tmp[17];
		a[j1+p09] = tmp[18];	a[j2+p09] = tmp[19];
		a[j1+p10] = tmp[20];	a[j2+p10] = tmp[21];
		a[j1+p11] = tmp[22];	a[j2+p11] = tmp[23];
		a[j1+p12] = tmp[24];	a[j2+p12] = tmp[25];
		a[j1+p13] = tmp[26];	a[j2+p13] = tmp[27];
		a[j1+p14] = tmp[28];	a[j2+p14] = tmp[29];
		a[j1+p15] = tmp[30];	a[j2+p15] = tmp[31];
		a[j1+p16] = tmp[32];	a[j2+p16] = tmp[33];
		a[j1+p17] = tmp[34];	a[j2+p17] = tmp[35];
		a[j1+p18] = tmp[36];	a[j2+p18] = tmp[37];
		a[j1+p19] = tmp[38];	a[j2+p19] = tmp[39];
		a[j1+p20] = tmp[40];	a[j2+p20] = tmp[41];
		a[j1+p21] = tmp[42];	a[j2+p21] = tmp[43];
		a[j1+p22] = tmp[44];	a[j2+p22] = tmp[45];
		a[j1+p23] = tmp[46];	a[j2+p23] = tmp[47];
		a[j1+p24] = tmp[48];	a[j2+p24] = tmp[49];
		a[j1+p25] = tmp[50];	a[j2+p25] = tmp[51];
		a[j1+p26] = tmp[52];	a[j2+p26] = tmp[53];
		a[j1+p27] = tmp[54];	a[j2+p27] = tmp[55];
		a[j1+p28] = tmp[56];	a[j2+p28] = tmp[57];
		a[j1+p29] = tmp[58];	a[j2+p29] = tmp[59];
		a[j1+p30] = tmp[60];	a[j2+p30] = tmp[61];
	#endif
	}
}

/*
DIT - Since the radix is prime - differs from DIF only in having the + and - of the final output-combinations swapped.
*/
void radix31_dit_pass1(double a[], int n)
{
	int j,j1,j2;
	static int n31,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30, first_entry=TRUE;

	if(!first_entry && (n/31) != n31)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n31=n/31;

/*   constant index offsets for array load/stores are here.	*/

		p01= n31;
		p02= p01+p01;
		p03= p02+p01;
		p04= p03+p01;
		p05= p04+p01;
		p06= p05+p01;
		p07= p06+p01;
		p08= p07+p01;
		p09= p08+p01;
		p10= p09+p01;
		p11= p10+p01;
		p12= p11+p01;
		p13= p12+p01;
		p14= p13+p01;
		p15= p14+p01;
		p16= p15+p01;
		p17= p16+p01;
		p18= p17+p01;
		p19= p18+p01;
		p20= p19+p01;
		p21= p20+p01;
		p22= p21+p01;
		p23= p22+p01;
		p24= p23+p01;
		p25= p24+p01;
		p26= p25+p01;
		p27= p26+p01;
		p28= p27+p01;
		p29= p28+p01;
		p30= p29+p01;

		p01= p01+ ( (p01>> DAT_BITS) << PAD_BITS );
		p02= p02+ ( (p02>> DAT_BITS) << PAD_BITS );
		p03= p03+ ( (p03>> DAT_BITS) << PAD_BITS );
		p04= p04+ ( (p04>> DAT_BITS) << PAD_BITS );
		p05= p05+ ( (p05>> DAT_BITS) << PAD_BITS );
		p06= p06+ ( (p06>> DAT_BITS) << PAD_BITS );
		p07= p07+ ( (p07>> DAT_BITS) << PAD_BITS );
		p08= p08+ ( (p08>> DAT_BITS) << PAD_BITS );
		p09= p09+ ( (p09>> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
		p15= p15+ ( (p15>> DAT_BITS) << PAD_BITS );
		p16= p16+ ( (p16>> DAT_BITS) << PAD_BITS );
		p17= p17+ ( (p17>> DAT_BITS) << PAD_BITS );
		p18= p18+ ( (p18>> DAT_BITS) << PAD_BITS );
		p19= p19+ ( (p19>> DAT_BITS) << PAD_BITS );
		p20= p20+ ( (p20>> DAT_BITS) << PAD_BITS );
		p21= p21+ ( (p21>> DAT_BITS) << PAD_BITS );
		p22= p22+ ( (p22>> DAT_BITS) << PAD_BITS );
		p23= p23+ ( (p23>> DAT_BITS) << PAD_BITS );
		p24= p24+ ( (p24>> DAT_BITS) << PAD_BITS );
		p25= p25+ ( (p25>> DAT_BITS) << PAD_BITS );
		p26= p26+ ( (p26>> DAT_BITS) << PAD_BITS );
		p27= p27+ ( (p27>> DAT_BITS) << PAD_BITS );
		p28= p28+ ( (p28>> DAT_BITS) << PAD_BITS );
		p29= p29+ ( (p29>> DAT_BITS) << PAD_BITS );
		p30= p30+ ( (p30>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-31 pass is here.	*/

	for(j=0; j < n31; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		RADIX_31_DFT(a+j1,0,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
	/* outputs: */	,a+j1,0,p30,p29,p28,p27,p26,p25,p24,p23,p22,p21,p20,p19,p18,p17,p16,p15,p14,p13,p12,p11,p10,p09,p08,p07,p06,p05,p04,p03,p02,p01);
	}
}

/***************/

void radix31x32_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-992 = 31x32 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2, *optr,*qptr;
	static int n992,first_entry=TRUE,p[992],q[992];
	static int o[992] = {
	   0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	 975,974,972,973,968,969,970,971,960,961,962,963,964,965,966,967,991,990,988,989,984,985,986,987,976,977,978,979,980,981,982,983,
	 951,950,948,949,944,945,946,947,959,958,956,957,952,953,954,955,943,942,940,941,936,937,938,939,928,929,930,931,932,933,934,935,
	 903,902,900,901,896,897,898,899,911,910,908,909,904,905,906,907,919,918,916,917,912,913,914,915,927,926,924,925,920,921,922,923,
	 891,890,888,889,895,894,892,893,887,886,884,885,880,881,882,883,871,870,868,869,864,865,866,867,879,878,876,877,872,873,874,875,
	 843,842,840,841,847,846,844,845,839,838,836,837,832,833,834,835,859,858,856,857,863,862,860,861,855,854,852,853,848,849,850,851,
	 819,818,816,817,823,822,820,821,827,826,824,825,831,830,828,829,811,810,808,809,815,814,812,813,807,806,804,805,800,801,802,803,
	 771,770,768,769,775,774,772,773,779,778,776,777,783,782,780,781,787,786,784,785,791,790,788,789,795,794,792,793,799,798,796,797,
	 765,764,767,766,763,762,760,761,755,754,752,753,759,758,756,757,739,738,736,737,743,742,740,741,747,746,744,745,751,750,748,749,
	 717,716,719,718,715,714,712,713,707,706,704,705,711,710,708,709,733,732,735,734,731,730,728,729,723,722,720,721,727,726,724,725,
	 693,692,695,694,691,690,688,689,701,700,703,702,699,698,696,697,685,684,687,686,683,682,680,681,675,674,672,673,679,678,676,677,
	 645,644,647,646,643,642,640,641,653,652,655,654,651,650,648,649,661,660,663,662,659,658,656,657,669,668,671,670,667,666,664,665,
	 633,632,635,634,637,636,639,638,629,628,631,630,627,626,624,625,613,612,615,614,611,610,608,609,621,620,623,622,619,618,616,617,
	 585,584,587,586,589,588,591,590,581,580,583,582,579,578,576,577,601,600,603,602,605,604,607,606,597,596,599,598,595,594,592,593,
	 561,560,563,562,565,564,567,566,569,568,571,570,573,572,575,574,553,552,555,554,557,556,559,558,549,548,551,550,547,546,544,545,
	 513,512,515,514,517,516,519,518,521,520,523,522,525,524,527,526,529,528,531,530,533,532,535,534,537,536,539,538,541,540,543,542,
	 510,511,509,508,505,504,507,506,497,496,499,498,501,500,503,502,481,480,483,482,485,484,487,486,489,488,491,490,493,492,495,494,
	 462,463,461,460,457,456,459,458,449,448,451,450,453,452,455,454,478,479,477,476,473,472,475,474,465,464,467,466,469,468,471,470,
	 438,439,437,436,433,432,435,434,446,447,445,444,441,440,443,442,430,431,429,428,425,424,427,426,417,416,419,418,421,420,423,422,
	 390,391,389,388,385,384,387,386,398,399,397,396,393,392,395,394,406,407,405,404,401,400,403,402,414,415,413,412,409,408,411,410,
	 378,379,377,376,382,383,381,380,374,375,373,372,369,368,371,370,358,359,357,356,353,352,355,354,366,367,365,364,361,360,363,362,
	 330,331,329,328,334,335,333,332,326,327,325,324,321,320,323,322,346,347,345,344,350,351,349,348,342,343,341,340,337,336,339,338,
	 306,307,305,304,310,311,309,308,314,315,313,312,318,319,317,316,298,299,297,296,302,303,301,300,294,295,293,292,289,288,291,290,
	 258,259,257,256,262,263,261,260,266,267,265,264,270,271,269,268,274,275,273,272,278,279,277,276,282,283,281,280,286,287,285,284,
	 252,253,254,255,250,251,249,248,242,243,241,240,246,247,245,244,226,227,225,224,230,231,229,228,234,235,233,232,238,239,237,236,
	 204,205,206,207,202,203,201,200,194,195,193,192,198,199,197,196,220,221,222,223,218,219,217,216,210,211,209,208,214,215,213,212,
	 180,181,182,183,178,179,177,176,188,189,190,191,186,187,185,184,172,173,174,175,170,171,169,168,162,163,161,160,166,167,165,164,
	 132,133,134,135,130,131,129,128,140,141,142,143,138,139,137,136,148,149,150,151,146,147,145,144,156,157,158,159,154,155,153,152,
	 120,121,122,123,124,125,126,127,116,117,118,119,114,115,113,112,100,101,102,103, 98, 99, 97, 96,108,109,110,111,106,107,105,104,
	  72, 73, 74, 75, 76, 77, 78, 79, 68, 69, 70, 71, 66, 67, 65, 64, 88, 89, 90, 91, 92, 93, 94, 95, 84, 85, 86, 87, 82, 83, 81, 80,
	  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 40, 41, 42, 43, 44, 45, 46, 47, 36, 37, 38, 39, 34, 35, 33, 32};
	int i1,i2,jj[32];
	double tmp[1984];
	static double    c     = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp(  i*twopi/16)	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32)	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
const int dat_bits = DAT_BITS+1, pad_bits = PAD_BITS+1;
	if(!first_entry && (n/992) != n992)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n992=n/992;
		/* Constant padded-array index offsets for array load/stores are here. */
		for(j = 0; j < 992; ++j)
		{
			p[j] = j*n992;
			p[j] = p[j] + ( (p[j] >> dat_bits) << pad_bits );
		}
		for(j = 0; j < 992; ++j)
		{
			o[j] = p[o[j]];
		}
	}
	/* Now that have o-indices, mimic the radix-31-related indexing in the loop to predefine a suitable permutation array for those: */
	jj[0] = 992;
	for(i1 = 1; i1 < 31; ++i1)
	{
		jj[i1] = jj[i1-1] - 32;
	}
	jj[0] = 0;
	for(i2 = 0; i2 < 32; ++i2)
	{
		/* Without the extra auxiliary q-array, the indexing here would be based on the p-array like so:
		RADIX_31_DIF(
			a+j1,p[jj[0]],p[jj[1]],p[jj[2]],p[jj[3]],p[jj[4]],p[jj[5]],p[jj[6]],p[jj[7]],p[jj[8]],p[jj[9]],p[jj[10]],p[jj[11]],p[jj[12]],p[jj[13]],p[jj[14]],p[jj[15]],p[jj[16]],p[jj[17]],p[jj[18]],p[jj[19]],p[jj[20]],p[jj[21]],p[jj[22]],p[jj[23]],p[jj[24]],p[jj[25]],p[jj[26]],p[jj[27]],p[jj[28]],p[jj[29]],p[jj[30]],
			tmp+i2*62
		);
		We want to replace the 31 disparate p[jj[0-30]] indices with a block of 31 contiguous memory locations, so do like so:
		*/
		for(i1 = 0; i1 < 31; ++i1)
		{
			q[i2*31 + i1] = p[jj[i1]];
			jj[i1] -= 31;
			jj[i1] += ( (-(int)((uint32)jj[i1] >> 31)) & 992);	/* (jj[i1]-62)%1984; */
		}
	}

/*...The radix-992 pass is here.	*/

    for(j=0; j < n992; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> dat_bits) << pad_bits );
	#else
		j1 = j + ( (j >> dat_bits) << pad_bits );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (992 64-bit complex, i.e 1984 64-bit reals) and do 32 radix-31 transforms...*/
	/*
	Twiddleless version arranges 32 sets of radix-31 DFT inputs as follows: 0 in upper left corner,
	decrement 32 horizontally and 31 vertically, all resulting indexing modulo 992:

		RADIX_31_DFT(000,960,928,896,864,832,800,768,736,704,672,640,608,576,544,512,480,448,416,384,352,320,288,256,224,192,160,128,096,064,032)
		RADIX_31_DFT(961,                                                    ...                                                            ,001)
		RADIX_31_DFT(930,                                                    ...                                                            ,962)
		RADIX_31_DFT(899,                                                    ...                                                            ,931)
		RADIX_31_DFT(868,                                                    ...                                                            ,900)
		RADIX_31_DFT(837,                                                    ...                                                            ,869)
		RADIX_31_DFT(806,                                                    ...                                                            ,838)
		RADIX_31_DFT(775,                                                    ...                                                            ,807)
		RADIX_31_DFT(744,                                                    ...                                                            ,776)
		RADIX_31_DFT(713,                                                    ...                                                            ,745)
		RADIX_31_DFT(682,                                                    ...                                                            ,714)
		RADIX_31_DFT(651,                                                    ...                                                            ,683)
		RADIX_31_DFT(620,                                                    ...                                                            ,652)
		RADIX_31_DFT(589,                                                    ...                                                            ,621)
		RADIX_31_DFT(558,                                                    ...                                                            ,590)
		RADIX_31_DFT(527,                                                    ...                                                            ,559)
		RADIX_31_DFT(496,                                                    ...                                                            ,528)
		RADIX_31_DFT(465,                                                    ...                                                            ,497)
		RADIX_31_DFT(434,                                                    ...                                                            ,466)
		RADIX_31_DFT(403,                                                    ...                                                            ,435)
		RADIX_31_DFT(372,                                                    ...                                                            ,404)
		RADIX_31_DFT(341,                                                    ...                                                            ,373)
		RADIX_31_DFT(310,                                                    ...                                                            ,342)
		RADIX_31_DFT(279,                                                    ...                                                            ,311)
		RADIX_31_DFT(248,                                                    ...                                                            ,280)
		RADIX_31_DFT(217,                                                    ...                                                            ,249)
		RADIX_31_DFT(186,                                                    ...                                                            ,218)
		RADIX_31_DFT(155,                                                    ...                                                            ,187)
		RADIX_31_DFT(124,                                                    ...                                                            ,156)
		RADIX_31_DFT(093,                                                    ...                                                            ,125)
		RADIX_31_DFT(062,                                                    ...                                                            ,094)
		RADIX_31_DFT(031,                                                    ...                                                            ,063)

	[The above indexing is w.r.to complex-paired data, so need to double all offsets in the actual DFT-macro calls, which are based on a real data array]

	We generate a length-992 table of p-offsets and use an assembly-line scheme to manage the 'wild dance of the indices' here.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-32 DFT outputs.
	*/
		for(i2 = 0; i2 < 32; ++i2)
		{
			qptr = q+i2*31;
			RADIX_31_DIF(
				a+j1,qptr,
				tmp+i2*62
			);
		}
		/*...and now do 31 radix-32 transforms:

		RADIX_32_DIF(tmp+ 0, 62,124,...,1922, a+j1,p[000],p[001],p[002],p[003],p[004],p[005],p[006],p[007],p[008],p[009],p[010],p[011],p[012],p[013],p[014],p[015],p[016],p[017],p[018],p[019],p[020],p[021],p[022],p[023],p[024],p[025],p[026],p[027],p[028],p[029],p[030],p[031]);
		RADIX_32_DIF(tmp+ 2, 64,126,...,1924, a+j1,p[031],...
		...                                        ...
		RADIX_32_DIF(tmp+60,122,184,...,1982, a+j1,p[960],p[961],p[962],p[963],p[964],p[965],p[966],p[967],p[968],p[969],p[970],p[971],p[972],p[973],p[974],p[975],p[976],p[977],p[978],p[979],p[980],p[981],p[982],p[983],p[984],p[985],p[986],p[987],p[988],p[989],p[990],p[991]);

		The required output permutation is stored in the o-array.

		In the implementation of the above, the jj-offsets here are raw real-array offsets, so need to be doubled relative to complex-indexing;
		The optr-offset is an index into the length-992 o-array, whose elements implicitly contain the needed doubling.
		*/
		jj[0] = 0;
		for(i2 = 1; i2 < 32; ++i2) {
			jj[i2] = jj[i2-1] + 62;	/* (i2*62) */
		}
		for(i1 = 0; i1 < 31; ++i1) {
			optr = o + (i1 << 5);
			RADIX_32_DIF(
				tmp,jj,
				a+j1,optr
			);
			for(i2 = 0; i2 < 32; ++i2) {
				jj[i2] = jj[i2] + 2;
			}
		}
		/* Totals: 32*radix31 + 31*radix32 = 32*(??? FADD, ??? FMUL) + 31*(??? FADD, ??? FMUL) = ????? FADD, ????? FMUL	*/
	}
}

/***************/

void radix31x32_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-992 = 31x32 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix992_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,j1,j2, *optr,*qptr;
	static int n992,first_entry=TRUE,p[992];
	static int q[992] = {0,1,3,2,7,6,5,4,15,14,13,12,11,10,9,8,31,30,29,28,27,26,25,24,
	23,22,21,20,19,18,17,16,975,974,973,972,971,970,969,968,967,966,965,964,963,962,961,960,983,982,981,
	980,979,978,977,976,987,986,985,984,989,988,990,991,951,950,949,948,947,946,945,944,955,954,953,952,
	957,956,958,959,935,934,933,932,931,930,929,928,939,938,937,936,941,940,942,943,903,902,901,900,899,
	898,897,896,907,906,905,904,909,908,910,911,923,922,921,920,925,924,926,927,915,914,913,912,917,916,
	918,919,891,890,889,888,893,892,894,895,883,882,881,880,885,884,886,887,875,874,873,872,877,876,878,
	879,867,866,865,864,869,868,870,871,843,842,841,840,845,844,846,847,835,834,833,832,837,836,838,839,
	851,850,849,848,853,852,854,855,861,860,862,863,857,856,858,859,819,818,817,816,821,820,822,823,829,
	828,830,831,825,824,826,827,803,802,801,800,805,804,806,807,813,812,814,815,809,808,810,811,771,770,
	769,768,773,772,774,775,781,780,782,783,777,776,778,779,797,796,798,799,793,792,794,795,789,788,790,
	791,785,784,786,787,765,764,766,767,761,760,762,763,757,756,758,759,753,752,754,755,749,748,750,751,
	745,744,746,747,741,740,742,743,737,736,738,739,717,716,718,719,713,712,714,715,709,708,710,711,705,
	704,706,707,725,724,726,727,721,720,722,723,729,728,730,731,734,735,732,733,693,692,694,695,689,688,
	690,691,697,696,698,699,702,703,700,701,677,676,678,679,673,672,674,675,681,680,682,683,686,687,684,
	685,645,644,646,647,641,640,642,643,649,648,650,651,654,655,652,653,665,664,666,667,670,671,668,669,
	657,656,658,659,662,663,660,661,633,632,634,635,638,639,636,637,625,624,626,627,630,631,628,629,617,
	616,618,619,622,623,620,621,609,608,610,611,614,615,612,613,585,584,586,587,590,591,588,589,577,576,
	578,579,582,583,580,581,593,592,594,595,598,599,596,597,606,607,604,605,602,603,600,601,561,560,562,
	563,566,567,564,565,574,575,572,573,570,571,568,569,545,544,546,547,550,551,548,549,558,559,556,557,
	554,555,552,553,513,512,514,515,518,519,516,517,526,527,524,525,522,523,520,521,542,543,540,541,538,
	539,536,537,534,535,532,533,530,531,528,529,510,511,508,509,506,507,504,505,502,503,500,501,498,499,
	496,497,494,495,492,493,490,491,488,489,486,487,484,485,482,483,480,481,462,463,460,461,458,459,456,
	457,454,455,452,453,450,451,448,449,470,471,468,469,466,467,464,465,474,475,472,473,476,477,479,478,
	438,439,436,437,434,435,432,433,442,443,440,441,444,445,447,446,422,423,420,421,418,419,416,417,426,
	427,424,425,428,429,431,430,390,391,388,389,386,387,384,385,394,395,392,393,396,397,399,398,410,411,
	408,409,412,413,415,414,402,403,400,401,404,405,407,406,378,379,376,377,380,381,383,382,370,371,368,
	369,372,373,375,374,362,363,360,361,364,365,367,366,354,355,352,353,356,357,359,358,330,331,328,329,
	332,333,335,334,322,323,320,321,324,325,327,326,338,339,336,337,340,341,343,342,348,349,351,350,344,
	345,347,346,306,307,304,305,308,309,311,310,316,317,319,318,312,313,315,314,290,291,288,289,292,293,
	295,294,300,301,303,302,296,297,299,298,258,259,256,257,260,261,263,262,268,269,271,270,264,265,267,
	266,284,285,287,286,280,281,283,282,276,277,279,278,272,273,275,274,252,253,255,254,248,249,251,250,
	244,245,247,246,240,241,243,242,236,237,239,238,232,233,235,234,228,229,231,230,224,225,227,226,204,
	205,207,206,200,201,203,202,196,197,199,198,192,193,195,194,212,213,215,214,208,209,211,210,216,217,
	219,218,223,222,221,220,180,181,183,182,176,177,179,178,184,185,187,186,191,190,189,188,164,165,167,
	166,160,161,163,162,168,169,171,170,175,174,173,172,132,133,135,134,128,129,131,130,136,137,139,138,
	143,142,141,140,152,153,155,154,159,158,157,156,144,145,147,146,151,150,149,148,120,121,123,122,127,
	126,125,124,112,113,115,114,119,118,117,116,104,105,107,106,111,110,109,108,96,97,99,98,103,102,101,
	100,72,73,75,74,79,78,77,76,64,65,67,66,71,70,69,68,80,81,83,82,87,86,85,84,95,94,93,92,91,90,89,88,
	48,49,51,50,55,54,53,52,63,62,61,60,59,58,57,56,32,33,35,34,39,38,37,36,47,46,45,44,43,42,41,40};
	static int o[992] = {
	 0, 32, 64, 96,128,160,192,224,256,288,320,352,384,416,448,480,512,544,576,608,640,672,704,736,768,800,832,864,896,928,960,
	31, 63, 95,127,159,191,223,255,287,319,351,383,415,447,479,511,543,575,607,639,671,703,735,767,799,831,863,895,927,959,991,
	62, 94,126,158,190,222,254,286,318,350,382,414,446,478,510,542,574,606,638,670,702,734,766,798,830,862,894,926,958,990, 30,
	93,125,157,189,221,253,285,317,349,381,413,445,477,509,541,573,605,637,669,701,733,765,797,829,861,893,925,957,989, 29, 61,
	124,156,188,220,252,284,316,348,380,412,444,476,508,540,572,604,636,668,700,732,764,796,828,860,892,924,956,988, 28, 60, 92,
	155,187,219,251,283,315,347,379,411,443,475,507,539,571,603,635,667,699,731,763,795,827,859,891,923,955,987, 27, 59, 91,123,
	186,218,250,282,314,346,378,410,442,474,506,538,570,602,634,666,698,730,762,794,826,858,890,922,954,986, 26, 58, 90,122,154,
	217,249,281,313,345,377,409,441,473,505,537,569,601,633,665,697,729,761,793,825,857,889,921,953,985, 25, 57, 89,121,153,185,
	248,280,312,344,376,408,440,472,504,536,568,600,632,664,696,728,760,792,824,856,888,920,952,984, 24, 56, 88,120,152,184,216,
	279,311,343,375,407,439,471,503,535,567,599,631,663,695,727,759,791,823,855,887,919,951,983, 23, 55, 87,119,151,183,215,247,
	310,342,374,406,438,470,502,534,566,598,630,662,694,726,758,790,822,854,886,918,950,982, 22, 54, 86,118,150,182,214,246,278,
	341,373,405,437,469,501,533,565,597,629,661,693,725,757,789,821,853,885,917,949,981, 21, 53, 85,117,149,181,213,245,277,309,
	372,404,436,468,500,532,564,596,628,660,692,724,756,788,820,852,884,916,948,980, 20, 52, 84,116,148,180,212,244,276,308,340,
	403,435,467,499,531,563,595,627,659,691,723,755,787,819,851,883,915,947,979, 19, 51, 83,115,147,179,211,243,275,307,339,371,
	434,466,498,530,562,594,626,658,690,722,754,786,818,850,882,914,946,978, 18, 50, 82,114,146,178,210,242,274,306,338,370,402,
	465,497,529,561,593,625,657,689,721,753,785,817,849,881,913,945,977, 17, 49, 81,113,145,177,209,241,273,305,337,369,401,433,
	496,528,560,592,624,656,688,720,752,784,816,848,880,912,944,976, 16, 48, 80,112,144,176,208,240,272,304,336,368,400,432,464,
	527,559,591,623,655,687,719,751,783,815,847,879,911,943,975, 15, 47, 79,111,143,175,207,239,271,303,335,367,399,431,463,495,
	558,590,622,654,686,718,750,782,814,846,878,910,942,974, 14, 46, 78,110,142,174,206,238,270,302,334,366,398,430,462,494,526,
	589,621,653,685,717,749,781,813,845,877,909,941,973, 13, 45, 77,109,141,173,205,237,269,301,333,365,397,429,461,493,525,557,
	620,652,684,716,748,780,812,844,876,908,940,972, 12, 44, 76,108,140,172,204,236,268,300,332,364,396,428,460,492,524,556,588,
	651,683,715,747,779,811,843,875,907,939,971, 11, 43, 75,107,139,171,203,235,267,299,331,363,395,427,459,491,523,555,587,619,
	682,714,746,778,810,842,874,906,938,970, 10, 42, 74,106,138,170,202,234,266,298,330,362,394,426,458,490,522,554,586,618,650,
	713,745,777,809,841,873,905,937,969,  9, 41, 73,105,137,169,201,233,265,297,329,361,393,425,457,489,521,553,585,617,649,681,
	744,776,808,840,872,904,936,968,  8, 40, 72,104,136,168,200,232,264,296,328,360,392,424,456,488,520,552,584,616,648,680,712,
	775,807,839,871,903,935,967,  7, 39, 71,103,135,167,199,231,263,295,327,359,391,423,455,487,519,551,583,615,647,679,711,743,
	806,838,870,902,934,966,  6, 38, 70,102,134,166,198,230,262,294,326,358,390,422,454,486,518,550,582,614,646,678,710,742,774,
	837,869,901,933,965,  5, 37, 69,101,133,165,197,229,261,293,325,357,389,421,453,485,517,549,581,613,645,677,709,741,773,805,
	868,900,932,964,  4, 36, 68,100,132,164,196,228,260,292,324,356,388,420,452,484,516,548,580,612,644,676,708,740,772,804,836,
	899,931,963,  3, 35, 67, 99,131,163,195,227,259,291,323,355,387,419,451,483,515,547,579,611,643,675,707,739,771,803,835,867,
	930,962,  2, 34, 66, 98,130,162,194,226,258,290,322,354,386,418,450,482,514,546,578,610,642,674,706,738,770,802,834,866,898,
	961,  1, 33, 65, 97,129,161,193,225,257,289,321,353,385,417,449,481,513,545,577,609,641,673,705,737,769,801,833,865,897,929};
	int i1,i2,jj[32];
	double tmp[1984];
	static double    c     = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp(  i*twopi/16)	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32)	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
const int dat_bits = DAT_BITS+1, pad_bits = PAD_BITS+1;
	if(!first_entry && (n/992) != n992)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n992=n/992;
		/* Constant padded-array index offsets for array load/stores are here. */
		for(j = 0; j < 992; ++j)
		{
			p[j] = j*n992;
			p[j] = p[j] + ( (p[j] >> dat_bits) << pad_bits );
		}
		for(j = 0; j < 992; ++j)
		{
			q[j] = p[q[j]];
			o[j] = p[o[j]];
		}
	}

/*...The radix-992 pass is here.	*/

    for(j=0; j < n992; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs...Combined DIT input-scramble array is stored in the q-array.
	*/
	/*...gather the needed data (992 64-bit complex, i.e. 1984 64-bit reals) and do 31 radix-32 transforms...*/
		/*
		The jj-offsets here are raw real-array offsets, so need to be doubled relative to complex-indexing;
		The qptr-offset is an index into the length-992 q-array, whose elements implicitly contain the needed doubling.
		*/
		jj[0] = 0;
		for(i2 = 1; i2 < 32; ++i2)
		{
			jj[i2] = jj[i2-1] + 62;	/* (I2*62) */
		}
		for(i1 = 0; i1 < 31; ++i1)
		{
			qptr = q + (i1 << 5);
			RADIX_32_DIT(
				a+j1,qptr,
				tmp,jj
			);
			for(i2 = 0; i2 < 32; ++i2)
			{
				jj[i2] += 2;
			}
		}
		for(i2 = 0; i2 < 32; ++i2)
		{
			optr = o+i2*31;
			RADIX_31_DIT(
				tmp+i2*62,
				a+j1,optr
			);
		}
		/*...and now do 32 radix-31 transforms.
		The required output permutation is stored in the o-array.
		*/
	}
}

