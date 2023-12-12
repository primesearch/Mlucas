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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef radix31_included
#define radix31_included

/*
Van Buskirk's original constants in 'JVB' commented blocks, EWM-derived ones under those initials.
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

	DC3,4,8,3: call auxiliary consts B_c3483 and D_c3483
	DC5,6,-4,5:    "          "      B_c5645 and D_c5645
	DC9,-10,11,-9: "          "      B_c9019 and D_c9019

	DC15,12,18,-15:  "        "      B_c5285 and D_c5285
	DC16,-13,19,16:  "        "      B_c6396 and D_c6396
	DC17,14,20,-17:  "        "      B_c7407 and D_c7407

	DC24,21,27,-24:  "        "      C_c4174 and E_c4174 (Use C,E to avoid name collision with _c5285 and _c6396 terms above)
	DC25,-22,28,25:  "        "      C_c5285 and E_c5285
	DC26,23,29,-26:  "        "      C_c6396 and E_c6396

...and there is an analogous set of terms based on the sine (DS*) components, named [BB,C,D,E]_s****.

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

If have 0 spare registers, best I have found so far requires 5 add/sub, 3 mul:

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

// JVB:
#define DC0    ((double)-0.16666666666666666666666666666666700)
#define DC1    ((double)-0.41848031916582736768981687098748646)
#define DC2    ((double) 1.03708005397391977997160126838915860)
#define DC3    ((double) 0.06899615305609372142805240007781954)
#define DC4    ((double)-0.71199788901815727304919706300135852)
#define DC5    ((double) 0.49457655717669416346035981441598679)
#define DC6    ((double)-0.96808390091782605854354886831133920)
#define DC7    ((double) 0.32508216495576250692240420538780022)
#define DC8    ((double) 0.16949439222093165653795560902818660)
#define DC9    ((double)-0.29637372110299413755460095857226913)
#define DC10   ((double) 0.04534684817389996223192362526889384)
#define DC11   ((double)-0.34172056927689409978652458384116289)
#define DC12   ((double) 0.07310633583025898405268087482211124)
#define DC13   ((double)-0.20167486904767495859952810665520680)
#define DC14   ((double)-0.84201778467603236410459766557180865)
#define DC15   ((double)-0.13041239332209095847634591191071144)
#define DC16   ((double)-0.45935884940861915763712570904324634)
#define DC17   ((double) 0.48562853299381401878246608937483484)
#define DC18   ((double) 0.05290024156536825355103499921807512)
#define DC19   ((double)-0.14593330617823442758015434230035912)
#define DC20   ((double)-0.60928979281795509222872426261479928)
#define DC21   ((double) 0.50714554918236722335933886804513386)
#define DC22   ((double)-0.09614247372823604317748361539732969)
#define DC23   ((double)-0.57652834602834558094251739609128330)
#define DC24   ((double) 0.06326947458204131524249727162375559)
#define DC25   ((double)-0.76932220806904896036979558584064568)
#define DC26   ((double)-0.41013896184380627585461186107586186)
#define DC27   ((double) 0.36697396683700721180660373150156373)
#define DC28   ((double)-0.06956934754225036504904405641510024)
#define DC29   ((double)-0.41717983028166295195804014713757581)
// EWM:
#define DC0_5th ((double)-0.033333333333333333333333333333333333)
/* 046,1[0234],2[123] */
#define B_c3483 ((double)-10.31938532050192055740019444)
#define D_c3483 ((double)-1.039447157476621531648472621)
#define nDC4    ((double) 0.71199788901815727304919706300135852)
#define B_c5645 ((double)-1.957399490271361553945056740)
#define D_c5645 ((double)-1.354874954742460424002602037)
#define B_c9019 ((double) 0.1530056308809554651383615683)
#define D_c9019 ((double)-6.668408724817368094978521568)
#define B_c5285 ((double)-0.5605781319395146459386195222)
#define D_c5285 ((double)-5.397694515310063808645584290)
#define B_c6396 ((double)-0.4390355585993655658779099885)
#define D_c6396 ((double)-8.169650232597136939166445099)
#define B_c7407 ((double)-1.733872141912958941825760722)
#define D_c7407 ((double)-1.459688060776070564000064511)
#define C_c4174 ((double) 8.015643444687584573541352865)
#define E_c4174 ((double)-1.021509017898617060353255097)
#define C_c5285 ((double)-0.1249703605587412988133367679)
#define E_c5285 ((double)-89.48778340022079724273433851)
#define C_c6396 ((double) 1.405690265164092306971350691)
#define E_c6396 ((double)-1.699387856677906977066853648)

// JVB:
#define DS0    ((double) 0.92796072713833698701991188315309154)
#define DS1    ((double) 0.75538558300353835543581057551249341)
#define DS2    ((double)-0.71798870119142816504365819783899618)
#define DS3    ((double)-1.14129841620797351615515640500369610)
#define DS4    ((double) 0.17619908884183581953599075099065958)
#define DS5    ((double)-0.44789045813036076345265894471834662)
#define DS6    ((double)-0.42330971501654535111149820716469997)
#define DS7    ((double)-0.54178961234959234550766744684833662)
#define DS8    ((double) 0.09389915421923158205500850212998998)
#define DS9    ((double) 0.10209749786491606368824206751661151)
#define DS10   ((double) 0.36010442196019251577804178188101286)
#define DS11   ((double)-0.25800692409527645208979971436440136)
#define DS12   ((double)-0.68362751420088060523602271315212364)
#define DS13   ((double) 1.42448046407360845685461718962915170)
#define DS14   ((double) 0.71705834223306922870900492103008074)
#define DS15   ((double)-0.43581585016667742322916851390840979)
#define DS16   ((double) 0.44515550282226738565114480689187057)
#define DS17   ((double) 0.52796329829794124367930049003364155)
#define DS18   ((double)-0.49467751640467748805281786601391635)
#define DS19   ((double) 1.03076374706570777841622291802646680)
#define DS20   ((double) 0.51886829082317972874576505841504815)
#define DS21   ((double) 0.82187006169806796277221650036735597)
#define DS22   ((double)-0.79391727157204814808670251421336608)
#define DS23   ((double)-0.75874209018537673384980278323642215)
#define DS24   ((double)-0.19833215631171570317346941647854694)
#define DS25   ((double)-0.43149551708674121082911820080227901)
#define DS26   ((double) 0.27843331019209505198390355841028348)
#define DS27   ((double) 0.59471076351191660168095301828858968)
#define DS28   ((double)-0.57448393456065017248053353113326010)
#define DS29   ((double)-0.54903093419716620563980211536607990)

// EWM:
#define DS0_5th ((double) 0.18559214542766739740398237663061831)
#define B_s3483 ((double)-0.1543847659293762455593100378)
#define D_s3483 ((double) 77.72870185137802409062726308)
#define nDS4    ((double)-0.17619908884183581953599075099065958)
#define B_s5645 ((double) 0.9451188506751776707508747134)
#define D_s5645 ((double) 1.689563031430332839157319596)
#define B_s9019 ((double)-3.527064124888175236105416334)
#define D_s9019 ((double)-1.112194193764683738807670520)
#define B_s5285 ((double) 1.568615537822747404332565815)
#define D_s5285 ((double)-1.561648155256061432564325550)
#define B_s6396 ((double)-3.199961485463981734674359059)
#define D_s6396 ((double)-1.134960866988575539262620301)
#define B_s7407 ((double) 1.358159448857025499908422489)
#define D_s7407 ((double)-1.749196678153277082987222794)
#define C_s4174 ((double)-4.143907256301630677016949990)
#define E_s4174 ((double)-1.080478024630318207542706141)
#define C_s5285 ((double)-1.839920092176650010161002721)
#define E_s5285 ((double)-1.408224849824341570791003234)
#define C_s6396 ((double)-2.725040655738747305153314734)
#define E_s6396 ((double)-1.186102149729676553233822995)

// JVB:
#define c50    ((double)-0.89442719099991587856366946749251043)
#define c51    ((double)-0.27639320225002103035908263312687237)
#define s52    ((double) 0.61803398874989484820458683436563806)
// EWM:
#define c51d50m1 ((double)-0.6909830056250525758977065828)	/* c51/c50 - 1 */

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
	ur1234 = xr0 + DC0*xr38d + DC1*xrsigma;\
	ui1234 = xi0 + DC0*xi38d + DC1*xisigma;\
	ur12 = ur1234 + DC2*xr27c5af;\
	ui12 = ui1234 + DC2*xi27c5af;\
	ur34 = ur1234 + DC7*xr16b49e;\
	ui34 = ui1234 + DC7*xi16b49e;\
	ur1 = ur12 + DC3*xr49e + DC4*xr5af;\
	ui1 = ui12 + DC3*xi49e + DC4*xi5af;\
	ur2 = ur12 + DC5*xr16b + DC6*xr27c;\
	ui2 = ui12 + DC5*xi16b + DC6*xi27c;\
	ur3 = ur34 + DC8*xr49e + DC3*xr5af;\
	ui3 = ui34 + DC8*xi49e + DC3*xi5af;\
	ur4 = ur34 - DC4*xr16b + DC5*xr27c;\
	ui4 = ui34 - DC4*xi16b + DC5*xi27c;\
	ur5 = DC0*ursigma - (ur1 + ur2 + ur3 + ur4) + 5*xr0;\
	ui5 = DC0*uisigma - (ui1 + ui2 + ui3 + ui4) + 5*xi0;\
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
	ur6 = DC9*pr6 + DC10*pr7;\
	ui6 = DC9*pi6 + DC10*pi7;\
	ur7 = DC11*pr6 + DC9*pr7;\
	ui7 = DC11*pi6 + DC9*pi7;\
	pr89 = pr8 + pr9;\
	pi89 = pi8 + pi9;\
	prcd = prc + prd;\
	picd = pic + pid;\
	ur89 = DC12*pr89 + DC15*prcd;\
	ui89 = DC12*pi89 + DC15*picd;\
	ur8 = ur89 + DC13*pr9 + DC16*prd;\
	ui8 = ui89 + DC13*pi9 + DC16*pid;\
	ur9 = ur89 + DC14*pr8 + DC17*prc;\
	ui9 = ui89 + DC14*pi8 + DC17*pic;\
	urcd = DC18*prcd - DC15*pr89;\
	uicd = DC18*picd - DC15*pi89;\
	urc = urcd - DC16*pr9 + DC19*prd;\
	uic = uicd - DC16*pi9 + DC19*pid;\
	urd = urcd - DC17*pr8 + DC20*prc;\
	uid = uicd - DC17*pi8 + DC20*pic;\
	prab = pra + prb;\
	piab = pia + pib;\
	pref = pre + prf;\
	pief = pie + pif;\
	urab = DC21*prab + DC24*pref;\
	uiab = DC21*piab + DC24*pief;\
	ura = urab + DC22*prb + DC25*prf;\
	uia = uiab + DC22*pib + DC25*pif;\
	urb = urab + DC23*pra + DC26*pre;\
	uib = uiab + DC23*pia + DC26*pie;\
	uref = DC27*pref - DC24*prab;\
	uief = DC27*pief - DC24*piab;\
	ure = uref - DC25*prb + DC28*prf;\
	uie = uief - DC25*pib + DC28*pif;\
	urf = uref - DC26*pra + DC29*pre;\
	uif = uief - DC26*pia + DC29*pie;\
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
	vr1234 = DS0*yr38d + DS1*yrsigma;\
	vi1234 = DS0*yi38d + DS1*yisigma;\
	vr12 = vr1234 + DS2*yr27c5af;\
	vi12 = vi1234 + DS2*yi27c5af;\
	vr34 = vr1234 + DS7*yr16b49e;\
	vi34 = vi1234 + DS7*yi16b49e;\
	vr1 = vr12 + DS3*yr49e + DS4*yr5af;\
	vi1 = vi12 + DS3*yi49e + DS4*yi5af;\
	vr2 = vr12 + DS5*yr16b + DS6*yr27c;\
	vi2 = vi12 + DS5*yi16b + DS6*yi27c;\
	vr3 = vr34 + DS8*yr49e + DS3*yr5af;\
	vi3 = vi34 + DS8*yi49e + DS3*yi5af;\
	vr4 = vr34 - DS4*yr16b + DS5*yr27c;\
	vi4 = vi34 - DS4*yi16b + DS5*yi27c;\
	vr5 = DS0*vrsigma - (vr1 + vr2 + vr3 + vr4);\
	vi5 = DS0*visigma - (vi1 + vi2 + vi3 + vi4);\
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
	vr6 = DS9*qr6 + DS10*qr7;\
	vi6 = DS9*qi6 + DS10*qi7;\
	vr7 = DS11*qr6 + DS9*qr7;\
	vi7 = DS11*qi6 + DS9*qi7;\
	qr89 = qr8 + qr9;\
	qi89 = qi8 + qi9;\
	qrcd = qrc + qrd;\
	qicd = qic + qid;\
	vr89 = DS12*qr89 + DS15*qrcd;\
	vi89 = DS12*qi89 + DS15*qicd;\
	vr8 = vr89 + DS13*qr9 + DS16*qrd;\
	vi8 = vi89 + DS13*qi9 + DS16*qid;\
	vr9 = vr89 + DS14*qr8 + DS17*qrc;\
	vi9 = vi89 + DS14*qi8 + DS17*qic;\
	vrcd = DS18*qrcd - DS15*qr89;\
	vicd = DS18*qicd - DS15*qi89;\
	vrc = vrcd - DS16*qr9 + DS19*qrd;\
	vic = vicd - DS16*qi9 + DS19*qid;\
	vrd = vrcd - DS17*qr8 + DS20*qrc;\
	vid = vicd - DS17*qi8 + DS20*qic;\
	qrab = qra + qrb;\
	qiab = qia + qib;\
	qref = qre + qrf;\
	qief = qie + qif;\
	vrab = DS21*qrab + DS24*qref;\
	viab = DS21*qiab + DS24*qief;\
	vra = vrab + DS22*qrb + DS25*qrf;\
	via = viab + DS22*qib + DS25*qif;\
	vrb = vrab + DS23*qra + DS26*qre;\
	vib = viab + DS23*qia + DS26*qie;\
	vref = DS27*qref - DS24*qrab;\
	vief = DS27*qief - DS24*qiab;\
	vre = vref - DS25*qrb + DS28*qrf;\
	vie = vief - DS25*qib + DS28*qif;\
	vrf = vref - DS26*qra + DS29*qre;\
	vif = vief - DS26*qia + DS29*qie;\
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
	u1234 = x0 + DC0*x38d + DC1*xsigma;\
	u12 = u1234 + DC2*x27c5af;\
	u34 = u1234 + DC7*x16b49e;\
	u1 = u12 + DC3*x49e + DC4*x5af;\
	u2 = u12 + DC5*x16b + DC6*x27c;\
	u3 = u34 + DC8*x49e + DC3*x5af;\
	u4 = u34 - DC4*x16b + DC5*x27c;\
	u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0;\
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
	u6 = DC9*p6 + DC10*p7;\
	u7 = DC11*p6 + DC9*p7;\
	p89 = p8 + p9;\
	pcd = pc + pd;\
	u89 = DC12*p89 + DC15*pcd;\
	u8 = u89 + DC13*p9 + DC16*pd;\
	u9 = u89 + DC14*p8 + DC17*pc;\
	ucd = DC18*pcd - DC15*p89;\
	uc = ucd - DC16*p9 + DC19*pd;\
	ud = ucd - DC17*p8 + DC20*pc;\
	pab = pa + pb;\
	pef = pe + pf;\
	uab = DC21*pab + DC24*pef;\
	ua = uab + DC22*pb + DC25*pf;\
	ub = uab + DC23*pa + DC26*pe;\
	uef = DC27*pef - DC24*pab;\
	ue = uef - DC25*pb + DC28*pf;\
	uf = uef - DC26*pa + DC29*pe;\
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
	v1234 = DS0*y38d + DS1*ysigma;\
	v12 = v1234 + DS2*y27c5af;\
	v34 = v1234 + DS7*y16b49e;\
	v1 = v12 + DS3*y49e + DS4*y5af;\
	v2 = v12 + DS5*y16b + DS6*y27c;\
	v3 = v34 + DS8*y49e + DS3*y5af;\
	v4 = v34 - DS4*y16b + DS5*y27c;\
	v5 = DS0*vsigma - (v1 + v2 + v3 + v4);\
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
	v6 = DS9*q6 + DS10*q7;\
	v7 = DS11*q6 + DS9*q7;\
	q89 = q8 + q9;\
	qcd = qc + qd;\
	v89 = DS12*q89 + DS15*qcd;\
	v8 = v89 + DS13*q9 + DS16*qd;\
	v9 = v89 + DS14*q8 + DS17*qc;\
	vcd = DS18*qcd - DS15*q89;\
	vc = vcd - DS16*q9 + DS19*qd;\
	vd = vcd - DS17*q8 + DS20*qc;\
	qab = qa + qb;\
	qef = qe + qf;\
	vab = DS21*qab + DS24*qef;\
	va = vab + DS22*qb + DS25*qf;\
	vb = vab + DS23*qa + DS26*qe;\
	vef = DS27*qef - DS24*qab;\
	ve = vef - DS25*qb + DS28*qf;\
	vf = vef - DS26*qa + DS29*qe;\
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
	u1234 = x0 + DC0*x38d + DC1*xsigma;\
	u12 = u1234 + DC2*x27c5af;\
	u34 = u1234 + DC7*x16b49e;\
	u1 = u12 + DC3*x49e + DC4*x5af;\
	u2 = u12 + DC5*x16b + DC6*x27c;\
	u3 = u34 + DC8*x49e + DC3*x5af;\
	u4 = u34 - DC4*x16b + DC5*x27c;\
	u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0;\
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
	u6 = DC9*p6 + DC10*p7;\
	u7 = DC11*p6 + DC9*p7;\
	p89 = p8 + p9;\
	pcd = pc + pd;\
	u89 = DC12*p89 + DC15*pcd;\
	u8 = u89 + DC13*p9 + DC16*pd;\
	u9 = u89 + DC14*p8 + DC17*pc;\
	ucd = DC18*pcd - DC15*p89;\
	uc = ucd - DC16*p9 + DC19*pd;\
	ud = ucd - DC17*p8 + DC20*pc;\
	pab = pa + pb;\
	pef = pe + pf;\
	uab = DC21*pab + DC24*pef;\
	ua = uab + DC22*pb + DC25*pf;\
	ub = uab + DC23*pa + DC26*pe;\
	uef = DC27*pef - DC24*pab;\
	ue = uef - DC25*pb + DC28*pf;\
	uf = uef - DC26*pa + DC29*pe;\
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
	v1234 = DS0*y38d + DS1*ysigma;\
	v12 = v1234 + DS2*y27c5af;\
	v34 = v1234 + DS7*y16b49e;\
	v1 = v12 + DS3*y49e + DS4*y5af;\
	v2 = v12 + DS5*y16b + DS6*y27c;\
	v3 = v34 + DS8*y49e + DS3*y5af;\
	v4 = v34 - DS4*y16b + DS5*y27c;\
	v5 = DS0*vsigma - (v1 + v2 + v3 + v4);\
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
	v6 = DS9*q6 + DS10*q7;\
	v7 = DS11*q6 + DS9*q7;\
	q89 = q8 + q9;\
	qcd = qc + qd;\
	v89 = DS12*q89 + DS15*qcd;\
	v8 = v89 + DS13*q9 + DS16*qd;\
	v9 = v89 + DS14*q8 + DS17*qc;\
	vcd = DS18*qcd - DS15*q89;\
	vc = vcd - DS16*q9 + DS19*qd;\
	vd = vcd - DS17*q8 + DS20*qc;\
	qab = qa + qb;\
	qef = qe + qf;\
	vab = DS21*qab + DS24*qef;\
	va = vab + DS22*qb + DS25*qf;\
	vb = vab + DS23*qa + DS26*qe;\
	vef = DS27*qef - DS24*qab;\
	ve = vef - DS25*qb + DS28*qf;\
	vf = vef - DS26*qa + DS29*qe;\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += __Br00 ;		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += __Br00;		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4);\
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
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4);\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += __Bi00 ;		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += __Bi00;		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4);\
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
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5);\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2);\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3);\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0);\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9);\
	xc += x8;\
	x9 += x0;\
	MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1);\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0);\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6);\
	xf += x8;\
	x6 += x0;\
	MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4);\
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
__A,__i00,__i01,__i02,__i03,__i04,__i05,__i06,__i07,__i08,__i09,__i10,__i11,__i12,__i13,__i14,__i15,__i16,__i17,__i18,__i19,__i20,__i21,__i22,__i23,__i24,__i25,__i26,__i27,__i28,__i29,__i30,/* Inputs: Base address plus 31 offsets */\
__B,__o00,__o01,__o02,__o03,__o04,__o05,__o06,__o07,__o08,__o09,__o10,__o11,__o12,__o13,__o14,__o15,__o16,__o17,__o18,__o19,__o20,__o21,__o22,__o23,__o24,__o25,__o26,__o27,__o28,__o29,__o30 /* Outputs: Base address plus 31 offsets */\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(__B+__o00);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(__B+__o00);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(__B+__o00+RE_IM_STRIDE);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(__B+__o00+RE_IM_STRIDE);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
__A,__idx, /* Inputs: Base address plus 31 (index) offsets */\
__out /* Outputs: Base address of a contiguous memblock corr. to 31 x complex = 62 x real data */\
)\
{\
	double TMP01,TMP02;	/* 2 double-tmps to mimic a pair of common spill locations */\
	/* 16 double-tmps to mimic 16 hardware-floating registers: */\
	double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf;\
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
	x5  = *(__A+__idx[0x0e]);\
	x6  = *(__A+__idx[0x05]);\
	x7  = *(__A+__idx[0x04]);\
	x8  = *(__A+__idx[0x1c]);\
	x9  = *(__A+__idx[0x0a]);\
	xa  = *(__A+__idx[0x08]);\
	xb  = *(__A+__idx[0x19]);\
	xc  = *(__A+__idx[0x14]);\
	xd  = *(__A+__idx[0x10]);\
	xe  = *(__A+__idx[0x13]);\
	xf  = *(__A+__idx[0x09]);		/* xr-terms: */\
	x1 += *(__A+__idx[0x1e]);		/* x1 = __Ar01 + __Ar30 */\
	x2 += *(__A+__idx[0x18]);		/* x2 = __Ar07 + __Ar24 */\
	x3 += *(__A+__idx[0x0d]);		/* x3 = __Ar18 + __Ar13 */\
	x4 += *(__A+__idx[0x1d]);		/* x4 = __Ar02 + __Ar29 */\
	x5 += *(__A+__idx[0x11]);		/* x5 = __Ar14 + __Ar17 */\
	x6 += *(__A+__idx[0x1a]);		/* x6 = __Ar05 + __Ar26 */\
	x7 += *(__A+__idx[0x1b]);		/* x7 = __Ar04 + __Ar27 */\
	x8 += *(__A+__idx[0x03]);		/* x8 = __Ar28 + __Ar03 */\
	x9 += *(__A+__idx[0x15]);		/* x9 = __Ar10 + __Ar21 */\
	xa += *(__A+__idx[0x17]);		/* xa = __Ar08 + __Ar23 */\
	xb += *(__A+__idx[0x06]);		/* xb = __Ar25 + __Ar06 */\
	xc += *(__A+__idx[0x0b]);		/* xc = __Ar20 + __Ar11 */\
	xd += *(__A+__idx[0x0f]);		/* xd = __Ar16 + __Ar15 */\
	xe += *(__A+__idx[0x0c]);		/* xe = __Ar19 + __Ar12 */\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(__out);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(__out);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(Aim+__idx[0x01]);\
	x2  = *(Aim+__idx[0x07]);\
	x3  = *(Aim+__idx[0x12]);\
	x4  = *(Aim+__idx[0x02]);\
	x5  = *(Aim+__idx[0x0e]);\
	x6  = *(Aim+__idx[0x05]);\
	x7  = *(Aim+__idx[0x04]);\
	x8  = *(Aim+__idx[0x1c]);\
	x9  = *(Aim+__idx[0x0a]);\
	xa  = *(Aim+__idx[0x08]);\
	xb  = *(Aim+__idx[0x19]);\
	xc  = *(Aim+__idx[0x14]);\
	xd  = *(Aim+__idx[0x10]);\
	xe  = *(Aim+__idx[0x13]);\
	xf  = *(Aim+__idx[0x09]);\
	x1 -= *(Aim+__idx[0x1e]);\
	x2 -= *(Aim+__idx[0x18]);\
	x3 -= *(Aim+__idx[0x0d]);\
	x4 -= *(Aim+__idx[0x1d]);\
	x5 -= *(Aim+__idx[0x11]);\
	x6 -= *(Aim+__idx[0x1a]);\
	x7 -= *(Aim+__idx[0x1b]);\
	x8 -= *(Aim+__idx[0x03]);\
	x9 -= *(Aim+__idx[0x15]);\
	xa -= *(Aim+__idx[0x17]);\
	xb -= *(Aim+__idx[0x06]);\
	xc -= *(Aim+__idx[0x0b]);\
	xd -= *(Aim+__idx[0x0f]);\
	xe -= *(Aim+__idx[0x0c]);\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
	x5  = *(Aim+__idx[0x0e]);\
	x6  = *(Aim+__idx[0x05]);\
	x7  = *(Aim+__idx[0x04]);\
	x8  = *(Aim+__idx[0x1c]);\
	x9  = *(Aim+__idx[0x0a]);\
	xa  = *(Aim+__idx[0x08]);\
	xb  = *(Aim+__idx[0x19]);\
	xc  = *(Aim+__idx[0x14]);\
	xd  = *(Aim+__idx[0x10]);\
	xe  = *(Aim+__idx[0x13]);\
	xf  = *(Aim+__idx[0x09]);		/* xi-terms: */\
	x1 += *(Aim+__idx[0x1e]);		/* x1 = __Ai01 + __Ai30 */\
	x2 += *(Aim+__idx[0x18]);		/* x2 = __Ai07 + __Ai24 */\
	x3 += *(Aim+__idx[0x0d]);		/* x3 = __Ai18 + __Ai13 */\
	x4 += *(Aim+__idx[0x1d]);		/* x4 = __Ai02 + __Ai29 */\
	x5 += *(Aim+__idx[0x11]);		/* x5 = __Ai14 + __Ai17 */\
	x6 += *(Aim+__idx[0x1a]);		/* x6 = __Ai05 + __Ai26 */\
	x7 += *(Aim+__idx[0x1b]);		/* x7 = __Ai04 + __Ai27 */\
	x8 += *(Aim+__idx[0x03]);		/* x8 = __Ai28 + __Ai03 */\
	x9 += *(Aim+__idx[0x15]);		/* x9 = __Ai10 + __Ai21 */\
	xa += *(Aim+__idx[0x17]);		/* xa = __Ai08 + __Ai23 */\
	xb += *(Aim+__idx[0x06]);		/* xb = __Ai25 + __Ai06 */\
	xc += *(Aim+__idx[0x0b]);		/* xc = __Ai20 + __Ai11 */\
	xd += *(Aim+__idx[0x0f]);		/* xd = __Ai16 + __Ai15 */\
	xe += *(Aim+__idx[0x0c]);		/* xe = __Ai19 + __Ai12 */\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(__out+1);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(__out+1);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
	/****************************************************************************************************************/\
	x1  = *(__A+__idx[0x01]);\
	x2  = *(__A+__idx[0x07]);\
	x3  = *(__A+__idx[0x12]);\
	x4  = *(__A+__idx[0x02]);\
	x5  = *(__A+__idx[0x0e]);\
	x6  = *(__A+__idx[0x05]);\
	x7  = *(__A+__idx[0x04]);\
	x8  = *(__A+__idx[0x1c]);\
	x9  = *(__A+__idx[0x0a]);\
	xa  = *(__A+__idx[0x08]);\
	xb  = *(__A+__idx[0x19]);\
	xc  = *(__A+__idx[0x14]);\
	xd  = *(__A+__idx[0x10]);\
	xe  = *(__A+__idx[0x13]);\
	xf  = *(__A+__idx[0x09]);\
	x1 -= *(__A+__idx[0x1e]);\
	x2 -= *(__A+__idx[0x18]);\
	x3 -= *(__A+__idx[0x0d]);\
	x4 -= *(__A+__idx[0x1d]);\
	x5 -= *(__A+__idx[0x11]);\
	x6 -= *(__A+__idx[0x1a]);\
	x7 -= *(__A+__idx[0x1b]);\
	x8 -= *(__A+__idx[0x03]);\
	x9 -= *(__A+__idx[0x15]);\
	xa -= *(__A+__idx[0x17]);\
	xb -= *(__A+__idx[0x06]);\
	xc -= *(__A+__idx[0x0b]);\
	xd -= *(__A+__idx[0x0f]);\
	xe -= *(__A+__idx[0x0c]);\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
	__in, /* Inputs: Base address of a contiguous memblock corr. to 31 x complex = 62 x real data */\
	__B,__odx	/* Outputs: Base address plus 31 (index) offsets */\
)\
{\
	double TMP01,TMP02;	/* 2 double-tmps to mimic a pair of common spill locations */\
	/* 16 double-tmps to mimic 16 hardware-floating registers: */\
	double x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf;\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(__B+__odx[0x0]);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(__B+__odx[0x0]);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	*(__B+__odx[0x0e]) = x4;	/* sb = xe + x4 - xc */\
	*(__B+__odx[0x09]) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(__B+__odx[0x02]) = x5;	/* sc = x5 - xd */\
	*(__B+__odx[0x0c]) = xd;	/* s2 = x5 + xd - x7 */\
	*(__B+__odx[0x0a]) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(__B+__odx[0x0f]) = x8;	/* s3 = x8 - x3 */\
	*(__B+__odx[0x03]) = x3;	/* s8 = x8 + x3 - x9 */\
	*(__B+__odx[0x0d]) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(__B+__odx[0x04]) = xb;	/* s9 = xb - x1 */\
	*(__B+__odx[0x07]) = x1;	/* se = xb + x1 - xf */\
	*(__B+__odx[0x0b]) = xf;	/* s4 = xb + xf */\
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
	/* Part 2: Compute the yi-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
 	*(__B+__odx[0x1e]) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(__B+__odx[0x02]);		x0 = *(__B+__odx[0x03]);\
 	x2 -= x5;	 					x0 += x3;	\
 	*(__B+__odx[0x02]) = x2;		*(__B+__odx[0x03]) = x0;\
 	x5 *= 2;	 					x3 *= 2;	\
 	x2 += x5;	 					x0 -= x3;	\
 	*(__B+__odx[0x1d]) = x2;		*(__B+__odx[0x1c]) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers to reduce serialization (but stop at 4, more likely no help): */\
 	x0 = *(__B+__odx[0x04]);		x2 = *(__B+__odx[0x05]);	x3 = *(__B+__odx[0x06]);	x5 = *(__B+__odx[0x07]);\
	x0 -= xb;						x2 -= x6;					x3 += xa;					x5 -= x1;	\
	*(__B+__odx[0x04]) = x0;		*(__B+__odx[0x05]) = x2;	*(__B+__odx[0x06]) = x3;	*(__B+__odx[0x07]) = x5;\
	xb *= 2;						x6 *= 2;					xa *= 2;					x1 *= 2;	\
	x0 += xb;						x2 += x6;					x3 -= xa;					x5 += x1;	\
	*(__B+__odx[0x1b]) = x0;		*(__B+__odx[0x1a]) = x2;	*(__B+__odx[0x19]) = x3;	*(__B+__odx[0x18]) = x5;\
\
 	x0 = *(__B+__odx[0x08]);		x2 = *(__B+__odx[0x09]);	x3 = *(__B+__odx[0x0a]);	x5 = *(__B+__odx[0x0b]);\
	x0 -= xe;						x2 -= xc;					x3 -= x7;					x5 += xf;	\
	*(__B+__odx[0x08]) = x0;		*(__B+__odx[0x09]) = x2;	*(__B+__odx[0x0a]) = x3;	*(__B+__odx[0x0b]) = x5;\
	xe *= 2;						xc *= 2;					x7 *= 2;					xf *= 2;	\
	x0 += xe;						x2 += xc;					x3 += x7;					x5 -= xf;	\
	*(__B+__odx[0x17]) = x0;		*(__B+__odx[0x16]) = x2;	*(__B+__odx[0x15]) = x3;	*(__B+__odx[0x14]) = x5;\
\
 	x0 = *(__B+__odx[0x0c]);		x2 = *(__B+__odx[0x0d]);	x3 = *(__B+__odx[0x0e]);	x5 = *(__B+__odx[0x0f]);\
	x0 += xd;						x2 += x9;					x3 -= x4;					x5 += x8;	\
	*(__B+__odx[0x0c]) = x0;		*(__B+__odx[0x0d]) = x2;	*(__B+__odx[0x0e]) = x3;	*(__B+__odx[0x0f]) = x5;\
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
	x8 *= DC0_5th;		/* DC0*x38d */\
	x0 *= DC1;			/* DC1*xsigma */\
	x8 += *(Bim+__odx[0x0]);		/* x0 + DC0*x38d */\
	x8 += x0;			/* u1234 = x0 + DC0*x38d + DC1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DC2;			/* DC2*x27c5af */\
	x8 *= DC7;			/* DC7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DC2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DC7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DC3,B_c3483,DC8,D_c3483, xe,x5) */\
	x5 *= B_c3483;\
	xe += x5;\
	x5 *= D_c3483;\
	x5 += xe;\
	xe *= DC3;\
	x5 *= DC8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DC5,B_c5645,-DC4,D_c5645, xb,x2) */\
	x2 *= B_c5645;\
	xb += x2;\
	x2 *= D_c5645;\
	x2 += xb;\
	xb *= DC5;\
	x2 *= nDC4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DC0_5th;		/* DC0*usigma/5 */\
	x0 += *(Bim+__odx[0x0]);		/* DC0*usigma/5 + x0 */\
	x0 *= 5;			/* DC0*usigma + 5*x0 */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DC0*usigma - (u1 + u2 + u3 + u4) + 5*x0 */\
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
	/* MAT_2x2_MUL(DC9,B_c9019,DC11,D_c9019, x7,x3) */\
	x3 *= B_c9019;\
	x7 += x3;\
	x3 *= D_c9019;\
	x3 += x7;\
	x7 *= DC9;\
	x3 *= DC11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DC15,B_c5285,DC18,D_c5285, xd,x0) */\
	x0 *= B_c5285;\
	xd += x0;\
	x0 *= D_c5285;\
	x0 += xd;\
	xd *= DC15;\
	x0 *= DC18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DC16,B_c6396,DC19,D_c6396, xc,x9) */\
	x9 *= B_c6396;\
	xc += x9;\
	x9 *= D_c6396;\
	x9 += xc;\
	xc *= DC16;\
	x9 *= DC19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DC17,B_c7407,DC20,D_c7407, xd,x1) */\
	x1 *= B_c7407;\
	xd += x1;\
	x1 *= D_c7407;\
	x1 += xd;\
	xd *= DC17;\
	x1 *= DC20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DC24,C_c4174,DC27,E_c4174, xa,x0) */\
	x0 *= C_c4174;\
	xa += x0;\
	x0 *= E_c4174;\
	x0 += xa;\
	xa *= DC24;\
	x0 *= DC27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DC25,C_c5285,DC28,E_c5285, xf,x6) */\
	x6 *= C_c5285;\
	xf += x6;\
	x6 *= E_c5285;\
	x6 += xf;\
	xf *= DC25;\
	x6 *= DC28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DC26,C_c6396,DC29,E_c6396, xa,x4) */\
	x4 *= C_c6396;\
	xa += x4;\
	x4 *= E_c6396;\
	x4 += xa;\
	xa *= DC26;\
	x4 *= DC29;\
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
	*(Bim+__odx[0x0e]) = x4;	/* sb = xe + x4 - xc */\
	*(Bim+__odx[0x09]) = xc;	/* s1 = xe + xc */\
/* abc -> 57d: */\
	x0  = x5;\
	x5 -= xd;	/* x5 - xd */\
	xd -= x7;\
	x7 += x0;	/* x5 + x7 */\
	xd += x0;	/* x5 - x7 + xd */\
\
	*(Bim+__odx[0x02]) = x5;	/* sc = x5 - xd */\
	*(Bim+__odx[0x0c]) = xd;	/* s2 = x5 + xd - x7 */\
	*(Bim+__odx[0x0a]) = x7;	/* s7 = x5 + x7 */\
/* abc -> 893: */\
	x0  = x8;\
	x8 -= x3;	/* x8 - x3 */\
	x3 -= x9;\
	x9 += x0;	/* x8 + x9 */\
	x3 += x0;	/* x8 - x9 + x3 */\
\
	*(Bim+__odx[0x0f]) = x8;	/* s3 = x8 - x3 */\
	*(Bim+__odx[0x03]) = x3;	/* s8 = x8 + x3 - x9 */\
	*(Bim+__odx[0x0d]) = x9;	/* sd = x8 + x9 */\
/* abc -> bf1: */\
	x0  = xb;\
	xb -= x1;	/* xb - x1 */\
	x1 -= xf;\
	xf += x0;	/* xb + xf */\
	x1 += x0;	/* xb - xf + x1 */\
\
	*(Bim+__odx[0x04]) = xb;	/* s9 = xb - x1 */\
	*(Bim+__odx[0x07]) = x1;	/* se = xb + x1 - xf */\
	*(Bim+__odx[0x0b]) = xf;	/* s4 = xb + xf */\
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
	/* Part 2: Compute the yr-terms. Same code sequence as for xr-terms, but all DC-consts --> DS, no DC-component: */\
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
	x8 *= DS0_5th;		/* DS0*x38d */\
	x0 *= DS1;			/* DS1*xsigma */\
	x8 += x0;			/* u1234 = DS0*x38d + DS1*xsigma */\
	TMP02 = x8;		/* spill u1234 */\
	x0 = x2;			/* x27c5af */\
	x8 = xb;			/* x16b49e */\
	x0 *= DS2;			/* DS2*x27c5af */\
	x8 *= DS7;			/* DS7*x16b49e */\
	x0 += TMP02;		/* u12 = u1234 + DS2*x27c5af */\
	x8 += TMP02;		/* u34 = u1234 + DS7*x16b49e */\
	xb -= xe;			/* x16b49e - x49e */\
	x2 -= x5;			/* x27c5af - x5af */\
/******** use 2x2 const-matrix x vector in-place mul algo twice: *********/\
	/* MAT_2x2_MUL(DS3,B_s3483,DS8,D_s3483, xe,x5) */\
	x5 *= B_s3483;\
	xe += x5;\
	x5 *= D_s3483;\
	x5 += xe;\
	xe *= DS3;\
	x5 *= DS8;\
\
	xe += x0;			/* u1 */\
	x5 += x8;			/* u3 */\
	/* MAT_2x2_MUL(DS5,B_s5645,-DS4,D_s5645, xb,x2) */\
	x2 *= B_s5645;\
	xb += x2;\
	x2 *= D_s5645;\
	x2 += xb;\
	xb *= DS5;\
	x2 *= nDS4;\
\
	xb += x0;			/* u2 */\
	x2 += x8;			/* u4 */\
	x0 = TMP01;		/* usigma */\
	x0 *= DS0;			/* DS0*usigma */\
	x0 -= xe;\
	x0 -= xb;\
	x0 -= x5;\
	x0 -= x2;\
	TMP02 = x0;		/* u5 = DS0*usigma - (u1 + u2 + u3 + u4) */\
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
	/* MAT_2x2_MUL(DS9,B_s9019,DS11,D_s9019, x7,x3) */\
	x3 *= B_s9019;\
	x7 += x3;\
	x3 *= D_s9019;\
	x3 += x7;\
	x7 *= DS9;\
	x3 *= DS11;\
\
	/* Store copies of x1,xd: */\
	x0 = x1;\
	TMP01 = xd;\
	x0 -= x9;	/* p89 */\
	xd += xc;	/* pcd */\
	/* MAT_2x2_MUL(DS15,B_s5285,DS18,D_s5285, xd,x0) */\
	x0 *= B_s5285;\
	xd += x0;\
	x0 *= D_s5285;\
	x0 += xd;\
	xd *= DS15;\
	x0 *= DS18;\
\
	x8 = xd;	/* u89; ucd in x0 */\
	xd = TMP01;\
	/* MAT_2x2_MUL(DS16,B_s6396,DS19,D_s6396, xc,x9) */\
	x9 *= B_s6396;\
	xc += x9;\
	x9 *= D_s6396;\
	x9 += xc;\
	xc *= DS16;\
	x9 *= DS19;\
\
	xc += x8;\
	x9 += x0;\
	/* MAT_2x2_MUL(DS17,B_s7407,DS20,D_s7407, xd,x1) */\
	x1 *= B_s7407;\
	xd += x1;\
	x1 *= D_s7407;\
	x1 += xd;\
	xd *= DS17;\
	x1 *= DS20;\
\
	xd += x8;\
	x1 += x0;\
	/* Same code sequence here, but all the DC-indices += 9 and replace x1,9,d,c --> x4,6,a,f: */\
	x0 = x4;\
	TMP01 = xa;\
	x0 -= x6;	/* pab */\
	xa += xf;	/* pef */\
	/* MAT_2x2_MUL(DS24,C_s4174,DS27,E_s4174, xa,x0) */\
	x0 *= C_s4174;\
	xa += x0;\
	x0 *= E_s4174;\
	x0 += xa;\
	xa *= DS24;\
	x0 *= DS27;\
\
	x8 = xa;	/* uab; uef in x0 */\
	xa = TMP01;\
	/* MAT_2x2_MUL(DS25,C_s5285,DS28,E_s5285, xf,x6) */\
	x6 *= C_s5285;\
	xf += x6;\
	x6 *= E_s5285;\
	x6 += xf;\
	xf *= DS25;\
	x6 *= DS28;\
\
	xf += x8;\
	x6 += x0;\
	/* MAT_2x2_MUL(DS26,C_s6396,DS29,E_s6396, xa,x4) */\
	x4 *= C_s6396;\
	xa += x4;\
	x4 *= E_s6396;\
	x4 += xa;\
	xa *= DS26;\
	x4 *= DS29;\
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
 	*(Bim+__odx[0x1e]) = x0;	/* x0,x2 */\
/* Now take advantage of 2 free registers to reduce serialization: */\
 	x2 = *(Bim+__odx[0x02]);	x0 = *(Bim+__odx[0x03]);\
 	x2 += x5;					x0 -= x3;	\
 	*(Bim+__odx[0x02]) = x2;	*(Bim+__odx[0x03]) = x0;\
 	x5 *= 2;	 				x3 *= 2;	\
 	x2 -= x5;	 				x0 += x3;	\
 	*(Bim+__odx[0x1d]) = x2;	*(Bim+__odx[0x1c]) = x0;	/* x0,2,3,5 */\
/* Now take advantage of 4 free registers toreduce serialization, but stop at 4, more likely no help: */\
 	x0 = *(Bim+__odx[0x04]);	x2 = *(Bim+__odx[0x05]);	x3 = *(Bim+__odx[0x06]);	x5 = *(Bim+__odx[0x07]);\
	x0 += xb;					x2 += x6;					x3 -= xa;					x5 += x1;	\
	*(Bim+__odx[0x04]) = x0;	*(Bim+__odx[0x05]) = x2;	*(Bim+__odx[0x06]) = x3;	*(Bim+__odx[0x07]) = x5;\
	xb *= 2;					x6 *= 2;					xa *= 2;					x1 *= 2;	\
	x0 -= xb;					x2 -= x6;					x3 += xa;					x5 -= x1;	\
	*(Bim+__odx[0x1b]) = x0;	*(Bim+__odx[0x1a]) = x2;	*(Bim+__odx[0x19]) = x3;	*(Bim+__odx[0x18]) = x5;\
\
 	x0 = *(Bim+__odx[0x08]);	x2 = *(Bim+__odx[0x09]);	x3 = *(Bim+__odx[0x0a]);	x5 = *(Bim+__odx[0x0b]);\
	x0 += xe;					x2 += xc;					x3 += x7;					x5 -= xf;	\
	*(Bim+__odx[0x08]) = x0;	*(Bim+__odx[0x09]) = x2;	*(Bim+__odx[0x0a]) = x3;	*(Bim+__odx[0x0b]) = x5;\
	xe *= 2;					xc *= 2;					x7 *= 2;					xf *= 2;	\
	x0 -= xe;					x2 -= xc;					x3 -= x7;					x5 += xf;	\
	*(Bim+__odx[0x17]) = x0;	*(Bim+__odx[0x16]) = x2;	*(Bim+__odx[0x15]) = x3;	*(Bim+__odx[0x14]) = x5;\
\
 	x0 = *(Bim+__odx[0x0c]);	x2 = *(Bim+__odx[0x0d]);	x3 = *(Bim+__odx[0x0e]);	x5 = *(Bim+__odx[0x0f]);\
	x0 -= xd;					x2 -= x9;					x3 += x4;					x5 -= x8;	\
	*(Bim+__odx[0x0c]) = x0;	*(Bim+__odx[0x0d]) = x2;	*(Bim+__odx[0x0e]) = x3;	*(Bim+__odx[0x0f]) = x5;\
	xd *= 2;					x9 *= 2;					x4 *= 2;					x8 *= 2;	\
	x0 += xd;					x2 += x9;					x3 -= x4;					x5 += x8;	\
	*(Bim+__odx[0x13]) = x0;	*(Bim+__odx[0x12]) = x2;	*(Bim+__odx[0x11]) = x3;	*(Bim+__odx[0x10]) = x5;\
\
/* Totals: 658 FADD, 234 FMUL. */\
}

#endif	/* #ifndef radix31_included */
