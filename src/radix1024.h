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
#ifndef radix1024_included
#define radix1024_included

#include "radix512.h"

	// 'bc -l' code for these: p2=8*a(1);d=p2/1024;t=-d;	t+=(d+d);c(t);s(t); [repeat 64 times]:
	// Of the odd-order 1024th roots, note that _4f,_53,_7f end up being unused by the radix-1024 DFT twiddles array:
	#define c1024_01 ((double)0.99998117528260114265)
	#define s1024_01 ((double)0.00613588464915447535)	/* exp(01*I*twopi/1024) */
	#define c1024_03 ((double)0.99983058179582342201)
	#define s1024_03 ((double)0.01840672990580482090)	/* exp(03*I*twopi/1024) */		
	#define c1024_05 ((double)0.99952941750109316308)
	#define s1024_05 ((double)0.03067480317663662588)	/* exp(05*I*twopi/1024) */
	#define c1024_07 ((double)0.99907772775264538289)
	#define s1024_07 ((double)0.04293825693494082301)	/* exp(07*I*twopi/1024) */		
	#define c1024_09 ((double)0.99847558057329475221)
	#define s1024_09 ((double)0.05519524434968993972)	/* exp(09*I*twopi/1024) */
	#define c1024_0b ((double)0.99772306664419160985)
	#define s1024_0b ((double)0.06744391956366405780)	/* exp(0b*I*twopi/1024) */		
	#define c1024_0d ((double)0.99682029929116571498)
	#define s1024_0d ((double)0.07968243797143012103)	/* exp(0d*I*twopi/1024) */
	#define c1024_0f ((double)0.99576741446765979399)
	#define s1024_0f ((double)0.09190895649713272849)	/* exp(0f*I*twopi/1024) */		
	#define c1024_11 ((double)0.99456457073425545213)
	#define s1024_11 ((double)0.10412163387205457897)	/* exp(11*I*twopi/1024) */
	#define c1024_13 ((double)0.99321194923479453312)
	#define s1024_13 ((double)0.11631863091190476708)	/* exp(13*I*twopi/1024) */		
	#define c1024_15 ((double)0.99170975366909952288)
	#define s1024_15 ((double)0.12849811079379317243)	/* exp(15*I*twopi/1024) */
	#define c1024_17 ((double)0.99005821026229710553)
	#define s1024_17 ((double)0.14065823933284923051)	/* exp(17*I*twopi/1024) */		
	#define c1024_19 ((double)0.98825756773074949143)
	#define s1024_19 ((double)0.15279718525844342750)	/* exp(19*I*twopi/1024) */
	#define c1024_1b ((double)0.98630809724459864790)
	#define s1024_1b ((double)0.16491312048996992118)	/* exp(1b*I*twopi/1024) */		
	#define c1024_1d ((double)0.98421009238692907323)
	#define s1024_1d ((double)0.17700422041214875594)	/* exp(1d*I*twopi/1024) */
	#define c1024_1f ((double)0.98196386910955526412)
	#define s1024_1f ((double)0.18906866414980621248)	/* exp(1f*I*twopi/1024) */		
	#define c1024_21 ((double)0.97956976568544053449)
	#define s1024_21 ((double)0.20110463484209191127)	/* exp(21*I*twopi/1024) */
	#define c1024_23 ((double)0.97702814265775435155)
	#define s1024_23 ((double)0.21311031991609137366)	/* exp(23*I*twopi/1024) */		
	#define c1024_25 ((double)0.97433938278557586059)
	#define s1024_25 ((double)0.22508391135979283567)	/* exp(25*I*twopi/1024) */
	#define c1024_27 ((double)0.97150389098625177561)
	#define s1024_27 ((double)0.23702360599436720653)	/* exp(27*I*twopi/1024) */		
	#define c1024_29 ((double)0.96852209427441731631)
	#define s1024_29 ((double)0.24892760574572016775)	/* exp(29*I*twopi/1024) */
	#define c1024_2b ((double)0.96539444169768937465)
	#define s1024_2b ((double)0.26079411791527551791)	/* exp(2b*I*twopi/1024) */		
	#define c1024_2d ((double)0.96212140426904159553)
	#define s1024_2d ((double)0.27262135544994898410)	/* exp(2d*I*twopi/1024) */
	#define c1024_2f ((double)0.95870347489587155549)
	#define s1024_2f ((double)0.28440753721127184321)	/* exp(2f*I*twopi/1024) */		
	#define c1024_31 ((double)0.95514116830577072162)
	#define s1024_31 ((double)0.29615088824362382370)	/* exp(31*I*twopi/1024) */
	#define c1024_33 ((double)0.95143502096900836968)
	#define s1024_33 ((double)0.30784964004153489325)	/* exp(33*I*twopi/1024) */		
	#define c1024_35 ((double)0.94758559101774113480)
	#define s1024_35 ((double)0.31950203081601567745)	/* exp(35*I*twopi/1024) */
	#define c1024_37 ((double)0.94359345816196036165)
	#define s1024_37 ((double)0.33110630575987640127)	/* exp(37*I*twopi/1024) */		
	#define c1024_39 ((double)0.93945922360218991213)
	#define s1024_39 ((double)0.34266071731199439711)	/* exp(39*I*twopi/1024) */
	#define c1024_3b ((double)0.93518350993894757782)
	#define s1024_3b ((double)0.35416352542049038186)	/* exp(3b*I*twopi/1024) */		
	#define c1024_3d ((double)0.93076696107898373214)
	#define s1024_3d ((double)0.36561299780477386950)	/* exp(3d*I*twopi/1024) */
	#define c1024_3f ((double)0.92621024213831134218)
	#define s1024_3f ((double)0.37700741021641825620)	/* exp(3f*I*twopi/1024) */
	#define c1024_41 ((double)0.92151403934204194368)
	#define s1024_41 ((double)0.38834504669882629109)	/* exp(41*I*twopi/1024) */
	#define c1024_43 ((double)0.91667905992104266335)
	#define s1024_43 ((double)0.39962419984564682799)	/* exp(43*I*twopi/1024) */		
	#define c1024_45 ((double)0.91170603200542985165)
	#define s1024_45 ((double)0.41084317105790394162)	/* exp(45*I*twopi/1024) */
	#define c1024_47 ((double)0.90659570451491536559)
	#define s1024_47 ((double)0.42200027079979968537)	/* exp(47*I*twopi/1024) */		
	#define c1024_49 ((double)0.90134884704602201485)
	#define s1024_49 ((double)0.43309381885315196790)	/* exp(49*I*twopi/1024) */
	#define c1024_4b ((double)0.89596624975618515621)
	#define s1024_4b ((double)0.44412214457042923104)	/* exp(4b*I*twopi/1024) */		
	#define c1024_4d ((double)0.89044872324475789026)
	#define s1024_4d ((double)0.45508358712634382292)	/* exp(4d*I*twopi/1024) */
	#define c1024_4f ((double)0.88479709843093778043)
	#define s1024_4f ((double)0.46597649576796617728)	/* exp(4f*I*twopi/1024) */		
	#define c1024_51 ((double)0.87901222642863347817)
	#define s1024_51 ((double)0.47679923006332213271)	/* exp(51*I*twopi/1024) */
	#define c1024_53 ((double)0.87309497841829009899)
	#define s1024_53 ((double)0.48755016014843595399)	/* exp(53*I*twopi/1024) */		
	#define c1024_55 ((double)0.86704624551569265185)
	#define s1024_55 ((double)0.49822766697278185175)	/* exp(55*I*twopi/1024) */
	#define c1024_57 ((double)0.86086693863776727973)
	#define s1024_57 ((double)0.50883014254310703626)	/* exp(57*I*twopi/1024) */		
	#define c1024_59 ((double)0.85455798836540052117)
	#define s1024_59 ((double)0.51935599016558958668)	/* exp(59*I*twopi/1024) */
	#define c1024_5b ((double)0.84812034480329725170)
	#define s1024_5b ((double)0.52980362468629466753)	/* exp(5b*I*twopi/1024) */		
	#define c1024_5d ((double)0.84155497743689841004)
	#define s1024_5d ((double)0.54017147272989288060)	/* exp(5d*I*twopi/1024) */
	#define c1024_5f ((double)0.83486287498638005676)
	#define s1024_5f ((double)0.55045797293660480227)	/* exp(5f*I*twopi/1024) */		
	#define c1024_61 ((double)0.82804504525775575255)
	#define s1024_61 ((double)0.56066157619733602312)	/* exp(61*I*twopi/1024) */
	#define c1024_63 ((double)0.82110251499110467956)
	#define s1024_63 ((double)0.57078074588696727951)	/* exp(63*I*twopi/1024) */		
	#define c1024_65 ((double)0.81403632970594836217)
	#define s1024_65 ((double)0.58081395809576454434)	/* exp(65*I*twopi/1024) */
	#define c1024_67 ((double)0.80684755354379927274)
	#define s1024_67 ((double)0.59075970185887422768)	/* exp(67*I*twopi/1024) */		
	#define c1024_69 ((double)0.79953726910790503405)
	#define s1024_69 ((double)0.60061647938386892590)	/* exp(69*I*twopi/1024) */
	#define c1024_6b ((double)0.79210657730021235236)
	#define s1024_6b ((double)0.61038280627630945196)	/* exp(6b*I*twopi/1024) */		
	#define c1024_6d ((double)0.78455659715557523362)
	#define s1024_6d ((double)0.62005721176328917788)	/* exp(6d*I*twopi/1024) */
	#define c1024_6f ((double)0.77688846567323245066)
	#define s1024_6f ((double)0.62963823891492702460)	/* exp(6f*I*twopi/1024) */		
	#define c1024_71 ((double)0.76910333764557963998)
	#define s1024_71 ((double)0.63912444486377574303)	/* exp(71*I*twopi/1024) */
	#define c1024_73 ((double)0.76120238548426181469)
	#define s1024_73 ((double)0.64851440102211244430)	/* exp(73*I*twopi/1024) */		
	#define c1024_75 ((double)0.75318679904361248316)
	#define s1024_75 ((double)0.65780669329707865614)	/* exp(75*I*twopi/1024) */
	#define c1024_77 ((double)0.74505778544146596311)
	#define s1024_77 ((double)0.66699992230363750586)	/* exp(77*I*twopi/1024) */		
	#define c1024_79 ((double)0.73681656887736987581)
	#define s1024_79 ((double)0.67609270357531595956)	/* exp(79*I*twopi/1024) */
	#define c1024_7b ((double)0.72846439044822519723)
	#define s1024_7b ((double)0.68508366777270038056)	/* exp(7b*I*twopi/1024) */		
	#define c1024_7d ((double)0.72000250796138162984)
	#define s1024_7d ((double)0.69397146088965400820)	/* exp(7d*I*twopi/1024) */
	#define c1024_7f ((double)0.71143219574521644231)
	#define s1024_7f ((double)0.70275474445722530165)	/* exp(7f*I*twopi/1024) */

#endif	/* #ifndef radix1024_included */
