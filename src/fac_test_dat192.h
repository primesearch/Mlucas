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
#ifndef fac_test_dat192_included
#define fac_test_dat192_included

	#include "types.h"

	struct testFac160{
		uint32 p;
		uint64 d2;
		uint64 d1;
		uint64 d0;
	};

	struct testFac192{
		uint32 p;
		uint64 d2;
		uint64 d1;
		uint64 d0;
	};

	/*******************************************/
	/*      Fermat-number test factors:        */
	/*******************************************/

	// Here interpret the above testFac struct as a minimalist [n,k]-pair format,
	// where Fn = 2^2^n+1 is the Fermat number and q = k.2^(n+2)+1 the factor:
	// To check any particular (alleged) factor q of Fn using Pari, use Mod(2,q)^(2^n)+1.

	// Testcases with factors < 2^192:
	static const struct testFac192 ffac192[] =
	{
		{ 86,0ull,0ull,	   20018578522347ull},		// 2012 M. Dangler & Rodenkirch
		{ 88,0ull,0ull,	     119942751127ull},		// 2001 T. Nohara & Durman
		{ 90,0ull,0ull,	     198922467387ull},		// 2001 P. Grobstich & Durman
		{ 91,0ull,0ull,	             1421ull},		// 1977 D. E. Shippee
		{ 93,0ull,0ull,2*	        92341ull},		// 1979 R. Baillie
		{ 94,0ull,0ull,2*	 482524552001ull},		// 2001 P. Grobstich & Durman
		{ 96,0ull,0ull,8*	3334131633063ull},		// 2008 M. Ptáček & Durman
		{107,0ull,0ull,4*	   1289179925ull},		// 1992 G. B. Gostin
		{116,0ull,0ull,4*	   3433149787ull},		// 1999 T. Taura
		{122,0ull,0ull,	          5234775ull},		// 1986 G. B. Gostin
		{125,0ull,0ull,	                5ull},		// 1956 R. M. Robinson
		{133,0ull,0ull,	      88075576149ull},		// 2001 P. Samidoost & Durman
		{142,0ull,0ull,2*	      8152599ull},		// 1986 G. B. Gostin
		{144,0ull,0ull,2*	           17ull},		// 1956 R. M. Robinson
		{146,0ull,0ull,	         37092477ull},		// 1987 G. B. Gostin
		{147,0ull,0ull,	             3125ull},		// 1979 G. B. Gostin & P. B. McLaughlin
		{147,0ull,0ull,	        124567335ull},		// 1990 G. B. Gostin
		{150,0ull,0ull,32*	         1575ull},		// 1956 R. M. Robinson
		{150,0ull,0ull,4*	         5439ull},		// 1980 G. B. Gostin & P. B. McLaughlin & H. Suyama
		{0,0ull,0ull,0ull}
	};

	/*******************************************/
	/*      Mersenne-number test factors:      */
	/*******************************************/

	/* Factors > 128 but <= 160 bits. If desired, we can construct more test factors
	by multiplying together a 64-bit factor q1 of M(p1) and a 96-bit factor q2 of M(p2)
	and checking whether q1*q2 divides M(p1*p2).*/
	static const struct testFac160 fac160[] =
	{
		{     629,       133ull,11545660419510266595ull,15875370168207932041ull},
		{     631,      1394ull,15571349859840161706ull,  509892144742137431ull},
		{     673,    121320ull, 4492854135134704005ull,14226674137430228263ull},
		{     695,2649519282ull,14842833464112563611ull,10174116463236461383ull},
		{     731, 655903171ull,17652352551621896287ull, 7660429456444636239ull},
		{     805,1083827012ull,18314245293386716597ull, 2219421057460140527ull},
		{     877,  13161208ull,18225246095436784582ull,12343089078196252631ull},
		{     957,      4730ull,14663183769241509326ull, 8097149896429635207ull},
		{     967,    215159ull,  881920578744577810ull,17184239148975426263ull},
		{    1017, 212724356ull, 9900144438119899815ull,17733134473107607967ull},
		{    1033,       261ull, 5238930328752646394ull, 2803405107698253561ull},
		{    1087,         1ull, 4415476118538293365ull,16346425147370540471ull},
		{    1087,     70130ull,11905462972019801043ull, 6167785434693019223ull},
		{    1131,   5800574ull,18429773635221665090ull,17951008765075981215ull},
		{    1157,  22381525ull,14500669099417213747ull,15903397166638806257ull},
		{    1283,        14ull, 3291757557782450881ull, 3893270457587058239ull},
		{    1319,      1552ull, 1390029428449091172ull,14288981644299514807ull},
		{    1483,      2674ull,14802171160149427175ull, 5085420234315110585ull},
		{    6659,       664ull,14291576310931480037ull, 4949688733053552967ull},
		{    8191,    617742ull, 6334326874596939334ull,11405337619840706193ull},
		{18031451,      2122ull, 5198971222801411122ull,12425019173815339143ull},	/* Note: composite factor! */
		{0,0ull,0ull,0ull}
	};

	/* Factors > 160 but <= 192 bits. We can construct more test factors by multiplying
	together smaller factors of M(p) with multiple factors, or for exponents p1, p2, p3, ...
	and corresponding factors q1, q2, q3, ... , checking whether q1*q2*q3*...
	divides M(p1*p2*p3*...). */
	static const struct testFac192 fac192[] =
	{
		{     677,     157590042578912ull,10558642444782195772ull,  329809049266961143ull},
		{     773,       9118322195022ull, 1933308633079010416ull,17814616685598394119ull},
		{     971,      70286054459973ull,17012949627558354271ull, 3547755741880899889ull},
		{     997,  492416983078691417ull, 8040689323464953445ull,16007877010440112335ull},
		{    1001,         59364131986ull, 9565712986615012496ull,10050950882119470361ull},
		{0,0ull,0ull,0ull}
	};

#endif	/* #ifndef fac_test_dat192_included */
