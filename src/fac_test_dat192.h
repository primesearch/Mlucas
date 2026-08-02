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

#if(defined(P3WORD) || defined(P4WORD))	// only the 192-bit modpow consumes these
	/* Constructed (p,k) vectors giving q = 2.k.p+1 > 2^128, for the q-pipelined twopmodq192_q4/_q8.

	Those routines take a *64-bit* k together with a 192-bit exponent p, and form q themselves. No
	catalogued Mersenne factor can reach q > 2^128 under that signature: q = 2.k.p+1 with k < 2^64
	forces p > 2^64, and no known Mersenne exponent is anywhere near that big. (That is also why the
	fac160[]/fac192[] loops in factor_test.h guard their q4 call on k < 2^64 and, correctly, skip
	every one of their entries.) So the only way to exercise the routines above 2^128 - the size
	range they exist for - is to construct vectors, which is what these are.

	Each entry is one exponent p plus the four k of a single q4/q8 call. For the lanes flagged in
	.res, p is prime, q = 2.k.p+1 is prime, and 2^p == 1 (mod q), i.e. q genuinely divides 2^p-1;
	the remaining lanes are verified non-factors, so .res is a non-trivial mask and lane crosstalk
	cannot hide. Several non-factor lanes carry k near its 64-bit ceiling, which drives q to just
	under 2^192 and exercises the full-width q = 2.k.p+1 formation.

	Only a lane whose q really is a factor can detect a wrong answer: the routines return one
	divisibility bit per lane, so on a non-factor q a corrupted residue still reads "not a factor".
	Measured on a build independently shown to mishandle 28 of 144 genuine factors above 2^128, a
	random-k differential test against an arbitrary-precision oracle flagged 5 of 16000 lanes; these
	vectors flag 19% of theirs. Hence genuine factors, not random q.

	Construction: pick a prime p of the wanted size, then the smallest k for which q = 2.k.p+1 is
	prime and 2^p == 1 (mod q). That last condition holds with probability 1/(2k), so the search is
	only tractable for small k - which is why every p here is within a few bits of its q, and why
	there is no entry with both a small p and a large k. Every lane was re-verified with GMP. */
	struct testFac192p {
		uint64 p2,p1,p0;	// exponent p, 192 bits, low word last
		uint64 k[4];		// the four trial factors of one q4 call: q_j = 2.k_j.p + 1
		uint32 res;			// expected q4 return: bit j set iff q_j divides 2^p-1
	};

	static const struct testFac192p fac192p[] =
	{
		{0x0000000000000000ull,0x00213a5cfb4ffaf1ull,0x7fc498a800cedc8full, {                8593ull,15214548644798751000ull,12775845238650594594ull,13075755774346653527ull}, 0x1},	// p 118 bits, q 132/182/182/182 bits
		{0x0000000000000000ull,0x008fe5d795b4c471ull,0x6d71034ef4524abfull, { 2763321911299042568ull,                3317ull,16679543456978578366ull,14224913552053546146ull}, 0x2},	// p 120 bits, q 182/132/185/184 bits
		{0x0000000000000000ull,0x0431c9000e2a4606ull,0xd8bfcb23808af8c1ull, {10707202604359032825ull, 4870861780386308287ull,                  72ull, 1324479603544283704ull}, 0x4},	// p 123 bits, q 187/186/130/184 bits
		{0x0000000000000000ull,0x1700d2cdc47834efull,0x882ea54142645801ull, {15783913230357963913ull,12847826680588602063ull,12904419274851956850ull,                   7ull}, 0x8},	// p 125 bits, q 190/190/190/129 bits
		{0x0000000000000000ull,0x518bbe5cd5cd507cull,0x692ede0d65b6241dull, {                1488ull, 5759559193453849237ull,15687860524863433765ull, 3693094640594910208ull}, 0x1},	// p 127 bits, q 138/190/192/190 bits
		{0x0000000000000000ull,0xc9fa05cfdebfd1afull,0x8ce952aa4839afbdull, { 3900782993366110328ull,                 395ull, 2206542630629832720ull, 9303172449162996934ull}, 0x2},	// p 128 bits, q 191/138/190/192 bits
		{0x0000000000000000ull,0xd2c0c7af7db202adull,0xa9748a1b0b753ee1ull, { 8505764787670883816ull, 8118806472317829564ull,                 267ull, 8060509614542063432ull}, 0x4},	// p 128 bits, q 192/192/137/192 bits
		{0x0000000000000003ull,0xe0ddc3f674698a89ull,0xd8cf1d4cefe6151full, { 2237421913063066265ull, 1340132015343656597ull,  930454171901463158ull,                 209ull}, 0x8},	// p 130 bits, q 192/192/191/139 bits
		{0x000000000000003full,0x9f86281a4026e6f7ull,0x8cb149419b2d9c2full, {                  37ull,   96135231250510595ull,   76867814451075485ull,   88800288255623438ull}, 0x1},	// p 134 bits, q 141/192/192/192 bits
		{0x00000000000003d4ull,0x29a4c79a1055a822ull,0x84d7046e5c9d59a9ull, {    1060076026885497ull,                  12ull,    3653361580340800ull,    8998231113037378ull}, 0x2},	// p 138 bits, q 189/143/191/192 bits
		{0x00000000000026e3ull,0x10ab923c90aaf476ull,0x2423c06e3f12c2cfull, {     475711864557646ull,     856067926826099ull,                   4ull,     591975508701151ull}, 0x4},	// p 142 bits, q 192/192/145/192 bits
		{0x0000000000075d26ull,0x20e30535c1a92b94ull,0x8ca3bd247109410dull, {      12957052590606ull,      15022157530405ull,      14212087166057ull,               12420ull}, 0x8},	// p 147 bits, q 192/192/192/162 bits
		{0x0000000000f49ea5ull,0xb1e494212cb76026ull,0x2ddc435b1584a115ull, {                   3ull,        339195661178ull,        345270076518ull,        389580843610ull}, 0x1},	// p 152 bits, q 155/192/192/192 bits
		{0x0000000010b355b6ull,0x93a9035cf075506full,0x6d238fa396db6159ull, {          4977490711ull,                1611ull,         18808043090ull,         30858213833ull}, 0x2},	// p 157 bits, q 190/168/192/192 bits
		{0x00000002ad400816ull,0x5cb3dab219806970ull,0xb380b4ca53beb057ull, {           654612475ull,           230478869ull,               29484ull,           490678292ull}, 0x4},	// p 162 bits, q 192/191/178/192 bits
		{0x00000044774cf897ull,0xc53d0423917bba5bull,0xb01c3286f1703141ull, {            24992275ull,            26461495ull,             1185738ull,                 911ull}, 0x8},	// p 167 bits, q 192/192/188/177 bits
		{0x00000e359624b691ull,0x01ce6d2edef1dc07ull,0x1f4c5667cf7e7e93ull, {               22624ull,              365880ull,              558918ull,              481967ull}, 0x1},	// p 172 bits, q 188/192/192/192 bits
		{0x0001ae2fbb528809ull,0xbdef6b6ed0d0739full,0x3e0ac4b267758513ull, {                2154ull,                   1ull,               10670ull,               19152ull}, 0x2},	// p 177 bits, q 189/178/192/192 bits
		{0x00139e67d602a053ull,0xec7efb8187dd52d0ull,0x61299dd9913dcf4bull, {                1651ull,                1030ull,                 117ull,                 602ull}, 0x4},	// p 181 bits, q 192/192/189/191 bits
		{0x01a68f7e96e1ab2eull,0xba58038b9102f259ull,0x4273ea77f1c98e57ull, {                  14ull,                  61ull,                  63ull,                   5ull}, 0x8},	// p 185 bits, q 190/192/192/189 bits
		{0x0ea3146a17417a61ull,0x0146d76d5e6c30a2ull,0xf48afcb93f22a1e5ull, {                   3ull,                   3ull,                   3ull,                   3ull}, 0xf},	// p 188 bits, q 191/191/191/191 bits
		{0x330ff598c0b907a2ull,0x0ac9447d0ec04cbdull,0xd75814e2362eb6cbull, {                   1ull,                   1ull,                   1ull,                   1ull}, 0xf},	// p 190 bits, q 191/191/191/191 bits
		{0x429cd78c5dfb44b0ull,0xbd3d719d6629dd42ull,0xbe7cedeeba47d14full, {                   1ull,                   1ull,                   1ull,                   1ull}, 0xf},	// p 191 bits, q 192/192/192/192 bits
		{0x6213754c091bd1adull,0x386dc79f173faeb4ull,0x2c491d550c8ca76bull, {                   1ull,                   1ull,                   1ull,                   1ull}, 0xf},	// p 191 bits, q 192/192/192/192 bits
		{0ull,0ull,0ull, {0ull,0ull,0ull,0ull}, 0}
	};

#endif	/* #if(defined(P3WORD) || defined(P4WORD)) */

#endif	/* #ifndef fac_test_dat192_included */
