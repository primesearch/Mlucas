/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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
#ifndef fac_test_dat128_included
#define fac_test_dat128_included

	#include "types.h"

	struct testFac128{
		uint32 p;
		uint64 d1;
		uint64 d0;
	};

	struct testFac128x2{
		uint64 phi;
		uint64 plo;
		uint64 d1;
		uint64 d0;
	};

	/*******************************************/
	/*      Fermat-number test factors:        */
	/*******************************************/

	// Here interpret the above testFac128 struct as a minimalist [n,k1,k0]-trio format,
	// where k = k1*2^64 + k0, Fn = 2^2^n+1 is the Fermat number and q = k.2^(n+2)+1 the factor.
	// To check any particular (alleged) factor q of Fn using Pari, use Mod(2,q)^(2^n)+1.

	// Testcases with factors < 2^128:
	static const struct testFac128 ffac128[] =
	{
		{ 7,        0ull,11141971095088142685ull},			// 1970 Morrison/Brillhart
		{11,0ull, 20506415569062558ull},{11,0ull,434673084282938711ull},	// 1988 Brent
		{13,      528ull,11889784005418124336ull},			// 1995 Brent
		{15, 69801146ull,13376893484985149619ull},			// 1997 Crandall/van Halewyn
		{16,       39ull, 1485176525646847404ull},			// 1996 Crandall/Dilcher
		{18,        0ull,   77509585098133576ull},			// 1999 Crandall/McIntosh/Tardif
		{19,971680160ull,15036099344401914846ull},			// 2009 Bessell/Woltman
		{22,208923546ull,18410182759848881875ull},			// 2010 Bessell/Woltman
		{31,        0ull,       5463561471303ull},			// 2001 Kruppa/Forbes
		{37,        0ull,          1275438465ull},			// 1991 Gostin
		{39,        0ull,    2864929972774011ull},			// 2012 Rajala/Woltman
		{42,        0ull,     222636358286122ull},			// 2011 Maznichenko/Rodenkirch
		{43,        0ull,        212675402445ull},			// 2000 Samidoost/Durman
		{48,        0ull,       2139543641769ull},			// 2001 Bodschwinna/Durman
		{52,0ull, 4119ull},{52,0ull, 21626655ull},{52,0ull,81909357657279ull},	// 1963 Wrathall, 1982 Keller, 2010 Vonck/Durman
		{58,        0ull,                 190ull},			// 1957 Robinson
		{61,        0ull,           439880504ull},			// 1986 Gostin
		{62,        0ull,                 697ull},			// 1977 Shippee
		{63,        0ull,                  36ull},			// 1956 Robinson
		{64,        0ull,            35707278ull},			// 1986 Gostin
		{65,        0ull,    2421791520862166ull},			// 2011 Maznichenko/Woltman
		{66,        0ull,               15102ull},			// 1977 Shippee
		{71,        0ull,                 683ull},			// 1977 Shippee
		{72,        0ull,            76432329ull},			// 1986 Gostin
		{73,        0ull,                   5ull},			// 1906 Morehead
		{75,        0ull,             3447431ull},			// 1982 Gostin
		{77,0ull, 425ull},{77,0ull,5940341195ull},			// 1957 Robinson/Selfridge, 1999 Taura
		{81,        0ull,                 542ull},			// 1957 Robinson/Selfridge
		{83,        0ull,       6383454640628ull},			// 2005 Danilov/Durman
		{88,        0ull,        119942751127ull},			// 2001 Nohara/Durman
		{99,        0ull,              129864ull},			// 1979 Gostin/McLaughlin/Suyama
		{117,       0ull,                  14ull},			// 1956 Robinson
		{0,0ull,0ull}
	};

	// Testcases with factors > 256 bits (i.e. larger than our current largest fixed-length
	// modpow support, hence necessitating the use of the arbitrary-length mi64-code modpow).
	// Here interpret the above testFac128 struct as a [n,pow2,k] trio, Fn = 2^2^n+1 is the
	// Fermat number and q = k.2^pow2 + 1 the factor. Limit ourselves to n < 10000 in order
	// to keep the resulting self-test timings reasonable.
	// To check any particular (alleged) factor q of Fn using Pari, use Mod(2,q)^(2^n)+1.
	static const struct testFac128 ffacBig[] =
	{
		{ 230, 232ull,372236097ull},
		{ 232, 236ull,70899775ull},
		{ 250, 252ull,403ull},
		{ 251, 254ull,85801657ull},
		{ 255, 257ull,629ull},
		{ 256, 258ull,36986355ull},
		{ 259, 262ull,36654265ull},
		{ 267, 271ull,177ull},
		{ 268, 276ull,21ull},
		{ 275, 279ull,22347ull},
		{ 284, 286ull,1061341513ull},
		{ 286, 288ull,78472588395ull},
		{ 287, 289ull,5915ull},
		{ 297, 301ull,72677552745ull},
		{ 298, 302ull,247ull},
		{ 299, 304ull,272392805475ull},
		{ 301, 304ull,7183437ull},
		{ 316, 320ull,7ull},
		{ 329, 333ull,1211ull},
		{ 334, 341ull,27609ull},
		{ 338, 342ull,27654487ull},
		{ 343, 345ull,4844391185ull},
		{ 353, 355ull,18908555ull},
		{ 370, 373ull,573230511ull},
		{ 375, 377ull,733251ull},
		{ 376, 378ull,810373ull},
		{ 380, 385ull,321116871ull},
		{ 398, 401ull,120845ull},
		{ 416, 419ull,38039ull},
		{ 417, 420ull,303472680883ull},
		{ 431, 434ull,5769285ull},
		{ 452, 455ull,27ull},
		{ 459, 465ull,5449229488169ull},
		{ 468, 471ull,27114089ull},
		{ 480, 484ull,5673968845ull},
		{ 517, 520ull,84977118993ull},
		{ 544, 547ull,225ull},
		{ 547, 550ull,77377ull},
		{ 556, 558ull,127ull},
		{ 569, 575ull,6616590375ull},
		{ 579, 581ull,63856313ull},
		{ 600, 605ull,6213186413ull},
		{ 620, 624ull,10084141ull},
		{ 635, 645ull,4258979ull},
		{ 637, 643ull,11969ull},
		{ 642, 644ull,52943971ull},
		{ 666, 668ull,217924552867ull},
		{ 667, 669ull,491628159ull},
		{ 692, 695ull,717ull},
		{ 723, 730ull,554815ull},
		{ 744, 747ull,17ull},
		{ 851, 859ull,497531ull},
		{ 885, 887ull,16578999ull},
		{ 906, 908ull,57063ull},
		{ 931, 933ull,1985ull},
		{ 943, 954ull,4785972759ull},
		{ 971, 976ull,541664191ull},
		{1069,1073ull,137883ull},
		{1082,1084ull,82165ull},
		{1114,1116ull,11618577ull},
		{1123,1125ull,25835ull},
		{1132,1136ull,10111717305ull},
		{1160,1162ull,2018719057ull},
		{1201,1203ull,837747239ull},
		{1225,1231ull,79707ull},
		{1229,1233ull,29139ull},
		{1394,1396ull,62705223ull},
		{1451,1454ull,13143ull},
		{1551,1553ull,291ull},
		{1598,1600ull,10923781ull},
		{1710,1719ull,351276975ull},
		{1722,1724ull,364182745ull},
		{1849,1851ull,98855ull},
		{1945,1947ull,5ull},
		{1990,1993ull,150863ull},
		{2023,2027ull,29ull},
		{2059,2063ull,591909ull},
		{2089,2099ull,431ull},
		{2420,2422ull,103257279ull},
		{2456,2458ull,85ull},
		{2606,2608ull,238451805ull},
		{3310,3313ull,5ull},
		{3314,3322ull,406860969ull},
		{3335,3337ull,43714055ull},
		{3506,3508ull,501ull},
		{3703,3706ull,262254673ull},
		{3723,3725ull,13308899ull},
		{4184,4189ull,465917283ull},
		{4250,4252ull,173373ull},
		{4258,4262ull,1435ull},
		{4260,4262ull,209161375ull},
		{4265,4269ull,72179955ull},
		{4332,4334ull,2466157ull},
		{4652,4654ull,143918649ull},
		{4724,4727ull,29ull},
		{5320,5323ull,21341ull},
		{5531,5533ull,1503975ull},
		{5792,5794ull,8872947ull},
		{5957,5960ull,421435ull},
		{6208,6210ull,763ull},
		{6355,6358ull,115185ull},
		{6390,6393ull,303ull},
		{6537,6539ull,17ull},
		{6835,6838ull,19ull},
		{6909,6912ull,6021ull},
		{7181,7187ull,168329ull},
		{7309,7312ull,145ull},
		{8239,8242ull,7473ull},
		{8269,8271ull,592131ull},
		{8298,8300ull,1054057ull},
		{8555,8557ull,645ull},
		{9322,9324ull,8247ull},
		{9428,9431ull,9ull},
		{9447,9449ull,5505161ull},
		{9448,9450ull,19ull},
		{9549,9551ull,1211ull},
		{9691,9693ull,260435ull},
		{9747,9749ull,44670651ull},
		{   0,   0ull,0ull}
	};

	/*******************************************/
	/*      Mersenne-number test factors:      */
	/*******************************************/

	/* Factors > 96 but <= 128 bits. If desired, we can construct more test factors
	by multiplying together a 63/64-bit factor q1 of M(p1) and a 65/64-bit factor q2 of M(p2)
	and checking whether q1*q2 divides M(p1*p2).*/
	static const struct testFac128 fac128[] =
	{
		{     695,   12240518780192025ull, 1654746039858251761ull},
		{     845, 2923447923687422893ull,  353773459776294223ull},
		{    1113,     128099917305337ull, 7733695761441692271ull},
		{    1145,   10811609908058563ull, 5349936413307099433ull},
		{    1149,        700245770430ull,  701890237964104231ull},
		{     733,  756146046438660814ull, 7804835620876695225ull},
		{     737,  106450884062962221ull,17050154159176967743ull},
		{     871, 7448657723978021346ull,15223106393317212577ull},
		{     947,     644719741813452ull,16055621295463638505ull},
		{     953,      44696312570505ull, 4961431981940124743ull},
		{     989,   99970972632587991ull, 1738540175825943815ull},
		{    1081,         67677680549ull,13887741953162944095ull},
		{    1091,    5287390011750720ull, 2894679571106043497ull},
		{    1097,      11129117045170ull,10375766809019373543ull},
		{    1099,       1551337752834ull, 8321741389535251703ull},
		{    1133,  133834206206032981ull, 6586095673132787791ull},
		{    1141,       5747037125100ull, 2460710484528304153ull},
		{    1181,      10824073357153ull, 7361144750966677159ull},
		{    1189,   32559650929209964ull, 8212830436061989903ull},
		{   27691,   94004235929829273ull, 4235226679561594903ull},
		{  319057,        103337218078ull, 8676403300852410079ull},
		{17363977,      62897895526806ull,14211535226588354713ull},
		{10624093,          5492917609ull,14854696485656401105ull},
		{10698673,          5799457823ull,10285356664749312993ull},
		{20799431,          4303087381ull,16578386512849109713ull},
		{33652757,          5202063708ull,18263664019678288919ull},
		{21823211,          7785579841ull, 7607475409566672241ull},
		{22330859,          7593776864ull, 5630449305759171207ull},
		{11808917,         20308449831ull, 9058738039473012457ull},
		{20090969,         15531431134ull, 5034609389988515233ull},
		{20313967,         18216394609ull, 8291172543411688687ull},
		{20544481,         16259503442ull,15859685870849762975ull},
		{22217387,         20551559047ull,11995354231649723881ull},
		{10207999,         28364424832ull,15122069645900159367ull},
		{19964723,         34441477586ull, 9636073161914837921ull},
		{21145199,         30977655046ull, 1304857345634219175ull},
		{22030163,         43144178324ull, 4788416424359163737ull},
		{33562153,         45963786472ull, 2258783450670948535ull},
		{33693587,         66325700032ull,15262466751214122975ull},
		{11865241,         57210216387ull, 3082735332820781609ull},
		{21801929,         80355238912ull,15689518004012743009ull},
		{19951201,        109346652057ull,10819675441336938065ull},
		{20616781,      17534809723250ull,10329047311584913071ull},
		{20648443,       1221873279710ull, 2595613477835803991ull},
		{21250771,      12549422209078ull, 8612165677489771129ull},
		{21547787,        112416184026ull, 9015544550402598895ull},
		{21675733,        142220976614ull,11385509628023387489ull},
		{15714269,      14320762091913ull, 2773697912020767049ull},
		{19687561,       1996508583829ull, 7515490546312285159ull},
		{20152333,        365842230851ull, 2388855518206098663ull},
		{20510053,        261078947686ull,  465403687781705377ull},
		{20759821,     199835753775288ull,17079803649869853575ull},
		{20989043,        202355339943ull,15105677628487752455ull},
		{33713123,      18738454648009ull,16692905930976531153ull},
		{20542751,        412571049040ull,18170931828058363183ull},
		{20812849,     534505286298455ull, 2216600112648316881ull},
		{0,0,0ull}
	};

	/* Factors > 96 but <= 128 bits, with p > 64 bits - most of these are from my Jan 2003 runs of primes near 2^89: */
	static const struct testFac128x2 fac128x2[] =
	{
		{33554431ull,18446744073709551175ull,      30899672449023ull,18446744073303442655ull},
		{33554431ull,18446744073709551295ull,      85098334519295ull,18446744072895454529ull},
		{33554431ull,18446744073709551513ull,  430360347665235967ull,18446742752658555745ull},
		{33554431ull,18446744073709551513ull,  259661119604391935ull,18446743276643598623ull},
		{33554431ull,18446744073709551513ull,  293843670505881599ull,18446743171715500217ull},
		{33554431ull,18446744073709551567ull,   12593691025735679ull,18446744055318810857ull},
		{33554431ull,18446744073709551595ull,            67108863ull,18446744073709551575ull},
		{33554431ull,18446744073709551595ull,        631158865919ull,18446744073709156607ull},
		{33554432ull,                  89ull,      16384159383552ull,            43457455ull},
		{33554432ull,                 705ull, 7006880245689090048ull,     147219019329871ull},
		{33554432ull,                 741ull,        882615779328ull,            19491265ull},
		{33554432ull,                 741ull,     220386851553280ull,          4866917641ull},
		{33554432ull,                 767ull,     460834488713216ull,         10533930447ull},
		{33554432ull,                 837ull,    3613005073874944ull,         90124763455ull},
		{33554432ull,                 837ull,    1504315717976064ull,         37524469375ull},
		{33554432ull,                1059ull,            67108864ull,                2119ull},
		{33554432ull,                1275ull,          4898947072ull,              186151ull},
		/*{33554432ull,                1275ull,    1096248325046272ull,         41655201151ull},*/
		{33554432ull,                1337ull,      16333760626688ull,           650830209ull},
		{33554432ull,                1337ull,       1481763717120ull,            59041921ull},
		{33554432ull,                1521ull,   19343955296518144ull,        876848578633ull},
		{33554432ull,                1547ull,          2147483648ull,               99009ull},
		{33554432ull,                1917ull,    1124806032359424ull,         64261351945ull},
		/*{33554432ull,                1917ull,    2192593957945344ull,        125265199465ull},*/
		{33554432ull,                1917ull,         22749904896ull,             1299727ull},
		{33554432ull,                2097ull,      34829164871680ull,          2176665031ull},
		{33554432ull,                2585ull,      37028188127232ull,          2852614711ull},
		/*{33554432ull,                2585ull,    1096267786616832ull,         84455377711ull},*/
		{33554432ull,                2661ull,          2080374784ull,              164983ull},
		{33554432ull,                2675ull,           805306368ull,               64201ull},
		{33554432ull,                2729ull,      42057326395392ull,          3420544975ull},
		{33554432ull,                2729ull,       4876868255744ull,           396638319ull},
		{33554432ull,                2907ull,         22615687168ull,             1959319ull},
		{33554432ull,                3045ull,    1301649700159488ull,        118122200281ull},
		{33554432ull,                3045ull,         99321118720ull,             9013201ull},
		{33554432ull,                3155ull,      26011932557312ull,          2445806481ull},
		{33554432ull,                3159ull,        413994582016ull,            38975743ull},
		{33554432ull,                3507ull,        250450280448ull,            26176249ull},
		{33554432ull,                4155ull,   62084915730055168ull,       7687891270471ull},
		{33554432ull,                4451ull,          7851737088ull,             1041535ull},
		{33554432ull,                4485ull,        434328567808ull,            58053841ull},
		{33554432ull,                4745ull,         82678120448ull,            11691681ull},
		{33554432ull,                4745ull,         71672266752ull,            10135321ull},
		{33554432ull,                5121ull,        143814295552ull,            21948607ull},
		{33554432ull,                5247ull,            67108864ull,               10495ull},
		{33554432ull,                5411ull,     544322411823104ull,         87777631593ull},
		{33554432ull,                5499ull,    1739631741632512ull,        285096017935ull},
		{33554432ull,                5735ull,      10563740499968ull,          1805515641ull},
		{33554432ull,                5837ull,       7317214986240ull,          1272874591ull},
		{0ull,0ull,0ull,0ull}
	};

#endif	/* #ifndef fac_test_dat128_included */
