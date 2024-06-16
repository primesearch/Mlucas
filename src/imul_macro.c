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

#include "imul_macro.h"
#include "factor.h"

/* Binary predicates for use of stdlib qsort(): Only support unsigned for multiword ints (for now): */

// Need extra () around deref-of-reinterpret-casts here; otherwise get "request for member ‘d1’ in something not a structure or union" errors:
int ncmp_uint128(const void * a, const void * b)
{
	if( CMPULT128( (*(uint128*)a) , (*(uint128*)b) ) ) {
		return -1;
	} else if CMPEQ128( (*(uint128*)a) , (*(uint128*)b) ) {
		return 0;
	} else {
		return +1;
	}
}

int ncmp_uint192(const void * a, const void * b)
{
	if( CMPULT192( (*(uint192*)a) , (*(uint192*)b) ) ) {
		return -1;
	} else if CMPEQ192( (*(uint192*)a) , (*(uint192*)b) ) {
		return 0;
	} else {
		return +1;
	}
}

int ncmp_uint256(const void * a, const void * b)
{
	if( CMPULT256( (*(uint256*)a) , (*(uint256*)b) ) ) {
		return -1;
	} else if CMPEQ256( (*(uint256*)a) , (*(uint256*)b) ) {
		return 0;
	} else {
		return +1;
	}
}

/*********************************************************************************/

/* Subroutine versions of the basic 64-bit integer MUL macros, for crappy compilers
that don't correctly inline the macro form of these.
*/
/* Had to move the subroutine-form definitions of these here to fix MSVC
"error LNK2005 already defined ... fatal error LNK1169 ... multiple symbols found" errors:
*/
#ifdef	MUL_LOHI64_SUBROUTINE

	/* 64x32=>96-bit product algorithm: represent the inputs as x = a + b*2^32, y = c ( < 2^32),
	then do 4 MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .
	COST: 2 MUL, 6 IOPs.

	Even though the high output for 64x32-bit is always < 2^32,
	declare y and hi here as 64-bit ints to allow flexibility for caller:
	*/
	void	MUL64x32(uint64 x, uint64 y, uint64 *lo, uint64 *hi)
	{
		uint64 t;

	/*	*lo = ((uint32)(x & 0x00000000ffffffff)) * y;*/	/* a*c, v1 */
	/* try this 32x64-bit form in hopes compiler can optimize it better than 64x64-bit: */
		*lo =  (uint32)(x) * y;							/* a*c, v2 */
		 t  = ((uint32)(x >> 32)) * y;					/* b*c */
		*hi = (t >> 32);
		 t <<= 32;
		*lo +=  t;
		*hi += (*lo < t);
	}

	/* Generic 128-bit product algorithm: represent the inputs as x = a + b*2^32, y = c + d*2^32, and
	then do 4 MULs and a bunch of add-with-carries to get x*y = b*d*2^64 + (a*d + b*c)*2^32 + a*c .
	COST: 4 MUL, 18 IOPs.
	*/
	void MUL_LOHI64(uint64 x, uint64 y, uint64 *lo, uint64 *hi)
	{
	#if 0
		uint64 a,b,c,d,ac,ad,bc;
		a = x & (uint64)0x00000000ffffffff;
		/*a = (x<<32)>>32;*/
		b =  x>>32;
	#else
	/* try this 32x64-bit form in hopes compiler can optimize it better than 64x64-bit: */
		uint32 a,b;
		uint64 c,d,ac,ad,bc;
		a = (uint32)(x);
		b = (uint32)(x>>32);
	#endif
		c = y & (uint64)0x00000000ffffffff;
		/*c = (y<<32)>>32;*/
		d = y>>32;
		/* Calculate 4 subproducts in order in which they are first used */
		 ac  = a*c;
		 bc  = b*c;
		*hi  = b*d;
		 ad  = a*d;
		*lo  =  ac;	/* use lo to store copy of ac */
		 ac +=       (bc<<32);	*hi += ( ac < *lo);
		*lo  =  ac + (ad<<32);	*hi += (*lo <  ac);
								*hi += (bc>>32) + (ad>>32);
	}

	/* Generic 128-bit squaring algorithm: represent the input as x = a + b*2^32, and
	then do 3 MULs and a bunch of add-with-carries to get x^2 = b^2*2^64 + a*b*2^33 + a^2 . */
	/*
	void SQR_LOHI64(uint64 x, uint64 *lo, uint64 *hi)
	{
		uint64 a = x & (uint64)0x00000000ffffffff;
		uint64 b = x>>32;
		uint64 aa = a*a;
		uint64 ab = a*b;
		uint64 bb = b*b;
		*hi =  bb;
		*lo =  aa + (ab<<33); if(*lo < aa) ++*hi;
		*hi = *hi + (ab>>31);
	}
	*/
	void SQR_LOHI64(uint64 x, uint64 *lo, uint64 *hi)
	{
	#if 0
		uint64 a,b,aa,ab;
		a = x & (uint64)0x00000000ffffffff;
		/*a = (x<<32)>>32;*/
		b = x>>32;
		 ab = a*b;
		 aa = a*a;
		*hi = b*b;
	#else
	/* try this 32x64-bit form in hopes compiler can optimize it better than 64x64-bit: */
		uint32 a,b;
		uint64 ab,aa;
		a = (uint32)(x);
		b = (uint32)(x>>32);
		 ab = a*(uint64)b;
		 aa = a*(uint64)a;
		*hi = b*(uint64)b;
	#endif
		*lo = aa + (ab<<33);
		*hi+= (ab>>31) + (*lo < aa);
	}

	/*
	uint32	__MULL32	(uint32 x32, uint32 y32)
	{
		return x32*y32;
	}

	uint32	__MULH32	(uint32 x32, uint32 y32)
	{
		return ((x32*(uint64)y32) >> 32);
	}
	*/

	uint64	__MULL64	(uint64 x, uint64 y)
	{
		return (uint64)x*y;
	}

	uint64	__MULH64	(uint64 x, uint64 y)
	{
		uint64 lo, hi;
		MUL_LOHI64(x, y, &lo, &hi);
		return hi;
	}

#endif

/************************/

/* Really dumb-but-reliable wide-mul emulation, using bitwise multiword add-and-accumulate: */
void mul_lohi64_via_bitwise_add(uint64 x, uint64 y, uint64 *lo, uint64 *hi)
{
	uint128 product = {0ull, 0ull}, yshift;
	yshift.d0 = y;
	yshift.d1 = 0ull;
	while(x)
	{
		if(x & 1ull)
			ADD128(product, yshift, product);
		x >>= 1;
		ADD128(yshift, yshift, yshift);
	}
	*lo = product.d0;
	*hi = product.d1;
}

/*
Various small factor self-tests.
*/
int test_mul()
{
	/* Some test inputs: */
	static const uint64 in64[] =
	{
		10500365774503436201ull,15641779111884102847ull,12915704346519849697ull, 9287694977263609817ull,10016045939325808271ull,
		12928465186953089863ull,10943084366908922353ull,12431505097071887473ull,15718414880628053417ull,17130259977555383783ull,
		17017539148249318519ull,17331582637129315543ull,17075071760459211223ull,10651829289881962249ull,13875348797185757897ull,
		 9612328233936311881ull,10983089768339907263ull,11011193468080521223ull,12625693644608762351ull,10754411617748801759ull,
		12013002667638145399ull,12825068851465148441ull,17170312062251360359ull,12031568180050170481ull,10566050601376248193ull,
		14859011986541991007ull,14989510186280279071ull,16214445839716767361ull,10598954534434896529ull,14057592431325128593ull,
		16560715604598816241ull,12811182984027587681ull,15539690406724276993ull, 9985641878418216991ull,11449566963261644983ull,
		13764994941016683383ull,16672552177382226647ull,11623423906484739527ull,17169614195768327983ull,17990551511959203257ull,
		12155150476665080887ull, 9789648295780602191ull,10985008068068673007ull,12658180920879996367ull,15044832325095723511ull,
		14046232068642146993ull,11179467653397253391ull,11418381227787294353ull, 9600311723751200999ull,10641782461606146511ull,
		 9238369956591864047ull,14284214035586927489ull,16132196293595357033ull,12132533847075902297ull,13711702151742952271ull,
		13553802287542899791ull, 9814223728786685249ull,11664104294072220647ull,11841941210452963607ull,13495717684024316977ull,
		11453374630895585671ull,12409236261334108927ull,13941936454603083961ull,16695225546321153647ull,15274150651947770057ull,
		10670964685324409993ull,13984615577478771737ull,16478002069230680743ull,11555285165538057473ull,13017176322370888999ull,
		10899353333633834737ull,13082852309834271431ull, 9469759966159211119ull,16562044284339342287ull,12857396115642922399ull,
		 9771683082091843559ull,13254207919455131647ull,15437851890296633887ull,16287562569554793953ull, 9275136914779640273ull,
		11025186727945055441ull, 9854261882174597857ull,15285818081668067519ull,11180960576914148191ull,14428852738137829687ull,
		11752043400810302297ull,15463714506600122897ull,12907674385704768017ull,11394039312173819999ull,10366044942826856209ull,
		17973037128357378841ull,18328139888296813361ull,10003907358728522969ull,15548404807307836393ull,16134040372770842257ull,
		17112902589623655359ull,13909553524167145703ull,11109231768403310369ull,10738982321410344871ull,15374699987470193449ull,
		13762137093959242159ull,14732963741405224223ull, 9989690110846448329ull,16284666399547421921ull,13620718892071208681ull,
		11262843094267061911ull,13846368038170428023ull,10504362932070630521ull,12772485980427473983ull,12631296696930554257ull,
		17931410955470787601ull,12702828058463407919ull,14627091909328461551ull,13932143274652230937ull,12646230872610358319ull,
		17458131644489262929ull,11268951794162923127ull,16341053115349924417ull,12007015401209325023ull,15918852021480399527ull,
		 9341752468258236431ull,10013306797089936121ull,16160582285842380959ull,15811238096354839487ull,12298642243571509687ull,
		12597566906366358017ull,18335410068788651849ull,15899391255100053209ull,16860276014915517791ull,16072139673628398673ull,
		12143859436446628321ull,16375574063043584023ull,16086636940765296889ull,17558635633759926961ull,14646809963723851673ull,
		17895543665490669193ull,16726174377018723457ull,14461463666320960511ull,13440101789246121791ull,10991970915886596511ull,
		11337011726859890543ull,16582660415820077063ull,15998585226804618071ull,14271788618292432257ull,10760230434599829143ull,
		10730782077481732817ull,14331950710727986663ull,10961703930971884903ull,15296117896318004719ull,17040351390832159609ull,
		14702865045484912391ull,11533383163345547927ull,12125342584232876329ull,15810596930722833329ull,16759807228908186473ull,
		14680112878087367873ull,10214864560585659127ull,18111383456669963153ull,11554056516029888111ull,14102697559575163223ull,
		10951036977490542007ull,12068660979042357751ull,13936597432397301953ull,17741531859834789367ull,13038897762986995439ull,
		12923341570853700409ull,12468959361471384367ull,12146036035731849143ull,17389980150980900167ull,13632984658968284321ull,
		12579003037333710127ull,18017349281491722863ull,16300711072026110783ull,12837482707778126729ull,14435333013673979377ull,
		10362937015011578407ull,14943059862662596097ull,13083808953378894799ull,10961371406875908433ull,15948970413144688847ull,
		11851486348839027743ull,13617531092155742599ull,16092742349020987151ull,10748446445825061761ull,13581925937697282097ull,
		17001973141862761433ull,18338048883043493801ull,11445451093632964279ull,11466372601071120199ull,13435999004029409279ull,
		18266526887928755321ull,11868117613510431529ull,14158058431764629761ull,16826194113484414487ull,14626742151491329369ull,
		14224052887860439639ull,13197622355391030409ull,13545314040804294319ull,14152269016994013583ull,16067045133341942353ull,
		15370375936060062839ull,15059986833339178897ull,17402936053714225153ull,15488725481110429159ull,12672311049187782367ull,
		11593765511293358047ull,11085431129434680449ull,10218100740052053263ull,14100978243368072207ull,12183869854186459801ull,
		11855169834318581681ull,17649706665176130473ull,16209706991529970753ull,10440895312670202409ull,14608723870843452887ull,
		16542024904771148857ull,13724050657344286793ull,18033822787942121599ull,13620902108656457471ull,13141770251775333463ull,
		 9593004739201446737ull,11207409667480506167ull, 9315521980684624207ull,12920971770507553729ull,16822696323471127217ull,
		13793486024804669009ull,13533131423401796863ull, 9952502816461400807ull,13236391347894754153ull,18194162210171783167ull,
		 9889326862010987471ull,15295542430741533281ull,14567378334050418241ull,10288195248504974393ull,10040647230145986479ull,
		12344397890070085393ull, 9942526727677974137ull,12597569013851067049ull,10288230960652346249ull,10607814077757516457ull,
		11028045969272017577ull, 9803688218926567111ull,10586347451730006521ull,12699616391325608273ull,10780482862846368239ull,
		11747242053545138911ull,10537253101166306807ull,10520779020741875777ull,10000450921799397911ull,12056468544232052593ull,
		10028147391258765553ull,11886145592975848529ull,11289182387241165217ull, 9523705169798429639ull,0ull
	};
	uint32 i,j,ntest;
	uint64 lo0,hi0,lo1,hi1;
	/* count 'em: */
	ntest = 0;
	while(in64[ntest++]){};

	/* Multiply matrix: */
	for(i=0;i<ntest;++i)
	{
		for(j=0;j<ntest;++j)
		{
			mul_lohi64_via_bitwise_add(in64[i],in64[j],&lo0,&hi0);
		#ifdef MUL_LOHI64_SUBROUTINE
			MUL_LOHI64(in64[i],in64[j],&lo1,&hi1);
		#else
			MUL_LOHI64(in64[i],in64[j], lo1, hi1);
		#endif
			ASSERT(lo1 == lo0, "test_mul() low-output mismatch!");
			ASSERT(hi1 == hi0, "test_mul() hi -output mismatch!");

		/* Squaring is a special case: */
		  if(i ==j)
		  {
		#ifdef MUL_LOHI64_SUBROUTINE
			SQR_LOHI64(in64[i],&lo1,&hi1);
		#else
			SQR_LOHI64(in64[i], lo1, hi1);
		#endif
			ASSERT(lo1 == lo0, "test_mul() low-output mismatch!");
			ASSERT(hi1 == hi0, "test_mul() hi -output mismatch!");
		  }
		}
	}

	return 0;
}

/************************/

/* 02/18/06: 53-bit-input error correction:

ERROR: MUL50X50 mismatch:
  a,b,c,d =     8374846658018695,     8389535175654807, 70261070628162943015706377060352,     2709987495556513
  xhi,xlo =        3808860270810, 10820356292439487904-
  yhi,ylo =        3808860270810, 10820356292439487905
  a,b mod 4, a^b = 11, 11, 00
ERROR: MUL50X50 mismatch:
  a,b,c,d =     7912389366102847,     7634031464243407, 60403429378174083930921552052224,    -3914508548372495
  xhi,xlo =        3274476467869, 10219256645582653424-
  yhi,ylo =        3274476467869, 10219256645582653425
  a,b mod 4, a^b = 11, 11, 00
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     7774576030739923,     5171619810855699, 40207151421578450634795078647808,     1044807672723369
  xhi,xlo =        2179634046036,  8296675321289177000-
  yhi,ylo =        2179634046036,  8296675321289177001
  a,b mod 4, a^b = 11, 11, 00
ERROR: MUL50X50 mismatch:
  a,b,c,d =     8340178993678043,     7722615070485185, 64407991987122015003084471664640,     4341871219628315
  xhi,xlo =        3491564241893,  6138244563698243868+
  yhi,ylo =        3491564241893,  6138244563698243867
  a,b mod 4, a^b = 11, 01, 10
ERROR: MUL50X50 mismatch:
  a,b,c,d =     6347818202160477,     7661859343595957, 48636090203671744523774971609088,    -4343683409217599
  xhi,xlo =        2636567733001,  1716031374246311872-
  yhi,ylo =        2636567733001,  1716031374246311873
  a,b mod 4, a^b = 01, 01, 00
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     4155331468855427,     8804079469637269, 36583868474487739815758028341248,     1893996263767615
  xhi,xlo =        1983215483898,  3289521724244229696+
  yhi,ylo =        1983215483898,  3289521724244229695
  a,b mod 4, a^b = 11, 01, 10
ERROR: MUL50X50 mismatch:
  a,b,c,d =     5314352150771271,     4096097024809011, 21768102033561572617061872435200,     -821416751712219
  xhi,xlo =        1180051175783, 12528192746593007652-
  yhi,ylo =        1180051175783, 12528192746593007653
  a,b mod 4, a^b = 11, 11, 00
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     7382608150332025,     7013861826885735, 51780593488469298141367048863744,    -4338555062700369
  xhi,xlo =        2807031597639,  6615952897171928752+
  yhi,ylo =        2807031597639,  6615952897171928751
  a,b mod 4, a^b = 01, 11, 10
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     6388136281807843,     7188581617218155, 45921639043688195935724999016448,     -348945078026783
  xhi,xlo =        2489417040763,  1756054909596466656-
  yhi,ylo =        2489417040763,  1756054909596466657
  a,b mod 4, a^b = 11, 11, 00
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     7742864853150641,     8838436707901963, 68434820942410571099029838495744,    -3472508339887461
  xhi,xlo =        3709859076970, 16632824515166724764+
  yhi,ylo =        3709859076970, 16632824515166724763
  a,b mod 4, a^b = 01, 11, 10
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     6795002443144771,     8301158853088421, 56406394687668667407574084091904,    -1549241117295313
  xhi,xlo =        3057796783122,  4781273563150171440+
  yhi,ylo =        3057796783122,  4781273563150171439
  a,b mod 4, a^b = 11, 01, 10
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     8798153817998953,     8297712816479939, 73004553696971821920736753221632,     -992798505717765
  xhi,xlo =        3957584785979, 10618495122833911804+
  yhi,ylo =        3957584785979, 10618495122833911803
  a,b mod 4, a^b = 01, 11, 10
BORROW
ERROR: MUL50X50 mismatch:
  a,b,c,d =     5958159304844551,     7127957962038291, 42469509056059243322829799686144,     3027877065016197
  xhi,xlo =        2302276699148, 11289048543255479172-
  yhi,ylo =        2302276699148, 11289048543255479173
  a,b mod 4, a^b = 11, 11, 00

...So we know if LSB of result is wrong if LSB(a&b) != LSB(a*b), but that doesn't tell us whether is the result is off by +1 or -1. The sign of the difference is given by the next-to-least-significant bit of a^b: if that = 0 the numeric product (x) is one lower than the exact (y); if it = 1, the numeric product is one higher. Sweet!

*/

