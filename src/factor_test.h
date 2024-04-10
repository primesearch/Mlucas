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

/*** This file is designed to be directly included by factor.c ***/

#undef	USE_FMADD	// Need to add 100-bit modpow routines before enabling this for build of this file

/*
Perform factor self-tests.
*/
int test_fac()
{
#ifndef __CUDA_ARCH__	// test_fac() only relevant in CPU-build mode
	/*...time-related stuff	*/
	clock_t clock1, clock2;
	double tdiff;

	#include "fac_test_dat64.h"		// 63,64 and 65-bit known-factor (p,q) pairs, with the 65-bit having hi bit of q implied
	#include "fac_test_dat96.h"
	#include "fac_test_dat128.h"
	#include "fac_test_dat192.h"	// 160 and 192-bit known-factor pairs
	#include "fac_test_dat256.h"

	uint32	ntest63,ntest64,ntest65,ntest96,ntest128,ntest128x2,ntest160,ntest192,ntest256;
	uint64	p64,q64,res64;
	uint96	q96,res96;
	uint128 p128,two_p128,pinv128,q128,x128,res128;
	uint192 p192,two_p192,q192,x192,y192,res192;
	uint256 p256,two_p256,q256,x256,res256;
	const uint64 two64mod60 = 16;
	uint64 two64modp, k, karr[64];	// max. of 64 k's per batch-modpow for now
	uint32 i,j,l;
	double dbl,rnd;
	uint32 pm60,km60;
	uint64 hi64,lo64;
#if defined(FAC_DEBUG) && (defined(P2WORD) || defined(P3WORD) || defined(P4WORD))
	uint32 i2,i3,i4,ii,jj;
#endif
	char cbuf0[STR_MAX_LEN], cbuf1[STR_MAX_LEN], cbuf2[STR_MAX_LEN];
	uint64 *p = 0x0, *q = 0x0, *q2 = 0x0, *two_p = 0x0, *u64_arr = 0x0;

/****12/12/05: check p%4==3 / q%60 correlations: *****/
	uint32 pqmod60arr[60][60] = {{0}};
/*****************************************************/

	/* 3/29/2006: by way of testing the 256-bit factoring routines, find all 70-digit base-2
	probable primes occurring in the first few thousand digits of Pi: */
#if TEST_256
	uint64 tmp64;
	/* Here are the first 10000 or so digits of Pi: */
	char Pi[10240] = "";
	const uint32 primelen = 64;

	/* To avoid compiler limits on character literal length, copy digits into Pi[] in 1Kbyte chunks: */
	strcpy(cbuf0, "314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912983367336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958537105079227968925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825334468503526193118817101000313783875288658753320838142061717766914730359825349042875546873115956286388235378759375195778185778053217122680661300192787661119590921642019893809525720106548586327");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "886593615338182796823030195203530185296899577362259941389124972177528347913151557485724245415069595082953311686172785588907509838175463746493931925506040092770167113900984882401285836160356370766010471018194295559619894676783744944825537977472684710404753464620804668425906949129331367702898915210475216205696602405803815019351125338243003558764024749647326391419927260426992279678235478163600934172164121992458631503028618297455570674983850549458858692699569092721079750930295532116534498720275596023648066549911988183479775356636980742654252786255181841757467289097777279380008164706001614524919217321721477235014144197356854816136115735255213347574184946843852332390739414333454776241686251898356948556209921922218427255025425688767179049460165346680498862723279178608578438382796797668145410095388378636095068006422512520511739298489608412848862694560424196528502221066118630674427862203919494504712371378696095636437191728746776465757396241389086583264599581339047802759009946576407895126946839835259570982582262052248");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "940772671947826848260147699090264013639443745530506820349625245174939965143142980919065925093722169646151570985838741059788595977297549893016175392846813826868386894277415599185592524595395943104997252468084598727364469584865383673622262609912460805124388439045124413654976278079771569143599770012961608944169486855584840635342207222582848864815845602850601684273945226746767889525213852254995466672782398645659611635488623057745649803559363456817432411251507606947945109659609402522887971089314566913686722874894056010150330861792868092087476091782493858900971490967598526136554978189312978482168299894872265880485756401427047755513237964145152374623436454285844479526586782105114135473573952311342716610213596953623144295248493718711014576540359027993440374200731057853906219838744780847848968332144571386875194350643021845319104848100537061468067491927819119793995206141966342875444064374512371819217999839101591956181467514269123974894090718649423196156794520809514655022523160388193014209376213785595663893778708303906");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "979207734672218256259966150142150306803844773454920260541466592520149744285073251866600213243408819071048633173464965145390579626856100550810665879699816357473638405257145910289706414011097120628043903975951567715770042033786993600723055876317635942187312514712053292819182618612586732157919841484882916447060957527069572209175671167229109816909152801735067127485832228718352093539657251210835791513698820914442100675103346711031412671113699086585163983150197016515116851714376576183515565088490998985998238734552833163550764791853589322618548963213293308985706420467525907091548141654985946163718027098199430992448895757128289059232332609729971208443357326548938239119325974636673058360414281388303203824903758985243744170291327656180937734440307074692112019130203303801976211011004492932151608424448596376698389522868478312355265821314495768572624334418930396864262434107732269780280731891544110104468232527162010526522721116603966655730925471105578537634668206531098965269186205647693125705863566201855810072936065987648");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "611791045334885034611365768675324944166803962657978771855608455296541266540853061434443185867697514566140680070023787765913440171274947042056223053899456131407112700040785473326993908145466464588079727082668306343285878569830523580893306575740679545716377525420211495576158140025012622859413021647155097925923099079654737612551765675135751782966645477917450112996148903046399471329621073404375189573596145890193897131117904297828564750320319869151402870808599048010941214722131794764777262241425485454033215718530614228813758504306332175182979866223717215916077166925474873898665494945011465406284336639379003976926567214638530673609657120918076383271664162748888007869256029022847210403172118608204190004229661711963779213375751149595015660496318629472654736425230817703675159067350235072835405670403867435136222247715891504953098444893330963408780769325993978054193414473774418426312986080998886874132604721569516239658645730216315981931951673538129741677294786724229246543668009806769282382806899640048243540370141631496");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "589794092432378969070697794223625082216889573837986230015937764716512289357860158816175578297352334460428151262720373431465319777741603199066554187639792933441952154134189948544473456738316249934191318148092777710386387734317720754565453220777092120190516609628049092636019759882816133231666365286193266863360627356763035447762803504507772355471058595487027908143562401451718062464362679456127531813407833033625423278394497538243720583531147711992606381334677687969597030983391307710987040859133746414428227726346594704745878477872019277152807317679077071572134447306057007334924369311383504931631284042512192565179806941135280131470130478164378851852909285452011658393419656213491434159562586586557055269049652098580338507224264829397285847831630577775606888764462482468579260395352773480304802900587607582510474709164396136267604492562742042083208566119062545433721315359584506877246029016187667952406163425225771954291629919306455377991403734043287526288896399587947572917464263574552540790914513571113694109119393251910");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "760208252026187985318877058429725916778131496990090192116971737278476847268608490033770242429165130050051683233643503895170298939223345172201381280696501178440874519601212285993716231301711444846409038906449544400619869075485160263275052983491874078668088183385102283345085048608250393021332197155184306354550076682829493041377655279397517546139539846833936383047461199665385815384205685338621867252334028308711232827892125077126294632295639898989358211674562701021835646220134967151881909730381198004973407239610368540664319395097901906996395524530054505806855019567302292191393391856803449039820595510022635353619204199474553859381023439554495977837790237421617271117236434354394782218185286240851400666044332588856986705431547069657474585503323233421073015459405165537906866273337995851156257843229882737231989875714159578111963583300594087306812160287649628674460477464915995054973742562690104903778198683593814657412680492564879855614537234786733039046883834363465537949864192705638729317487233208376011230299113679386");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "270894387993620162951541337142489283072201269014754668476535761647737946752004907571555278196536213239264061601363581559074220202031872776052772190055614842555187925303435139844253223415762336106425063904975008656271095359194658975141310348227693062474353632569160781547818115284366795706110861533150445212747392454494542368288606134084148637767009612071512491404302725386076482363414334623518975766452164137679690314950191085759844239198629164219399490723623464684411739403265918404437805133389452574239950829659122850855582157250310712570126683024029295252201187267675622041542051618416348475651699981161410100299607838690929160302884002691041407928862150784245167090870006992821206604183718065355672525325675328612910424877618258297651579598470356222629348600341587229805349896502262917487882027342092222453398562647669149055628425039127577102840279980663658254889264880254566101729670266407655904290994568150652653053718294127033693137851786090407086671149655834343476933857817113864558736781230145876871266034891390956");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "200993936103102916161528813843790990423174733639480457593149314052976347574811935670911013775172100803155902485309066920376719220332290943346768514221447737939375170344366199104033751117354719185504644902636551281622882446257591633303910722538374218214088350865739177150968288747826569959957449066175834413752239709683408005355984917541738188399944697486762655165827658483588453142775687900290951702835297163445621296404352311760066510124120065975585127617858382920419748442360800719304576189323492292796501987518721272675079812554709589045563579212210333466974992356302549478024901141952123828153091140790738602515227429958180724716259166854513331239480494707911915326734302824418604142636395480004480026704962482017928964766975831832713142517029692348896276684403232609275249603579964692565049368183609003238092934595889706953653494060340216654437558900456328822505452556405644824651518754711962184439658253375438856909411303150952617937800297412076651479394259029896959469955657612186561967337862362561252163208628692221");
	strcat(Pi, cbuf0);
	strcpy(cbuf0, "032748892186543648022967807057656151446320469279068212073883778142335628236089632080682224680122482611771858963814091839036736722208883215137556003727983940041529700287830766709444745601345564172543709069793961225714298946715435784687886144458123145935719849225284716050492212424701412147805734551050080190869960330276347870810817545011930714122339086639383395294257869050764310063835198343893415961318543475464955697810382930971646514384070070736041123735998434522516105070270562352660127648483084076118301305279320542746286540360367453286510570658748822569815793678976697422057505968344086973502014102067235850200724522563265134105592401902742162484391403599895353945909440704691209140938700126456001623742880210927645793106579229552498872758461012648369998922569596881592056001016552563756785667227966198857827948488558343975187445455129656344348039664205579829368043522027709842942325330225763418070394769941597915945300697521482933665556615678736400536665641654732170439035213295435291694145990416087532018683793702348");
	strcat(Pi, cbuf0);

//	printf("First %d digits of Pi = %s\n", strlen(Pi), Pi);

	strcpy(cbuf0,"\0");
	for(i = 0; i < 10240; i++)
	{
		// If less than (primelen) digits remaining, quit:
		if(strlen(Pi+i) < primelen)
		{
			break;
		}

		strncpy(cbuf0, Pi+i, primelen);
		// Make sure surrent digit string ends after (primelen) digits:
		cbuf0[primelen] = '\0';
		p256 = convert_base10_char_uint256(cbuf0);
		/* Only test odds: */
		// Workaround for an MSVC compiler bug:
		//printf("p.d0 = %s\n", &cbuf1[convert_uint64_base10_char(cbuf1, p256.d0)]);
		tmp64 = p256.d0 & (uint64)0x0000000000000001ull;
		if(tmp64 == 0)
		{
			continue;
		}
		SUB256(p256,ONE256,q256);
		x256 = twopmodq256(q256,p256);
		if(CMPEQ256(x256,ONE256))
			printf("%s is a base-2 probable prime\n", cbuf0);
	}
#endif

#ifdef FACTOR_STANDALONE
	printf("Mfactor build flags:\n");

	/* TRYQ: */
	#ifndef TRYQ
		/* This flag is required: */
		ASSERT(0,"TRYQ not defined!");
	#else
		i = TRYQ;
		printf("TRYQ = %u\n", i);
	#endif

	/* THREE_OP128: */
	#ifndef THREE_OP128
	//	printf("THREE_OP128 not defined\n");
	#elif(THREE_OP128 == 0)
		printf("THREE_OP128 = FALSE\n");
	#else
		i = THREE_OP128;
		printf("THREE_OP128 = %u\n", i);
		/* iF NONZERO, Must = 1 : */
		ASSERT((THREE_OP128 == 1),"THREE_OP128 Must = 0 or 1!");
		/* Only relevant for TRYQ = 4 or 8: */
		#if(TRYQ != 4 && TRYQ != 8)
			#error	THREE_OP128 Only relevant for TRYQ = 4 or 8!
		#endif
		/* Only relevant for factoring up to 128 bits: */
		#if(defined(P3WORD) || defined(P4WORD))
			#error	THREE_OP128 Only relevant for factoring up to 128 bits!
		#endif
		/* Only relevant if using fully 128-bit modmul routines: */
		#if USE_128x96 != 0
			#error	THREE_OP128 Only relevant if using fully 128-bit modmul routines - undef USE_128x96 or set = 0!
		#endif
	#endif

	/* NUM_SIEVING_PRIME: */
	#ifndef NUM_SIEVING_PRIME
		/* This flag is required: */
		ASSERT(0,"NUM_SIEVING_PRIME not defined!");
	#else
		i = NUM_SIEVING_PRIME;
		printf("NUM_SIEVING_PRIME = %u\n", i);
	#endif

	/* TF_CLASSES: */
	#ifndef TF_CLASSES
		/* This flag is required: */
		ASSERT(0,"TF_CLASSES not defined!");
	#else
		i = TF_CLASSES;
		printf("TF_CLASSES = %u\n", i);
	#endif

	/* MUL_LOHI64_SUBROUTINE: */
	#ifndef MUL_LOHI64_SUBROUTINE
	//	printf("MUL_LOHI64_SUBROUTINE not defined\n");
	#else
		printf("MUL_LOHI64_SUBROUTINE = true\n");
	#endif

	/* MULH64_FAST: */
	#ifndef MULH64_FAST
	//	printf("MULH64_FAST not defined\n");
	#else
		printf("MULH64_FAST = true\n");
	#endif

	/* USE_FLOAT: */
	#ifndef USE_FLOAT
	//	printf("USE_FLOAT not defined\n");
	#else
		printf("USE_FLOAT = true\n");
	#endif

	/* USE_FMADD: */
	#ifndef USE_FMADD
	//	printf("USE_FMADD not defined\n");
	#else
		printf("USE_FMADD = true\n");
	#endif

	/* FACTOR_STANDALONE: */
	#ifndef FACTOR_STANDALONE
	//	printf("FACTOR_STANDALONE not defined\n");
	#else
		printf("FACTOR_STANDALONE = true\n");
	#endif

	/* FAC_DEBUG: */
	#ifndef FAC_DEBUG
	//	printf("FAC_DEBUG not defined\n");
	#else
		printf("FAC_DEBUG = true\n");
	#endif

	/* DBG_SIEVE: */
	#ifndef DBG_SIEVE
	//	printf("DBG_SIEVE not defined\n");
	#else
		i = DBG_SIEVE;
		printf("DBG_SIEVE = true\n");
	#endif

	/* NOBRANCH: */
	#ifndef NOBRANCH
	//	printf("NOBRANCH not defined\n");
	#else
		printf("NOBRANCH = true\n");
	#endif

	/* QUIT_WHEN_FACTOR_FOUND: */
	#ifndef QUIT_WHEN_FACTOR_FOUND
	//	printf("QUIT_WHEN_FACTOR_FOUND not defined\n");
	#else
		printf("QUIT_WHEN_FACTOR_FOUND = true\n");
	#endif

	/* USE_65BIT: */
	#ifndef USE_65BIT
	//	printf("USE_65BIT not defined\n");
	#else
		printf("USE_65BIT = true\n");
	#endif

	/* USE_128x96: */
	#ifndef USE_128x96
	//	printf("USE_128x96 not defined\n");
	#else
		i = USE_128x96;
		printf("USE_128x96 = %u\n", i);
		ASSERT(i <= 2,"Only USE_128x96 = 0-2 are recognized values!\n");
		/* Only relevant for factoring up to 128 bits: */
		#if(defined(P3WORD) || defined(P4WORD))
			#warning USE_128x96 Only relevant for factoring up to 128 bits!
		#endif
	#endif

	/* P2WORD: */
	#ifndef P2WORD
	//	printf("P2WORD not defined\n");
	#else
		printf("P2WORD = true\n");
		#ifndef USE_128x96
			printf("    USE_128x96 not defined\n");
		#else
			i = USE_128x96;
			printf("    USE_128x96 = %u\n", i);
			ASSERT(i <= 2,"Only USE_128x96 = 0-2 are recognized values!\n");
		#endif

	#endif

	/* P3WORD: */
	#ifndef P3WORD
	//	printf("P3WORD not defined\n");
	#else
		printf("P3WORD = true\n");

		#ifndef PIPELINE_MUL192
			printf("    PIPELINE_MUL192 not defined\n");
		#else
			printf("    PIPELINE_MUL192 = %u\n", PIPELINE_MUL192);
		#endif

	#endif

	/* P4WORD: */
	#ifndef P4WORD
	//	printf("P4WORD not defined\n");
	#else
		printf("P4WORD = true\n");
	#endif
#endif

#ifdef FACTOR_STANDALONE
	printf("Mfactor self-tests:\n");
#endif

	/***********************************************/
	/* Basic tests of the mi64 vector-int package: */
	/***********************************************/

	j = 607;
	l = (j + 63)>>6;	// #64-bit words needed
	q     = (uint64 *)calloc(l, sizeof(uint64));
	q2    = (uint64 *)calloc(l, sizeof(uint64));
	mi64_nega(q,q,l);
	ASSERT(mi64_iszero(q,l), "mi64 -0 == 0 check fails!");
	q[0] = 1;	mi64_nega(q,q,l);
	mi64_add_scalar(q,1,q,l);
	ASSERT(mi64_iszero(q,l), "mi64 -1 + 1 == 0 check fails!");

	// Sep 2015 Bugfix: Hit case with len = 3 and these addends, which give a ripple carry into the top word:
	q[0] =  6216518070457578443ull;	q2[0] = 12230226003251973173ull;
	q[1] = 16881888488052985758ull;	q2[1] =  1564855585656565857ull;
	q[2] =       65307107850795ull;	q2[2] =           2051081684ull;
	mi64_add(q,q2,q,3);
	ASSERT(q[0] == 0ull && q[1] == 0ull && q[2] == 65309158932480ull, "Sep 2015 mi64_add bugfix test fails!");

	/* Init the RNG: */
	rng_isaac_init(TRUE);
	for(i = 0; i < l; i++)
	{
		q [i] = rng_isaac_rand();
		q2[i] = q[i];
	}
	mi64_nega(q,q,l);
	mi64_negl(q2,q2,l);
	mi64_add_scalar(q2,1,q2,l);
	ASSERT(mi64_cmp_eq(q,q2,l), "mi64 -q == ~q+1 check fails!");
	free((void*)q);	free((void*)q2);
	q = q2 = 0x0;

	/* 06/18/2012: Test out new, streamlined version of mi64_mul_vector_hi_qmmp routine with
	Ex.: q = 2.k.M(127) + 1 with k = 7143819210136784550, i.e.
	q = 2430915709680614116949754105299803650411408301848040235701 ;
	y =  915005412744957807408012591600653057424688130286064771258 = y0 + 2^64*y1 + 2^128*y2,
	with y0 = 2294959606785646778; y1 = 10167084567166165345; y2 = 2688959234133783535 .

	Exact result:

		UMULH_192(q,y) = 354351598245602020483095922210514413558224553895064094733 = u0 + 2^64*u1 + 2^128*u2,
		with u0 = 141525868296128525, u1 = 4269430960237156763, u2 = 1041345754856384950 .
	*/
	k = 7143819210136784550ull;	p64 = 127;
	p192.d0 = 2294959606785646778ull; p192.d1 = 10167084567166165345ull; p192.d2 = 2688959234133783535ull;
	mi64_mul_vector_hi_qmmp((uint64*)&p192, p64, k, (uint64*)&q192, 192);
	ASSERT(q192.d0 == 141525868296128525ull && q192.d1 == 4269430960237156763ull && q192.d2 == 1041345754856384950ull, "mi64_mul_vector_hi_qmmp test fails!");

	/* 09/30/2015: Adapt above to test Fermat-factor analog of above, mi64_mul_vector_hi_qferm:
	Ex.: q = 2.k.2^128 + 1; k = 3571909605068392275, i.e.
	q = 2430915709680614116949754105299803650425695940268313804801 ;
	y =  915005412744957807408012591600653057424688130286064771258 = y0 + 2^64*y1 + 2^128*y2,
	with y0 = 2294959606785646778; y1 = 10167084567166165345; y2 = 2688959234133783535 .

	Exact result:

		UMULH_192(q,y) = 354351598245602020483095922210514413560307245404776864634 = u0 + 2^64*u1 + 2^128*u2,
		with u0 = 4269430960237156763, u1 = 1041345754856384950, u2 = 1041345754856384950 .
	*/
	k = 3571909605068392275ull;	p64 = 128;
	p192.d0 = 2294959606785646778ull; p192.d1 = 10167084567166165345ull; p192.d2 = 2688959234133783535ull;
	mi64_mul_vector_hi_qferm((uint64*)&p192, p64, k, (uint64*)&q192, 192);
	ASSERT(q192.d0 == 2224217378008898426ull && q192.d1 == 4269430960237156763ull && q192.d2 == 1041345754856384950ull, "mi64_mul_vector_hi_qferm test fails!");

	// Apr 2015: mi64_div bug debug - 0-pad both inputs to yield a length-4 mi64 array:
	// Use 2^256 as a template for our 0-padding, but use 1 less leading 0 because convert_base10_char_mi64
	// assumes worst-case value of 9999... for the leading digits in precomuting the #words for the calloc:
	//                    2^256 = 115792089237316195423570985008687907853269984665640564039457584007913129639936:
	// Feb 2020: Chnages to length-setting logic in convert_base10_char_mi64 mean we must init i,j = 0 prior to calling that function:
	i = 0; p = convert_base10_char_mi64( "00000000000000000000000000000000000000364131549958466711308970009901738230041", &i);
	ASSERT(mi64_getlen(p, i) == 3 && i == 4,"Bad p-length(s) in Apr2015 mi64_div test!");
	j = 0; q = convert_base10_char_mi64( "00000000000000000000000000000000000000000000000000000000019437941122649628431", &j);
	ASSERT(mi64_getlen(q, j) == 2 && j == 4,"Bad q-length(s) in Apr2015 mi64_div test!");
	q2      = (uint64 *)calloc(4, sizeof(uint64));	// for quotient
	u64_arr = (uint64 *)calloc(4, sizeof(uint64));	// for remainder
	mi64_div(p,q,i,i,q2,u64_arr);
	ASSERT(mi64_getlen(     q2, i) == 2 && q2[1] == 1 && q2[0] ==   286286737571717471ull, "bad quotient!");
	ASSERT(mi64_getlen(u64_arr, i) == 1 &&          u64_arr[0] ==   618006351061617544ull, "bad remainder!");
	fprintf(stderr,"Apr2015 mi64_div quicktest passes.\n");
	free((void*)p); free((void*)q); free((void*)q2); free((void*)u64_arr);
	p = 0x0; q = 0x0; q2 = 0x0; u64_arr = 0x0;

	/* 01/09/2008: mi64_div bug debug: */
	i = 0; p = convert_base10_char_mi64("531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127", &i);
	two_p = (uint64 *)calloc(i, sizeof(uint64));
	mi64_add(p,p,two_p,i);
	j = 0; q = convert_base10_char_mi64("4969289881134175801642878989330437804491760137935869781219375395913301677808943323410612629818326630668131744420258226244511522022525093242408710254941677603671849301746980479735516135243111", &j);
	ASSERT(i==j,"0");
	q2      = (uint64 *)calloc(i, sizeof(uint64));
	u64_arr = (uint64 *)calloc(i, sizeof(uint64));
	mi64_div(q,two_p,i,i,q2,u64_arr);
	ASSERT(mi64_getlen(q2, i) == 1 , "k must be 64-bit!");
	ASSERT(q2[0] == 4677965, "k != expected value of 9355930!");
	if(!mi64_cmp_eq_scalar(u64_arr, 1ull, i)) {		// Remainder = 1
		fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_mi64_base10_char(cbuf0, p, i, 0)],
					&cbuf1[convert_mi64_base10_char(cbuf1, q, i, 0)],
					&cbuf2[convert_mi64_base10_char(cbuf2, u64_arr, i, 0)]);
		ASSERT(0,"0");
	} else {
		fprintf(stderr,"mi64_div quicktest passes.\n");
	}
	free((void*)p); free((void*)q); free((void*)q2); free((void*)two_p); free((void*)u64_arr);
	p = 0x0; q = 0x0; q2 = 0x0; two_p = 0x0; u64_arr = 0x0;

	// 4/19/2012: mi64 gives bad result for e.g. the following M607 printout when built optimized in GCC...after adding
	//            this test code, however, even the opt-build passes - weird.
	// 11/2013: Under Debian/gcc4.6 this barfs at any opt-level > 0 ... traced it to the x86_64 inline asm version of mi64_add,
	//          but the asm is not at fault. No choice but to revert to C version of mi64_add() until further notice, though.
	// 4/21/2015: Added test of known-factor-of-MM31, after use of mi64_twopmodq to check if 2^-p == 1 (mod q) in a TF run failed.
// Base-2 PRP test of M127:
	clock1 = clock();
	j = convert_base10_char_uint64("127");
	i = (j + 63)>>6;	// #64-bit words needed
	p     = (uint64 *)calloc(i, sizeof(uint64));
	q     = (uint64 *)calloc(i, sizeof(uint64));
	p[0] = 1;	mi64_shl(p,p,j,i);	// 2^n
	mi64_sub_scalar(p,1,p,i);	// p = 2^n - 1;
	convert_mi64_base10_char(cbuf0, p, i, 0);
	ASSERT(STREQ(cbuf0, "170141183460469231731687303715884105727"), "M127 string-conversion test failed!");
	mi64_set_eq    (q, p, i);
	mi64_sub_scalar(q ,1ull,q ,i);	// q = p-1
	j = mi64_twopmodq(q, i, 0, p, i, 0x0);
	ASSERT(j == 1, "M127 base-2 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-2 PRP test of M127 passed: Time =%s\n",get_time_str(tdiff));

// Test of 3rd & 4th known factors-of-MM31 use same array dims as MM127 PRP test:
/* Same case using 2-word-specialized routines:
	p128.d1 =  0;	p128.d0 = 2147483647;
	q128.d1 = 13;	q128.d0 = 2749942686469094193ull;
	res128 = twopmodq128(p128, q128);
	exit(0);
*/
	uint32 lenP = 1, lenQ = 2;
	q2 = (uint64 *)calloc(lenQ, sizeof(uint64));	// Will store 2^-p (mod q) residue
	p[0] = 2147483647;
	q[1] = mi64_mul_scalar( p, 2*56474845800ull, q, lenP);
	q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
	if(mi64_twopmodq(p, lenP, 56474845800ull, q, lenQ, q2) != 1) {
		printf("ERROR: res = %s != 1\n", &cbuf[convert_mi64_base10_char(cbuf0, q2, lenQ, 0)]);
		ASSERT(0, "MM31 known-factor (k = 56474845800) test failed!");
	}
	q[1] = mi64_mul_scalar( p, 2*41448832329225ull, q, lenP);
	q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
	if(mi64_twopmodq(p, lenP, 41448832329225ull, q, lenQ, q2) != 1) {
		printf("ERROR: res = %s != 1\n", &cbuf[convert_mi64_base10_char(cbuf0, q2, lenQ, 0)]);
		ASSERT(0, "MM31 known-factor (k = 41448832329225) test failed!");
	}
	free((void*)p);	free((void*)q);	free((void*)q2);	p = q = q2 = 0x0;

// Base-2 PRP test of M607:
	j = convert_base10_char_uint64("607");
	i = (j + 63)>>6;	// #64-bit words needed
	p     = (uint64 *)calloc(i, sizeof(uint64));
	q     = (uint64 *)calloc(i, sizeof(uint64));
	p[0] = 1;	mi64_shl(p,p,j,i);	// 2^n
	mi64_sub_scalar(p,1,p,i);	// p = 2^n - 1;
	convert_mi64_base10_char(cbuf0, p, i, 0);
	ASSERT(STREQ(cbuf0, "531137992816767098689588206552468627329593117727031923199444138200403559860852242739162502265229285668889329486246501015346579337652707239409519978766587351943831270835393219031728127"), "M607 string-conversion test failed!");
	mi64_set_eq    (q, p, i);
	mi64_sub_scalar(q ,1ull,q ,i);	// q = p-1
	clock1 = clock();
	j = mi64_twopmodq(q, i, 0, p, i, 0x0);
	ASSERT(j == 1, "M607 base-2 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-2 PRP test of M607 passed: Time =%s\n",get_time_str(tdiff));
	// Try the general-base PRP routine on the same number:
	clock1 = clock();
	j = mi64_pprimeF(p, 3, i);
	ASSERT(j == 1, "M607 base-3 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-3 PRP test of M607 passed: Time =%s\n",get_time_str(tdiff));
	free((void*)p);	free((void*)q);	p = q = 0x0;

// Base-2 PRP test of M4423:
	j = convert_base10_char_uint64("4423");
	i = (j + 63)>>6;	// #64-bit words needed
	p     = (uint64 *)calloc(i, sizeof(uint64));
	q     = (uint64 *)calloc(i, sizeof(uint64));
	p[0] = 1;	mi64_shl(p,p,j,i);	// 2^n
	mi64_sub_scalar(p,1,p,i);	// p = 2^n - 1;
	mi64_set_eq    (q, p, i);
	mi64_sub_scalar(q ,1ull,q ,i);	// q = p-1
	clock1 = clock();
	j = mi64_twopmodq(q, i, 0, p, i, 0x0);
	ASSERT(j == 1, "M4423 base-2 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-2 PRP test of M4423 passed: Time =%s\n",get_time_str(tdiff));
	// Try the general-base PRP routine on the same number:
	clock1 = clock();
	j = mi64_pprimeF(p, 3, i);
	ASSERT(j == 1, "M4423 base-3 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-3 PRP test of M4423 passed: Time =%s\n",get_time_str(tdiff));
	free((void*)p);	free((void*)q);	p = q = 0x0;

// Base-3 PRP test of the Mersenne cofactor M7331/458072843161:
#if 0
	j = convert_base10_char_uint64("7331");
	i = (j + 63)>>6;	// #64-bit words needed
	p     = (uint64 *)calloc(i, sizeof(uint64));
	q     = (uint64 *)calloc(i, sizeof(uint64));
	p[0] = 1;	mi64_shl(p,p,j,i);	// 2^n
	mi64_sub_scalar(p,1,p,i);	// p = 2^n - 1; next we p /= 458072843161 :
ASSERT(0 == mi64_div_by_scalar64(p, 458072843161ull, i, p), "M7331/458072843161 divisibility test fails!");
	mi64_set_eq    (q, p, i);
	mi64_sub_scalar(q ,1ull,q ,i);	// q = p-1
	clock1 = clock();
	j = mi64_twopmodq(q, i, 0, p, i, 0x0);
	ASSERT(j == 1, "M7331 cofactor base-2 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-2 PRP test of M7331 cofactor passed: Time =%s\n",get_time_str(tdiff));
	// Try the general-base PRP routine on the same number:
	clock1 = clock();
	j = mi64_pprimeF(p, 3, i);
	ASSERT(j == 1, "M7331 cofactor base-3 PRP test failed!");
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-3 PRP test of M7331 cofactor passed: Time =%s\n",get_time_str(tdiff));
	free((void*)p);	free((void*)q);	p = q = 0x0;
#endif

// Base-2 PRP test of M11213:
#if 0
	j = convert_base10_char_uint64("11213");
	i = (j + 63)>>6;	// #64-bit words needed
	p     = (uint64 *)calloc(i, sizeof(uint64));
	q     = (uint64 *)calloc(i, sizeof(uint64));
	p[0] = 1;	mi64_shl(p,p,j,i);	// 2^n
	mi64_sub_scalar(p,1,p,i);	// p = 2^n - 1;
	mi64_set_eq    (q, p, i);
	mi64_sub_scalar(q ,1ull,q ,i);	// q = p-1
	clock1 = clock();
	j = mi64_twopmodq(q, i, 0, p, i, 0x0);
	ASSERT(j == 1, "M11213 base-2 PRP test failed!");
	free((void*)p);	free((void*)q);	p = q = 0x0;
	clock2 = clock();	tdiff = (double)(clock2 - clock1);	clock1 = clock2;
	printf	("Base-2 PRP test of M11213 passed: Time =%s\n",get_time_str(tdiff));
#endif
	/* 4/20/2012: Did a test TF run of mm607 to k = 3e7, printing debug data about sample factor candidates
	for every [2^14]th q which passes the first-10000-small-primes sieve (primes <= 104743). All 96 such
	sample-q are flagged by my mi64 modpow routine as composite (via base-2 Fermat PRP test). Rechecked
	using Pari/GP to both test for small factors (none showed any factors < sievelimit, but increasing the
	the sieving depth to 2^32 factors 50 of the 96 candidates) ... more importantly Pari indicates one of the
	96 candidates, with k = 28115877, is in fact prime. That warrants getting added to the self-test suite:
	*/
	j = 607;
	i = (j + 63)>>6;	// #64-bit words needed
	q     = (uint64 *)calloc(i, sizeof(uint64));
	q2    = (uint64 *)calloc(i, sizeof(uint64));
	q[0] = 1;	mi64_shl(q,q,j,i);	// 2^607
	mi64_sub_scalar(q,1,q,i);		// p = 2^607 - 1;
	// Mul by any scalar < 2^33 should have no carry out of the 10th 64-bit word
	ASSERT(0 == mi64_mul_scalar(q,2*28115877,q,i), "2.k.M607 (k = 28115877) illegal carryout on scalar-mul!");
	mi64_set_eq    (q2, q, i);		// q2 = q-1
	mi64_add_scalar(q ,1ull,q ,i);	// q = 2.k.p + 1
	convert_mi64_base10_char(cbuf0, q, i, 0);
	ASSERT(STREQ(cbuf0, "29866820952126214568806646392159603944715357116119498255498035716027095678819717544056871993402815945328710228895559628455719074056369970920495232704087963394016941839123205985860254232344759"), "q = 2.k.M607+1 (k = 28115877) string-conversion test failed!");
	ASSERT(mi64_twopmodq(q2, i, 0, q, i, 0x0) == 1, "q = 2.k.M607+1 (k = 28115877) base-2 PRP test failed!");
	free((void*)q);	free((void*)q2);
	q = q2 = 0x0;

	/* Test fast MM(p)-factor-candidate-specific hi-mul algorithm.
	NB: p = 1231 not an M-prime exponent, but 1279 too close to 1280 = 64*20 for our purposes.
	*/
	// 2nd multiplicand is just leading digits of Pi, sans decimal point:
	j = 0; q2 = convert_base10_char_mi64("3141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521", &j);
	ASSERT(j == 20, "vector lengths should be 20!");
	q     = (uint64 *)calloc(j, sizeof(uint64));	// output array
	mi64_mul_vector_hi_qmmp(q2, 1231, 60773088284ull, q, (j<<6));	// q = 2.k.M(p) + 1 with k = 60773088284
	convert_mi64_base10_char(cbuf0, q, j, 0);
	ASSERT(STREQ(cbuf0, "678299328487875406787553667584424766193319571425229812042632483796223090743976740829512533956144441574815272835626612961160454952708658437402700559999225654073147100413573556498251710301510338504761109128343850675314104893353603303495634850631971760134667616782442458276408663375682004856646999060481786800862572039635523841600325205075025327991817191734342347965082117753555537"), "mi64_mul_vector_hi_qmmp test failed!");
	free((void*)q);	free((void*)q2);
	q = q2 = 0x0;

	// 7/21/2005 bugfix test:
	q192.d2=506560280167ull; q192.d1=18446744073709551615ull; q192.d0=18446743060588991281ull;
	p192 = q192; p192.d0 -= 1;
	x192 = twopmodq192(p192,q192);
	ASSERT(CMPEQ192(x192, ONE192),"Bad twopmodq192 output");

#if 0
	/* 12/23/2008: Use this to help debug the mi64 powering routine: */
	j = mi64_twopmodq(&p192.d0, 3, 0, &q192.d0, 3, 0x0);
	if(j != 1) {
		printf("12/23/2008 mi64_twopmodq Test failed!\n");
	//	ASSERT(j == 1, "mi64_twopmodq != 1");
	//	exit(0);
	}
#endif

	/* 1/23/2004 bugfix test (cf. def. of MULH192 in imul_macro.h for details):

	X =                   14*2^128 +  3291757557782450881*2^64 +  3893270457587058239;
	Y = 14334090244500356821*2^128 + 14839649155872891980*2^64 + 12743646208794279936;
	! and w3 = 2846153632803221901, rather than the true w2 = 0 and w3 = 2846153632803221902.
	!
	! However, this kind of thing is only likely to happen in the one q*(qinv << shift) multiply
	! that we do just prior to entering the modmul-based powering loop. So we could define two
	! versions of the MULH routine: an exact one which will be used only for the initial mul,
	! and one which approximates the (presumably quasirandom) lower bits of the mul so as to get
	! the proper carry into the upper half with very high probability, which we use everywhere else.
	! ****TODO???****

	*/
	p192.d2=                  14ull; p192.d1= 3291757557782450881ull; p192.d0= 3893270457587058239ull;
	q192.d2=14334090244500356821ull; q192.d1=14839649155872891980ull; q192.d0=12743646208794279936ull;
	MULH192(p192,q192,x192);
	MULH192_TRUNC(p192,q192,0ull,y192);	// Expected value of 64-bit carry layer at top of low-half product = 0
	/* Reference value to compare to: */
	q192.d2=                  11ull; q192.d1=  320947345442520101ull; q192.d0= 2846153632803221902ull;
	ASSERT(CMPEQ192(x192, q192),"MULH192       fails!");
	ASSERT(CMPEQ192(y192, q192),"MULH192_TRUNC fails!");

	/* Count the # of test q's of the various sizes: */
	for(ntest63    = 0; fac63   [ntest63   ].p          != 0; ++ntest63   ){}
	for(ntest64    = 0; fac64   [ntest64   ].p          != 0; ++ntest64   ){}
	for(ntest65    = 0; fac65   [ntest65   ].p          != 0; ++ntest65   ){}
	for(ntest96    = 0; fac96   [ntest96   ].p          != 0; ++ntest96   ){}
	for(ntest128   = 0; fac128  [ntest128  ].p          != 0; ++ntest128  ){}
	for(ntest128x2 = 0; fac128x2[ntest128x2].plo        != 0; ++ntest128x2){}
	for(ntest160   = 0; fac160  [ntest160  ].p          != 0; ++ntest160  ){}
	for(ntest192   = 0; fac192  [ntest192  ].p          != 0; ++ntest192  ){}
	for(ntest256   = 0; STRNEQ(fac256[ntest256].p, "")		; ++ntest256  ){}

	// Sep 2015: Test 64-bit Fermat factors using the 64-bit modpow routines:
#ifdef FACTOR_STANDALONE
	printf("Testing 64-bit Fermat factors...\n");
#endif
	for(i = 0; ffac64[i].p != 0; i++)
	{
		// testFac struct uses [n,k]-pair format, where Fn = 2^2^n+1 is the Fermat number and q = k.2^(n+2)+1 the factor.
		// Here recast the factor as q = 2.k.2^n + 1 to leverage the sieve infrastructure developed
		// for Mersenne-number trial-factoring, where factors of Mp are of form q = 2.k.p + 1.
		// In the Fermat cae we let 2^n play the role of the Mersenne exponent p and generalize from there.
		p64 = 1ull << ffac64[i].p; k = ffac64[i].q << 1;	// Factors of Fn have form q = k.2^(n+2) + 1; n stored in .p, k in .q
		q64 = 2*k*p64 + 1;	// p64 now stores 2^n
		ASSERT(q64%(p64<<2)==1, "test_fac : q64 % 2^(n+2) != 1 !");
		pm60 = p64%60;
		km60 = k  %60;
		if(!CHECK_PKMOD60(&p64,1, km60, 0x0)) {
			fprintf(stderr,"Illegal (p,k) mod 60 pair: p,p mod 60, k,k mod 60 = %llu %4u %llu %4u\n",p64,pm60,k,km60);
			ASSERT(0,"0");
		}
		pm60 = p64%4620;
		km60 = k  %4620;
		if(!CHECK_PKMOD4620(&p64,1, km60, 0x0)) {
			fprintf(stderr,"Illegal (p,k) mod 4620 pair: p,p mod 4620, k,k mod 4620 = %llu %4u %llu %4u\n",p64,pm60,k,km60);
			ASSERT(0,"0");
		}
		res64 = twopmodq64(p64, q64);
		if(res64 != q64-1ull) {	// Nov 2021: fiddled twopmodq64() to return true-mod
			fprintf(stderr,"ERROR: twopmodq64(F%u, k = %llu) returns non-unity result %u\n",(uint32)ffac64[i].p,k, (uint32)res64);
			ASSERT(0,"0");
		}
	}

	// Test 128-bit factors using the 128-bit (for both exponent and k-multiplier) modpow routines:
#ifdef FACTOR_STANDALONE
	printf("Testing 128-bit Fermat factors...\n");
#endif
	for(i = 0; ffac128[i].p != 0; i++)
	{
		p64 = ffac128[i].p;	// Fermat index n
		p128.d1 = 0ull; p128.d0 = 1ull;	LSHIFT128(p128,p64,p128);	// p128 holds powering exponent 2^n
		q128.d1 = (uint64)ffac128[i].d1; q128.d0 = ffac128[i].d0;	// q128 holds k...
		LSHIFT128(q128,(p64+2),q128);	q128.d0 += 1ull;	// ...and now q128 holds q = k.2^(n+2) + 1 .

		/* This uses the generic 128-bit mod function to calculate q%(2*p): */
		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(res128,ONE128)) {
			fprintf(stderr,"ERROR: twopmodq128(F%u, %s ) returns non-unity result %s\n",(uint32)p64,
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, res128)]);
			ASSERT(0,"0");
		}
	}

	// Test 192-bit factors using the 192-bit (for both exponent and k-multiplier) modpow routines:
#ifdef FACTOR_STANDALONE
	printf("Testing 192-bit Fermat factors...\n");
#endif
	for(i = 0; ffac192[i].p != 0; i++)
	{
		p64 = ffac192[i].p;	// Fermat index n
		p192.d2 = p192.d1 = 0ull; p192.d0 = 1ull;	LSHIFT192(p192,p64,p192);	// p192 holds powering exponent 2^n
		q192.d2 = q192.d1 = 0ull; q192.d0 = ffac192[i].d0;	// q192 holds k (As of 2015, no k > 64-bit known)
		LSHIFT192(q192,(p64+2),q192);	q192.d0 += 1ull;	// ...and now q192 holds q = k.2^(n+2) + 1 .

		/* This uses the generic 192-bit mod function to calculate q%(2*p): */
		res192 = twopmodq192(p192, q192);
		if(!CMPEQ192(res192,ONE192)) {
			fprintf(stderr,"ERROR: twopmodq192(F%u, %s ) returns non-unity result %s\n",(uint32)p64,
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}
	}

	// Test 256-bit factors using the 256-bit (for both exponent and k-multiplier) modpow routines:
#ifdef FACTOR_STANDALONE
	printf("Testing 256-bit Fermat factors...\n");
#endif
	for(i = 0; ffac256[i].n != 0; i++)
	{
		p64 = ffac256[i].n;	// Fermat index n
		p256.d2 = p256.d1 = 0ull; p256.d0 = 1ull;	LSHIFT256(p256,p64,p256);	// p256 holds powering exponent 2^n
		q256.d2 = q256.d1 = 0ull; q256.d0 = ffac256[i].k;	// q256 holds k (As of 2015, no k > 64-bit known)
		LSHIFT256(q256,(p64+2),q256);	q256.d0 += 1ull;	// ...and now q256 holds q = k.2^(n+2) + 1 .

		/* This uses the generic 256-bit mod function to calculate q%(2*p): */
		res256 = twopmodq256(p256, q256);
		if(!CMPEQ256(res256,ONE256)) {
			fprintf(stderr,"ERROR: twopmodq256(F%u, %s ) returns non-unity result %s\n",(uint32)p64,
					&cbuf1[convert_uint256_base10_char(cbuf1, q256)],
					&cbuf2[convert_uint256_base10_char(cbuf2, res256)]);
			ASSERT(0,"0");
		}
	}

	// Test > 256-bit factors using the mi64 modpow routines:
#ifdef FACTOR_STANDALONE
	printf("Testing > 256-bit Fermat factors...\n");
#endif
	p  = (uint64 *)calloc(640, sizeof(uint64));	// ffacBig Fermat-factor array limited to Fn with n < 10000
	q  = (uint64 *)calloc(640, sizeof(uint64));
	q2 = (uint64 *)calloc(640, sizeof(uint64));
	for(i = 0; ffacBig[i].p != 0; i++)
	{
		j =         ffacBig[i].p;	// Fermat index n
		if(j > 1000) break;			// Tune this as desired to skip larger time-consuming cases
		l = (uint32)ffacBig[i].d1;	// Power of 2 appearing in factor q = k*2^l + 1
		ASSERT(l >= (j+2), "Power of 2 appearing in factor of Fn must be >= [n+2]!");
		k =         ffacBig[i].d0;	// Factor k; must be odd in this schema
		ASSERT(1ull == (k & 1ull), "k must be odd!");
		lenP = (j+63)>>6;	// Assume Fermat index increases as we traverse ffacBig array, thus this overwrites previous
		p[0] = 1ull;	p[lenP] = mi64_shl(p,p,j,lenP);	lenP += (p[lenP] != 0ull);	// case's p = (1 << j) array elements.
		lenQ = (l+63)>>6;
		q[0] = k;		q[lenQ] = mi64_shl(q,q,l,lenQ);	lenQ += (q[lenQ] != 0ull);
		q[0] += 1;	// q = 2.k.p + 1; No need to check for carry since 2.k.p even
	//printf("Testing F%u, q = %llu * 2^%u + 1, lenQ = %u...\n",j,k,l,lenQ);
		uint32 res1 = mi64_twopmodq(p, lenP, k << (l-j-1), q, lenQ, q2);	// Fiddle k to put q in Mersenne-like form = 2.k'.2^j + 1
			//	res1 = mi64_twopmodq_qferm(j, k << (l-j), q2);
		if(res1 != 1) {
			fprintf(stderr,"ERROR: mi64_twopmodq(F%u, q = %llu * 2^%u + 1 = %s) returns non-unity result %s\n",j,k,l,
					&cbuf1[convert_mi64_base10_char(cbuf1, q, lenQ, 0)],
					&cbuf2[convert_mi64_base10_char(cbuf2,q2, lenQ, 0)]);
			ASSERT(0,"0");
		}
	}

	/* Test 63-bit factors using the 63, 64 and 96-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 63-bit factors...\n");
#endif
	for(i = 0; fac63[i].p != 0; i++)
	{
		p64 = fac63[i].p; q64 = fac63[i].q;
		/* Make sure the MSB = 0: */
		ASSERT(( int64)p64 > 0, "test_fac : ( int64)p64 > 0");
		ASSERT(q64%(2*p64) ==1, "test_fac : q64%(2*p64) ==1");
		k = (q64-1)/(2*p64);	for(j = 0; j < 64; j++) { karr[j] = k; }
		pm60 = p64%60;
		km60 = k  %60;
		/* Since we know q%60 != 0, use the zero column to store the total count of q's for each p%60 value */
		++pqmod60arr[pm60][0];
		++pqmod60arr[pm60][q64%60];

		/* This property only applies for prime exponents, so use a quick base-2 Fermat
		compositeness test as an exponent filter: */
		if(twopmodq64(p64-1, p64) == 1ull && !CHECK_PKMOD60(&p64,1, km60, 0x0))
		{
			fprintf(stderr,"Illegal (p,k) mod 60 pair: p,p mod 60, k,k mod 60 = %llu %4u %llu %4u\n",p64,pm60,k,km60);
			ASSERT(0,"0");
		}

		if((res64 = twopmodq63(p64, q64)) != 1ull)
		{
			fprintf(stderr,"ERROR: twopmodq63(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		if((res64 = twopmodq64(p64, q64)) != 1ull)
		{
			fprintf(stderr,"ERROR: twopmodq64(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#ifdef USE_FLOAT
		k = (q64-1)/(2*p64);
		res64 = twopmodq78_3WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	/* this is currently sse2/msvc only :
		p192.d0 = p64; p192.d1 = p192.d2 = 0;
		x256 = twopmodq200_8WORD_DOUBLE((uint64*)&p192, k);	res64 = !x256.d3 && (uint64)CMPEQ192(x256, ONE192);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq200_8WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	*/
	#endif

		res96 = twopmodq96(p64, k);
		if(!CMPEQ96(ONE96,res96))
		{
			fprintf(stderr,"ERROR: twopmodq96(%llu, k = %llu) returns non-unity result %s\n",p64,k,
					&cbuf2[convert_uint96_base10_char(cbuf2, res96)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128_96(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq128_96(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
		}

	#ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif

		/* Also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 2)
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q2(p64,k,k, 0,0);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q2( %llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q2(p64,k,k);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q2( %llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	#elif(TRYQ == 4)
		res64 = twopmodq63_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq63_q4( %llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q4(p64, k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q4( %llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q4( %llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif

		res64 = twopmodq96_q4(p64,k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq96_q4( %llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q4( %llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 8)
		res64 = twopmodq63_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq63_q8( %llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #if defined(USE_FLOAT) && defined(USE_SSE2) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q8(p64, karr, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q8( %llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
		res64 = twopmodq96_q8(p64,k,k,k,k,k,k,k,k, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq96_q8( %llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q8( %llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 16)
	  #if defined(USE_FLOAT) && defined(USE_AVX)&& defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q16(p64 ,karr, 0,0);
		if(res64 != 0xffff)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q16( %llu, k = %llu x 16) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #else
		#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
	  #endif
	#elif(TRYQ >= 32)
		res64 = twopmodq78_3WORD_DOUBLE_q32(p64 ,karr, 0,0);
		if(res64 != 0xffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q32( %llu, k = %llu x 32) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
//	#elif(TRYQ == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q64(p64 ,karr, 0,0);
		if(res64 != 0xffffffffffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q64( %llu, k = %llu x 64) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif
	}

	/* Test 64-bit factors using the 64 and 96-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 64-bit factors...\n");
#endif
	for(i = 0; fac64[i].p != 0; i++)
	{
		p64 = fac64[i].p; q64 = fac64[i].q;

		ASSERT(q64%(2*p64)==1, "test_fac : q64%(2*p64)==1");

		k = (q64-1)/(2*p64);	for(j = 0; j < 64; j++) { karr[j] = k; }
		pm60 = p64%60;
		km60 = k  %60;
		/* Since we know q%60 != 0, use the zero column to store the total count of q's for each p%60 value */
		++pqmod60arr[pm60][0];
		++pqmod60arr[pm60][q64%60];

		/* This property only applies for prime exponents, so use a quick base-2 Fermat
		compositeness test as an exponent filter: */
		if(twopmodq64(p64-1, p64) == 1ull && !CHECK_PKMOD60(&p64,1, km60, 0x0))
		{
			fprintf(stderr,"Illegal (p,k) mod 60 pair: p,p mod 60, k,k mod 60 = %llu %4u %llu %4u\n",p64,pm60,k,km60);
			ASSERT(0,"0");
		}

		if(q64%(2*p64) != 1)
		{
			fprintf(stderr,"ERROR : (p, q) = ( %llu, %llu ) : q mod (2p) = %llu != 1!\n",p64,q64, q64%(2*p64));
			ASSERT(0,"0");
		}

		if((res64 = twopmodq64(p64, q64)) != 1ull)
		{
			fprintf(stderr,"ERROR: twopmodq64(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#ifdef USE_FLOAT
		k = (q64-1)/(2*p64);
		res64 = twopmodq78_3WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif

		res96 = twopmodq96(p64, k);
		if(!CMPEQ96(ONE96,res96))
		{
			fprintf(stderr,"ERROR: twopmodq96(%llu, k = %llu) returns non-unity result %s\n",p64,k,
					&cbuf2[convert_uint96_base10_char(cbuf2, res96)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128_96(p64,k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq128_96(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE(p64,k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#endif

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 2)
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q2(p64, k,k, 0,0);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q2( %llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q2(p64, k,k);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q2(%llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	#elif(TRYQ == 4)
		res64 = twopmodq64_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq64_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q4(p64, k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
		res64 = twopmodq96_q4(p64,k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 8)
		res64 = twopmodq64_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq64_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #if defined(USE_FLOAT) && defined(USE_SSE2) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q8(p64, karr, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q8(%llu, k = %llu x 4 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
		res64 = twopmodq96_q8(p64,k,k,k,k,k,k,k,k, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 16)
	  #if defined(USE_FLOAT) && defined(USE_AVX)&& defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q16(p64 ,karr, 0,0);
		if(res64 != 0xffff)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q16( %llu, k = %llu x 16) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #else
		#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
	  #endif
	#elif(TRYQ >= 32)
		res64 = twopmodq78_3WORD_DOUBLE_q32(p64 ,karr, 0,0);
		if(res64 != 0xffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q32( %llu, k = %llu x 32) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
//	#elif(TRYQ == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q64(p64 ,karr, 0,0);
		if(res64 != 0xffffffffffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q64( %llu, k = %llu x 64) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif
	}

	/* Test 65-bit factors using the 65 and 96-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 65-bit factors...\n");
#endif
	for(i = 0; fac65[i].p != 0; i++)
	{
		p64 = fac65[i].p; q64 = fac65[i].q;
		q128.d1 = (uint64)1; q128.d0 = q64;

		/* Modify this so it'll work with 65-bit q's: */
		ASSERT(((q64-1)/2 + 0x8000000000000000ull)%p64==0, "test_fac : ((q64-1)/2 + 0x8000000000000000ull)%p64==0");

		k = ((q64-1)/2 + 0x8000000000000000ull)/p64;	for(j = 0; j < 64; j++) { karr[j] = k; }
		pm60 = p64%60;
		km60 = k  %60;

		/* This property only applies for prime exponents, so use a quick base-2 Fermat
		compositeness test as an exponent filter: */
		if(twopmodq64(p64-1, p64) == 1ull && !CHECK_PKMOD60(&p64,1, km60, 0x0))
		{
			fprintf(stderr,"Illegal (p,k) mod 60 pair: p,p mod 60, k,k mod 60 = %llu %4u %llu %4u\n",p64,pm60,k,km60);
			ASSERT(0,"0");
		}
		if((res64 = twopmodq65(p64,k)) != 1)
		{
			fprintf(stderr,"ERROR: twopmodq65(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif

		res96 = twopmodq96(p64, k);
		if(!CMPEQ96(ONE96,res96))
		{
			fprintf(stderr,"ERROR: twopmodq96(%llu, k = %llu) returns non-unity result %s\n",p64,k,
					&cbuf2[convert_uint96_base10_char(cbuf2, res96)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128_96(p64,k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq128_96(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}

	#ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 2)
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q2(p64, k,k, 0,0);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q2(%llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q2(p64, k,k);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q2(%llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	#elif(TRYQ == 4)
		res64 = twopmodq65_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq65_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #ifdef USE_FLOAT
		res64 = twopmodq78_3WORD_DOUBLE_q4(p64, k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
		res64 = twopmodq96_q4(p64,k,k,k,k, 0,0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q4(p64,k,k,k,k);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 8)
		res64 = twopmodq65_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq65_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #if defined(USE_FLOAT) && defined(USE_SSE2) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q8(p64, karr, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q8(%llu, k = %llu x 4 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif
		res64 = twopmodq96_q8(p64,k,k,k,k,k,k,k,k, 0,0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
		res64 = twopmodq128_96_q8(p64,k,k,k,k,k,k,k,k);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq128_96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#elif(TRYQ == 16)
	  #if defined(USE_FLOAT) && defined(USE_AVX)&& defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q16(p64 ,karr, 0,0);
		if(res64 != 0xffff)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q16( %llu, k = %llu x 16) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #else
		#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
	  #endif
	#elif(TRYQ >= 32)
		res64 = twopmodq78_3WORD_DOUBLE_q32(p64 ,karr, 0,0);
		if(res64 != 0xffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q32( %llu, k = %llu x 32) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
//	#elif(TRYQ == 64)
		res64 = twopmodq78_3WORD_DOUBLE_q64(p64 ,karr, 0,0);
		if(res64 != 0xffffffffffffffff) {
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q64( %llu, k = %llu x 64) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif
	}

	/* Test 96-bit factors using the 96 and 128-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 96-bit factors...\n");
#endif
	for(i = 0; fac96[i].p != 0; i++)
	{
		p64 = fac96[i].p;
		q128.d1 = (uint64)fac96[i].d1; q128.d0 = fac96[i].d0;
		q128.d0 -= 1ull;	// Sub off the 1 in q = 2.k.p+1 and check that result == 0 (mod 2.p)
		res64 = mi64_div_binary((uint64*)&q128, &p64, 2,1, 0x0,0x0,0x0);
		if(res64 != 1) {
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint64_base10_char (cbuf2, res64)]);
			ASSERT(0,"0");
		}
		q128.d0 += 1ull;

		/* To find the quotient k = (q-1)/(2*p), which may be > 64 bits, use mod-inverse with base 2^128 arithmetic.
		Since the Newtonian mod-inverse algorithm only works for odd inputs, instead of finding (q-1)/(2*p), we find ((q-1)/2)/p.
		First, find inverse (mod 2^128) of p in preparation for modular multiply. See twopmodq128 for an explanation of this:
		*/
		pinv128.d0 = (p64 + p64 + p64) ^ (uint64)2;	pinv128.d1 = (uint64)0;
		for(j = 0; j < 4; j++)
		{
			hi64 = p64*pinv128.d0;
			pinv128.d0 = pinv128.d0*((uint64)2 - hi64);
		}
		/* pinv128 has 128 bits, but only the upper 64 get modified here. */
	#ifdef MUL_LOHI64_SUBROUTINE
		pinv128.d1 = -pinv128.d0*__MULH64(p64, pinv128.d0);
	#else
		MULH64(p64, pinv128.d0, hi64);
		pinv128.d1 = -pinv128.d0*hi64;
	#endif
		/* k is simply the bottom 128 bits of ((q-1)/2)*pinv128: */
		x128.d0	= ((q128.d0-1) >> 1) + (q128.d1 << 63);	x128.d1	= (q128.d1 >> 1);	/* (q-1)/2. */
		MULL128(x128, pinv128, x128);
	#if 0
		fprintf(stderr,"(p, q) = ( %s, %s ) : k = %s\n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)]);
	#endif
		pm60 = p64%60;
		km60 = (x128.d1*two64mod60 + x128.d0%60)%60;
/*
if((q128.d1 >> 14) == 0) {
	q96.d1 = (uint64)fac96[i].d1; q96.d0 = fac96[i].d0;
	dbl = (double)q96.d0 + (double)q96.d1*TWO64FLOAT;
	rnd = log(dbl)/log(2.0);
	if(rnd > 77)
		printf("p = %10llu, p,k (mod 60) = %2u, %2u, lg(q) = %10.5f\n",p64,pm60,km60,rnd);
}
*/
		/* This property only applies for prime exponents, so use a quick base-2 Fermat
		compositeness test as an exponent filter: */
		if(twopmodq64(p64-1, p64) == 1ull && !CHECK_PKMOD60(&p64,1, km60, 0x0))
		{
			fprintf(stderr,"Illegal (p,k) mod 60 pair: p,p mod 60, k,k mod 60 = %llu %4u %s %4u\n",p64,pm60,
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)],km60);
			ASSERT(0,"0");
		}

	/* Here use full 96-bit q in both floating and 96-bit modmul, so compute for both: */
		q96.d1 = (uint64)fac96[i].d1; q96.d0 = fac96[i].d0;
		/* For twopmodq*() versions taking p and k-args, need to ensure high 64 bits of k (stored in x128.d1) zero: */
		k = x128.d0;	for(j = 0; j < 64; j++) { karr[j] = k; }

	#ifdef USE_FLOAT
	  if((q96.d1 >> 14) == 0)
	  {
		/* Integer-truncation-on-store should obviate the need to subtract 1 from q, and (double)q is only accurate to 53 bits to begin with): */
		ASSERT(x128.d1 == 0, "High half of exactly-computed k nonzero!");
		dbl = (double)q96.d0 + (double)q96.d1*TWO64FLOAT;
		dbl /= (2.0*p64);
		rnd = DNINT(dbl);
		k = (uint64)rnd;
		ASSERT(x128.d0 == k, "Approx and exactly-computed k differ!");
		res64 = twopmodq78_3WORD_DOUBLE(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  }
	#endif
	if(x128.d1 == 0) {	// x128 holds k
		res96 = twopmodq96(p64, k);
		if(!CMPEQ96(ONE96,res96))
		{
			fprintf(stderr,"ERROR: twopmodq96(%llu, k = %llu) returns non-unity result %s\n",p64,k,
					&cbuf2[convert_uint96_base10_char(cbuf2, res96)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128_96(p64, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq128_96(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	}
		p128.d0 = p64; p128.d1 = 0;
		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR: twopmodq128(%s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, res128)]);
			ASSERT(0,"0");
		}

	if(x128.d1 == 0) {
		res64 = twopmodq128x2((uint64 *)&p128, k);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq128x2(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	}

	#ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE(p64, q128);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	#endif

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 2)
	  #ifdef USE_FLOAT
		if((q96.d1 >> 14) == 0 && (x128.d1 == 0))
		{
			res64 = twopmodq78_3WORD_DOUBLE_q2(p64, k,k, 0,0);
			if(res64 != 3)
			{
				fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q2(%llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
		}
	  #endif
	  #ifdef USE_FMADD
		/* Test any FMADD-based modmul routines, if def'd: */
		res64 = twopmodq100_2WORD_DOUBLE_q2(p64, k,k);
		if(res64 != 3)
		{
			fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q2(%llu, k = %llu x 2 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
			ASSERT(0,"0");
		}
	  #endif

	#elif(TRYQ == 4)

		if(x128.d1 == 0)	// k must be 64-bit
		{
		#ifdef USE_FLOAT
			if((q96.d1 >> 14) == 0)
			{
				res64 = twopmodq78_3WORD_DOUBLE_q4(p64, k,k,k,k, 0,0);
				if(res64 != 15)
				{
					fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
					ASSERT(0,"0");
				}
			}
		#endif
		#ifdef USE_FMADD
			/* Test any FMADD-based modmul routines, if def'd: */
			res64 = twopmodq100_2WORD_DOUBLE_q4(p64,k,k,k,k);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq100_2WORD_DOUBLE_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
		#endif
			res64 = twopmodq96_q4(p64,k,k,k,k, 0,0);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
			res64 = twopmodq128_96_q4(p64,k,k,k,k);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq128_96_q4(%llu, k = %llu x 4 ) failed to find factor, res = 0x%1X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
		}	// k must be 64-bit

	#elif(TRYQ == 8)

		if(x128.d1 == 0)	// k must be 64-bit
		{
		#if defined(USE_FLOAT) && defined(USE_SSE2) && (OS_BITS == 64)
			if((q96.d1 >> 14) == 0)
			{
				res64 = twopmodq78_3WORD_DOUBLE_q8(p64, karr, 0,0);
				if(res64 != 255)
				{
					fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q8(%llu, k = %llu x 4 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
					ASSERT(0,"0");
				}
			}
		#endif
			res64 = twopmodq96_q8(p64,k,k,k,k,k,k,k,k, 0,0);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
			res64 = twopmodq128_96_q8(p64,k,k,k,k,k,k,k,k);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq128_96_q8(%llu, k = %llu x 8 ) failed to find factor, res = 0x%2X.\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}
			res64 = twopmodq128_q8((uint64 *)&p128,k,k,k,k,k,k,k,k);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq128_q8( %s, %s x 8 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint128_base10_char(cbuf0,p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1,q128)], (uint32)res64);
				ASSERT(0,"0");
			}
		}	// k must be 64-bit

	#elif(TRYQ == 16)
		if(x128.d1 == 0)	// k must be 64-bit
		{
		#if defined(USE_FLOAT) && defined(USE_AVX)&& defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
			if((q96.d1 >> 14) == 0)
			{
				res64 = twopmodq78_3WORD_DOUBLE_q16(p64 ,karr, 0,0);
				if(res64 != 0xffff)
				{
					fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q16( %llu, k = %llu x 16) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
					ASSERT(0,"0");
				}
			}
		#else
			#error (TRYQ == 16) only supported for 64-bit/P1WORD/GCC/AVX builds!
		#endif
		}	// k must be 64-bit
	#elif(TRYQ >= 32)
		if((x128.d0 >> 52) == 0 && x128.d1 == 0) {	// k must be 64-bit, but 32-operand powmod also needs k < 2^52
			if((q96.d1 >> 14) == 0) {
				res64 = twopmodq78_3WORD_DOUBLE_q32(p64 ,karr, 0,0);
				if(res64 != 0xffffffff) {
					fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q32( %llu, k = %llu x 32) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
					ASSERT(0,"0");
				}
			}
		}	// k must be 52-bit or less
//	#elif(TRYQ == 64)
		if((x128.d0 >> 52) == 0 && x128.d1 == 0) {	// k must be 64-bit, but 64-operand powmod also needs k < 2^52
			if((q96.d1 >> 14) == 0) {
				res64 = twopmodq78_3WORD_DOUBLE_q64(p64 ,karr, 0,0);
				if(res64 != 0xffffffffffffffff) {
					fprintf(stderr,"ERROR: twopmodq78_3WORD_DOUBLE_q64( %llu, k = %llu x 64) failed to find factor, res = 0x%4X.\n",p64,k, (uint32)res64);
					ASSERT(0,"0");
				}
			}
		}	// k must be 52-bit or less
	#endif	// TRYQ
	}

#if(defined(P2WORD) || defined(P3WORD) || defined(P4WORD))

	/* Test 128-bit factors using the 128-bit modmul routines. */
  #ifdef FACTOR_STANDALONE
	printf("Testing 128-bit factors...\n");
  #endif
	/* 128_A: Factors > 96 but <= 128 bits. */
	for(i = 0; fac128[i].p != 0; i++)
	{
/* Comment out the left-justified continue's below to enable factor reconstruction code. */
		p64 = fac128[i].p;
		q128.d1 = fac128[i].d1; q128.d0 = fac128[i].d0;
		q128.d0 -= 1ull;	// Sub off the 1 in q = 2.k.p+1 and check that result == 0 (mod 2.p)
		res64 = mi64_div_binary((uint64*)&q128, &p64, 2,1, 0x0,0x0,0x0);
		if(res64 != 1) {
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint64_base10_char (cbuf2, res64)]);
			ASSERT(0,"0");
		}
		q128.d0 += 1ull;

		/* This uses the generic 128-bit mod function to calculate q%(2*p): */
		p128.d0 = p64; p128.d1 = 0;
		ADD128(p128, p128, two_p128);
		x128 = xmody128(q128, two_p128);
		if(!CMPEQ128(x128, ONE128))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint128_base10_char(cbuf0, p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)]);
			ASSERT(0,"0");
		}

		/* To find the quotient k = (q-1)/(2*p), which may be > 64 bits, use mod-inverse with base 2^128 arithmetic.
		Since the Newtonian mod-inverse algorithm only works for odd inputs, instead of finding (q-1)/(2*p), we find ((q-1)/2)/p.
		First, find inverse (mod 2^128) of p in preparation for modularmultiply. See twopmodq128 for an explanation of this:
		*/
		pinv128.d0 = (p64 + p64 + p64) ^ (uint64)2;  pinv128.d1 = (uint64)0;
		for(j = 0; j < 4; j++)
		{
			hi64 = p64*pinv128.d0;
			pinv128.d0 = pinv128.d0*((uint64)2 - hi64);
		}
		/* pinv128 has 128 bits, but only the upper 64 get modified here. */
		#ifdef MUL_LOHI64_SUBROUTINE
			pinv128.d1 = -pinv128.d0*__MULH64(p64, pinv128.d0);
		#else
			MULH64(p64, pinv128.d0, hi64);
			pinv128.d1 = -pinv128.d0*hi64;
		#endif
		/* k is simply the bottom 128 bits of ((q-1)/2)*pinv128: */
		x128.d0 = ((q128.d0-1) >> 1) + (q128.d1 << 63); x128.d1 = (q128.d1 >> 1);       /* (q-1)/2. */
		MULL128(x128, pinv128, x128);
		pm60 = p64%60;
		/* For twopmodq*() versions taking p and k-args, need to ensure high 64 bits of k (stored in x128.d1) zero: */
		k = x128.d0;
		km60 = (x128.d1*two64mod60 + x128.d0%60)%60;

		/* This property only applies for prime exponents, so use a quick base-2 Fermat
		compositeness test as an exponent filter: */
		if(twopmodq64(p64-1, p64) == 1ull && !CHECK_PKMOD60(&p64,1, km60, 0x0))
		{
			fprintf(stderr,"ERROR: Illegal (p,k) mod 60 pair: p, p mod 60, q128, k mod 60 = %s %4u %s %4u\n",
					&cbuf0[convert_uint64_base10_char (cbuf0,  p64)], pm60,
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)], km60);
			ASSERT(0,"0");
		}

		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR: twopmodq128(%s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, res128)]);
			ASSERT(0,"0");

		#if 0
			/* 10^31 in binary form - need this to reconstruct large factors that were truncated at 30 digits in the PrimeNet report printout: */
			const uint128 ten31 = {542101086242ull, 13875954555633532928ull};
			/* In case this is a factor sufficiently large (> 10^31) that it was truncated (low-order digits lost)
			at 30 digits in the PrimeNet report printout, attempt to reconstruct the lower portion of the factor
			by repeatedly multiplying the retained 30-upper-digit-portion by 10 and filling in the vacated lower decimal digits:
			*/
			printf("Attempting to reconstruct lower digits of the factor: 10^\n");

		    k = 1;
		    for(idx = 1; idx < 10; idx++)       /* Reconstruct up to 9 lower decimal digits */
		    {
				printf("...%d\n", idx);
				k *= 10;

				/* Multiply the reported factor by 10^idx, unless that would cause the result to overflow: */
				if(k > ((uint32)1 << leadz64(q128.d1)))
					goto CONT128;

				MUL_LOHI64(q128.d0, k, x128.d0, x128.d1);
				x128.d1 += k*q128.d1;
				/* Subtract whatever small constant is needed to make the initial candidate == 1 mod p: */
				j = x128.d1%p64;
				j = (j*two64modp)%p64;
				j += x128.d0%p64;
				j = j%p64 - 1;
				res64 = x128.d0 - j;
				x128.d1 -= (res64 > x128.d0);     /* Borrow from high word is checked here. */
				x128.d0 = res64;
				/*printf("   Subtracting %s  gives %s \n",
					&cbuf0[convert_uint64_base10_char (cbuf0, j)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)]);*/

				kmax = k/p64;
				for(j = 0; j <= kmax; j++)
				{
					x128.d0 += p64;	x128.d1 += (x128.d0 < p64);

					/* As a quick filter, only try the candidate if it == +- 1 (mod 8): */
					res64 = x128.d0 & (uint64)7;
					if(res64 == 3 || res64 == 5)
						continue;
					/*
					// 2) Is not divisible by 3,5,7,11,13,17,19:
					// Use that 2^64 == 1 mod  3, i.e. x128 %  3 = (   x128.d1+x128.d0)% 3:
					//	  2^64 == 1 mod  5, i.e. x128 %  5 = (   x128.d1+x128.d0)% 5:
					//	  2^64 == 2 mod  7, i.e. x128 %  7 = ( 2*x128.d1+x128.d0)% 7:
					//	  2^64 == 5 mod 11, i.e. x128 % 11 = ( 5*x128.d1+x128.d0)%11:
					//	  2^64 == 3 mod 13, i.e. x128 % 13 = ( 3*x128.d1+x128.d0)%13:
					//	  2^64 == 1 mod 17, i.e. x128 % 17 = (   x128.d1+x128.d0)%17:
					//	  2^64 ==17 mod 19, i.e. x128 % 19 = (17*x128.d1+x128.d0)%19:
					res64 = x128.d1+x128.d0;
					if(res64%3 == 0 || res64%5 == 0 || res64%17 == 0)
						continue;

					res64 += x128.d1;
					if(res64%7 == 0)
						continue;

					res64 += x128.d1;
					if(res64%13 == 0)
						continue;

					res64 += (x128.d1 << 1);
					if(res64%11 == 0)
						continue;

					res64 += (x128.d1 << 1);
					if(res64%11 == 0)
						continue;
					*/
					/* 3) Satisfies the (p,k)%60 rules: */
					lo.d0   = ((x128.d0-1) >> 1) + (x128.d1 << 63);	lo.d1   = (x128.d1 >> 1);       /* (q-1)/2. */
					MULL128(lo, pinv128, lo);
					km60 = (lo.d1*two64mod60 + lo.d0%60)%60;
					if(!CHECK_PKMOD60(&p64,1, km60, 0x0))
						continue;

					if(twopmodq128(p64, x128) == 1)
					{
						q128 = x128;
						printf("***Reconstructed factor: Q =  %s \n",
					&cbuf0[convert_uint128_base10_char(cbuf0, q128)]);
						goto CONT128;
					}
				}
		    }

			/* Attempt to reconstruct the upper portion by repeatedly adding 10^31 and retrying: */
			printf("Attempting to reconstruct upper digits of the factor...\n");
			x128 = q128;
			for(j = 0; j < 0x10000000; j++)
			{
				ADD128(ten31, x128, x128);

				/* As a quick filter, only try the candidate if it satisfies the (p,k)%60 rules: */
				lo.d0   = ((x128.d0-1) >> 1) + (x128.d1 << 63);lo.d1   = (x128.d1 >> 1);       /* (q-1)/2. */
				MULL128(lo, pinv128, lo);
				km60 = (lo.d1*two64mod60 + lo.d0%60)%60;

				if(!CHECK_PKMOD60(&p64,1, km60, 0x0))
				{
					continue;
				}

				if(twopmodq128(p64, x128) == 1)
				{
					q128 = x128;
					printf("***Reconstructed factor: Q =  %s \n",
					&cbuf0[convert_uint128_base10_char(cbuf0, q128)]);
					break;
				}
			}
		CONT128:
			continue;

		#endif	/* if(0) */
		}

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
		if(x128.d1 == 0)	// k must be 64-bit
		{
			p128.d0 = p64; p128.d1 = 0;
			res64 = twopmodq128x2((uint64 *)&p128, k);
			if(res64 != 1)
			{
				fprintf(stderr,"ERROR: twopmodq128x2(%llu, k = %llu) returns non-unity result %u\n",p64,k, (uint32)res64);
				ASSERT(0,"0");
			}

		#if(TRYQ == 4)
			res64 = twopmodq128_q4((uint64 *)&p128,k,k,k,k);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq128_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint128_base10_char(cbuf0,p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1,q128)], (uint32)res64);
				ASSERT(0,"0");
			}
		#elif(TRYQ == 8)
			res64 = twopmodq128_q8((uint64 *)&p128,k,k,k,k,k,k,k,k);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq128_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
					&cbuf0[convert_uint128_base10_char(cbuf0,p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1,q128)], (uint32)res64);
				ASSERT(0,"0");
			}
		#endif
		}	// k must be 64-bit
	}

#endif	/* #endif(defined(P2WORD) || defined(P3WORD) || defined(P4WORD)) */

#if(defined(P2WORD))

	/*** Only do this time-consuming series of 128-bit factor tests in debug mode: ***/
  #ifdef FAC_DEBUG
	/* 128_B: Construct more 128-bit test factors by multiplying together
	a 63-bit factor q1 of M(p1) and a 64-bit factor q2 of M(p2)
	and checking whether q1*q2 divides M(p1*p2).
	*/
   #ifdef FACTOR_STANDALONE
	printf("Testing 63*64-bit factors...");
   #endif
	/*for(i = 0; fac63[i].p != 0; i++)*/
	for(i = 0; i < 100; i++)
	{
	  for(i2 = 0; fac64[i2].p != 0; i2++)
	  {
		if(fac63[i].p == fac64[i2].p)
			continue;

		p64 = (uint64)fac63[i].p * (uint64)fac64[i2].p;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(fac63[i].q, fac64[i2].q,&q128.d0,&q128.d1);
	#else
		MUL_LOHI64(fac63[i].q, fac64[i2].q, q128.d0, q128.d1);
	#endif

		/* Skip the q%(2*p) == 1 and (p%60,q%60) checks, as they don't apply
		to composite factors which are a product of prime factors of
		different-exponent M(p)'s. */

		p128.d0 = p64; p128.d1 = 0;
		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR_63x64[%u][%u]: twopmodq128(%u*%u, %s*%s ) returns non-unity result %s\n",
					i,i2, fac63[i].p, fac64[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac64[i2].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, q128)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128x2B((uint64*)&p128, q128);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR63x64[%u][%u]: twopmodq128x2(%u*%u, %s*%s ) returns non-unity result %u\n",
					i,i2, fac63[i].p, fac64[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac64[i2].q)], (uint32)res64);
			ASSERT(0,"0");
		}
	  }
	}
   #ifdef FACTOR_STANDALONE
	printf("\n");
   #endif

	/* 128_C: Construct more 128-bit test factors by multiplying together
	a 64-bit factor q1 of M(p1) and a 64-bit factor q2 of M(p2) (p1 != p2)
	and checking whether q1*q2 divides M(p1*p2).
	*/
   #ifdef FACTOR_STANDALONE
	printf("Testing 64*64-bit factors...");
   #endif
	/*for(i = 0; fac64[i].p != 0; i++)*/
	for(i = 0; i < 100; i++)
	{
	  for(i2 = 0; fac64[i2].p != 0; i2++)
	  {
		if(fac64[i].p == fac64[i2].p)
			continue;

		p64 = (uint64)fac64[i].p * (uint64)fac64[i2].p;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(fac64[i].q, fac64[i2].q,&q128.d0,&q128.d1);
	#else
		MUL_LOHI64(fac64[i].q, fac64[i2].q, q128.d0, q128.d1);
	#endif

		/* Skip the q%(2*p) == 1 and (p%60,q%60) checks, as they don't apply
		to composite factors which are a product of prime factors of
		different-exponent M(p)'s. */

		p128.d0 = p64; p128.d1 = 0;
		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR_64x64[%u][%u]: twopmodq128(%u*%u, %s*%s ) returns non-unity result %s\n",
					i,i2, fac64[i].p, fac64[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac64[i].q)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac64[i2].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, q128)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128x2B((uint64*)&p128, q128);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR64x64[%u][%u]: twopmodq128x2(%u*%u, %s*%s ) returns non-unity result %u\n",
					i,i2, fac64[i].p, fac64[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac64[i].q)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac64[i2].q)], (uint32)res64);
			ASSERT(0,"0");
		}
	  }
	}
   #ifdef FACTOR_STANDALONE
	printf("\n");
   #endif

	/* 128_D: Construct more 128-bit test factors by multiplying together
	a 63-bit factor q1 of M(p1) and a 65-bit factor q2 of M(p2)
	and checking whether q1*q2 divides M(p1*p2).
	*/
   #ifdef FACTOR_STANDALONE
	printf("Testing 63*65-bit factors...");
   #endif
	/*for(i = 0; fac63[i].p != 0; i++)*/
	for(i = 0; i < 100; i++)
	{
	  for(i2 = 0; fac65[i2].p != 0; i2++)
	  {
		if(fac63[i].p == fac65[i2].p)
			continue;

		p64 = (uint64)fac63[i].p * (uint64)fac65[i2].p;
		x128.d0 = fac65[i2].q;	x128.d1 = 1;	/* Store full q65 in a 128-bit temp for printing purposes */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(fac63[i].q, fac65[i2].q,&q128.d0,&q128.d1);
	#else
		MUL_LOHI64(fac63[i].q, fac65[i2].q, q128.d0, q128.d1);
	#endif
		/* fac65.q's assumed to have (hidden) 65th bit = 1, so need
		to add 2^64*fac63.q to the output of MUL_LOHI64 here: */
		q128.d1 += fac63[i].q;
		if(q128.d1 <= fac63[i].q)
		{
			fprintf(stderr,"ERROR_63x65[%u][%u]: (p1*p2, q1*q2) = (%u*%u, %s*%s )\n",
					i,i2, fac63[i].p, fac65[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)]);
			fprintf(stderr," q128.d1 += fac63[i].q overflows!\n");
			ASSERT(q128.d1 > fac63[i].q,"q128.d1 > fac63[i].q");	/* Make sure sum didn't overflow */
		}

		/* Skip the q%(2*p) == 1 and (p%60,q%60) checks, as they don't apply
		to composite factors which are a product of prime factors of
		different-exponent M(p)'s. */

		p128.d0 = p64; p128.d1 = 0;
		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR_63x65[%u][%u]: twopmodq128(%u*%u, %s*%s ) returns non-unity result %s\n",
					i,i2, fac63[i].p, fac65[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, res128)]);
			ASSERT(0,"0");
		}

		res64 = twopmodq128x2B((uint64*)&p128, q128);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR63x65[%u][%u]: twopmodq128x2(%u*%u, %s*%s ) returns non-unity result %u\n",
					i,i2, fac63[i].p, fac65[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)], (uint32)res64);
			ASSERT(0,"0");
		}
	  }
	}
   #ifdef FACTOR_STANDALONE
	printf("\n");
   #endif
  #endif	/* #ifdef FAC_DEBUG */

#ifdef FACTOR_STANDALONE
	printf("Testing 128x2-bit factors...\n");
#endif
	for(i = 0; fac128x2[i].plo != 0; i++)
	{
		p128.d1 = fac128x2[i].phi;	p128.d0 = fac128x2[i].plo;
		q128.d1 = fac128x2[i].d1;	q128.d0 = fac128x2[i].d0;

		/* This uses the generic -bit mod function to calculate q%(2*p): */
		ADD128(p128, p128, two_p128);
		x128 = xmody128(q128, two_p128);
		if(!CMPEQ128(x128, ONE128))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint128_base10_char(cbuf0, p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)]);
			ASSERT(0,"0");
		}

		res128 = twopmodq128(p128, q128);
		if(!CMPEQ128(ONE128,res128))
		{
			fprintf(stderr,"ERROR: twopmodq128( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint128_base10_char(cbuf0, p128)],
					&cbuf1[convert_uint128_base10_char(cbuf1, q128)],
					&cbuf2[convert_uint128_base10_char(cbuf2, res128)]);
			ASSERT(0,"0");
		}
	}

#endif	/* #endif(defined(P2WORD)) */

#if(defined(P3WORD) || defined(P4WORD))

	/* Test 128-bit factors using the 160 and 192-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 128x2-bit factors using 160 and 192-bit modmul routines...\n");
#endif
	for(i = 0; fac128x2[i].plo != 0; i++)
	{
		p192.d2 = 0; p192.d1 = fac128x2[i].phi;	p192.d0 = fac128x2[i].plo;
		q192.d2 = 0; q192.d1 = fac128x2[i].d1;	q192.d0 = fac128x2[i].d0;

		/* This uses the generic -bit mod function to calculate q%(2*p): */
		ADD192(p192, p192, two_p192);
		x192 = xmody192(q192, two_p192, 0x0);
		if(!CMPEQ192(x192, ONE192))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, x192)]);
			ASSERT(0,"0");
		}

		// Now compute k = (q-1)/2p, while verifying that q%2p = 1:
		mi64_div((uint64*)&q192, (uint64*)&two_p192, 3,3, (uint64*)&x192, (uint64*)&res192);	// x192 contains k
		ASSERT(x192.d2 == 0 && x192.d1 == 0,"k > 2^64!");
		if(!CMPEQ192(res192, ONE192))
		{
			fprintf(stderr,"ERROR: twopmodq192( %s, %s ) returns non-unity result!\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)]);
			ASSERT(0,"0");
		}

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 4)
/*
		res64 = twopmodq160_q4(p192,q192,q192,q192,q192);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq160_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
			ASSERT(0,"0");
		}
*/
		res64 = twopmodq192_q4((uint64*)&p192,x192.d0,x192.d0,x192.d0,x192.d0);
		if(res64 != 15)
		{
			fprintf(stderr,"ERROR: twopmodq192_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
			ASSERT(0,"0");	// *** disable this to allow fast-UMULH192 timing-testing ***
		}
	#elif(TRYQ == 8)
/*
		res64 = twopmodq160_q8(p192,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq160_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
			ASSERT(0,"0");
		}
*/
		res64 = twopmodq192_q8(p192,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0);
		if(res64 != 255)
		{
			fprintf(stderr,"ERROR: twopmodq192_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
			ASSERT(0,"0");
		}
	#endif
	}

	/* Test 160-bit factors using the 160 and 192-bit modmul routines */
#ifdef FACTOR_STANDALONE
	printf("Testing 160-bit factors using 160 and 192-bit modmul routines...\n");
#endif
	for(i = 0; fac160[i].p != 0; i++)
	{
		p64 = fac160[i].p;
		q192.d2 = fac160[i].d2; q192.d1 = fac160[i].d1; q192.d0 = fac160[i].d0;
		p192.d0 = p64;	p192.d1 = p192.d2 = 0;		ADD192(p192, p192, two_p192);

		two64modp = 0x8000000000000000ull%p64;
		two64modp = (two64modp + two64modp)%p64;
		/* Really large factors may have high part sufficiently large that q.d2*two64modp
		overflows 64 bits, so stick an extra mod-p in there: */
		if((((q192.d2%p64)*two64modp + q192.d1%p64)*two64modp + q192.d0%p64)%p64 != 1)
		{
			fprintf(stderr,"ERROR: q != 1 modulo p for M( %s ), q = %s \n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)]);
			ASSERT(0,"0");
		}

		// Now compute k = (q-1)/2p, while verifying that q%2p = 1:
		mi64_div((uint64*)&q192, (uint64*)&two_p192, 3,3, (uint64*)&x192, (uint64*)&res192);	// x192 contains k
		if(!CMPEQ192(res192, ONE192))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}
/*
		res192 = twopmodq160(p192, q192);
		if(!CMPEQ192(ONE192,res192))
		{
			fprintf(stderr,"ERROR: twopmodq160( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}
*/
		res192 = twopmodq192(p192, q192);
		if(!CMPEQ192(res192, ONE192))
		{
			fprintf(stderr,"ERROR: twopmodq192( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 4)
		if(x192.d2 == 0 && x192.d1 == 0)	// k must be 64-bit for these
		{
		/*
			res64 = twopmodq160_q4(p192,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq160_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		*/
			res64 = twopmodq192_q4((uint64*)&p192,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq192_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		}
	#elif(TRYQ == 8)
		if(x192.d2 == 0 && x192.d1 == 0)	// k must be 64-bit for these
		{
		/*
			res64 = twopmodq160_q8(p192,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq160_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		*/
			res64 = twopmodq192_q8(p192,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq192_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		}
	#endif
	}

#endif	/* #if(defined(P3WORD) || defined(P4WORD)) */

#if(defined(P3WORD) || defined(P4WORD))

	/* Test 192-bit factors using the 192 and 256-bit modmul routines */
	/* 192_A: prime factors < 2^192: */
#ifdef FACTOR_STANDALONE
	printf("Testing 192-bit factors using the 192 and 256-bit modmul routines...\n");
#endif
	for(i = 0; fac192[i].p != 0; i++)
	{
		p64 = fac192[i].p;
		q192.d2 = fac192[i].d2; q192.d1 = fac192[i].d1; q192.d0 = fac192[i].d0;
		p192.d0 = p64;	p192.d1 = p192.d2 = 0;		ADD192(p192, p192, two_p192);

		two64modp = 0x8000000000000000ull%p64;
		two64modp = (two64modp + two64modp)%p64;
		/* Really large factors may have high part sufficiently large that q.d2*two64modp
		overflows 64 bits, so stick an extra mod-p in there: */
		if((((q192.d2%p64)*two64modp + q192.d1%p64)*two64modp + q192.d0%p64)%p64 != 1)
		{
			fprintf(stderr,"ERROR: q != 1 modulo p for M( %s ), q = %s \n",
					&cbuf0[convert_uint64_base10_char (cbuf0, p64)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)]);
			ASSERT(0,"0");
		}

		// Now compute k = (q-1)/2p, while verifying that q%2p = 1:
		mi64_div((uint64*)&q192, (uint64*)&two_p192, 3,3, (uint64*)&x192, (uint64*)&res192);	// x192 contains k
		if(!CMPEQ192(res192, ONE192))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}

		res192 = twopmodq192(p192, q192);
		if(!CMPEQ192(res192, ONE192))
		{
			fprintf(stderr,"ERROR: twopmodq192( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)],
					&cbuf2[convert_uint192_base10_char(cbuf2, res192)]);
			ASSERT(0,"0");
		}

	/* this is currently sse2/msvc only :
		// Floating-double-based impl. allows q up to 200 bits:
		if(x192.d1 == 0 && x192.d2 == 0) {
			x256 = twopmodq200_8WORD_DOUBLE((uint64*)&p192, x192.d0);	res64 = !x256.d3 && (uint64)CMPEQ192(x256, ONE192);
			if(res64 != 1)
			{
				fprintf(stderr,"ERROR: twopmodq200( %s, %s ) returns non-unity result %llu\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1, q192)], res64);
				ASSERT(0,"0");
			}
		}
	*/
		/* This uses the generic 256-bit mod function to calculate q%(2*p): */
		p256.d0 = p64;	p256.d1 = p256.d2 = p256.d3 = 0;
		q256.d0 = q192.d0;
		q256.d1 = q192.d1;
		q256.d2 = q192.d2;
		q256.d3 =       0;
		ADD256(p256, p256, two_p256);
		x256 = xmody256(q256, two_p256, 0x0);
		if(!CMPEQ256(x256, ONE256))
		{
			fprintf(stderr,"ERROR : (p, q) = ( %s, %s ) : q mod (2p) = %s != 1!\n",
					&cbuf0[convert_uint256_base10_char(cbuf0, p256)],
					&cbuf1[convert_uint256_base10_char(cbuf1, q256)],
					&cbuf2[convert_uint256_base10_char(cbuf2, x256)]);
			ASSERT(0,"0");
		}

		res256 = twopmodq256(p256, q256);
		if(!CMPEQ256(ONE256,res256))
		{
			fprintf(stderr,"ERROR: twopmodq256( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint256_base10_char(cbuf0, p256)],
					&cbuf1[convert_uint256_base10_char(cbuf1, q256)],
					&cbuf2[convert_uint256_base10_char(cbuf2, res256)]);
			ASSERT(0,"0");
		}

		/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
	#if(TRYQ == 4)
		if(x192.d2 == 0 && x192.d1 == 0)	// k must be 64-bit for these
		{
			res64 = twopmodq192_q4((uint64*)&p192,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR: twopmodq192_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		}
	#elif(TRYQ == 8)
		if(x192.d2 == 0 && x192.d1 == 0)	// k must be 64-bit for these
		{
			res64 = twopmodq192_q8(p192,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0,x192.d0);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR: twopmodq192_q8( %s, %s x 8 ) failed to find factor, res = 0x%2X.\n",
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint192_base10_char(cbuf1, q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		}
	#endif
	}

#endif

#if(0)//defined(P3WORD))

	/*** Only do this time-consuming series of 192-bit factor tests in debug mode: ***/
  #ifdef FAC_DEBUG
	/* 192_B: Construct more 192-bit test factors by multiplying together
	a 63-bit factor q1 of M(p1), a 64-bit factor q2 of M(p2) and a
	65-bit factor q3 of M(p3) and checking whether q1*q2*q3 divides M(p1*p2*p3).
	*/
   #ifdef FACTOR_STANDALONE
	printf("Testing 63*64*65-bit factors...");
   #endif
	/*for(i = 0; fac63[i].p != 0; i++)*/
	for(i = 0; i < 5; i++)
	{
	  /*for(i2 = 0; fac65[i2].p != 0; i2++)*/
	  for(i2 = 0; i2 < 100; i2++)
	  {
		if(fac63[i].p == fac65[i2].p)
			continue;

		p64 = (uint64)fac63[i].p * (uint64)fac65[i2].p;
		x128.d0 = fac65[i2].q;	x128.d1 = 1;	/* Store full q65 in a 128-bit temp for printing purposes */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(fac63[i].q, fac65[i2].q,&q128.d0,&q128.d1);
	#else
		MUL_LOHI64(fac63[i].q, fac65[i2].q, q128.d0, q128.d1);
	#endif
		/* fac65.q's assumed to have (hidden) 65th bit = 1, so need
		to add 2^64*fac63.q to the output of MUL_LOHI64 here: */
		q128.d1 += fac63[i].q;
		if(q128.d1 <= fac63[i].q)
		{
			fprintf(stderr,"ERROR192_63x65[%u][%u]: (p1*p2, q1*q2) = (%u*%u, %s*%s )\n", i,i2, fac63[i].p, fac65[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)]);
			fprintf(stderr," q128.d1 += fac63[i].q overflows!\n");
			ASSERT(q128.d1 > fac63[i].q,"q128.d1 > fac63[i].q");	/* Make sure sum didn't overflow */
		}

		/* Now multiply the 128-bit 63x65-bit factor product by each 64-bit test factor in turn. */
		/*for(i3 = 0; fac64[i3].p != 0; i3++)*/
		for(i3 = 0; i3 < 100; i3++)
		{
			if(fac64[i3].p == fac63[i].p || fac64[i3].p == fac65[i2].p)
				continue;

			/* Since the product of three test exponents will generally
			overflow 64-bits, store that in the lower 2 words of the p192 variable:
			*/
		#ifdef MUL_LOHI64_SUBROUTINE
			/* Multiply to get 3-exponent product: */
			MUL_LOHI64(p64, fac64[i3].p,&p192.d0,&p192.d1);	p192.d2 = 0;
			/* Low  128 bits of the 192-bit 3-factor product: */
			MUL_LOHI64(q128.d0,fac64[i3].q,&q192.d0,&q192.d1);
			/* High 128 bits of the 192-bit 3-factor product: */
			MUL_LOHI64(q128.d1,fac64[i3].q,  &tmp64,&q192.d2);	q192.d1 += tmp64;	q192.d2 += (q192.d1 < tmp64);
		#else
			/* Multiply to get 3-exponent product: */
			MUL_LOHI64(p64, fac64[i3].p, p192.d0, p192.d1);	p192.d2 = 0;
			/* Low  128 bits of the 192-bit 3-factor product: */
			MUL_LOHI64(q128.d0,fac64[i3].q, q192.d0, q192.d1);
			/* High 128 bits of the 192-bit 3-factor product: */
			MUL_LOHI64(q128.d1,fac64[i3].q,   tmp64, q192.d2);	q192.d1 += tmp64;	q192.d2 += (q192.d1 < tmp64);
		#endif

			/* Skip the q%(2*p) == 1 and (p%60,q%60) checks, as they don't apply
			to composite factors which are a product of prime factors of
			different-exponent M(p)'s. */

			res192 = twopmodq192(p192, q192);
			if(!CMPEQ192(res192, ONE192))
			{
				fprintf(stderr,"ERROR_63x65x64[%u][%u][%u]: twopmodq192(%u*%u*%u = %s, %s*%s*%s = %s) returns non-unity result %s\n", i,i2,i3, fac63[i].p, fac65[i2].p, fac64[i3].p,
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac63[i].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)],
					&cbuf3[convert_uint64_base10_char (cbuf3, fac64[i3].q)],
					&cbuf4[convert_uint192_base10_char(cbuf4, q192)]);
				ASSERT(0,"0");
			}

			p256.d0 = p192.d0;	q256.d0 = q192.d0;
			p256.d1 = p192.d1;	q256.d1 = q192.d1;
			p256.d2 = p192.d2;	q256.d2 = q192.d2;
			p256.d3 =       0;	q256.d3 =       0;
			res256 = twopmodq256(p256, q256);
			if(!CMPEQ256(ONE256,res256))
			{
				fprintf(stderr,"ERROR_63x65x64[%u][%u][%u]: twopmodq256(%u*%u*%u = %s, %s*%s*%s = %s) returns non-unity result %s\n", i,i2,i3, fac63[i].p, fac65[i2].p, fac64[i3].p,
						&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
						&cbuf1[convert_uint64_base10_char (cbuf1, fac63[i].q)],
						&cbuf2[convert_uint128_base10_char(cbuf2, x128)],
						&cbuf3[convert_uint64_base10_char (cbuf3, fac64[i3].q)],
						&cbuf4[convert_uint256_base10_char(cbuf4, q256)],
						&cbuf5[convert_uint256_base10_char(cbuf5, res256)]);
				ASSERT(0,"0");
			}

			/* In debug mode, also test the multiple-q versions of the modular exponentiation routines: */
		#if(TRYQ == 4)
			res64 = twopmodq192_q4(p192,q192,q192,q192,q192);
			if(res64 != 15)
			{
				fprintf(stderr,"ERROR_63x65x64[%u][%u][%u]: (p1*p2*p3, q1*q2*q3) = (%u*%u*%u, %s*%s*%s )\n", i,i2,i3, fac63[i].p, fac65[i2].p, fac64[i3].p,
					&cbuf1[convert_uint64_base10_char (cbuf1, fac63[i].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)],
					&cbuf3[convert_uint64_base10_char (cbuf3, fac64[i3].q)]);
				fprintf(stderr,"ERROR: twopmodq192_q4( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1,q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		#elif(TRYQ == 8)
			res64 = twopmodq192_q8(p192,q192,q192,q192,q192,q192,q192,q192,q192);
			if(res64 != 255)
			{
				fprintf(stderr,"ERROR_63x65x64[%u][%u][%u]: (p1*p2*p3, q1*q2*q3) = (%u*%u*%u, %s*%s*%s )\n", i,i2,i3, fac63[i].p, fac65[i2].p, fac64[i3].p,
					&cbuf1[convert_uint64_base10_char (cbuf1, fac63[i].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)],
					&cbuf3[convert_uint64_base10_char (cbuf3, fac64[i3].q)]);
				fprintf(stderr,"ERROR: twopmodq192_q8( %s, %s x 4 ) failed to find factor, res = 0x%1X.\n",
					&cbuf0[convert_uint192_base10_char(cbuf0, p192)],
					&cbuf1[convert_uint192_base10_char(cbuf1,q192)], (uint32)res64);
				ASSERT(0,"0");
			}
		#endif
		}
	  }
	}
   #ifdef FACTOR_STANDALONE
	printf("\n");
   #endif
  #endif	/* #ifdef FAC_DEBUG */

#endif	/* #if(defined(P3WORD)) */

#if(defined(P4WORD))

  #ifdef FACTOR_STANDALONE
	printf("Testing 256-bit factors...\n");
  #endif
	for(i = 0; i < ntest256; i++)
	{
		p256 = convert_base10_char_uint256(fac256[i].p);	ADD256(p256,p256,two_p256);
		q256 = convert_base10_char_uint256(fac256[i].q);
		ASSERT(CMPEQ256(xmody256(q256, two_p256, &x256), ONE256), "ERROR: q%(2p) != 1");
		res256 = twopmodq256(p256, q256);
		if(!CMPEQ256(res256, ONE256))
		{
			fprintf(stderr,"ERROR: twopmodq256( %s, %s ) returns non-unity result %s\n",
					&cbuf0[convert_uint256_base10_char(cbuf0, p256)],
					&cbuf1[convert_uint256_base10_char(cbuf1, q256)],
					&cbuf2[convert_uint256_base10_char(cbuf2, res256)]);
			ASSERT(0,"0");
		}
	#if 0	/************* need to use k-based for FP200! **********/
	/* this is currently sse2/msvc only :
	  // Floating-double-based impl. allows p up to 128 bits, k up to 64 bits:
	  if((q256.d3 == 0 && q256.d2 == 0 && q256.d1 == 0) && (x256.d3 == 0 && x256.d2 == 0 && x256.d1 == 0))
	  {
		p128.d0 = p192.d0;
		p128.d1 = p192.d1;
	printf("twopmodq200, p = %s, k = %llu\n", fac256->p, x256.d0);
		x256 = twopmodq200_8WORD_DOUBLE(p128, x256.d0);	res64 = !x256.d3 && (uint64)CMPEQ192(x256, ONE192);
		if(res64 != 1)
		{
			fprintf(stderr,"ERROR: twopmodq200( %s, %s ) returns non-unity result %llu\n",
					&cbuf0[convert_uint256_base10_char(cbuf0, p256)],
					&cbuf1[convert_uint256_base10_char(cbuf1, q256)], res64);
			ASSERT(0,"0");
		}
	  }
	*/
	#endif
	}

	/*** Only do this time-consuming series of 256-bit factor tests in debug mode: ***/
  #ifdef FAC_DEBUG
	/* Construct 256-bit test factors by multiplying together a randomly chosen quartet
	consisting of a 63-bit factor q1 of M(p1), 2 (distinct) 64-bit factors q2 of M(p2) and q3 of M(p3)
	and a 65-bit factor q4 of M(p3) and checking whether q1*q2*q3*q4 divides M(p1*p2*p3*p4).
	*/
	#define NTEST256	1000000

   #ifdef FACTOR_STANDALONE
	printf("Testing 63*64*64*65-bit factors...");
   #endif
	/*for(i = 0; fac63[i].p != 0; i++)*/
	for(i = 0; i < 5; i++)
	{
	  /*for(i2 = 0; fac65[i2].p != 0; i2++)*/
	  for(i2 = 0; i2 < 100; i2++)
	  {
		if(fac63[i].p == fac65[i2].p)
			continue;

		p63 = (uint64)fac63[i].p * (uint64)fac65[i2].p;
		x128.d0 = fac65[i2].q;	x128.d1 = 1;	/* Store full q65 in a 128-bit temp for printing purposes */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(fac63[i].q, fac65[i2].q,&q128.d0,&q128.d1);
	#else
		MUL_LOHI64(fac63[i].q, fac65[i2].q, q128.d0, q128.d1);
	#endif
		/* fac65.q's assumed to have (hidden) 65th bit = 1, so need
		to add 2^64*fac63.q to the output of MUL_LOHI64 here: */
		q128.d1 += fac63[i].q;
		if(q128.d1 <= fac63[i].q)
		{
			fprintf(stderr,"ERROR128_63x65[%u][%u]: (p1*p2, q1*q2) = (%u*%u, %s*%s )\n", i,i2, fac63[i].p, fac65[i2].p,
					&cbuf0[convert_uint64_base10_char (cbuf0, fac63[i].q)],
					&cbuf1[convert_uint128_base10_char(cbuf1, x128)]);
			fprintf(stderr," q128.d1 += fac63[i].q overflows!\n");
			ASSERT(q128.d1 > fac63[i].q,"q128.d1 > fac63[i].q");	/* Make sure sum didn't overflow */
		}

		/* Now multiply the 128-bit 63x65-bit factor product by the product of each pair of 64-bit test factors in turn. */
		srand(1);	/* Init stdlib RNG */
		/*for(i3 = 0; fac64[i3].p != 0; i3++)*/
		for(i3 = 0; i3 < 30; i3++)
		{
			ii = rand()%100;
			if(fac64[ii].p == fac63[i].p || fac64[ii].p == fac65[i2].p)
				continue;

			for(i4 = 0; i4 < 30; i4++)
			{
				jj = rand()%100;
				if(fac64[jj].p == fac64[ii].p)
					continue;

				p64 = (uint64)fac64[ii].p * (uint64)fac64[jj].p;
				/* Multiply the two 2-exponent products (p63 and p64) to get 4-exponent product: */
				MUL_LOHI64(p63,p64,&p256.d0,&p256.d1);	p256.d3 = p256.d2 = 0;

				/* Now do a pair of 64-bit scalar*vector MULs to build the 256-bit factor product: */
				q256.d2 = mi64_mul_scalar(&q128.d0, fac64[ii].q, &q256.d0, 2);
				q256.d3 = mi64_mul_scalar(&q256.d0, fac64[jj].q, &q256.d0, 3);

				res256 = twopmodq256(p256, q256);
				if(!CMPEQ256(ONE256,res256))
				{
					fprintf(stderr,"ERROR_63x65x64x64[%u][%u][%u][%u]: twopmodq256(%u*%u*%u*%u = %s, %s*%s*%s*%s = %s) returns non-unity result %s\n", i,i2,i3,i4, fac63[i].p, fac65[i2].p, fac64[ii].p, fac64[jj].p,
					&cbuf0[convert_uint256_base10_char(cbuf0, p256)],
					&cbuf1[convert_uint64_base10_char (cbuf1, fac63[i].q)],
					&cbuf2[convert_uint128_base10_char(cbuf2, x128)],
					&cbuf3[convert_uint64_base10_char (cbuf3, fac64[ii].q)],
					&cbuf4[convert_uint64_base10_char (cbuf4, fac64[jj].q)],
					&cbuf5[convert_uint256_base10_char(cbuf5, q256)],
					&cbuf6[convert_uint256_base10_char(cbuf6, res256)]);
					ASSERT(0,"0");
				}
			}
		}
	  }
	}
   #ifdef FACTOR_STANDALONE
	printf("\n");
   #endif
  #endif	/* #ifdef FAC_DEBUG */

#endif	/* #if(defined(P4WORD)) */

  #ifdef FACTOR_STANDALONE
	printf("Factoring self-tests completed successfully.\n");
  #endif
#endif	// __CUDA_ARCH__ ?
	return 0;
}

