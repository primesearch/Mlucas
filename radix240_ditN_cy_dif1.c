/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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
#include "radix16.h"

#define RADIX 240	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 15	// ODD_RADIX = [radix >> trailz(radix)]

#define USE_COMPACT_OBJ_CODE	1

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

/* Use for toggling higher-accuracy version of the twiddles computation */
//#define HIACC 0	<*** prefer to set via compile-time flag; default is FALSE [= LOACC]

#define EPS 1e-10

// See the radix28 version of this routine for details about the
// small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.

#ifdef USE_SSE2

  #ifdef USE_AVX
	#define OFF1	0x3c0
	#define OFF2	0x780
	#define OFF3	0xb40
	#define OFF4	0xf00
  #else
	#define OFF1	0x1e0
	#define OFF2	0x3c0
	#define OFF3	0x5a0
	#define OFF4	0x780
  #endif

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // For Fermat-mod in AVX mode we need RADIX*4 = 960 [if HIACC] or 12 [if not] vec_dbl slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(68,960) = 960 if AVX+HIACC, max(68,12) = 68 if AVX+LOACC, 20 if SSE2
  // to (half_arr_offset + RADIX) to get required value of radix240_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset240 = 0x480;	// + RADIX = 0x570;  Used for thread local-storage-integrity checking
   #if HIACC	// Need extra 4*RADIX vec_dbl slots in this mode
	const int radix240_creals_in_local_store = 0x940;	// AVX+HIACC: 0x570 + 0x3c0, round up to next even multiple of 0x10
   #else
	const int radix240_creals_in_local_store = 0x5c0;	// AVX+LOACC: 0x570 +  0x44, round up to next even multiple of 0x10
   #endif
	const uint64 radix240_avx_negadwt_consts[RADIX] = {	// 8 entries per line ==> RADIX/8 lines:
		0x3FF0000000000000ull,0x3FEFFFD3151E5533ull,0x3FEFFF4C54F76E1Cull,0x3FEFFE6BC105954Eull,0x3FEFFD315BBF4275ull,0x3FEFFB9D2897136Eull,0x3FEFF9AF2BFBC297ull,0x3FEFF7676B581A63ull,
		0x3FEFF4C5ED12E61Dull,0x3FEFF1CAB88EDFF7ull,0x3FEFEE75D62A9C46ull,0x3FEFEAC74F40720Cull,0x3FEFE6BF2E2660AFull,0x3FEFE25D7E2DF2F8ull,0x3FEFDDA24BA41F4Dull,0x3FEFD88DA3D12526ull,
		0x3FEFD31F94F867C6ull,0x3FEFCD582E584632ull,0x3FEFC7378029F05Full,0x3FEFC0BD9BA139AEull,0x3FEFB9EA92EC689Bull,0x3FEFB2BE793403B9ull,0x3FEFAB39629A9BE3ull,0x3FEFA35B643C93B9ull,
		0x3FEF9B24942FE45Cull,0x3FEF92950983DF6Bull,0x3FEF89ACDC40EE4Bull,0x3FEF806C25684EA8ull,0x3FEF76D2FEF3CC4Bull,0x3FEF6CE183D57825ull,0x3FEF6297CFF75CB0ull,0x3FEF57F6003B2F91ull,
		0x3FEF4CFC327A0080ull,0x3FEF41AA8583E57Eull,0x3FEF3601191FA459ull,0x3FEF2A000E0A5970ull,0x3FEF1DA785F71BCEull,0x3FEF10F7A38E9E90ull,0x3FEF03F08A6ECF94ull,0x3FEEF6925F2A7380ull,
		0x3FEEE8DD4748BF15ull,0x3FEEDAD16944EDD0ull,0x3FEECC6EEC8DD5E9ull,0x3FEEBDB5F9857999ull,0x3FEEAEA6B98095C0ull,0x3FEE9F4156C62DDAull,0x3FEE8F85FC8F1553ull,0x3FEE7F74D705762Bull,
		0x3FEE6F0E134454FFull,0x3FEE5E51DF571265ull,0x3FEE4D406A38E9ABull,0x3FEE3BD9E3D46CEFull,0x3FEE2A1E7D02FE9Full,0x3FEE180E678C4853ull,0x3FEE05A9D625AF0Full,0x3FEDF2F0FC71C4E5ull,
		0x3FEDDFE40EFFB805ull,0x3FEDCC83434ABF29ull,0x3FEDB8CECFB98376ull,0x3FEDA4C6EB9D87C2ull,0x3FED906BCF328D46ull,0x3FED7BBDB39DF5C3ull,0x3FED66BCD2EE2313ull,0x3FED51696819D42Bull,
		0x3FED3BC3AEFF7F95ull,0x3FED25CBE464AB60ull,0x3FED0F8245F5427Full,0x3FECF8E71242E7ABull,0x3FECE1FA88C445BBull,0x3FECCABCE9D45D78ull,0x3FECB32E76B1D0F4ull,0x3FEC9B4F717E2C63ull,
		0x3FEC83201D3D2C6Dull,0x3FEC6AA0BDD40210ull,0x3FEC51D198089406ull,0x3FEC38B2F180BDB1ull,0x3FEC1F4510C18B95ull,0x3FEC05883D2E7560ull,0x3FEBEB7CBF08957Dull,0x3FEBD122DF6DDE43ull,
		0x3FEBB67AE8584CAAull,0x3FEB9B85249D18A2ull,0x3FEB8041DFEBE2FBull,0x3FEB64B166CDE0EEull,0x3FEB48D406A50540ull,0x3FEB2CAA0DAB2702ull,0x3FEB1033CAF125F6ull,0x3FEAF3718E5E0C9Cull,
		0x3FEAD663A8AE2FDCull,0x3FEAB90A6B724C62ull,0x3FEA9B66290EA1A3ull,0x3FEA7D7734BA0A8Eull,0x3FEA5F3DE27D13F2ull,0x3FEA40BA87311090ull,0x3FEA21ED787F2AEFull,0x3FEA02D70CDF74DBull,
		0x3FE9E3779B97F4A8ull,0x3FE9C3CF7CBBB030ull,0x3FE9A3DF0929B594ull,0x3FE983A69A8C21B8ull,0x3FE963268B572492ull,0x3FE9425F36C80335ull,0x3FE92150F8E417B1ull,0x3FE8FFFC2E77CEBAull,
		0x3FE8DE613515A328ull,0x3FE8BC806B151741ull,0x3FE89A5A2F91ABE4ull,0x3FE877EEE269D586ull,0x3FE8553EE43DEF13ull,0x3FE8324A966F2AA5ull,0x3FE80F125B1E8028ull,0x3FE7EB96952B99DCull,
		0x3FE7C7D7A833BEC2ull,0x3FE7A3D5F890BAF9ull,0x3FE77F91EB57C602ull,0x3FE75B0BE65866FBull,0x3FE73644501B56CDull,0x3FE7113B8FE16056ull,0x3FE6EBF20DA23E86ull,0x3FE6C668320B7884ull,
		0x3FE6A09E667F3BCDull,0x3FE67A951513345Cull,0x3FE6544CA88F62DBull,0x3FE62DC58C6CF0DBull,0x3FE607002CD5031Dull,0x3FE5DFFCF69F89EDull,0x3FE5B8BC57520F97ull,0x3FE5913EBD1E84E7ull,
		0x3FE5698496E20BD8ull,0x3FE5418E5423C050ull,0x3FE5195C65137F0Cull,0x3FE4F0EF3A88AAADull,0x3FE4C8474600EEEEull,0x3FE49F64F99F0207ull,0x3FE47648C8296447ull,0x3FE44CF325091DD6ull,
		0x3FE4236484487ABEull,0x3FE3F99D5A91C51Full,0x3FE3CF9E1D2DFDB2ull,0x3FE3A56742039280ull,0x3FE37AF93F9513EAull,0x3FE350548CFFE7F2ull,0x3FE32579A1FAFBDAull,0x3FE2FA68F6D5740Cull,
		0x3FE2CF2304755A5Eull,0x3FE2A3A844564AA5ull,0x3FE277F930881DAFull,0x3FE24C1643AD9295ull,0x3FE21FFFF8FAF674ull,0x3FE1F3B6CC34CA8Bull,0x3FE1C73B39AE68C8ull,0x3FE19A8DBE48A6C1ull,
		0x3FE16DAED770771Dull,0x3FE1409F031D897Eull,0x3FE1135EBFD0E8D7ull,0x3FE0E5EE8C939850ull,0x3FE0B84EE8F52E9Dull,0x3FE08A80550A6FE5ull,0x3FE05C83516BE635ull,0x3FE02E585F347876ull,
		0x3FE0000000000000ull,0x3FDFA2F56BD3B979ull,0x3FDF459207170FCEull,0x3FDEE7D6D7F64AD2ull,0x3FDE89C4E59427B1ull,0x3FDE2B5D3806F63Bull,0x3FDDCCA0D855B380ull,0x3FDD6D90D07521CBull,
		0x3FDD0E2E2B44DE01ull,0x3FDCAE79F48C726Cull,0x3FDC4E7538F866FCull,0x3FDBEE2106174F02ull,0x3FDB8D7E6A56D476ull,0x3FDB2C8E7500C0C6ull,0x3FDACB523638033Bull,0x3FDA69CABEF5B501ull,
		0x3FDA07F921061AD1ull,0x3FD9A5DE6F05A44Bull,0x3FD9437BBC5DE90Aull,0x3FD8E0D21D42A377ull,0x3FD87DE2A6AEA963ull,0x3FD81AAE6E60E271ull,0x3FD7B7368AD93C61ull,0x3FD7537C13559D33ull,
		0x3FD6EF801FCED33Cull,0x3FD68B43C8F5832Aull,0x3FD626C8282F1408ull,0x3FD5C20E57929942ull,0x3FD55D1771E5BAB9ull,0x3FD4F7E492999AEEull,0x3FD49276D5C7BB48ull,0x3FD42CCF582EDE82ull,
		0x3FD3C6EF372FE950ull,0x3FD360D790CAC12Eull,0x3FD2FA89839B2985ull,0x3FD294062ED59F06ull,0x3FD22D4EB2443163ull,0x3FD1C6642E435B69ull,0x3FD15F47C3BED971ull,0x3FD0F7FA942E7E48ull,
		0x3FD0907DC1930690ull,0x3FD028D26E72EA99ull,0x3FCF81F37BAE5D8Cull,0x3FCEB1E9A690650Eull,0x3FCDE189A594FBCCull,0x3FCD10D5C1B71B7Full,0x3FCC3FD044DD3D45ull,0x3FCB6E7B79D2ECC9ull,
		0x3FCA9CD9AC4258F6ull,0x3FC9CAED28ADE228ull,0x3FC8F8B83C69A60Bull,0x3FC8263D35950926ull,0x3FC7537E63143E2Eull,0x3FC6807E1489CB33ull,0x3FC5AD3E9A500CADull,0x3FC4D9C24572B693ull,
		0x3FC4060B67A85375ull,0x3FC3321C534BC1BBull,0x3FC25DF75B55AF15ull,0x3FC1899ED3561233ull,0x3FC0B5150F6DA2D1ull,0x3FBFC0B8C88EA05Aull,0x3FBE16EE4E236BF8ull,0x3FBC6CCF5AF11FD1ull,
		0x3FBAC2609B3C576Cull,0x3FB917A6BC29B42Cull,0x3FB76CA66BB0BC83ull,0x3FB5C164588EB8DAull,0x3FB415E532398E49ull,0x3FB26A2DA8D2974Aull,0x3FB0BE426D197A8Bull,0x3FAE245060BE0012ull,
		0x3FAACBC748EFC90Eull,0x3FA772F2F75F573Cull,0x3FA419DCD176E0F7ull,0x3FA0C08E3D596AEEull,0x3F9ACE214390CA91ull,0x3F941ADACC128E22ull,0x3F8ACEB7C72CA0A8ull,0x3F7ACEDD6862D0D7ull
	};
  #else
	const int half_arr_offset240 = 0x4f8;	// + RADIX = 0x5e8; Used for thread local-storage-integrity checking
	const int radix240_creals_in_local_store = 0x1000;	// (half_arr_offset + RADIX) + 0x14(=20), round up to next even multiple of 0x10
  #endif

	#include "sse2_macro.h"
	#include "radix15_sse_macro.h"

#endif	// USE_SSE2

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data:
		int iter;
		int tid;
		int ndivr;

		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;
		int wts_idx_inc2;
		int icycle[ODD_RADIX];
	#ifdef USE_SSE2
		int jcycle[ODD_RADIX];
	  #ifdef USE_AVX
		int kcycle[ODD_RADIX];
		int lcycle[ODD_RADIX];
	  #endif
	#endif

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *rn0;
		struct complex *rn1;
	#ifdef USE_SSE2
		vec_dbl *r00;
		vec_dbl *half_arr;
	#else
		double *r00;
		double *half_arr;
	#endif
		uint64*sm_ptr;
		uint32 bjmodnini;
		int bjmodn0;
	// For large radix0 use thread-local arrays for DWT indices/carries - only caveat is these must be SIMD-aligned:
	// Since GCC uses __BIGGEST_ALIGNMENT__ = 16 on x86, which is too small to be useful for avx data,
	// we are forced to resort to fugly hackage - add pad slots to a garbage-named struct-internal array along with
	// a pointer-to-be-inited-at-runtime, when we set ptr to the lowest-index array element having the desired alginment:
		double *cy_r,*cy_i;
		double cy_dat[2*RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	};

#endif

/***************/

int radix240_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-240 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-240 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix240_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl), sz_vd_m1 = sz_vd-1;
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#else
	const int sz_vd = sizeof(double), sz_vd_m1 = sz_vd-1;
#endif

	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
#endif
	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;

	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0;
	static int poff[RADIX>>2], p_out_hi[ODD_RADIX];
#ifndef MULTITHREAD
	static int plo[16], p_in_hi[ODD_RADIX], po_br[16];
	int k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	uint64 i64;
	const uint64 dif_perm16[ODD_RADIX] = {
		0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x76450123fecd89abull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,0x32017645ba89fecdull,
		0xefdc98ba10325476ull,0x45672310cdefab98ull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull};
	const uint64 dit_perm16[ODD_RADIX] = {
		0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x98abefcd10236745ull,
		0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,
		0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x89bafedc01327654ull};
#endif
	static double radix_inv, n2inv;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *addr, *addi;
	struct complex t[RADIX], *tptr;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
	static int ii[ODD_RADIX] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		,icycle[ODD_RADIX],ic;
#ifdef USE_SSE2
	static int jcycle[ODD_RADIX],jc;
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX];	// NB: kc already declared as part of k0-f set above
	static int lcycle[ODD_RADIX],lc;
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7, *add8, *add9, *adda, *addb, *addc, *addd, *adde, *addf;
  #endif	// MULTITHREAD

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2, *cc0, *ss0, *sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,
		*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,*r0f,
		*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,*r1f,
		*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,*r2f,
		*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,*r3f,
		*r40,*r41,*r42,*r43,*r44,*r45,*r46,*r47,*r48,*r49,*r4a,*r4b,*r4c,*r4d,*r4e,*r4f,
		*r50,*r51,*r52,*r53,*r54,*r55,*r56,*r57,*r58,*r59,*r5a,*r5b,*r5c,*r5d,*r5e,*r5f,
		*r60,*r61,*r62,*r63,*r64,*r65,*r66,*r67,*r68,*r69,*r6a,*r6b,*r6c,*r6d,*r6e,*r6f,
		*r70,*r71,*r72,*r73,*r74,*r75,*r76,*r77,*r78,*r79,*r7a,*r7b,*r7c,*r7d,*r7e,*r7f,
		*r80,*r81,*r82,*r83,*r84,*r85,*r86,*r87,*r88,*r89,*r8a,*r8b,*r8c,*r8d,*r8e,*r8f,
		*r90,*r91,*r92,*r93,*r94,*r95,*r96,*r97,*r98,*r99,*r9a,*r9b,*r9c,*r9d,*r9e,*r9f,
		*ra0,*ra1,*ra2,*ra3,*ra4,*ra5,*ra6,*ra7,*ra8,*ra9,*raa,*rab,*rac,*rad,*rae,*raf,
		*rb0,*rb1,*rb2,*rb3,*rb4,*rb5,*rb6,*rb7,*rb8,*rb9,*rba,*rbb,*rbc,*rbd,*rbe,*rbf,
		*rc0,*rc1,*rc2,*rc3,*rc4,*rc5,*rc6,*rc7,*rc8,*rc9,*rca,*rcb,*rcc,*rcd,*rce,*rcf,
		*rd0,*rd1,*rd2,*rd3,*rd4,*rd5,*rd6,*rd7,*rd8,*rd9,*rda,*rdb,*rdc,*rdd,*rde,*rdf,
		*re0,*re1,*re2,*re3,*re4,*re5,*re6,*re7,*re8,*re9,*rea,*reb,*rec,*red,*ree,*ref,
		*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,*s1p0f,
		*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*s1p18,*s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,*s1p1f,
		*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,*s1p2f,
		*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*s1p38,*s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,*s1p3f,
		*s1p40,*s1p41,*s1p42,*s1p43,*s1p44,*s1p45,*s1p46,*s1p47,*s1p48,*s1p49,*s1p4a,*s1p4b,*s1p4c,*s1p4d,*s1p4e,*s1p4f,
		*s1p50,*s1p51,*s1p52,*s1p53,*s1p54,*s1p55,*s1p56,*s1p57,*s1p58,*s1p59,*s1p5a,*s1p5b,*s1p5c,*s1p5d,*s1p5e,*s1p5f,
		*s1p60,*s1p61,*s1p62,*s1p63,*s1p64,*s1p65,*s1p66,*s1p67,*s1p68,*s1p69,*s1p6a,*s1p6b,*s1p6c,*s1p6d,*s1p6e,*s1p6f,
		*s1p70,*s1p71,*s1p72,*s1p73,*s1p74,*s1p75,*s1p76,*s1p77,*s1p78,*s1p79,*s1p7a,*s1p7b,*s1p7c,*s1p7d,*s1p7e,*s1p7f,
		*s1p80,*s1p81,*s1p82,*s1p83,*s1p84,*s1p85,*s1p86,*s1p87,*s1p88,*s1p89,*s1p8a,*s1p8b,*s1p8c,*s1p8d,*s1p8e,*s1p8f,
		*s1p90,*s1p91,*s1p92,*s1p93,*s1p94,*s1p95,*s1p96,*s1p97,*s1p98,*s1p99,*s1p9a,*s1p9b,*s1p9c,*s1p9d,*s1p9e,*s1p9f,
		*s1pa0,*s1pa1,*s1pa2,*s1pa3,*s1pa4,*s1pa5,*s1pa6,*s1pa7,*s1pa8,*s1pa9,*s1paa,*s1pab,*s1pac,*s1pad,*s1pae,*s1paf,
		*s1pb0,*s1pb1,*s1pb2,*s1pb3,*s1pb4,*s1pb5,*s1pb6,*s1pb7,*s1pb8,*s1pb9,*s1pba,*s1pbb,*s1pbc,*s1pbd,*s1pbe,*s1pbf,
		*s1pc0,*s1pc1,*s1pc2,*s1pc3,*s1pc4,*s1pc5,*s1pc6,*s1pc7,*s1pc8,*s1pc9,*s1pca,*s1pcb,*s1pcc,*s1pcd,*s1pce,*s1pcf,
		*s1pd0,*s1pd1,*s1pd2,*s1pd3,*s1pd4,*s1pd5,*s1pd6,*s1pd7,*s1pd8,*s1pd9,*s1pda,*s1pdb,*s1pdc,*s1pdd,*s1pde,*s1pdf,
		*s1pe0,*s1pe1,*s1pe2,*s1pe3,*s1pe4,*s1pe5,*s1pe6,*s1pe7,*s1pe8,*s1pe9,*s1pea,*s1peb,*s1pec,*s1ped,*s1pee,*s1pef,
		*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
		*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
		*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer

#endif	// USE_SSE2

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy240_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
  #if PFETCH
	double *addp;
  #endif
	int bjmodn[RADIX];
	double temp,frac,
		cy_r[RADIX],cy_i[RADIX];

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static int *_bjmodnini = 0x0, *_bjmodn[RADIX];
	static double *_maxerr = 0x0, *_cy_r[RADIX],*_cy_i[RADIX];
	if(!_maxerr) {
		_cy_r[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

	foo_array[0] = base[0];
	foo_array[1] = base[1];
	foo_array[2] = baseinv[0];
	foo_array[3] = baseinv[1];
	wt_arr    = foo_array + 4;
	wtinv_arr = wt_arr    + ODD_RADIX;
	bs_arr    = wtinv_arr + ODD_RADIX;
	bsinv_arr = bs_arr    + ODD_RADIX;

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/%d in %s.\n", iter,RADIX,func);
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

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		if(CY_THREADS < NTHREADS)	{ WARN(HERE, "CY_THREADS < NTHREADS", "", 1); return(ERR_ASSERT); }
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, j);

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#if 0//def OS_TYPE_MACOSX

			if(CY_THREADS > 1) {
				main_work_units = CY_THREADS/2;
				pool_work_units = CY_THREADS - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
			} else {
				main_work_units = 1;
				printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
			}

		#else

			main_work_units = 0;
			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
			tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < 4; l++) {
				if( ((uint32)&tdat[ithread].cy_dat[l] & sz_vd_m1) == 0 ) {
					tdat[ithread].cy_r = &tdat[ithread].cy_dat[l];
					tdat[ithread].cy_i = tdat[ithread].cy_r + RADIX;
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",sz_vd,l);
					break;
				}
			}
			ASSERT(HERE, l < 4, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix240_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix240_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr + 0;		tm2 = tmp + 0x100;
		r00 = tmp + 0x00;		r80 = tm2 + 0x00;
		r01 = tmp + 0x02;		r81 = tm2 + 0x02;
		r02 = tmp + 0x04;		r82 = tm2 + 0x04;
		r03 = tmp + 0x06;		r83 = tm2 + 0x06;
		r04 = tmp + 0x08;		r84 = tm2 + 0x08;
		r05 = tmp + 0x0a;		r85 = tm2 + 0x0a;
		r06 = tmp + 0x0c;		r86 = tm2 + 0x0c;
		r07 = tmp + 0x0e;		r87 = tm2 + 0x0e;
		r08 = tmp + 0x10;		r88 = tm2 + 0x10;
		r09 = tmp + 0x12;		r89 = tm2 + 0x12;
		r0a = tmp + 0x14;		r8a = tm2 + 0x14;
		r0b = tmp + 0x16;		r8b = tm2 + 0x16;
		r0c = tmp + 0x18;		r8c = tm2 + 0x18;
		r0d = tmp + 0x1a;		r8d = tm2 + 0x1a;
		r0e = tmp + 0x1c;		r8e = tm2 + 0x1c;
		r0f = tmp + 0x1e;		r8f = tm2 + 0x1e;
		r10 = tmp + 0x20;		r90 = tm2 + 0x20;
		r11 = tmp + 0x22;		r91 = tm2 + 0x22;
		r12 = tmp + 0x24;		r92 = tm2 + 0x24;
		r13 = tmp + 0x26;		r93 = tm2 + 0x26;
		r14 = tmp + 0x28;		r94 = tm2 + 0x28;
		r15 = tmp + 0x2a;		r95 = tm2 + 0x2a;
		r16 = tmp + 0x2c;		r96 = tm2 + 0x2c;
		r17 = tmp + 0x2e;		r97 = tm2 + 0x2e;
		r18 = tmp + 0x30;		r98 = tm2 + 0x30;
		r19 = tmp + 0x32;		r99 = tm2 + 0x32;
		r1a = tmp + 0x34;		r9a = tm2 + 0x34;
		r1b = tmp + 0x36;		r9b = tm2 + 0x36;
		r1c = tmp + 0x38;		r9c = tm2 + 0x38;
		r1d = tmp + 0x3a;		r9d = tm2 + 0x3a;
		r1e = tmp + 0x3c;		r9e = tm2 + 0x3c;
		r1f = tmp + 0x3e;		r9f = tm2 + 0x3e;
		r20 = tmp + 0x40;		ra0 = tm2 + 0x40;
		r21 = tmp + 0x42;		ra1 = tm2 + 0x42;
		r22 = tmp + 0x44;		ra2 = tm2 + 0x44;
		r23 = tmp + 0x46;		ra3 = tm2 + 0x46;
		r24 = tmp + 0x48;		ra4 = tm2 + 0x48;
		r25 = tmp + 0x4a;		ra5 = tm2 + 0x4a;
		r26 = tmp + 0x4c;		ra6 = tm2 + 0x4c;
		r27 = tmp + 0x4e;		ra7 = tm2 + 0x4e;
		r28 = tmp + 0x50;		ra8 = tm2 + 0x50;
		r29 = tmp + 0x52;		ra9 = tm2 + 0x52;
		r2a = tmp + 0x54;		raa = tm2 + 0x54;
		r2b = tmp + 0x56;		rab = tm2 + 0x56;
		r2c = tmp + 0x58;		rac = tm2 + 0x58;
		r2d = tmp + 0x5a;		rad = tm2 + 0x5a;
		r2e = tmp + 0x5c;		rae = tm2 + 0x5c;
		r2f = tmp + 0x5e;		raf = tm2 + 0x5e;
		r30 = tmp + 0x60;		rb0 = tm2 + 0x60;
		r31 = tmp + 0x62;		rb1 = tm2 + 0x62;
		r32 = tmp + 0x64;		rb2 = tm2 + 0x64;
		r33 = tmp + 0x66;		rb3 = tm2 + 0x66;
		r34 = tmp + 0x68;		rb4 = tm2 + 0x68;
		r35 = tmp + 0x6a;		rb5 = tm2 + 0x6a;
		r36 = tmp + 0x6c;		rb6 = tm2 + 0x6c;
		r37 = tmp + 0x6e;		rb7 = tm2 + 0x6e;
		r38 = tmp + 0x70;		rb8 = tm2 + 0x70;
		r39 = tmp + 0x72;		rb9 = tm2 + 0x72;
		r3a = tmp + 0x74;		rba = tm2 + 0x74;
		r3b = tmp + 0x76;		rbb = tm2 + 0x76;
		r3c = tmp + 0x78;		rbc = tm2 + 0x78;
		r3d = tmp + 0x7a;		rbd = tm2 + 0x7a;
		r3e = tmp + 0x7c;		rbe = tm2 + 0x7c;
		r3f = tmp + 0x7e;		rbf = tm2 + 0x7e;
		r40 = tmp + 0x80;		rc0 = tm2 + 0x80;
		r41 = tmp + 0x82;		rc1 = tm2 + 0x82;
		r42 = tmp + 0x84;		rc2 = tm2 + 0x84;
		r43 = tmp + 0x86;		rc3 = tm2 + 0x86;
		r44 = tmp + 0x88;		rc4 = tm2 + 0x88;
		r45 = tmp + 0x8a;		rc5 = tm2 + 0x8a;
		r46 = tmp + 0x8c;		rc6 = tm2 + 0x8c;
		r47 = tmp + 0x8e;		rc7 = tm2 + 0x8e;
		r48 = tmp + 0x90;		rc8 = tm2 + 0x90;
		r49 = tmp + 0x92;		rc9 = tm2 + 0x92;
		r4a = tmp + 0x94;		rca = tm2 + 0x94;
		r4b = tmp + 0x96;		rcb = tm2 + 0x96;
		r4c = tmp + 0x98;		rcc = tm2 + 0x98;
		r4d = tmp + 0x9a;		rcd = tm2 + 0x9a;
		r4e = tmp + 0x9c;		rce = tm2 + 0x9c;
		r4f = tmp + 0x9e;		rcf = tm2 + 0x9e;
		r50 = tmp + 0xa0;		rd0 = tm2 + 0xa0;
		r51 = tmp + 0xa2;		rd1 = tm2 + 0xa2;
		r52 = tmp + 0xa4;		rd2 = tm2 + 0xa4;
		r53 = tmp + 0xa6;		rd3 = tm2 + 0xa6;
		r54 = tmp + 0xa8;		rd4 = tm2 + 0xa8;
		r55 = tmp + 0xaa;		rd5 = tm2 + 0xaa;
		r56 = tmp + 0xac;		rd6 = tm2 + 0xac;
		r57 = tmp + 0xae;		rd7 = tm2 + 0xae;
		r58 = tmp + 0xb0;		rd8 = tm2 + 0xb0;
		r59 = tmp + 0xb2;		rd9 = tm2 + 0xb2;
		r5a = tmp + 0xb4;		rda = tm2 + 0xb4;
		r5b = tmp + 0xb6;		rdb = tm2 + 0xb6;
		r5c = tmp + 0xb8;		rdc = tm2 + 0xb8;
		r5d = tmp + 0xba;		rdd = tm2 + 0xba;
		r5e = tmp + 0xbc;		rde = tm2 + 0xbc;
		r5f = tmp + 0xbe;		rdf = tm2 + 0xbe;
		r60 = tmp + 0xc0;		re0 = tm2 + 0xc0;
		r61 = tmp + 0xc2;		re1 = tm2 + 0xc2;
		r62 = tmp + 0xc4;		re2 = tm2 + 0xc4;
		r63 = tmp + 0xc6;		re3 = tm2 + 0xc6;
		r64 = tmp + 0xc8;		re4 = tm2 + 0xc8;
		r65 = tmp + 0xca;		re5 = tm2 + 0xca;
		r66 = tmp + 0xcc;		re6 = tm2 + 0xcc;
		r67 = tmp + 0xce;		re7 = tm2 + 0xce;
		r68 = tmp + 0xd0;		re8 = tm2 + 0xd0;
		r69 = tmp + 0xd2;		re9 = tm2 + 0xd2;
		r6a = tmp + 0xd4;		rea = tm2 + 0xd4;
		r6b = tmp + 0xd6;		reb = tm2 + 0xd6;
		r6c = tmp + 0xd8;		rec = tm2 + 0xd8;
		r6d = tmp + 0xda;		red = tm2 + 0xda;
		r6e = tmp + 0xdc;		ree = tm2 + 0xdc;
		r6f = tmp + 0xde;		ref = tm2 + 0xde;
		r70 = tmp + 0xe0;
		r71 = tmp + 0xe2;
		r72 = tmp + 0xe4;
		r73 = tmp + 0xe6;
		r74 = tmp + 0xe8;
		r75 = tmp + 0xea;
		r76 = tmp + 0xec;
		r77 = tmp + 0xee;
		r78 = tmp + 0xf0;
		r79 = tmp + 0xf2;
		r7a = tmp + 0xf4;
		r7b = tmp + 0xf6;
		r7c = tmp + 0xf8;
		r7d = tmp + 0xfa;
		r7e = tmp + 0xfc;
		r7f = tmp + 0xfe;
		tmp += 0x1e0;			tm2 += 0x1e0;
		s1p00 = tmp + 0x00;		s1p80 = tm2 + 0x00;
		s1p01 = tmp + 0x02;		s1p81 = tm2 + 0x02;
		s1p02 = tmp + 0x04;		s1p82 = tm2 + 0x04;
		s1p03 = tmp + 0x06;		s1p83 = tm2 + 0x06;
		s1p04 = tmp + 0x08;		s1p84 = tm2 + 0x08;
		s1p05 = tmp + 0x0a;		s1p85 = tm2 + 0x0a;
		s1p06 = tmp + 0x0c;		s1p86 = tm2 + 0x0c;
		s1p07 = tmp + 0x0e;		s1p87 = tm2 + 0x0e;
		s1p08 = tmp + 0x10;		s1p88 = tm2 + 0x10;
		s1p09 = tmp + 0x12;		s1p89 = tm2 + 0x12;
		s1p0a = tmp + 0x14;		s1p8a = tm2 + 0x14;
		s1p0b = tmp + 0x16;		s1p8b = tm2 + 0x16;
		s1p0c = tmp + 0x18;		s1p8c = tm2 + 0x18;
		s1p0d = tmp + 0x1a;		s1p8d = tm2 + 0x1a;
		s1p0e = tmp + 0x1c;		s1p8e = tm2 + 0x1c;
		s1p0f = tmp + 0x1e;		s1p8f = tm2 + 0x1e;
		s1p10 = tmp + 0x20;		s1p90 = tm2 + 0x20;
		s1p11 = tmp + 0x22;		s1p91 = tm2 + 0x22;
		s1p12 = tmp + 0x24;		s1p92 = tm2 + 0x24;
		s1p13 = tmp + 0x26;		s1p93 = tm2 + 0x26;
		s1p14 = tmp + 0x28;		s1p94 = tm2 + 0x28;
		s1p15 = tmp + 0x2a;		s1p95 = tm2 + 0x2a;
		s1p16 = tmp + 0x2c;		s1p96 = tm2 + 0x2c;
		s1p17 = tmp + 0x2e;		s1p97 = tm2 + 0x2e;
		s1p18 = tmp + 0x30;		s1p98 = tm2 + 0x30;
		s1p19 = tmp + 0x32;		s1p99 = tm2 + 0x32;
		s1p1a = tmp + 0x34;		s1p9a = tm2 + 0x34;
		s1p1b = tmp + 0x36;		s1p9b = tm2 + 0x36;
		s1p1c = tmp + 0x38;		s1p9c = tm2 + 0x38;
		s1p1d = tmp + 0x3a;		s1p9d = tm2 + 0x3a;
		s1p1e = tmp + 0x3c;		s1p9e = tm2 + 0x3c;
		s1p1f = tmp + 0x3e;		s1p9f = tm2 + 0x3e;
		s1p20 = tmp + 0x40;		s1pa0 = tm2 + 0x40;
		s1p21 = tmp + 0x42;		s1pa1 = tm2 + 0x42;
		s1p22 = tmp + 0x44;		s1pa2 = tm2 + 0x44;
		s1p23 = tmp + 0x46;		s1pa3 = tm2 + 0x46;
		s1p24 = tmp + 0x48;		s1pa4 = tm2 + 0x48;
		s1p25 = tmp + 0x4a;		s1pa5 = tm2 + 0x4a;
		s1p26 = tmp + 0x4c;		s1pa6 = tm2 + 0x4c;
		s1p27 = tmp + 0x4e;		s1pa7 = tm2 + 0x4e;
		s1p28 = tmp + 0x50;		s1pa8 = tm2 + 0x50;
		s1p29 = tmp + 0x52;		s1pa9 = tm2 + 0x52;
		s1p2a = tmp + 0x54;		s1paa = tm2 + 0x54;
		s1p2b = tmp + 0x56;		s1pab = tm2 + 0x56;
		s1p2c = tmp + 0x58;		s1pac = tm2 + 0x58;
		s1p2d = tmp + 0x5a;		s1pad = tm2 + 0x5a;
		s1p2e = tmp + 0x5c;		s1pae = tm2 + 0x5c;
		s1p2f = tmp + 0x5e;		s1paf = tm2 + 0x5e;
		s1p30 = tmp + 0x60;		s1pb0 = tm2 + 0x60;
		s1p31 = tmp + 0x62;		s1pb1 = tm2 + 0x62;
		s1p32 = tmp + 0x64;		s1pb2 = tm2 + 0x64;
		s1p33 = tmp + 0x66;		s1pb3 = tm2 + 0x66;
		s1p34 = tmp + 0x68;		s1pb4 = tm2 + 0x68;
		s1p35 = tmp + 0x6a;		s1pb5 = tm2 + 0x6a;
		s1p36 = tmp + 0x6c;		s1pb6 = tm2 + 0x6c;
		s1p37 = tmp + 0x6e;		s1pb7 = tm2 + 0x6e;
		s1p38 = tmp + 0x70;		s1pb8 = tm2 + 0x70;
		s1p39 = tmp + 0x72;		s1pb9 = tm2 + 0x72;
		s1p3a = tmp + 0x74;		s1pba = tm2 + 0x74;
		s1p3b = tmp + 0x76;		s1pbb = tm2 + 0x76;
		s1p3c = tmp + 0x78;		s1pbc = tm2 + 0x78;
		s1p3d = tmp + 0x7a;		s1pbd = tm2 + 0x7a;
		s1p3e = tmp + 0x7c;		s1pbe = tm2 + 0x7c;
		s1p3f = tmp + 0x7e;		s1pbf = tm2 + 0x7e;
		s1p40 = tmp + 0x80;		s1pc0 = tm2 + 0x80;
		s1p41 = tmp + 0x82;		s1pc1 = tm2 + 0x82;
		s1p42 = tmp + 0x84;		s1pc2 = tm2 + 0x84;
		s1p43 = tmp + 0x86;		s1pc3 = tm2 + 0x86;
		s1p44 = tmp + 0x88;		s1pc4 = tm2 + 0x88;
		s1p45 = tmp + 0x8a;		s1pc5 = tm2 + 0x8a;
		s1p46 = tmp + 0x8c;		s1pc6 = tm2 + 0x8c;
		s1p47 = tmp + 0x8e;		s1pc7 = tm2 + 0x8e;
		s1p48 = tmp + 0x90;		s1pc8 = tm2 + 0x90;
		s1p49 = tmp + 0x92;		s1pc9 = tm2 + 0x92;
		s1p4a = tmp + 0x94;		s1pca = tm2 + 0x94;
		s1p4b = tmp + 0x96;		s1pcb = tm2 + 0x96;
		s1p4c = tmp + 0x98;		s1pcc = tm2 + 0x98;
		s1p4d = tmp + 0x9a;		s1pcd = tm2 + 0x9a;
		s1p4e = tmp + 0x9c;		s1pce = tm2 + 0x9c;
		s1p4f = tmp + 0x9e;		s1pcf = tm2 + 0x9e;
		s1p50 = tmp + 0xa0;		s1pd0 = tm2 + 0xa0;
		s1p51 = tmp + 0xa2;		s1pd1 = tm2 + 0xa2;
		s1p52 = tmp + 0xa4;		s1pd2 = tm2 + 0xa4;
		s1p53 = tmp + 0xa6;		s1pd3 = tm2 + 0xa6;
		s1p54 = tmp + 0xa8;		s1pd4 = tm2 + 0xa8;
		s1p55 = tmp + 0xaa;		s1pd5 = tm2 + 0xaa;
		s1p56 = tmp + 0xac;		s1pd6 = tm2 + 0xac;
		s1p57 = tmp + 0xae;		s1pd7 = tm2 + 0xae;
		s1p58 = tmp + 0xb0;		s1pd8 = tm2 + 0xb0;
		s1p59 = tmp + 0xb2;		s1pd9 = tm2 + 0xb2;
		s1p5a = tmp + 0xb4;		s1pda = tm2 + 0xb4;
		s1p5b = tmp + 0xb6;		s1pdb = tm2 + 0xb6;
		s1p5c = tmp + 0xb8;		s1pdc = tm2 + 0xb8;
		s1p5d = tmp + 0xba;		s1pdd = tm2 + 0xba;
		s1p5e = tmp + 0xbc;		s1pde = tm2 + 0xbc;
		s1p5f = tmp + 0xbe;		s1pdf = tm2 + 0xbe;
		s1p60 = tmp + 0xc0;		s1pe0 = tm2 + 0xc0;
		s1p61 = tmp + 0xc2;		s1pe1 = tm2 + 0xc2;
		s1p62 = tmp + 0xc4;		s1pe2 = tm2 + 0xc4;
		s1p63 = tmp + 0xc6;		s1pe3 = tm2 + 0xc6;
		s1p64 = tmp + 0xc8;		s1pe4 = tm2 + 0xc8;
		s1p65 = tmp + 0xca;		s1pe5 = tm2 + 0xca;
		s1p66 = tmp + 0xcc;		s1pe6 = tm2 + 0xcc;
		s1p67 = tmp + 0xce;		s1pe7 = tm2 + 0xce;
		s1p68 = tmp + 0xd0;		s1pe8 = tm2 + 0xd0;
		s1p69 = tmp + 0xd2;		s1pe9 = tm2 + 0xd2;
		s1p6a = tmp + 0xd4;		s1pea = tm2 + 0xd4;
		s1p6b = tmp + 0xd6;		s1peb = tm2 + 0xd6;
		s1p6c = tmp + 0xd8;		s1pec = tm2 + 0xd8;
		s1p6d = tmp + 0xda;		s1ped = tm2 + 0xda;
		s1p6e = tmp + 0xdc;		s1pee = tm2 + 0xdc;
		s1p6f = tmp + 0xde;		s1pef = tm2 + 0xde;
		s1p70 = tmp + 0xe0;
		s1p71 = tmp + 0xe2;
		s1p72 = tmp + 0xe4;
		s1p73 = tmp + 0xe6;
		s1p74 = tmp + 0xe8;
		s1p75 = tmp + 0xea;
		s1p76 = tmp + 0xec;
		s1p77 = tmp + 0xee;
		s1p78 = tmp + 0xf0;
		s1p79 = tmp + 0xf2;
		s1p7a = tmp + 0xf4;
		s1p7b = tmp + 0xf6;
		s1p7c = tmp + 0xf8;
		s1p7d = tmp + 0xfa;
		s1p7e = tmp + 0xfc;
		s1p7f = tmp + 0xfe;
		tmp += 0x1e0;	// sc_ptr += 0x3c0
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;
		tmp += 0x1e;	// += 0x3c => sc_ptr + 0x3fc
		// DFT-roots:
		isrt2     = tmp + 0x00;
		cc0       = tmp + 0x01;
		ss0       = tmp + 0x02;
		sse2_c3m1 = tmp + 0x03;
		sse2_s    = tmp + 0x04;
		sse2_cn1  = tmp + 0x05;
		sse2_cn2  = tmp + 0x06;
		sse2_ss3  = tmp + 0x07;
		sse2_sn1  = tmp + 0x08;
		sse2_sn2  = tmp + 0x09;
		tmp += 0x0a;	// += 0xa => sc_ptr + 0x406
	  #ifdef USE_AVX
		cy_r = tmp;	cy_i = tmp+0x3c;	tmp += 2*0x3c;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x78 + 2 => sc_ptr += 0x480
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x78;	tmp += 2*0x78;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xf0 + 2 => sc_ptr += 0x4f8
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(HERE, half_arr_offset240 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT(HERE, (radix240_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)r00) + (20 << l2_sz_vd), "radix240_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(isrt2,ISRT2);
		VEC_DBL_INIT(cc0  ,  c16);
		VEC_DBL_INIT(ss0  ,  s16);
		VEC_DBL_INIT(sse2_c3m1, c3m1);
		VEC_DBL_INIT(sse2_s   , s   );
		VEC_DBL_INIT(sse2_cn1 , cn1 );
		VEC_DBL_INIT(sse2_cn2 , cn2 );
		VEC_DBL_INIT(sse2_ss3 , ss3 );
		VEC_DBL_INIT(sse2_sn1 , sn1 );
		VEC_DBL_INIT(sse2_sn2 , sn2 );

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)cy_r - (int)isrt2;	// #bytes in 1st of above block of consts
		tmp = isrt2;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = sz_vd;	// sse2_rnd is a solo (in the SIMD-vector) datum
		tmp = sse2_rnd;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.

		In 4-way SIMD (AVX) mode, we expand this from 2^2 2-vector table entries to 2^4 4-vector entries.
		*/
		tmp = half_arr;

	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

	  #if HIACC
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		// Simple qfloat-based loop to crank out the roots which populate the radix240_avx_negadwt_consts table:
			struct qfloat qc,qs,qx,qy,qt,qn,qmul;
			qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
			qc = qfcos(qt);	qs = qfsin(qt);
			qx = QONE;		qy = QZRO;
			for(j = 0; j < RADIX; j++) {
				printf("j = %3u: cos = 0x%16llX\n",j,qfdbl_as_uint64(qx));
				// Up-multiply the complex exponential:
				qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
				qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
			}
			exit(0);
		#endif

		tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
		// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
		tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = radix240_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix240_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix240_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix240_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix240_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix240_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix240_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*sz_vd/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix240_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix240_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix240_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix240_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix240_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix240_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix240_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix240_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 8;	// reset to point to start of above block
		nbytes = 4*sz_vd;	// 2 AVX-register-sized complex data

	  #endif	// HIACC toggle

		// Propagate the above consts to the remaining threads:
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// AVX?
	}
	else
	{
	#ifdef USE_AVX
		/* Forward-weight multipliers: */
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;

		nbytes = 64 << l2_sz_vd;

	#elif defined(USE_SSE2)

		ctmp = (struct complex *)tmp;
		/* Forward-weight multipliers: */
		ctmp->re = 1.0;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = .50;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = 1.0;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		ctmp->re = .25;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .25;	++ctmp;
		ctmp->re = .25;	ctmp->im = .25;	++ctmp;
		/* Forward-base[] multipliers: */
		ctmp->re = base   [0];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [0];	ctmp->im = base   [1];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [1];	++ctmp;
		/* Inverse-base[] multipliers: */
		ctmp->re = baseinv[0];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[0];	ctmp->im = baseinv[1];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[1];	++ctmp;

		nbytes = 16 << l2_sz_vd;
	#endif

		// Propagate the above consts to the remaining threads:
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	}

	/* Floating-point sign mask used for FABS on packed doubles: */
	sign_mask = sm_ptr;
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sign_mask+i) = (uint64)0x7FFFFFFFFFFFFFFFull;
	}
	// Set up the SIMD-tupled-32-bit-int SSE constants used by the carry macros:
	sse_bw  = sm_ptr + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
	tmp64 = (uint64)bw;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_bw+i) = tmp64;
	}

	sse_sw  = sse_bw + RE_IM_STRIDE;
	tmp64 = (uint64)sw;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_sw+i) = tmp64;
	}

	sse_n   = sse_sw + RE_IM_STRIDE;
	tmp64 = (uint64)n;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_n +i) = tmp64;
	}

	nbytes = 4 << l2_sz_vd;

#ifdef USE_AVX
	n_minus_sil   = (struct uint32x4 *)sse_n + 1;
	n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
	sinwt         = (struct uint32x4 *)sse_n + 3;
	sinwtm1       = (struct uint32x4 *)sse_n + 4;
	nbytes += 64;;
#endif

	// Propagate the above consts to the remaining threads:
	tmp = (vec_dbl *)sm_ptr;
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, nbytes);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

	// For large radices, array-access to bjmodn means only init base-ptr here:
	#ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	#endif

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
	#ifndef MULTITHREAD
		// Needed for the length-16 loops which manage the DFT16 and carry-macro calls in the compact-object-code build:
		po_br[0x0] = 0; po_br[0x1] = p8; po_br[0x2] = p4; po_br[0x3] = pc; po_br[0x4] = p2; po_br[0x5] = pa; po_br[0x6] = p6; po_br[0x7] = pe; po_br[0x8] = p1; po_br[0x9] = p9; po_br[0xa] = p5; po_br[0xb] = pd; po_br[0xc] = p3; po_br[0xd] = pb; po_br[0xe] = p7; po_br[0xf] = pf;
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3; plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7; plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb; plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
	#endif
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );

		poff[     0] =   0; poff[     1] =     p4; poff[     2] =     p8; poff[     3] =     pc;
		poff[0x04+0] = p10; poff[0x04+1] = p10+p4; poff[0x04+2] = p10+p8; poff[0x04+3] = p10+pc;
		poff[0x08+0] = p20; poff[0x08+1] = p20+p4; poff[0x08+2] = p20+p8; poff[0x08+3] = p20+pc;
		poff[0x0c+0] = p30; poff[0x0c+1] = p30+p4; poff[0x0c+2] = p30+p8; poff[0x0c+3] = p30+pc;
		poff[0x10+0] = p40; poff[0x10+1] = p40+p4; poff[0x10+2] = p40+p8; poff[0x10+3] = p40+pc;
		poff[0x14+0] = p50; poff[0x14+1] = p50+p4; poff[0x14+2] = p50+p8; poff[0x14+3] = p50+pc;
		poff[0x18+0] = p60; poff[0x18+1] = p60+p4; poff[0x18+2] = p60+p8; poff[0x18+3] = p60+pc;
		poff[0x1c+0] = p70; poff[0x1c+1] = p70+p4; poff[0x1c+2] = p70+p8; poff[0x1c+3] = p70+pc;
		poff[0x20+0] = p80; poff[0x20+1] = p80+p4; poff[0x20+2] = p80+p8; poff[0x20+3] = p80+pc;
		poff[0x24+0] = p90; poff[0x24+1] = p90+p4; poff[0x24+2] = p90+p8; poff[0x24+3] = p90+pc;
		poff[0x28+0] = pa0; poff[0x28+1] = pa0+p4; poff[0x28+2] = pa0+p8; poff[0x28+3] = pa0+pc;
		poff[0x2c+0] = pb0; poff[0x2c+1] = pb0+p4; poff[0x2c+2] = pb0+p8; poff[0x2c+3] = pb0+pc;
		poff[0x30+0] = pc0; poff[0x30+1] = pc0+p4; poff[0x30+2] = pc0+p8; poff[0x30+3] = pc0+pc;
		poff[0x34+0] = pd0; poff[0x34+1] = pd0+p4; poff[0x34+2] = pd0+p8; poff[0x34+3] = pd0+pc;
		poff[0x38+0] = pe0; poff[0x38+1] = pe0+p4; poff[0x38+2] = pe0+p8; poff[0x38+3] = pe0+pc;

	#ifndef MULTITHREAD
		p_in_hi [0x0] =   0; p_in_hi [0x1] = p20; p_in_hi [0x2] = p10; p_in_hi [0x3] = pe0; p_in_hi [0x4] = pd0; p_in_hi [0x5] = pc0; p_in_hi [0x6] = pb0; p_in_hi [0x7] = pa0; p_in_hi [0x8] = p90; p_in_hi [0x9] = p80; p_in_hi [0xa] = p70; p_in_hi [0xb] = p60; p_in_hi [0xc] = p50; p_in_hi [0xd] = p40; p_in_hi [0xe] = p30;
		p_out_hi[0x0] =   0; p_out_hi[0x1] = p10; p_out_hi[0x2] = p20; p_out_hi[0x3] = p30; p_out_hi[0x4] = p40; p_out_hi[0x5] = p50; p_out_hi[0x6] = p60; p_out_hi[0x7] = p70; p_out_hi[0x8] = p80; p_out_hi[0x9] = p90; p_out_hi[0xa] = pa0; p_out_hi[0xb] = pb0; p_out_hi[0xc] = pc0; p_out_hi[0xd] = pd0; p_out_hi[0xe] = pe0;
	#endif

		if(_cy_r[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;
			for(i = 0; i < RADIX; i++) {
				free((void *)_bjmodn[i]); _bjmodn[i] = 0x0;
				free((void *)  _cy_r[i]);   _cy_r[i] = 0x0;
				free((void *)  _cy_i[i]);   _cy_i[i] = 0x0;
			}
			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		for(i = 0; i < RADIX; i++) {
			_bjmodn[i]	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn[i]== 0x0);
		}
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		for(i = 0; i < RADIX; i++) {
			_cy_r[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r[i]== 0x0);
			_cy_i[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i[i]== 0x0);
		}
		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays!");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in %s.\n", func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = NDIVR/CY_THREADS;
		}
		else
		{
			jhi = NDIVR/CY_THREADS/2;
		}

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

		// In non-power-of-2-runlength case, both Mersenne and Fermat-mod share these next 2 loops:
		if(CY_THREADS > 1)
		{
			for(ithread = 1; ithread < CY_THREADS; ithread++)
			{
				_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
			}
		}

		// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
		j = _bjmodnini[CY_THREADS];
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_bjmodn[0][ithread] = _bjmodnini[ithread];
			for(i = 1; i < RADIX; i++) {
				MOD_ADD32(_bjmodn[i-1][ithread], j, n, _bjmodn[i][ithread]);
			}

			// Every (ODD_RADIX)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			{
				/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
				fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
				*/
				for(i = 0; i < RADIX; i += ODD_RADIX) {
					_bjmodn[i][ithread] = n;
				}
			}
		}

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
			so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
			*/
			/* indices into IBDWT weights arrays (mod NWT) is here: */
			ii[0]= 0;
			ii[1]= (SW_DIV_N*(NDIVR >> 1)) % nwt;	// nwt *not* a power of 2, must use library-mod!
			for(i = 2; i < ODD_RADIX; i++) {
				MOD_ADD32(ii[i-1], ii[1], nwt, ii[i]);
			}

			/* Find the circular-index-shift (cf. the head-of-file comments of radix28_ditN_cy_dif1.c) by searching bjmodn01 ... bjmodn[nwt] for the one == bw: */
			for(i = 1; i < ODD_RADIX; i++) {
				if( _bjmodn[i][0] == bw ) {
					wts_idx_incr = i;
					break;
				};
			}
			ASSERT(HERE, wts_idx_incr != 0, "wts_idx_incr init failed!");

		#ifdef USE_SSE2
			wts_idx_inc2 = wts_idx_incr << (2*l2_sz_vd - 3);	/* In the SIMD version, use icycle0-6 as actual address
							offsets, so wts_idx_incr includes a *sizeof(vec_dbl) for the array-of-vector-doubles indexing, and another
							doubling|quadrupling|... to reflect the fact that the SIMD version of the loop is equivalent to 2|4|... scalar
							loop executions, i.e. corresponds to [#doubles in each vec_dbl] scalar-code increments of the icycle indices. */
			wts_idx_inc2 %= nwt16;	/* Need an extra mod since [2|4|...]*wts_idx_incr may be >= nwt */
		#endif
			/* Subtract nwt from the increments to ease fast-mod */
			wts_idx_incr -= nwt;
		#ifdef USE_SSE2
			wts_idx_inc2 -= nwt16;
		#endif

			for(i = 0; i < ODD_RADIX; i++) {
				/* Need this both in scalar mode and to ease the SSE2-array init */
				j = _bjmodn[i][0] > sw;	bs_arr[i] = base[j];	bsinv_arr[i] = baseinv[j];
				wt_arr[i] = wt0[ii[i]];	// inverse wts must be reinited on each pass, since these have a *scale multiplier
				/* Give icycle indices their proper starting values: */
				icycle[i] = i;
			}

		#ifdef USE_SSE2
			tmp = half_arr;
			for(i = 0; i < ODD_RADIX; i++, tmp++) {
				tmp->d0 = wt_arr[icycle[i]];
			/* Now set the imaginary parts to the values corresponding to the 2nd of each pair of scalar-mode loop passes.
			Use this sequence for mod-add, as it is faster than general-mod '% nwt': */
				jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
				tmp->d1 = wt_arr[jcycle[i]];
		  #ifdef USE_AVX
				kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
				tmp->d2 = wt_arr[kcycle[i]];
				lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
				tmp->d3 = wt_arr[lcycle[i]];
		  #endif
			}

			// Propagate the above wts-consts to the remaining threads:
			nbytes = ODD_RADIX*sz_vd;
			tmp = half_arr;
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

			tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */

		  #ifdef USE_AVX

			// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-4 index subsequences, which will be accessed via their head element
			// by the [ijkl]cycle* index quartets in the respective carry-macro call:
			tmp->d0 = bs_arr[0x0];	tmp->d1 = bs_arr[0x1];	tmp->d2 = bs_arr[0x2];	tmp->d3 = bs_arr[0x3];	++tmp;
			tmp->d0 = bs_arr[0x1];	tmp->d1 = bs_arr[0x2];	tmp->d2 = bs_arr[0x3];	tmp->d3 = bs_arr[0x4];	++tmp;
			tmp->d0 = bs_arr[0x2];	tmp->d1 = bs_arr[0x3];	tmp->d2 = bs_arr[0x4];	tmp->d3 = bs_arr[0x5];	++tmp;
			tmp->d0 = bs_arr[0x3];	tmp->d1 = bs_arr[0x4];	tmp->d2 = bs_arr[0x5];	tmp->d3 = bs_arr[0x6];	++tmp;
			tmp->d0 = bs_arr[0x4];	tmp->d1 = bs_arr[0x5];	tmp->d2 = bs_arr[0x6];	tmp->d3 = bs_arr[0x7];	++tmp;
			tmp->d0 = bs_arr[0x5];	tmp->d1 = bs_arr[0x6];	tmp->d2 = bs_arr[0x7];	tmp->d3 = bs_arr[0x8];	++tmp;
			tmp->d0 = bs_arr[0x6];	tmp->d1 = bs_arr[0x7];	tmp->d2 = bs_arr[0x8];	tmp->d3 = bs_arr[0x9];	++tmp;
			tmp->d0 = bs_arr[0x7];	tmp->d1 = bs_arr[0x8];	tmp->d2 = bs_arr[0x9];	tmp->d3 = bs_arr[0xa];	++tmp;
			tmp->d0 = bs_arr[0x8];	tmp->d1 = bs_arr[0x9];	tmp->d2 = bs_arr[0xa];	tmp->d3 = bs_arr[0xb];	++tmp;
			tmp->d0 = bs_arr[0x9];	tmp->d1 = bs_arr[0xa];	tmp->d2 = bs_arr[0xb];	tmp->d3 = bs_arr[0xc];	++tmp;
			tmp->d0 = bs_arr[0xa];	tmp->d1 = bs_arr[0xb];	tmp->d2 = bs_arr[0xc];	tmp->d3 = bs_arr[0xd];	++tmp;
			tmp->d0 = bs_arr[0xb];	tmp->d1 = bs_arr[0xc];	tmp->d2 = bs_arr[0xd];	tmp->d3 = bs_arr[0xe];	++tmp;
			tmp->d0 = bs_arr[0xc];	tmp->d1 = bs_arr[0xd];	tmp->d2 = bs_arr[0xe];	tmp->d3 = bs_arr[0x0];	++tmp;
			tmp->d0 = bs_arr[0xd];	tmp->d1 = bs_arr[0xe];	tmp->d2 = bs_arr[0x0];	tmp->d3 = bs_arr[0x1];	++tmp;
			tmp->d0 = bs_arr[0xe];	tmp->d1 = bs_arr[0x0];	tmp->d2 = bs_arr[0x1];	tmp->d3 = bs_arr[0x2];	++tmp;

			tmp->d0 = bsinv_arr[0x0];	tmp->d1 = bsinv_arr[0x1];	tmp->d2 = bsinv_arr[0x2];	tmp->d3 = bsinv_arr[0x3];	++tmp;
			tmp->d0 = bsinv_arr[0x1];	tmp->d1 = bsinv_arr[0x2];	tmp->d2 = bsinv_arr[0x3];	tmp->d3 = bsinv_arr[0x4];	++tmp;
			tmp->d0 = bsinv_arr[0x2];	tmp->d1 = bsinv_arr[0x3];	tmp->d2 = bsinv_arr[0x4];	tmp->d3 = bsinv_arr[0x5];	++tmp;
			tmp->d0 = bsinv_arr[0x3];	tmp->d1 = bsinv_arr[0x4];	tmp->d2 = bsinv_arr[0x5];	tmp->d3 = bsinv_arr[0x6];	++tmp;
			tmp->d0 = bsinv_arr[0x4];	tmp->d1 = bsinv_arr[0x5];	tmp->d2 = bsinv_arr[0x6];	tmp->d3 = bsinv_arr[0x7];	++tmp;
			tmp->d0 = bsinv_arr[0x5];	tmp->d1 = bsinv_arr[0x6];	tmp->d2 = bsinv_arr[0x7];	tmp->d3 = bsinv_arr[0x8];	++tmp;
			tmp->d0 = bsinv_arr[0x6];	tmp->d1 = bsinv_arr[0x7];	tmp->d2 = bsinv_arr[0x8];	tmp->d3 = bsinv_arr[0x9];	++tmp;
			tmp->d0 = bsinv_arr[0x7];	tmp->d1 = bsinv_arr[0x8];	tmp->d2 = bsinv_arr[0x9];	tmp->d3 = bsinv_arr[0xa];	++tmp;
			tmp->d0 = bsinv_arr[0x8];	tmp->d1 = bsinv_arr[0x9];	tmp->d2 = bsinv_arr[0xa];	tmp->d3 = bsinv_arr[0xb];	++tmp;
			tmp->d0 = bsinv_arr[0x9];	tmp->d1 = bsinv_arr[0xa];	tmp->d2 = bsinv_arr[0xb];	tmp->d3 = bsinv_arr[0xc];	++tmp;
			tmp->d0 = bsinv_arr[0xa];	tmp->d1 = bsinv_arr[0xb];	tmp->d2 = bsinv_arr[0xc];	tmp->d3 = bsinv_arr[0xd];	++tmp;
			tmp->d0 = bsinv_arr[0xb];	tmp->d1 = bsinv_arr[0xc];	tmp->d2 = bsinv_arr[0xd];	tmp->d3 = bsinv_arr[0xe];	++tmp;
			tmp->d0 = bsinv_arr[0xc];	tmp->d1 = bsinv_arr[0xd];	tmp->d2 = bsinv_arr[0xe];	tmp->d3 = bsinv_arr[0x0];	++tmp;
			tmp->d0 = bsinv_arr[0xd];	tmp->d1 = bsinv_arr[0xe];	tmp->d2 = bsinv_arr[0x0];	tmp->d3 = bsinv_arr[0x1];	++tmp;
			tmp->d0 = bsinv_arr[0xe];	tmp->d1 = bsinv_arr[0x0];	tmp->d2 = bsinv_arr[0x1];	tmp->d3 = bsinv_arr[0x2];	++tmp;

		  #else

			/* In SSE2 mode, because we apply doubled weights to data arranged as [a.re,b.re,...],[a.im,b.im,...] but apply
			doubled base multipliers to shuffled data [a.re,a.im],[b.re,b.im],... (i.e. shuffled to yield same data layout as
			in the scalar case), the weights need to have disparate real and imag parts, but the base/baseinv terms do not: */
			for(i = 0; i < ODD_RADIX; i++) {
				VEC_DBL_INIT(tmp, bs_arr[i]);	++tmp;
			}
			for(i = 0; i < ODD_RADIX; i++) {
				VEC_DBL_INIT(tmp, bsinv_arr[i]);	++tmp;
			}

		  #endif

			// Propagate the above consts to the remaining threads:
			nbytes <<= 1;	// [base+binv] ==> 2x as many consts as [wts], since the wtinv data done each pass of outer-loop
			tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] <<= l2_sz_vd;		jcycle[i] <<= l2_sz_vd;
			#ifdef USE_AVX
				kcycle[i] <<= l2_sz_vd;		lcycle[i] <<= l2_sz_vd;
			#endif
			}

		#endif	// USE_SSE2 ?
		}	// MODULUS_TYPE_FERMAT ?

	#ifdef USE_PTHREAD

		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r00      = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].r00      = (double *)foo_array;
			tdat[ithread].half_arr = (double *)&wts_idx_incr;
		#endif	// USE_SSE2
		}

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			// These inits must occur just once, in loop-simulated full-pass mode,
			// in order to get the array-index-offset values of the icycle/jcycle indices right:
			for(i = 0; i < ODD_RADIX; i++) {
				tdat[0].icycle[i] = icycle[i];
			#ifdef USE_SSE2
				tdat[0].wts_idx_inc2 = wts_idx_inc2;
				tdat[0].jcycle[i] = jcycle[i];
			  #ifdef USE_AVX
				tdat[0].kcycle[i] = kcycle[i];
				tdat[0].lcycle[i] = lcycle[i];
			  #endif
			//	printf("Thread 0 idx pair[%2d] = [%2d,%2d]\n",i,icycle[i],jcycle[i]);
			#endif
			}
			// For remaining threads, simulate the loop-evolution of the above indices.
			// Note that the non-thread-associated *cycle[] arry data will get changed fom their above-inited
			// values in the loop here, but that's OK because in || mode only the thread-associated values matter:
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				jstart = 0;
				jhi = NDIVR/CY_THREADS;	// The earlier setting = NDIVR/CY_THREADS/2 was for simulating bjmodn evolution, must double that here
				// khi = 1 for Fermat-mod, thus no outer loop needed here
				for(j = jstart; j < jhi; j += stride)
				{
					for(i = 0; i < ODD_RADIX; i++) {
					#ifndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:
						icycle[i] += wts_idx_incr;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt);
					#else
						icycle[i] += wts_idx_inc2;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt16);
						jcycle[i] += wts_idx_inc2;		jcycle[i] += ( (-(int)((uint32)jcycle[i] >> 31)) & nwt16);
					  #ifdef USE_AVX
						kcycle[i] += wts_idx_inc2;		kcycle[i] += ( (-(int)((uint32)kcycle[i] >> 31)) & nwt16);
						lcycle[i] += wts_idx_inc2;		lcycle[i] += ( (-(int)((uint32)lcycle[i] >> 31)) & nwt16);
					  #endif
					#endif
					}
				}
				for(i = 0; i < ODD_RADIX; i++) {
					tdat[ithread].icycle[i] = icycle[i];
				#ifdef USE_SSE2
					tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
					tdat[ithread].jcycle[i] = jcycle[i];
				  #ifdef USE_AVX
					tdat[ithread].kcycle[i] = kcycle[i];
					tdat[ithread].lcycle[i] = lcycle[i];
				  #endif
			//	printf("Thread %d idx pair[%2d] = [%2d,%2d]\n",ithread,i,icycle[i],jcycle[i]);
				#endif
				}
			}
			// Restore the original loop-start values of the cycle arrays, since we use these for init of inv-wts below:
			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] = tdat[0].icycle[i];
			#ifdef USE_SSE2
				jcycle[i] = tdat[0].jcycle[i];
			  #ifdef USE_AVX
				kcycle[i] = tdat[0].kcycle[i];
				lcycle[i] = tdat[0].lcycle[i];
			  #endif
			#endif
			}
		}

	#endif	// USE_PTHREAD

	}	/* endif(first_entry) */

/*...The radix-240 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy_r[i][ithread] = 0;
			_cy_i[i][ithread] = 0;
		}
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r[0][0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/
		khi = 1;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}

		// Now that full_pass-dependent scale factor known, init inverse weights tiny-table used for Fermat-mod
		for(i = 0; i < ODD_RADIX; i++) {
			wtinv_arr[i] = scale*wt1[ii[i]];
		}

	// In threaded mode, the master *cycle[] values are unmodified during main loop exec; only the thread-associated
	// copies of these index arrays get modified. In non-threaded mode we must separately store copies of the masters
	// in order to so;ve the save/restore issue. We start from the (static, unmodified during loop) ii[]-index values:
	#ifndef MULTITHREAD
		for(i = 0; i < ODD_RADIX; i++) {
			/* Reinit *cycle indices their proper starting values - recall in SIMD mode these all are ( << 4): */
			icycle[i] = i;
		#ifdef USE_SSE2
			jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
		  #ifdef USE_AVX
			kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
			lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
			kcycle[i] <<= l2_sz_vd;		lcycle[i] <<= l2_sz_vd;
		  #endif
			icycle[i] <<= l2_sz_vd;		jcycle[i] <<= l2_sz_vd;
		#endif
		}
	#endif

	#ifdef USE_SSE2
		// Remember: *cycle[] entries all << l2_sz_vd here - must left-shift-on-the-fly before using:
		tm2 = half_arr + ODD_RADIX;
		for(i = 0; i < ODD_RADIX; i++, tm2++) {
			tm2->d0 = wtinv_arr[icycle[i] >> l2_sz_vd];
			tm2->d1 = wtinv_arr[jcycle[i] >> l2_sz_vd];
		#ifdef USE_AVX
			tm2->d2 = wtinv_arr[kcycle[i] >> l2_sz_vd];
			tm2->d3 = wtinv_arr[lcycle[i] >> l2_sz_vd];
		#endif
		}

		// Propagate the above inv-wts to the remaining threads - surrounding consts are unchanged:
		nbytes = ODD_RADIX*sz_vd;
		tmp = half_arr + ODD_RADIX;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// USE_SSE2 ?
	}	// MODULUS_TYPE_FERMAT ?

#if defined(USE_SSE2) && defined(USE_PTHREAD)

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, sz_vd);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

#endif	// USE_PTHREAD

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		tdat[ithread].iter = iter;
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
				tdat[ithread].cy_i[i] = _cy_i[i][ithread];
			}
		}
	}
#endif

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
	#ifdef USE_SSE2
	//	VEC_DBL_INIT(max_err, 0.0);	*** must do this in conjunction with thread-local-data-copy
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			col = _col[ithread];
			co2 = _co2[ithread];
			co3 = _co3[ithread];

			for(l = 0; l < RADIX; l++) {
				bjmodn[l] = _bjmodn[l][ithread];
			}
			/* init carries	*/
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];
			}
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
			}
		#else
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];
			}
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				tmp->d0 = _cy_r[l  ][ithread];		tm2->d0 = _cy_i[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];		tm2->d1 = _cy_i[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];		tm2->d2 = _cy_i[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];		tm2->d3 = _cy_i[l+3][ithread];
			}
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				tmp->d0 = _cy_r[l][ithread];		tmp->d1 = _cy_i[l][ithread];
			}
		#else
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];		cy_i[l] = _cy_i[l][ithread];
			}
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix240_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
				_cy_r[l+2][ithread] = tmp->d2;
				_cy_r[l+3][ithread] = tmp->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];
			}
		#endif
		}
		else
		{
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				_cy_r[l  ][ithread] = tmp->d0;		_cy_i[l  ][ithread] = tm2->d0;
				_cy_r[l+1][ithread] = tmp->d1;		_cy_i[l+1][ithread] = tm2->d1;
				_cy_r[l+2][ithread] = tmp->d2;		_cy_i[l+2][ithread] = tm2->d2;
				_cy_r[l+3][ithread] = tmp->d3;		_cy_i[l+3][ithread] = tm2->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				_cy_r[l][ithread] = tmp->d0;		_cy_i[l][ithread] = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];		_cy_i[l][ithread] = cy_i[l];
			}
		#endif
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
		ASSERT(HERE, 0x0 == cy240_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = tdat[ithread].cy_r[l];
			}
		}
		else
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = tdat[ithread].cy_r[l];
				_cy_i[l][ithread] = tdat[ithread].cy_i[l];
			}
		}
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/RADIX set of contiguous data into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the forward DIF FFT of the first block of RADIX complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the RADIX outputs of (1);
	(3) Reweight and perform a forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next RADIX elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		for(l = 0; l < RADIX; l++) {
			t[l].re = _cy_r[l][CY_THREADS - 1];
		}
		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = _cy_r[l][ithread-1];
			}
		}
		_cy_r[0][0] =+t[RADIX-1].re;	/* ...The wraparound carry is here: */
		for(l = 1; l < RADIX; l++) {
			_cy_r[l][0] = t[l-1].re;
		}
	}
	else
	{
		j = CY_THREADS - 1;
		for(l = 0; l < RADIX; l++) {
			t[l].re = _cy_r[l][j];		t[l].im = _cy_i[l][j];
		}
		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = _cy_r[l][ithread-1];	_cy_i[l][ithread] = _cy_i[l][ithread-1];
			}
		}
		_cy_r[0][0] =-t[RADIX-1].im;	_cy_i[0][0] =+t[RADIX-1].re;	/* ...The 2 Mo"bius carries are here: */
		for(l = 1; l < RADIX; l++) {
			_cy_r[l][0] = t[l-1].re;	_cy_i[l][0] = t[l-1].im;
		}
	}

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + poff[l];	// poff[] = p0,4,8,...
				a[jt   ] *= radix_inv;
				a[jt+p1] *= radix_inv;
				a[jt+p2] *= radix_inv;
				a[jt+p3] *= radix_inv;
			}
		}
	}
}	/* endfor(outer) */

	dtmp = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(l = 0; l < RADIX; l++) {
			dtmp += fabs(_cy_r[l][ithread]) + fabs(_cy_i[l][ithread]);
		}
		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}
	if(dtmp != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
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

/****************/

void radix240_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-240 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int i, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0;
	static int plo[16], p_in_hi[ODD_RADIX], p_out_hi[ODD_RADIX];
	uint64 i64;
	const uint64 perm16[ODD_RADIX] = {
		0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x76450123fecd89abull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,0x32017645ba89fecdull,
		0xefdc98ba10325476ull,0x45672310cdefab98ull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull};
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in-place DFT macros:
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT(HERE, (double *)t == &(t[0x00].re), "Unexpected value for Tmp-array-start pointer!");
		first_entry=FALSE;
		NDIVR=n/RADIX;

		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3; plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7; plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb; plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p_in_hi [0x0] =   0; p_in_hi [0x1] = p20; p_in_hi [0x2] = p10; p_in_hi [0x3] = pe0; p_in_hi [0x4] = pd0; p_in_hi [0x5] = pc0; p_in_hi [0x6] = pb0; p_in_hi [0x7] = pa0; p_in_hi [0x8] = p90; p_in_hi [0x9] = p80; p_in_hi [0xa] = p70; p_in_hi [0xb] = p60; p_in_hi [0xc] = p50; p_in_hi [0xd] = p40; p_in_hi [0xe] = p30;
		p_out_hi[0x0] =   0; p_out_hi[0x1] = p10; p_out_hi[0x2] = p20; p_out_hi[0x3] = p30; p_out_hi[0x4] = p40; p_out_hi[0x5] = p50; p_out_hi[0x6] = p60; p_out_hi[0x7] = p70; p_out_hi[0x8] = p80; p_out_hi[0x9] = p90; p_out_hi[0xa] = pa0; p_out_hi[0xb] = pb0; p_out_hi[0xc] = pc0; p_out_hi[0xd] = pd0; p_out_hi[0xe] = pe0;
	}

/*...The radix-240 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (240 64-bit complex) and do 6 radix-15 transforms...
	Twiddleless version arranges 16 sets of radix-15 DFT inputs as follows:
	0 in upper left corner, decrement 16 horizontally and 15 vertically, indexing modulo 240
	(we can auto-generate these by compiling test_fft_radix.c with -DTTYPE=0 -DRADIX=240, running
	the resulting executable and snarfing the first set of index-outputs, "DIT input-scramble array"):

		RADIX_15_DFT(00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10)
		RADIX_15_DFT(e1,d1,c1,b1,a1,91,81,71,61,51,41,31,21,11,01)
		RADIX_15_DFT(d2,c2,b2,a2,92,82,72,62,52,42,32,22,12,02,e2)
		RADIX_15_DFT(c3,b3,a3,93,83,73,63,53,43,33,23,13,03,e3,d3)
		RADIX_15_DFT(b4,a4,94,84,74,64,54,44,34,24,14,04,e4,d4,c4)
		RADIX_15_DFT(a5,95,85,75,65,55,45,35,25,15,05,e5,d5,c5,b5)
		RADIX_15_DFT(96,86,76,66,56,46,36,26,16,06,e6,d6,c6,b6,a6)
		RADIX_15_DFT(87,77,67,57,47,37,27,17,07,e7,d7,c7,b7,a7,97)
		RADIX_15_DFT(78,68,58,48,38,28,18,08,e8,d8,c8,b8,a8,98,88)
		RADIX_15_DFT(69,59,49,39,29,19,09,e9,d9,c9,b9,a9,99,89,79)
		RADIX_15_DFT(5a,4a,3a,2a,1a,0a,ea,da,ca,ba,aa,9a,8a,7a,6a)
		RADIX_15_DFT(4b,3b,2b,1b,0b,eb,db,cb,bb,ab,9b,8b,7b,6b,5b)
		RADIX_15_DFT(3c,2c,1c,0c,ec,dc,cc,bc,ac,9c,8c,7c,6c,5c,4c)
		RADIX_15_DFT(2d,1d,0d,ed,dd,cd,bd,ad,9d,8d,7d,6d,5d,4d,3d)
		RADIX_15_DFT(1e,0e,ee,de,ce,be,ae,9e,8e,7e,6e,5e,4e,3e,2e)
		RADIX_15_DFT(0f,ef,df,cf,bf,af,9f,8f,7f,6f,5f,4f,3f,2f,1f)

	If we subtract costant offsets 1,2,3,...,e,f from the last 15 rows as we do in the implementation to reduce the number
	of index offsets needing to be stored, we decrement 16 horizontally and 16 vertically:

			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + 0
			e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00 + 1
			d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0 + 2
			c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0 + 3
			b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0 + 4
			a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0 + 5
			90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0 + 6
			80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90 + 7
			70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80 + 8
			60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70 + 9
			50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60 + a
			40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50 + b
			30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40 + c
			20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30 + d
			10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20 + e
			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + f

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in computable-index fashion. RHS is simply i = 0,..,f.
	For the LHS we need 2 ingredients:
		[1] A formula yielding the lead part of the leftmost term of each row: 0,e,d,...,1,0. Note (0xf-i) = f,e,d,...,1,0,
			so now just need apply a bitmask which == 0 if i==0 and all-ones otherwise, thus idx = (0xf-i) & (-(i > 0)) .
		[2] An efficient decrement (mod 15) scheme to yield the remaining leading parts of each row's elements:
			idx--; idx += (-(idx < 0))&15;
	*/
		tptr = t;
		for(i = 0; i < 16; i++) {
		// [1] here:
			k0 = plo[i];	// p0,..,f
			jt = j1 + k0; jp = j2 + k0;
		// [2] here:
			k0 = (15-i)&(-(i>0)); /* 0,e,d,...,0 */	// Now get the resulting p* offsets:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = p_out_hi[k0];
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = p_out_hi[k1];
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = p_out_hi[k2];
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = p_out_hi[k3];
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = p_out_hi[k4];
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = p_out_hi[k5];
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = p_out_hi[k6];
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = p_out_hi[k7];
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = p_out_hi[k8];
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = p_out_hi[k9];
			kb = ka-1; kb += (-(kb < 0))&15;		ka = p_out_hi[ka];
			kc = kb-1; kc += (-(kc < 0))&15;		kb = p_out_hi[kb];
			kd = kc-1; kd += (-(kd < 0))&15;		kc = p_out_hi[kc];
			ke = kd-1; ke += (-(ke < 0))&15;		kd = p_out_hi[kd];
													ke = p_out_hi[ke];
			RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
			);	tptr += 0xf;
		}

		/*...and now do 15 radix-16 transforms.
		The required output permutation is:
			00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f   0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f + p00
			25,24,27,26,23,22,20,21,2d,2c,2f,2e,2b,2a,28,29   5,4,7,6,3,2,0,1,d,c,f,e,b,a,8,9 + p20
			1a,1b,19,18,1e,1f,1d,1c,16,17,15,14,11,10,13,12   a,b,9,8,e,f,d,c,6,7,5,4,1,0,3,2 + p10
			e7,e6,e4,e5,e0,e1,e2,e3,ef,ee,ec,ed,e8,e9,ea,eb   7,6,4,5,0,1,2,3,f,e,c,d,8,9,a,b + pe0
			d9,d8,db,da,dd,dc,df,de,d5,d4,d7,d6,d3,d2,d0,d1   9,8,b,a,d,c,f,e,5,4,7,6,3,2,0,1 + pd0
			c2,c3,c1,c0,c6,c7,c5,c4,ca,cb,c9,c8,ce,cf,cd,cc   2,3,1,0,6,7,5,4,a,b,9,8,e,f,d,c + pc0
			bb,ba,b8,b9,bf,be,bc,bd,b7,b6,b4,b5,b0,b1,b2,b3   b,a,8,9,f,e,c,d,7,6,4,5,0,1,2,3 + pb0
			a1,a0,a3,a2,a5,a4,a7,a6,a9,a8,ab,aa,ad,ac,af,ae = 1,0,3,2,5,4,7,6,9,8,b,a,d,c,f,e + pa0
			9c,9d,9e,9f,9a,9b,99,98,92,93,91,90,96,97,95,94   c,d,e,f,a,b,9,8,2,3,1,0,6,7,5,4 + p90
			83,82,80,81,87,86,84,85,8b,8a,88,89,8f,8e,8c,8d   3,2,0,1,7,6,4,5,b,a,8,9,f,e,c,d + p80
			7e,7f,7d,7c,79,78,7b,7a,71,70,73,72,75,74,77,76   e,f,d,c,9,8,b,a,1,0,3,2,5,4,7,6 + p70
			64,65,66,67,62,63,61,60,6c,6d,6e,6f,6a,6b,69,68   4,5,6,7,2,3,1,0,c,d,e,f,a,b,9,8 + p60
			5d,5c,5f,5e,5b,5a,58,59,53,52,50,51,57,56,54,55   d,c,f,e,b,a,8,9,3,2,0,1,7,6,4,5 + p50
			46,47,45,44,41,40,43,42,4e,4f,4d,4c,49,48,4b,4a   6,7,5,4,1,0,3,2,e,f,d,c,9,8,b,a + p40
			38,39,3a,3b,3c,3d,3e,3f,34,35,36,37,32,33,31,30   8,9,a,b,c,d,e,f,4,5,6,7,2,3,1,0 + p30

		As in the [better-documented] radix240_dit_pass1, simply encode each 16-perm as a hex-char string.
		This means we must extract each p-offset in little-endian fashion, e.g. low 4 bits have rightmost p-offset above.
		*/
		tptr = t;
		for(i = 0; i < 15; i++) {
			i64 = perm16[i];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jt = j1 + p_in_hi[i]; jp = j2 + p_in_hi[i];	// p_in_hi[] = p0,p20,...,p30
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr++;
		}
	}
}

/**************/

void radix240_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-240 complex inverse DIT FFT pass on the data in the length-N real vector A.
*/
	int i, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0;
#if USE_COMPACT_OBJ_CODE
	static int plo[16], p_in_hi[ODD_RADIX], p_out_hi[ODD_RADIX];
	uint64 i64;
	const uint64 perm16[ODD_RADIX] = {
		0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x98abefcd10236745ull,
		0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,
		0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x89bafedc01327654ull};
#endif
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in-place DFT macros:
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT(HERE, (double *)t == &(t[0x00].re), "Unexpected value for Tmp-array-start pointer!");
		first_entry=FALSE;
		NDIVR=n/RADIX;

		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3; plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7; plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb; plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p_in_hi [0x0] =   0; p_in_hi [0x1] = p20; p_in_hi [0x2] = p10; p_in_hi [0x3] = pe0; p_in_hi [0x4] = pd0; p_in_hi [0x5] = pc0; p_in_hi [0x6] = pb0; p_in_hi [0x7] = pa0; p_in_hi [0x8] = p90; p_in_hi [0x9] = p80; p_in_hi [0xa] = p70; p_in_hi [0xb] = p60; p_in_hi [0xc] = p50; p_in_hi [0xd] = p40; p_in_hi [0xe] = p30;
		p_out_hi[0x0] =   0; p_out_hi[0x1] = p10; p_out_hi[0x2] = p20; p_out_hi[0x3] = p30; p_out_hi[0x4] = p40; p_out_hi[0x5] = p50; p_out_hi[0x6] = p60; p_out_hi[0x7] = p70; p_out_hi[0x8] = p80; p_out_hi[0x9] = p90; p_out_hi[0xa] = pa0; p_out_hi[0xb] = pb0; p_out_hi[0xc] = pc0; p_out_hi[0xd] = pd0; p_out_hi[0xe] = pe0;
	}

/*...The radix-240 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;
	/*
	Gather the needed data (240 64-bit complex) and do 15 radix-16 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array =
		00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08   0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00
		25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d   5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20
		1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16   a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p10
		e7,e6,e5,e4,e3,e2,e1,e0,eb,ea,e9,e8,ed,ec,ee,ef   7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + pe0
		d9,d8,da,db,de,df,dc,dd,d1,d0,d2,d3,d6,d7,d4,d5   9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + pd0
		c2,c3,c0,c1,c4,c5,c7,c6,cc,cd,cf,ce,c8,c9,cb,ca   2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + pc0
		bb,ba,b9,b8,bd,bc,be,bf,b3,b2,b1,b0,b5,b4,b6,b7   b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + pb0
		a1,a0,a2,a3,a6,a7,a4,a5,ae,af,ac,ad,aa,ab,a8,a9 = 1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0
		9c,9d,9f,9e,98,99,9b,9a,94,95,97,96,90,91,93,92   c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p90
		83,82,81,80,85,84,86,87,8d,8c,8e,8f,89,88,8a,8b   3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p80
		7e,7f,7c,7d,7a,7b,78,79,76,77,74,75,72,73,70,71   e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + p70
		64,65,67,66,60,61,63,62,68,69,6b,6a,6f,6e,6d,6c   4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p60
		5d,5c,5e,5f,59,58,5a,5b,55,54,56,57,51,50,52,53   d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p50
		46,47,44,45,42,43,40,41,4a,4b,48,49,4c,4d,4f,4e   6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p40
		38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34   8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p30

	...which has the same pattern of p*0 offsets as the outputs of the DIF, but different row-internal permutation patterns.

	In order to support a compact-object-code version of the 240-DIT we need to wrap the 15 radix-16
	subtransforms in a loop with some auxiliary-index-vector scheme telling each subDFT which of the
	above p-offset sets to use. First break each 16-perm into four 4-perms and see what pattern emerges.
	We index each of the 4 sub-4-perm quartet members here according to its leading index:

					16-perm					[0-3]-perm	[4-7]-perm	[8-b]-perm	[c-f]-perm	4-perms assembled:
		-------------------------------		----------	----------	----------	----------	------
		0,1,3,2 7,6,5,4 f,e,d,c b,a,9,8		0,1,3,2[0]	7,6,5,4[7]	b,a,9,8[b]	f,e,d,c[d]	[07fb]
		5,4,6,7 1,0,2,3 9,8,a,b e,f,c,d		1,0,2,3[1]	5,4,6,7[5]	9,8,a,b[9]	e,f,c,d[c]	[519e]
		a,b,8,9 c,d,f,e 2,3,0,1 4,5,7,6		2,3,0,1[2]	4,5,7,6[4]	a,b,8,9[a]	c,d,f,e[a]	[ac24]
		7,6,5,4 3,2,1,0 b,a,9,8 d,c,e,f		3,2,1,0[3]	7,6,5,4[7]	b,a,9,8[b]	d,c,e,f[b]	[73bd]
		9,8,a,b e,f,c,d 1,0,2,3 6,7,4,5		1,0,2,3[1]	6,7,4,5[6]	9,8,a,b[9]	e,f,c,d[c]	[9e16]
		2,3,0,1 4,5,7,6 c,d,f,e 8,9,b,a		2,3,0,1[2]	4,5,7,6[4]	8,9,b,a[8]	c,d,f,e[a]	[24c8]
		b,a,9,8 d,c,e,f 3,2,1,0 5,4,6,7		3,2,1,0[3]	5,4,6,7[5]	b,a,9,8[b]	d,c,e,f[b]	[bd35]
		1,0,2,3 6,7,4,5 e,f,c,d a,b,8,9		1,0,2,3[1]	6,7,4,5[6]	a,b,8,9[a]	e,f,c,d[c]	[16ea]
		c,d,f,e 8,9,b,a 4,5,7,6 0,1,3,2		0,1,3,2[0]	4,5,7,6[4]	8,9,b,a[8]	c,d,f,e[a]	[c840]
		3,2,1,0 5,4,6,7 d,c,e,f 9,8,a,b		3,2,1,0[3]	5,4,6,7[5]	9,8,a,b[9]	d,c,e,f[b]	[35d9]
		e,f,c,d a,b,8,9 6,7,4,5 2,3,0,1		2,3,0,1[2]	6,7,4,5[6]	a,b,8,9[a]	e,f,c,d[c]	[ea62]
		4,5,7,6 0,1,3,2 8,9,b,a f,e,d,c		0,1,3,2[0]	4,5,7,6[4]	8,9,b,a[8]	f,e,d,c[d]	[408f]
		d,c,e,f 9,8,a,b 5,4,6,7 1,0,2,3		1,0,2,3[1]	5,4,6,7[5]	9,8,a,b[9]	d,c,e,f[b]	[d951]
		6,7,4,5 2,3,0,1 a,b,8,9 c,d,f,e		2,3,0,1[2]	6,7,4,5[6]	a,b,8,9[a]	c,d,f,e[a]	[62ac]
		8,9,b,a f,e,d,c 0,1,3,2 7,6,5,4		0,1,3,2[0]	7,6,5,4[7]	8,9,b,a[8]	f,e,d,c[d]	[8f07]

	Thus each 16-perm needs four x 4-bits = 16 bits to specify for a total of 240 bits,
	rather more than we would like, but I see no good way to reduce that further, to say <= 64 bits.
	In fact the bookkeeping needed even for this index-pattern compression is nontrivial, so just
	specify each 16-perm directly via 16 x 4-bits, with each 4-bit subfield indexing directly into
	a static p[0-f] offsets array. We simply encode each 16-perm as a hex-char string: *NOTE* this
	means we must extract each p-offset in little-endian fashion, e.g. low 4 bits have rightmost p-offset above.
	*/
	#if USE_COMPACT_OBJ_CODE
		tptr = t;
		for(i = 0; i < 15; i++) {
			i64 = perm16[i];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jt = j1 + p_in_hi[i]; jp = j2 + p_in_hi[i];	// p_in_hi[] = p0,p20,...,p30
			RADIX_16_DIT(a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
				c16,s16);	tptr++;
		}
	#else
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_16_DIT(a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_16_DIT(a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p10; jp = j2+p10;	RADIX_16_DIT(a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+pe0; jp = j2+pe0;	RADIX_16_DIT(a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+pd0; jp = j2+pd0;	RADIX_16_DIT(a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+pc0; jp = j2+pc0;	RADIX_16_DIT(a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+pb0; jp = j2+pb0;	RADIX_16_DIT(a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+pa0; jp = j2+pa0;	RADIX_16_DIT(a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p90; jp = j2+p90;	RADIX_16_DIT(a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p80; jp = j2+p80;	RADIX_16_DIT(a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p70; jp = j2+p70;	RADIX_16_DIT(a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p60; jp = j2+p60;	RADIX_16_DIT(a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p50; jp = j2+p50;	RADIX_16_DIT(a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],a[jt+p9],a[jp+p9],a[jt+p8],a[jp+p8],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p1],a[jp+p1],a[jt   ],a[jp   ],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_16_DIT(a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
		jt = j1+p30; jp = j2+p30;	RADIX_16_DIT(a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pb],a[jp+pb],a[jt+pa],a[jp+pa],a[jt+pf],a[jp+pf],a[jt+pe],a[jp+pe],a[jt+pd],a[jp+pd],a[jt+pc],a[jp+pc],a[jt   ],a[jp   ],a[jt+p1],a[jp+p1],a[jt+p3],a[jp+p3],a[jt+p2],a[jp+p2],a[jt+p7],a[jp+p7],a[jt+p6],a[jp+p6],a[jt+p5],a[jp+p5],a[jt+p4],a[jp+p4],
												tptr->re,tptr->im,(tptr+0x0f)->re,(tptr+0x0f)->im,(tptr+0x1e)->re,(tptr+0x1e)->im,(tptr+0x2d)->re,(tptr+0x2d)->im,(tptr+0x3c)->re,(tptr+0x3c)->im,(tptr+0x4b)->re,(tptr+0x4b)->im,(tptr+0x5a)->re,(tptr+0x5a)->im,(tptr+0x69)->re,(tptr+0x69)->im,(tptr+0x78)->re,(tptr+0x78)->im,(tptr+0x87)->re,(tptr+0x87)->im,(tptr+0x96)->re,(tptr+0x96)->im,(tptr+0xa5)->re,(tptr+0xa5)->im,(tptr+0xb4)->re,(tptr+0xb4)->im,(tptr+0xc3)->re,(tptr+0xc3)->im,(tptr+0xd2)->re,(tptr+0xd2)->im,(tptr+0xe1)->re,(tptr+0xe1)->im,
												c16,s16);	tptr++;
	#endif
		/*...and now do 16 radix-15 transforms.
		Using the initial out-index patterning:

			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 1
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 2
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 3
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 4
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 5
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 6
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 7
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 8
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + 9
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + a
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + b
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + c
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + d
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + e
			00,10,20,30,40,50,60,70,80,90,a0,b0,c0,d0,e0 + f

		...the required output permutation is:

			00,0f,1e,2d,3c,4b,5a,69,78,87,96,a5,b4,c3,d2
			e1,e0,ef,0e,1d,2c,3b,4a,59,68,77,86,95,a4,b3
			c2,d1,d0,df,ee,0d,1c,2b,3a,49,58,67,76,85,94
			a3,b2,c1,c0,cf,de,ed,0c,1b,2a,39,48,57,66,75
			84,93,a2,b1,b0,bf,ce,dd,ec,0b,1a,29,38,47,56
			65,74,83,92,a1,a0,af,be,cd,dc,eb,0a,19,28,37
			46,55,64,73,82,91,90,9f,ae,bd,cc,db,ea,09,18
			27,36,45,54,63,72,81,80,8f,9e,ad,bc,cb,da,e9
			08,17,26,35,44,53,62,71,70,7f,8e,9d,ac,bb,ca
			d9,e8,07,16,25,34,43,52,61,60,6f,7e,8d,9c,ab
			ba,c9,d8,e7,06,15,24,33,42,51,50,5f,6e,7d,8c
			9b,aa,b9,c8,d7,e6,05,14,23,32,41,40,4f,5e,6d
			7c,8b,9a,a9,b8,c7,d6,e5,04,13,22,31,30,3f,4e
			5d,6c,7b,8a,99,a8,b7,c6,d5,e4,03,12,21,20,2f
			3e,4d,5c,6b,7a,89,98,a7,b6,c5,d4,e3,02,11,10
			1f,2e,3d,4c,5b,6a,79,88,97,a6,b5,c4,d3,e2,01 ,

		but our initial patterning has unit strides going downward in each column, thus the real o-perm
		results from casting the above from a 16 x 15 matrix to a 16 x 15 [first convert to linear indexing,
		chop into 16-element rows rather than 15-element ones]:

			00,0f,1e,2d,3c,4b,5a,69,78,87,96,a5,b4,c3,d2,e1,
			e0,ef,0e,1d,2c,3b,4a,59,68,77,86,95,a4,b3,c2,d1,
			d0,df,ee,0d,1c,2b,3a,49,58,67,76,85,94,a3,b2,c1,
			c0,cf,de,ed,0c,1b,2a,39,48,57,66,75,84,93,a2,b1,
			b0,bf,ce,dd,ec,0b,1a,29,38,47,56,65,74,83,92,a1,
			a0,af,be,cd,dc,eb,0a,19,28,37,46,55,64,73,82,91,
			90,9f,ae,bd,cc,db,ea,09,18,27,36,45,54,63,72,81,
			80,8f,9e,ad,bc,cb,da,e9,08,17,26,35,44,53,62,71,
			70,7f,8e,9d,ac,bb,ca,d9,e8,07,16,25,34,43,52,61,
			60,6f,7e,8d,9c,ab,ba,c9,d8,e7,06,15,24,33,42,51,
			50,5f,6e,7d,8c,9b,aa,b9,c8,d7,e6,05,14,23,32,41,
			40,4f,5e,6d,7c,8b,9a,a9,b8,c7,d6,e5,04,13,22,31,
			30,3f,4e,5d,6c,7b,8a,99,a8,b7,c6,d5,e4,03,12,21,
			20,2f,3e,4d,5c,6b,7a,89,98,a7,b6,c5,d4,e3,02,11,
			10,1f,2e,3d,4c,5b,6a,79,88,97,a6,b5,c4,d3,e2,01,

		...and transposing to yield a 16 x 15 perm-matrix with a simpler index pattern, to boot:

			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + 0
			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + f
			10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20 + e
			20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30 + d
			30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40 + c
			40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50 + b
			50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60 + a
			60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70 + 9
			70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80 + 8
			80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90 + 7
			90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0 + 6
			a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0 + 5
			b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0 + 4
			c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0 + 3
			d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0 + 2
			e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00 + 1 .

		In order to make the above amenable to loop-based execution, we need to encode the indices
		to the left & right of the + in computable-index fashion. RHS is simple: (16 - i)&0xf for i = 0,..,f.
		For the LHS we need 2 ingredients:
			[1] A formula yielding the lead parts: 0,0,1,2,...,e. Note (i-1) = -1,0,1,2,...,e, so now
				just need apply a bitmask which == 0 if i==0 and all-ones otherwise, thus idx = (i-1) & (-(i > 0)) .
			[2] An efficient decrement (mod 15) scheme to yield the remaining leading parts of each row's elements:
				idx--; idx += (-(idx < 0))&15;
		*/
	#if USE_COMPACT_OBJ_CODE
		tptr = t;
		for(i = 0; i < 16; i++) {
		// [1] here:
			k0 = plo[(16 - i)&0xf];	// p0,f,e,...,2,1
			jt = j1 + k0; jp = j2 + k0;
		// [2] here:
			k0 = (i-1)&(-(i>0)); /* 0,0,1,2...e */	// Now get the resulting p* offsets:
			k1 = k0-1; k1 += (-(k1 < 0))&15;		k0 = p_out_hi[k0];
			k2 = k1-1; k2 += (-(k2 < 0))&15;		k1 = p_out_hi[k1];
			k3 = k2-1; k3 += (-(k3 < 0))&15;		k2 = p_out_hi[k2];
			k4 = k3-1; k4 += (-(k4 < 0))&15;		k3 = p_out_hi[k3];
			k5 = k4-1; k5 += (-(k5 < 0))&15;		k4 = p_out_hi[k4];
			k6 = k5-1; k6 += (-(k6 < 0))&15;		k5 = p_out_hi[k5];
			k7 = k6-1; k7 += (-(k7 < 0))&15;		k6 = p_out_hi[k6];
			k8 = k7-1; k8 += (-(k8 < 0))&15;		k7 = p_out_hi[k7];
			k9 = k8-1; k9 += (-(k9 < 0))&15;		k8 = p_out_hi[k8];
			ka = k9-1; ka += (-(ka < 0))&15;		k9 = p_out_hi[k9];
			kb = ka-1; kb += (-(kb < 0))&15;		ka = p_out_hi[ka];
			kc = kb-1; kc += (-(kc < 0))&15;		kb = p_out_hi[kb];
			kd = kc-1; kd += (-(kd < 0))&15;		kc = p_out_hi[kc];
			ke = kd-1; ke += (-(ke < 0))&15;		kd = p_out_hi[kd];
													ke = p_out_hi[ke];
			RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke]
			);	tptr += 0xf;
		}
	#else
		tptr = t;
		jt = j1   ; jp = j2   ;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10]);	tptr += 0xf;
		jt = j1+pf; jp = j2+pf;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10]);	tptr += 0xf;
		jt = j1+pe; jp = j2+pe;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20]);	tptr += 0xf;
		jt = j1+pd; jp = j2+pd;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30]);	tptr += 0xf;
		jt = j1+pc; jp = j2+pc;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40]);	tptr += 0xf;
		jt = j1+pb; jp = j2+pb;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50]);	tptr += 0xf;
		jt = j1+pa; jp = j2+pa;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60]);	tptr += 0xf;
		jt = j1+p9; jp = j2+p9;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70]);	tptr += 0xf;
		jt = j1+p8; jp = j2+p8;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80]);	tptr += 0xf;
		jt = j1+p7; jp = j2+p7;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90]);	tptr += 0xf;
		jt = j1+p6; jp = j2+p6;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0]);	tptr += 0xf;
		jt = j1+p5; jp = j2+p5;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0]);	tptr += 0xf;
		jt = j1+p4; jp = j2+p4;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0]);	tptr += 0xf;
		jt = j1+p3; jp = j2+p3;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0]);	tptr += 0xf;
		jt = j1+p2; jp = j2+p2;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ],a[jt+pe0],a[jp+pe0]);	tptr += 0xf;
		jt = j1+p1; jp = j2+p1;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
										tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
										a[jt+pe0],a[jp+pe0],a[jt+pd0],a[jp+pd0],a[jt+pc0],a[jp+pc0],a[jt+pb0],a[jp+pb0],a[jt+pa0],a[jp+pa0],a[jt+p90],a[jp+p90],a[jt+p80],a[jp+p80],a[jt+p70],a[jp+p70],a[jt+p60],a[jp+p60],a[jt+p50],a[jp+p50],a[jt+p40],a[jp+p40],a[jt+p30],a[jp+p30],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt    ],a[jp    ]);
	#endif
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy240_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
		struct complex *tptr;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0;
		int poff[RADIX>>2], plo[16], p_in_hi[ODD_RADIX], p_out_hi[ODD_RADIX], po_br[16];
		int k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
		uint64 i64;
		const uint64 dif_perm16[ODD_RADIX] = {
			0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x76450123fecd89abull,0x98badcfe54763201ull,
			0x23106754ab98efdcull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,0x32017645ba89fecdull,
			0xefdc98ba10325476ull,0x45672310cdefab98ull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull};
		const uint64 dit_perm16[ODD_RADIX] = {
			0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x98abefcd10236745ull,
			0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,
			0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x89bafedc01327654ull};
		int j,j1,j2,jt,jp,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */

	#ifdef USE_SSE2

	  #ifdef USE_AVX
		const int l2_sz_vd = 5;
	  #else
		const int l2_sz_vd = 4;
	  #endif
		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7, *add8, *add9, *adda, *addb, *addc, *addd, *adde, *addf;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2;	// utility ptrs
		int *itmp;			// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2, *cc0, *ss0, *sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,
			*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,*r0f,
			*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,*r1f,
			*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,*r2f,
			*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,*r3f,
			*r40,*r41,*r42,*r43,*r44,*r45,*r46,*r47,*r48,*r49,*r4a,*r4b,*r4c,*r4d,*r4e,*r4f,
			*r50,*r51,*r52,*r53,*r54,*r55,*r56,*r57,*r58,*r59,*r5a,*r5b,*r5c,*r5d,*r5e,*r5f,
			*r60,*r61,*r62,*r63,*r64,*r65,*r66,*r67,*r68,*r69,*r6a,*r6b,*r6c,*r6d,*r6e,*r6f,
			*r70,*r71,*r72,*r73,*r74,*r75,*r76,*r77,*r78,*r79,*r7a,*r7b,*r7c,*r7d,*r7e,*r7f,
			*r80,*r81,*r82,*r83,*r84,*r85,*r86,*r87,*r88,*r89,*r8a,*r8b,*r8c,*r8d,*r8e,*r8f,
			*r90,*r91,*r92,*r93,*r94,*r95,*r96,*r97,*r98,*r99,*r9a,*r9b,*r9c,*r9d,*r9e,*r9f,
			*ra0,*ra1,*ra2,*ra3,*ra4,*ra5,*ra6,*ra7,*ra8,*ra9,*raa,*rab,*rac,*rad,*rae,*raf,
			*rb0,*rb1,*rb2,*rb3,*rb4,*rb5,*rb6,*rb7,*rb8,*rb9,*rba,*rbb,*rbc,*rbd,*rbe,*rbf,
			*rc0,*rc1,*rc2,*rc3,*rc4,*rc5,*rc6,*rc7,*rc8,*rc9,*rca,*rcb,*rcc,*rcd,*rce,*rcf,
			*rd0,*rd1,*rd2,*rd3,*rd4,*rd5,*rd6,*rd7,*rd8,*rd9,*rda,*rdb,*rdc,*rdd,*rde,*rdf,
			*re0,*re1,*re2,*re3,*re4,*re5,*re6,*re7,*re8,*re9,*rea,*reb,*rec,*red,*ree,*ref,
			*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,*s1p0f,
			*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*s1p18,*s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,*s1p1f,
			*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,*s1p2f,
			*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*s1p38,*s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,*s1p3f,
			*s1p40,*s1p41,*s1p42,*s1p43,*s1p44,*s1p45,*s1p46,*s1p47,*s1p48,*s1p49,*s1p4a,*s1p4b,*s1p4c,*s1p4d,*s1p4e,*s1p4f,
			*s1p50,*s1p51,*s1p52,*s1p53,*s1p54,*s1p55,*s1p56,*s1p57,*s1p58,*s1p59,*s1p5a,*s1p5b,*s1p5c,*s1p5d,*s1p5e,*s1p5f,
			*s1p60,*s1p61,*s1p62,*s1p63,*s1p64,*s1p65,*s1p66,*s1p67,*s1p68,*s1p69,*s1p6a,*s1p6b,*s1p6c,*s1p6d,*s1p6e,*s1p6f,
			*s1p70,*s1p71,*s1p72,*s1p73,*s1p74,*s1p75,*s1p76,*s1p77,*s1p78,*s1p79,*s1p7a,*s1p7b,*s1p7c,*s1p7d,*s1p7e,*s1p7f,
			*s1p80,*s1p81,*s1p82,*s1p83,*s1p84,*s1p85,*s1p86,*s1p87,*s1p88,*s1p89,*s1p8a,*s1p8b,*s1p8c,*s1p8d,*s1p8e,*s1p8f,
			*s1p90,*s1p91,*s1p92,*s1p93,*s1p94,*s1p95,*s1p96,*s1p97,*s1p98,*s1p99,*s1p9a,*s1p9b,*s1p9c,*s1p9d,*s1p9e,*s1p9f,
			*s1pa0,*s1pa1,*s1pa2,*s1pa3,*s1pa4,*s1pa5,*s1pa6,*s1pa7,*s1pa8,*s1pa9,*s1paa,*s1pab,*s1pac,*s1pad,*s1pae,*s1paf,
			*s1pb0,*s1pb1,*s1pb2,*s1pb3,*s1pb4,*s1pb5,*s1pb6,*s1pb7,*s1pb8,*s1pb9,*s1pba,*s1pbb,*s1pbc,*s1pbd,*s1pbe,*s1pbf,
			*s1pc0,*s1pc1,*s1pc2,*s1pc3,*s1pc4,*s1pc5,*s1pc6,*s1pc7,*s1pc8,*s1pc9,*s1pca,*s1pcb,*s1pcc,*s1pcd,*s1pce,*s1pcf,
			*s1pd0,*s1pd1,*s1pd2,*s1pd3,*s1pd4,*s1pd5,*s1pd6,*s1pd7,*s1pd8,*s1pd9,*s1pda,*s1pdb,*s1pdc,*s1pdd,*s1pde,*s1pdf,
			*s1pe0,*s1pe1,*s1pe2,*s1pe3,*s1pe4,*s1pe5,*s1pe6,*s1pe7,*s1pe8,*s1pe9,*s1pea,*s1peb,*s1pec,*s1ped,*s1pee,*s1pef,
			*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
			*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
			*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
						cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
						cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
						ss3 =  0.95105651629515357210,	/*  sin(u) */
						sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
						sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy_r = thread_arg->cy_r,*cy_i = thread_arg->cy_i, temp,frac;
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		struct complex t[RADIX];
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;
	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;
	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3; plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7; plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb; plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );

		poff[     0] =   0; poff[     1] =     p4; poff[     2] =     p8; poff[     3] =     pc;
		poff[0x04+0] = p10; poff[0x04+1] = p10+p4; poff[0x04+2] = p10+p8; poff[0x04+3] = p10+pc;
		poff[0x08+0] = p20; poff[0x08+1] = p20+p4; poff[0x08+2] = p20+p8; poff[0x08+3] = p20+pc;
		poff[0x0c+0] = p30; poff[0x0c+1] = p30+p4; poff[0x0c+2] = p30+p8; poff[0x0c+3] = p30+pc;
		poff[0x10+0] = p40; poff[0x10+1] = p40+p4; poff[0x10+2] = p40+p8; poff[0x10+3] = p40+pc;
		poff[0x14+0] = p50; poff[0x14+1] = p50+p4; poff[0x14+2] = p50+p8; poff[0x14+3] = p50+pc;
		poff[0x18+0] = p60; poff[0x18+1] = p60+p4; poff[0x18+2] = p60+p8; poff[0x18+3] = p60+pc;
		poff[0x1c+0] = p70; poff[0x1c+1] = p70+p4; poff[0x1c+2] = p70+p8; poff[0x1c+3] = p70+pc;
		poff[0x20+0] = p80; poff[0x20+1] = p80+p4; poff[0x20+2] = p80+p8; poff[0x20+3] = p80+pc;
		poff[0x24+0] = p90; poff[0x24+1] = p90+p4; poff[0x24+2] = p90+p8; poff[0x24+3] = p90+pc;
		poff[0x28+0] = pa0; poff[0x28+1] = pa0+p4; poff[0x28+2] = pa0+p8; poff[0x28+3] = pa0+pc;
		poff[0x2c+0] = pb0; poff[0x2c+1] = pb0+p4; poff[0x2c+2] = pb0+p8; poff[0x2c+3] = pb0+pc;
		poff[0x30+0] = pc0; poff[0x30+1] = pc0+p4; poff[0x30+2] = pc0+p8; poff[0x30+3] = pc0+pc;
		poff[0x34+0] = pd0; poff[0x34+1] = pd0+p4; poff[0x34+2] = pd0+p8; poff[0x34+3] = pd0+pc;
		poff[0x38+0] = pe0; poff[0x38+1] = pe0+p4; poff[0x38+2] = pe0+p8; poff[0x38+3] = pe0+pc;

		p_in_hi [0x0] =   0; p_in_hi [0x1] = p20; p_in_hi [0x2] = p10; p_in_hi [0x3] = pe0; p_in_hi [0x4] = pd0; p_in_hi [0x5] = pc0; p_in_hi [0x6] = pb0; p_in_hi [0x7] = pa0; p_in_hi [0x8] = p90; p_in_hi [0x9] = p80; p_in_hi [0xa] = p70; p_in_hi [0xb] = p60; p_in_hi [0xc] = p50; p_in_hi [0xd] = p40; p_in_hi [0xe] = p30;
		p_out_hi[0x0] =   0; p_out_hi[0x1] = p10; p_out_hi[0x2] = p20; p_out_hi[0x3] = p30; p_out_hi[0x4] = p40; p_out_hi[0x5] = p50; p_out_hi[0x6] = p60; p_out_hi[0x7] = p70; p_out_hi[0x8] = p80; p_out_hi[0x9] = p90; p_out_hi[0xa] = pa0; p_out_hi[0xb] = pb0; p_out_hi[0xc] = pc0; p_out_hi[0xd] = pd0; p_out_hi[0xe] = pe0;
	#if !USE_SCALAR_DFT_MACRO || !defined(USE_SSE2)	// SIMD build uses these for radix-256 DFTs; Scalar-double build uses these for carry-macro loops
		// Needed for the length-16 loops which manage the DFT16 and carry-macro calls in the compact-object-code build:
		po_br[0x0] = 0; po_br[0x1] = p8; po_br[0x2] = p4; po_br[0x3] = pc; po_br[0x4] = p2; po_br[0x5] = pa; po_br[0x6] = p6; po_br[0x7] = pe; po_br[0x8] = p1; po_br[0x9] = p9; po_br[0xa] = p5; po_br[0xb] = pd; po_br[0xc] = p3; po_br[0xd] = pb; po_br[0xe] = p7; po_br[0xf] = pf;
	#endif

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		r00 = thread_arg->r00;
		tmp = r00 + 0;			tm2 = tmp + 0x100;
		r00 = tmp + 0x00;		r80 = tm2 + 0x00;
		r01 = tmp + 0x02;		r81 = tm2 + 0x02;
		r02 = tmp + 0x04;		r82 = tm2 + 0x04;
		r03 = tmp + 0x06;		r83 = tm2 + 0x06;
		r04 = tmp + 0x08;		r84 = tm2 + 0x08;
		r05 = tmp + 0x0a;		r85 = tm2 + 0x0a;
		r06 = tmp + 0x0c;		r86 = tm2 + 0x0c;
		r07 = tmp + 0x0e;		r87 = tm2 + 0x0e;
		r08 = tmp + 0x10;		r88 = tm2 + 0x10;
		r09 = tmp + 0x12;		r89 = tm2 + 0x12;
		r0a = tmp + 0x14;		r8a = tm2 + 0x14;
		r0b = tmp + 0x16;		r8b = tm2 + 0x16;
		r0c = tmp + 0x18;		r8c = tm2 + 0x18;
		r0d = tmp + 0x1a;		r8d = tm2 + 0x1a;
		r0e = tmp + 0x1c;		r8e = tm2 + 0x1c;
		r0f = tmp + 0x1e;		r8f = tm2 + 0x1e;
		r10 = tmp + 0x20;		r90 = tm2 + 0x20;
		r11 = tmp + 0x22;		r91 = tm2 + 0x22;
		r12 = tmp + 0x24;		r92 = tm2 + 0x24;
		r13 = tmp + 0x26;		r93 = tm2 + 0x26;
		r14 = tmp + 0x28;		r94 = tm2 + 0x28;
		r15 = tmp + 0x2a;		r95 = tm2 + 0x2a;
		r16 = tmp + 0x2c;		r96 = tm2 + 0x2c;
		r17 = tmp + 0x2e;		r97 = tm2 + 0x2e;
		r18 = tmp + 0x30;		r98 = tm2 + 0x30;
		r19 = tmp + 0x32;		r99 = tm2 + 0x32;
		r1a = tmp + 0x34;		r9a = tm2 + 0x34;
		r1b = tmp + 0x36;		r9b = tm2 + 0x36;
		r1c = tmp + 0x38;		r9c = tm2 + 0x38;
		r1d = tmp + 0x3a;		r9d = tm2 + 0x3a;
		r1e = tmp + 0x3c;		r9e = tm2 + 0x3c;
		r1f = tmp + 0x3e;		r9f = tm2 + 0x3e;
		r20 = tmp + 0x40;		ra0 = tm2 + 0x40;
		r21 = tmp + 0x42;		ra1 = tm2 + 0x42;
		r22 = tmp + 0x44;		ra2 = tm2 + 0x44;
		r23 = tmp + 0x46;		ra3 = tm2 + 0x46;
		r24 = tmp + 0x48;		ra4 = tm2 + 0x48;
		r25 = tmp + 0x4a;		ra5 = tm2 + 0x4a;
		r26 = tmp + 0x4c;		ra6 = tm2 + 0x4c;
		r27 = tmp + 0x4e;		ra7 = tm2 + 0x4e;
		r28 = tmp + 0x50;		ra8 = tm2 + 0x50;
		r29 = tmp + 0x52;		ra9 = tm2 + 0x52;
		r2a = tmp + 0x54;		raa = tm2 + 0x54;
		r2b = tmp + 0x56;		rab = tm2 + 0x56;
		r2c = tmp + 0x58;		rac = tm2 + 0x58;
		r2d = tmp + 0x5a;		rad = tm2 + 0x5a;
		r2e = tmp + 0x5c;		rae = tm2 + 0x5c;
		r2f = tmp + 0x5e;		raf = tm2 + 0x5e;
		r30 = tmp + 0x60;		rb0 = tm2 + 0x60;
		r31 = tmp + 0x62;		rb1 = tm2 + 0x62;
		r32 = tmp + 0x64;		rb2 = tm2 + 0x64;
		r33 = tmp + 0x66;		rb3 = tm2 + 0x66;
		r34 = tmp + 0x68;		rb4 = tm2 + 0x68;
		r35 = tmp + 0x6a;		rb5 = tm2 + 0x6a;
		r36 = tmp + 0x6c;		rb6 = tm2 + 0x6c;
		r37 = tmp + 0x6e;		rb7 = tm2 + 0x6e;
		r38 = tmp + 0x70;		rb8 = tm2 + 0x70;
		r39 = tmp + 0x72;		rb9 = tm2 + 0x72;
		r3a = tmp + 0x74;		rba = tm2 + 0x74;
		r3b = tmp + 0x76;		rbb = tm2 + 0x76;
		r3c = tmp + 0x78;		rbc = tm2 + 0x78;
		r3d = tmp + 0x7a;		rbd = tm2 + 0x7a;
		r3e = tmp + 0x7c;		rbe = tm2 + 0x7c;
		r3f = tmp + 0x7e;		rbf = tm2 + 0x7e;
		r40 = tmp + 0x80;		rc0 = tm2 + 0x80;
		r41 = tmp + 0x82;		rc1 = tm2 + 0x82;
		r42 = tmp + 0x84;		rc2 = tm2 + 0x84;
		r43 = tmp + 0x86;		rc3 = tm2 + 0x86;
		r44 = tmp + 0x88;		rc4 = tm2 + 0x88;
		r45 = tmp + 0x8a;		rc5 = tm2 + 0x8a;
		r46 = tmp + 0x8c;		rc6 = tm2 + 0x8c;
		r47 = tmp + 0x8e;		rc7 = tm2 + 0x8e;
		r48 = tmp + 0x90;		rc8 = tm2 + 0x90;
		r49 = tmp + 0x92;		rc9 = tm2 + 0x92;
		r4a = tmp + 0x94;		rca = tm2 + 0x94;
		r4b = tmp + 0x96;		rcb = tm2 + 0x96;
		r4c = tmp + 0x98;		rcc = tm2 + 0x98;
		r4d = tmp + 0x9a;		rcd = tm2 + 0x9a;
		r4e = tmp + 0x9c;		rce = tm2 + 0x9c;
		r4f = tmp + 0x9e;		rcf = tm2 + 0x9e;
		r50 = tmp + 0xa0;		rd0 = tm2 + 0xa0;
		r51 = tmp + 0xa2;		rd1 = tm2 + 0xa2;
		r52 = tmp + 0xa4;		rd2 = tm2 + 0xa4;
		r53 = tmp + 0xa6;		rd3 = tm2 + 0xa6;
		r54 = tmp + 0xa8;		rd4 = tm2 + 0xa8;
		r55 = tmp + 0xaa;		rd5 = tm2 + 0xaa;
		r56 = tmp + 0xac;		rd6 = tm2 + 0xac;
		r57 = tmp + 0xae;		rd7 = tm2 + 0xae;
		r58 = tmp + 0xb0;		rd8 = tm2 + 0xb0;
		r59 = tmp + 0xb2;		rd9 = tm2 + 0xb2;
		r5a = tmp + 0xb4;		rda = tm2 + 0xb4;
		r5b = tmp + 0xb6;		rdb = tm2 + 0xb6;
		r5c = tmp + 0xb8;		rdc = tm2 + 0xb8;
		r5d = tmp + 0xba;		rdd = tm2 + 0xba;
		r5e = tmp + 0xbc;		rde = tm2 + 0xbc;
		r5f = tmp + 0xbe;		rdf = tm2 + 0xbe;
		r60 = tmp + 0xc0;		re0 = tm2 + 0xc0;
		r61 = tmp + 0xc2;		re1 = tm2 + 0xc2;
		r62 = tmp + 0xc4;		re2 = tm2 + 0xc4;
		r63 = tmp + 0xc6;		re3 = tm2 + 0xc6;
		r64 = tmp + 0xc8;		re4 = tm2 + 0xc8;
		r65 = tmp + 0xca;		re5 = tm2 + 0xca;
		r66 = tmp + 0xcc;		re6 = tm2 + 0xcc;
		r67 = tmp + 0xce;		re7 = tm2 + 0xce;
		r68 = tmp + 0xd0;		re8 = tm2 + 0xd0;
		r69 = tmp + 0xd2;		re9 = tm2 + 0xd2;
		r6a = tmp + 0xd4;		rea = tm2 + 0xd4;
		r6b = tmp + 0xd6;		reb = tm2 + 0xd6;
		r6c = tmp + 0xd8;		rec = tm2 + 0xd8;
		r6d = tmp + 0xda;		red = tm2 + 0xda;
		r6e = tmp + 0xdc;		ree = tm2 + 0xdc;
		r6f = tmp + 0xde;		ref = tm2 + 0xde;
		r70 = tmp + 0xe0;
		r71 = tmp + 0xe2;
		r72 = tmp + 0xe4;
		r73 = tmp + 0xe6;
		r74 = tmp + 0xe8;
		r75 = tmp + 0xea;
		r76 = tmp + 0xec;
		r77 = tmp + 0xee;
		r78 = tmp + 0xf0;
		r79 = tmp + 0xf2;
		r7a = tmp + 0xf4;
		r7b = tmp + 0xf6;
		r7c = tmp + 0xf8;
		r7d = tmp + 0xfa;
		r7e = tmp + 0xfc;
		r7f = tmp + 0xfe;
		tmp += 0x1e0;			tm2 += 0x1e0;
		s1p00 = tmp + 0x00;		s1p80 = tm2 + 0x00;
		s1p01 = tmp + 0x02;		s1p81 = tm2 + 0x02;
		s1p02 = tmp + 0x04;		s1p82 = tm2 + 0x04;
		s1p03 = tmp + 0x06;		s1p83 = tm2 + 0x06;
		s1p04 = tmp + 0x08;		s1p84 = tm2 + 0x08;
		s1p05 = tmp + 0x0a;		s1p85 = tm2 + 0x0a;
		s1p06 = tmp + 0x0c;		s1p86 = tm2 + 0x0c;
		s1p07 = tmp + 0x0e;		s1p87 = tm2 + 0x0e;
		s1p08 = tmp + 0x10;		s1p88 = tm2 + 0x10;
		s1p09 = tmp + 0x12;		s1p89 = tm2 + 0x12;
		s1p0a = tmp + 0x14;		s1p8a = tm2 + 0x14;
		s1p0b = tmp + 0x16;		s1p8b = tm2 + 0x16;
		s1p0c = tmp + 0x18;		s1p8c = tm2 + 0x18;
		s1p0d = tmp + 0x1a;		s1p8d = tm2 + 0x1a;
		s1p0e = tmp + 0x1c;		s1p8e = tm2 + 0x1c;
		s1p0f = tmp + 0x1e;		s1p8f = tm2 + 0x1e;
		s1p10 = tmp + 0x20;		s1p90 = tm2 + 0x20;
		s1p11 = tmp + 0x22;		s1p91 = tm2 + 0x22;
		s1p12 = tmp + 0x24;		s1p92 = tm2 + 0x24;
		s1p13 = tmp + 0x26;		s1p93 = tm2 + 0x26;
		s1p14 = tmp + 0x28;		s1p94 = tm2 + 0x28;
		s1p15 = tmp + 0x2a;		s1p95 = tm2 + 0x2a;
		s1p16 = tmp + 0x2c;		s1p96 = tm2 + 0x2c;
		s1p17 = tmp + 0x2e;		s1p97 = tm2 + 0x2e;
		s1p18 = tmp + 0x30;		s1p98 = tm2 + 0x30;
		s1p19 = tmp + 0x32;		s1p99 = tm2 + 0x32;
		s1p1a = tmp + 0x34;		s1p9a = tm2 + 0x34;
		s1p1b = tmp + 0x36;		s1p9b = tm2 + 0x36;
		s1p1c = tmp + 0x38;		s1p9c = tm2 + 0x38;
		s1p1d = tmp + 0x3a;		s1p9d = tm2 + 0x3a;
		s1p1e = tmp + 0x3c;		s1p9e = tm2 + 0x3c;
		s1p1f = tmp + 0x3e;		s1p9f = tm2 + 0x3e;
		s1p20 = tmp + 0x40;		s1pa0 = tm2 + 0x40;
		s1p21 = tmp + 0x42;		s1pa1 = tm2 + 0x42;
		s1p22 = tmp + 0x44;		s1pa2 = tm2 + 0x44;
		s1p23 = tmp + 0x46;		s1pa3 = tm2 + 0x46;
		s1p24 = tmp + 0x48;		s1pa4 = tm2 + 0x48;
		s1p25 = tmp + 0x4a;		s1pa5 = tm2 + 0x4a;
		s1p26 = tmp + 0x4c;		s1pa6 = tm2 + 0x4c;
		s1p27 = tmp + 0x4e;		s1pa7 = tm2 + 0x4e;
		s1p28 = tmp + 0x50;		s1pa8 = tm2 + 0x50;
		s1p29 = tmp + 0x52;		s1pa9 = tm2 + 0x52;
		s1p2a = tmp + 0x54;		s1paa = tm2 + 0x54;
		s1p2b = tmp + 0x56;		s1pab = tm2 + 0x56;
		s1p2c = tmp + 0x58;		s1pac = tm2 + 0x58;
		s1p2d = tmp + 0x5a;		s1pad = tm2 + 0x5a;
		s1p2e = tmp + 0x5c;		s1pae = tm2 + 0x5c;
		s1p2f = tmp + 0x5e;		s1paf = tm2 + 0x5e;
		s1p30 = tmp + 0x60;		s1pb0 = tm2 + 0x60;
		s1p31 = tmp + 0x62;		s1pb1 = tm2 + 0x62;
		s1p32 = tmp + 0x64;		s1pb2 = tm2 + 0x64;
		s1p33 = tmp + 0x66;		s1pb3 = tm2 + 0x66;
		s1p34 = tmp + 0x68;		s1pb4 = tm2 + 0x68;
		s1p35 = tmp + 0x6a;		s1pb5 = tm2 + 0x6a;
		s1p36 = tmp + 0x6c;		s1pb6 = tm2 + 0x6c;
		s1p37 = tmp + 0x6e;		s1pb7 = tm2 + 0x6e;
		s1p38 = tmp + 0x70;		s1pb8 = tm2 + 0x70;
		s1p39 = tmp + 0x72;		s1pb9 = tm2 + 0x72;
		s1p3a = tmp + 0x74;		s1pba = tm2 + 0x74;
		s1p3b = tmp + 0x76;		s1pbb = tm2 + 0x76;
		s1p3c = tmp + 0x78;		s1pbc = tm2 + 0x78;
		s1p3d = tmp + 0x7a;		s1pbd = tm2 + 0x7a;
		s1p3e = tmp + 0x7c;		s1pbe = tm2 + 0x7c;
		s1p3f = tmp + 0x7e;		s1pbf = tm2 + 0x7e;
		s1p40 = tmp + 0x80;		s1pc0 = tm2 + 0x80;
		s1p41 = tmp + 0x82;		s1pc1 = tm2 + 0x82;
		s1p42 = tmp + 0x84;		s1pc2 = tm2 + 0x84;
		s1p43 = tmp + 0x86;		s1pc3 = tm2 + 0x86;
		s1p44 = tmp + 0x88;		s1pc4 = tm2 + 0x88;
		s1p45 = tmp + 0x8a;		s1pc5 = tm2 + 0x8a;
		s1p46 = tmp + 0x8c;		s1pc6 = tm2 + 0x8c;
		s1p47 = tmp + 0x8e;		s1pc7 = tm2 + 0x8e;
		s1p48 = tmp + 0x90;		s1pc8 = tm2 + 0x90;
		s1p49 = tmp + 0x92;		s1pc9 = tm2 + 0x92;
		s1p4a = tmp + 0x94;		s1pca = tm2 + 0x94;
		s1p4b = tmp + 0x96;		s1pcb = tm2 + 0x96;
		s1p4c = tmp + 0x98;		s1pcc = tm2 + 0x98;
		s1p4d = tmp + 0x9a;		s1pcd = tm2 + 0x9a;
		s1p4e = tmp + 0x9c;		s1pce = tm2 + 0x9c;
		s1p4f = tmp + 0x9e;		s1pcf = tm2 + 0x9e;
		s1p50 = tmp + 0xa0;		s1pd0 = tm2 + 0xa0;
		s1p51 = tmp + 0xa2;		s1pd1 = tm2 + 0xa2;
		s1p52 = tmp + 0xa4;		s1pd2 = tm2 + 0xa4;
		s1p53 = tmp + 0xa6;		s1pd3 = tm2 + 0xa6;
		s1p54 = tmp + 0xa8;		s1pd4 = tm2 + 0xa8;
		s1p55 = tmp + 0xaa;		s1pd5 = tm2 + 0xaa;
		s1p56 = tmp + 0xac;		s1pd6 = tm2 + 0xac;
		s1p57 = tmp + 0xae;		s1pd7 = tm2 + 0xae;
		s1p58 = tmp + 0xb0;		s1pd8 = tm2 + 0xb0;
		s1p59 = tmp + 0xb2;		s1pd9 = tm2 + 0xb2;
		s1p5a = tmp + 0xb4;		s1pda = tm2 + 0xb4;
		s1p5b = tmp + 0xb6;		s1pdb = tm2 + 0xb6;
		s1p5c = tmp + 0xb8;		s1pdc = tm2 + 0xb8;
		s1p5d = tmp + 0xba;		s1pdd = tm2 + 0xba;
		s1p5e = tmp + 0xbc;		s1pde = tm2 + 0xbc;
		s1p5f = tmp + 0xbe;		s1pdf = tm2 + 0xbe;
		s1p60 = tmp + 0xc0;		s1pe0 = tm2 + 0xc0;
		s1p61 = tmp + 0xc2;		s1pe1 = tm2 + 0xc2;
		s1p62 = tmp + 0xc4;		s1pe2 = tm2 + 0xc4;
		s1p63 = tmp + 0xc6;		s1pe3 = tm2 + 0xc6;
		s1p64 = tmp + 0xc8;		s1pe4 = tm2 + 0xc8;
		s1p65 = tmp + 0xca;		s1pe5 = tm2 + 0xca;
		s1p66 = tmp + 0xcc;		s1pe6 = tm2 + 0xcc;
		s1p67 = tmp + 0xce;		s1pe7 = tm2 + 0xce;
		s1p68 = tmp + 0xd0;		s1pe8 = tm2 + 0xd0;
		s1p69 = tmp + 0xd2;		s1pe9 = tm2 + 0xd2;
		s1p6a = tmp + 0xd4;		s1pea = tm2 + 0xd4;
		s1p6b = tmp + 0xd6;		s1peb = tm2 + 0xd6;
		s1p6c = tmp + 0xd8;		s1pec = tm2 + 0xd8;
		s1p6d = tmp + 0xda;		s1ped = tm2 + 0xda;
		s1p6e = tmp + 0xdc;		s1pee = tm2 + 0xdc;
		s1p6f = tmp + 0xde;		s1pef = tm2 + 0xde;
		s1p70 = tmp + 0xe0;
		s1p71 = tmp + 0xe2;
		s1p72 = tmp + 0xe4;
		s1p73 = tmp + 0xe6;
		s1p74 = tmp + 0xe8;
		s1p75 = tmp + 0xea;
		s1p76 = tmp + 0xec;
		s1p77 = tmp + 0xee;
		s1p78 = tmp + 0xf0;
		s1p79 = tmp + 0xf2;
		s1p7a = tmp + 0xf4;
		s1p7b = tmp + 0xf6;
		s1p7c = tmp + 0xf8;
		s1p7d = tmp + 0xfa;
		s1p7e = tmp + 0xfc;
		s1p7f = tmp + 0xfe;
		tmp += 0x1e0;	// sc_ptr += 0x3c0
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;
		tmp += 0x1e;	// += 0x3c => sc_ptr + 0x3fc
		ASSERT(HERE, (tmp->d0 == tmp->d1) && (tmp->d0 == ISRT2), "thread-local memcheck failed!");
		// DFT-roots:
		isrt2     = tmp + 0x00;
		cc0       = tmp + 0x01;
		ss0       = tmp + 0x02;
		sse2_c3m1 = tmp + 0x03;
		sse2_s    = tmp + 0x04;
		sse2_cn1  = tmp + 0x05;
		sse2_cn2  = tmp + 0x06;
		sse2_ss3  = tmp + 0x07;
		sse2_sn1  = tmp + 0x08;
		sse2_sn2  = tmp + 0x09;
		tmp += 0x0a;	// += 0xa => sc_ptr + 0x406
	  #ifdef USE_AVX
		cy_r = tmp;	cy_i = tmp+0x3c;	tmp += 2*0x3c;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x78 + 2 => sc_ptr += 0x480
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #else
		cy_r = tmp;	cy_i = tmp+0x78;	tmp += 2*0x78;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0xf0 + 2 => sc_ptr += 0x4f8
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif

		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix240_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;

		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
		wts_idx_incr = *(int *)thread_arg->half_arr;
		base      = (double *)thread_arg->r00;
		baseinv   = base + 2;
		wt_arr    = base + 4;
		wtinv_arr = wt_arr    + ODD_RADIX;
		bs_arr    = wtinv_arr + ODD_RADIX;
		bsinv_arr = bs_arr    + ODD_RADIX;

	#endif	// USE_SSE2 ?

		// Can't simply use thread-associated values of these *cycle index arrays here, since
		// thread values must be ***read-only*** so as to retain the proper first-init values
		// on each entry to this thread-task. Instead use the bjmodn data storage block - which
		// is otherwise unused in Fermat-Mod mode - for local storage of these cycle tables:
		int *icycle = bjmodn,ic;
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int *jcycle = icycle + ODD_RADIX,jc;
	  #ifdef USE_AVX
		int *kcycle = jcycle + ODD_RADIX;	// NB: kc already declared as part of k0-f set above
		int *lcycle = kcycle + ODD_RADIX,lc;
	  #endif
	#endif
		for(j = 0; j < ODD_RADIX; j++) {
			icycle[j] = thread_arg->icycle[j];
		#ifdef USE_SSE2
			jcycle[j] = thread_arg->jcycle[j];
		  #ifdef USE_AVX
			kcycle[j] = thread_arg->kcycle[j];
			lcycle[j] = thread_arg->lcycle[j];
		  #endif
		#endif
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */
			uint32 bjmodnini = thread_arg->bjmodnini;
			bjmodn[0] = thread_arg->bjmodn0;
			for(l = 1; l < RADIX; l++) {	// must use e.g. l for loop idx here as i is used for dwt indexing
				MOD_ADD32(bjmodn[l-1], bjmodnini, n, bjmodn[l]);
			}

			/* init carries	*/
			addr = thread_arg->cy_r;
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				tmp->d0 = *(addr+l  );
				tmp->d1 = *(addr+l+1);
				tmp->d2 = *(addr+l+2);
				tmp->d3 = *(addr+l+3);
			}
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				tmp->d0 = *(addr+l  );
				tmp->d1 = *(addr+l+1);
			}
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = *(addr+l);
			}
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				tmp->d0 = *(addr+l  );		tm2->d0 = *(addi+l  );
				tmp->d1 = *(addr+l+1);		tm2->d1 = *(addi+l+1);
				tmp->d2 = *(addr+l+2);		tm2->d2 = *(addi+l+2);
				tmp->d3 = *(addr+l+3);		tm2->d3 = *(addi+l+3);
			}
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				tmp->d0 = *(addr+l  );		tmp->d1 = *(addi+l  );
			}
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = *(addr+l);		cy_i[l] = *(addi+l);
			}
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix240_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			addr = thread_arg->cy_r;
		#ifdef USE_AVX
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				*(addr+l  ) = tmp->d0;
				*(addr+l+1) = tmp->d1;
				*(addr+l+2) = tmp->d2;
				*(addr+l+3) = tmp->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				*(addr+l  ) = tmp->d0;
				*(addr+l+1) = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				*(addr+l) = cy_r[l];
			}
		#endif
		}
		else
		{
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#ifdef USE_AVX
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tm2->d0;
				*(addr+l+1) = tmp->d1;		*(addi+l+1) = tm2->d1;
				*(addr+l+2) = tmp->d2;		*(addi+l+2) = tm2->d2;
				*(addr+l+3) = tmp->d3;		*(addi+l+3) = tm2->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				*(addr+l) = cy_r[l];		*(addi+l) = cy_i[l];
			}
		#endif
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}
		return 0x0;
	}
#endif

#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

	#undef OFF1
	#undef OFF2
	#undef OFF3
	#undef OFF4

#endif

#undef RADIX
#undef ODD_RADIX
