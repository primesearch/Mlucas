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
#include "radix32.h"

#define RADIX 224	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 7	// ODD_RADIX = [radix >> trailz(radix)]

#define EPS 1e-10

#ifndef PFETCH_DIST
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

// See the radix28 version of this routine for details about the
// small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.

#ifdef USE_SSE2

	#include "sse2_macro_gcc64.h"

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]),
  // For Fermat-mod in AVX mode we need RADIX*4 = 224*4 [if HIACC] or 12 [if not] slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(132,224) = 224 if AVX, 40 if SSE2
  // to (half_arr_offset224 + RADIX) to get required value of radix224_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/4 = 0x38 fewer carry slots than AVX:
	const int half_arr_offset224 = 0x3d0;	// + RADIX = + 0xe0 = 0x4b0; Used for thread local-storage-integrity checking
	const int radix224_creals_in_local_store = 0x850;	// (half_arr_offset + RADIX*4) = 0x4b0 + 0x3a0 = 0x850 and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset224 = 0x408;	// + RADIX = + 0xe0 = 0x4e8; Used for thread local-storage-integrity checking
	const int radix224_creals_in_local_store = 0x888;	// (half_arr_offset + RADIX*4) = 0x4e8 + 0x3a0 = 0x888 and round up to nearest multiple of 4
  #else
	const int half_arr_offset224 = 0x478;	// + RADIX = + 0xe0 = 0x558; Used for thread local-storage-integrity checking
	const int radix224_creals_in_local_store = 0x5a0;	// (half_arr_offset + RADIX) + 40 = 0x578 + 0x28 and round up to nearest multiple of 8
  #endif

  #ifdef USE_AVX
	const uint64 radix224_avx_negadwt_consts[RADIX] = {	// 8 entries per line ==> RADIX/8 lines:
		0x3FF0000000000000ull,0x3FEFFFCC6FC2A580ull,0x3FEFFF31BFB0C29Bull,0x3FEFFE2FF1BCDB0Cull,0x3FEFFCC70925C367ull,0x3FEFFAF70A7696A3ull,0x3FEFF8BFFB86A774ull,0x3FEFF621E3796D7Eull,
		0x3FEFF31CCABE6E4Cull,0x3FEFEFB0BB112225ull,0x3FEFEBDDBF78D4AEull,0x3FEFE7A3E448815Full,0x3FEFE303371EABCBull,0x3FEFDDFBC6E533BCull,0x3FEFD88DA3D12526ull,0x3FEFD2B8DF6283E5ull,
		0x3FEFCC7D8C64135Full,0x3FEFC5DBBEEB19EFull,0x3FEFBED38C57202Eull,0x3FEFB7650B51AC12ull,0x3FEFAF9053CDF7E6ull,0x3FEFA7557F08A517ull,0x3FEF9EB4A7876AE4ull,0x3FEF95ADE918C0E4ull,
		0x3FEF8C4160D38565ull,0x3FEF826F2D169FB1ull,0x3FEF78376D889E2Eull,0x3FEF6D9A4317505Aull,0x3FEF6297CFF75CB0ull,0x3FEF573037A3D269ull,0x3FEF4B639EDDB726ull,0x3FEF3F322BAB907Cull,
		0x3FEF329C0558E969ull,0x3FEF25A15475D3B1ull,0x3FEF184242D66526ull,0x3FEF0A7EFB9230D7ull,0x3FEEFC57AB03BC37ull,0x3FEEEDCC7EC7F026ull,0x3FEEDEDDA5BD85F6ull,0x3FEECF8B5004705Cull,
		0x3FEEBFD5AEFD405Aull,0x3FEEAFBCF548861Full,0x3FEE9F4156C62DDAull,0x3FEE8E630894D893ull,0x3FEE7D22411130F5ull,0x3FEE6B7F37D53C1Full,0x3FEE597A25B7A677ull,0x3FEE471344CB0C75ull,
		0x3FEE344AD05D3F86ull,0x3FEE212104F686E5ull,0x3FEE0D962058DC8Cull,0x3FEDF9AA617F262Bull,0x3FEDE55E089C6A31ull,0x3FEDD0B1571B00E8ull,0x3FEDBBA48F9BC19Full,0x3FEDA637F5F52BF9ull,
		0x3FED906BCF328D46ull,0x3FED7A4061932205ull,0x3FED63B5F4893380ull,0x3FED4CCCD0B9318Bull,0x3FED35853FF8C869ull,0x3FED1DDF8D4DF2DCull,0x3FED05DC04EE0859ull,0x3FECED7AF43CC773ull,
		0x3FECD4BCA9CB5C71ull,0x3FECBBA17557641Cull,0x3FECA229A7C9EAC3ull,0x3FEC88559336677Cull,0x3FEC6E258AD9B3A1ull,0x3FEC5399E318FE91ull,0x3FEC38B2F180BDB1ull,0x3FEC1D710CC398BEull,
		0x3FEC01D48CB95263ull,0x3FEBE5DDCA5DAD22ull,0x3FEBC98D1FCF4C8Full,0x3FEBACE2E84E92E3ull,0x3FEB8FDF803C7AE7ull,0x3FEB728345196E3Eull,0x3FEB54CE95841811ull,0x3FEB36C1D1383422ull,
		0x3FEB185D590D5A44ull,0x3FEAF9A18EF5C645ull,0x3FEADA8ED5FD1C43ull,0x3FEABB2592472980ull,0x3FEA9B66290EA1A3ull,0x3FEA7B5100A3D882ull,0x3FEA5AE6806B7862ull,0x3FEA3A2710DD34C6ull,
		0x3FEA19131B8279C4ull,0x3FE9F7AB0AF517E3ull,0x3FE9D5EF4ADDEC96ull,0x3FE9B3E047F38741ull,0x3FE9917E6FF8CAE3ull,0x3FE96ECA31BB8C5Bull,0x3FE94BC3FD132D51ull,0x3FE9286C42DF33C3ull,
		0x3FE904C37505DE4Bull,0x3FE8E0CA0672B507ull,0x3FE8BC806B151741ull,0x3FE897E717DEC5D0ull,0x3FE872FE82C26A36ull,0x3FE84DC722B21A82ull,0x3FE828416F9DD9FFull,0x3FE8026DE27216ABull,
		0x3FE7DC4CF5162385ull,0x3FE7B5DF226AAFAFull,0x3FE78F24E6483A6Full,0x3FE7681EBD7D8412ull,0x3FE740CD25CDFBB2ull,0x3FE719309DF029E7ull,0x3FE6F149A58C1872ull,0x3FE6C918BD39B6CFull,
		0x3FE6A09E667F3BCDull,0x3FE677DB23CF8422ull,0x3FE64ECF78886E06ull,0x3FE6257BE8F131D5ull,0x3FE5FBE0FA38B7C6ull,0x3FE5D1FF3273EAB5ull,0x3FE5A7D7189C0802ull,0x3FE57D69348CECA0ull,
		0x3FE552B60F035F34ull,0x3FE527BE319B5775ull,0x3FE4FC8226CE42AAull,0x3FE4D10279F1456Cull,0x3FE4A53FB7337A9Bull,0x3FE4793A6B9C2F9Cull,0x3FE44CF325091DD6ull,0x3FE4206A722CA183ull,
		0x3FE3F3A0E28BEDD1ull,0x3FE3C697067D3E5Bull,0x3FE3994D6F260600ull,0x3FE36BC4AE791B21ull,0x3FE33DFD5734E146ull,0x3FE30FF7FCE17035ull,0x3FE2E1B533CEB87Eull,0x3FE2B3359112A581ull,
		0x3FE28479AA873CFEull,0x3FE2558216C8BC20ull,0x3FE2264F6D33B221ull,0x3FE1F6E245E3187Full,0x3FE1C73B39AE68C8ull,0x3FE1975AE227B00Eull,0x3FE16741D9999FF8ull,0x3FE136F0BB059D86ull,
		0x3FE106682221CD8Aull,0x3FE0D5A8AB571ED6ull,0x3FE0A4B2F3BF522Bull,0x3FE073879922FFEEull,0x3FE0422739F79BACull,0x3FE01092755D7570ull,0x3FDFBD93D63B71DCull,0x3FDF599C7750D525ull,
		0x3FDEF5401024C4F4ull,0x3FDE907FE4268A9Cull,0x3FDE2B5D3806F63Bull,0x3FDDC5D951B44853ull,0x3FDD5FF578561769ull,0x3FDCF9B2F44931B3ull,0x3FDC93130F1B7ADEull,0x3FDC2C171387C600ull,
		0x3FDBC4C04D71ABC1ull,0x3FDB5D1009E15CC0ull,0x3FDAF50796FF7054ull,0x3FDA8CA84410AFA4ull,0x3FDA23F36171DD2Dull,0x3FD9BAEA409378C3ull,0x3FD9518E33F58017ull,0x3FD8E7E08F232BD5ull,
		0x3FD87DE2A6AEA963ull,0x3FD81395D02CD150ull,0x3FD7A8FB6230DA82ull,0x3FD73E14B4480A32ull,0x3FD6D2E31EF560C0ull,0x3FD66767FBAD436Full,0x3FD5FBA4A4D12317ull,0x3FD58F9A75AB1FDDull,
		0x3FD5234ACA69A9FEull,0x3FD4B6B7001B1FB3ull,0x3FD449E074A9684Eull,0x3FD3DCC886D58C89ull,0x3FD36F7096334C28ull,0x3FD301DA0324B0F0ull,0x3FD294062ED59F06ull,0x3FD225F67B3762C0ull,
		0x3FD1B7AC4AFC3C02ull,0x3FD149290192E725ull,0x3FD0DA6E03222383ull,0x3FD06B7CB48437B2ull,0x3FCFF8ACF684E6EDull,0x3FCF19F97B215F1Bull,0x3FCE3AE1C4919687ull,0x3FCD5B68A1CC5001ull,
		0x3FCC7B90E3024582ull,0x3FCB9B5D5995173Eull,0x3FCABAD0D80E36CEull,0x3FC9D9EE3215CEABull,0x3FC8F8B83C69A60Bull,0x3FC81731CCD4013Dull,0x3FC7355DBA227EA6ull,0x3FC6533EDC1CF07Bull,
		0x3FC570D80B7C3350ull,0x3FC48E2C21E101A4ull,0x3FC3AB3DF9CAC47Dull,0x3FC2C8106E8E613Aull,0x3FC1E4A65C4D04B5ull,0x3FC101029FEAEBCBull,0x3FC01D281706297Eull,0x3FBE72333FDAD55Full,
		0x3FBCA9B4332D6F61ull,0x3FBAE0D8C72C6757ull,0x3FB917A6BC29B42Cull,0x3FB74E23D38E7485ull,0x3FB58455CFC8625Dull,0x3FB3BA427437435Dull,0x3FB1EFEF851A5624ull,0x3FB02562C77DBCC2ull,
		0x3FACB544024FC940ull,0x3FA91F65F10DD814ull,0x3FA58936E93C0B44ull,0x3FA1F2C279E5B3E0ull,0x3F9CB82865EB9D38ull,0x3F958A6F4A2381B0ull,0x3F8CB8E185D04D38ull,0x3F7CB90FCE395705ull
	};
  #endif

	#include "radix32_ditN_cy_dif1_asm.h"

#endif	// USE_SSE2

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data:
		int iter;
		int tid;
		int ndivr;
		int target_idx, target_set;	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
		double target_cy;

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
	#endif
	#ifdef USE_AVX
		int kcycle[ODD_RADIX];
		int lcycle[ODD_RADIX];
	#endif
	#ifdef USE_AVX512
		int mcycle[ODD_RADIX];
		int ncycle[ODD_RADIX];
		int ocycle[ODD_RADIX];
		int pcycle[ODD_RADIX];
	#endif

	// double data:
		double maxerr;
		double scale;
		double prp_mult;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		double *wts_mult, *inv_mult;
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
	#ifdef USE_AVX512
		double cy_dat[2*RADIX+8] __attribute__ ((__aligned__(8)));
	#else
		double cy_dat[2*RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	#endif
	};

#endif

/***************/

int radix224_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-224 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-224 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix224_ditN_cy_dif1";
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs mers-value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 14|28|56 for avx512,avx,sse, respectively.
	// But fixed-incr too restrictive here, so 'divide 14|28|56 into pieces' via increment-array whose elts sum to 14|28|56:
	const int *incr,*inc_arr = 0x0;	// Fermat-mod only uses inc_arr in AVX+ mode, init = 0 to quiet 'uninit' warnings.
  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
	const int incr_long[] = {14}, incr_med[] = {7,7}, incr_short[] = {3,4,3,4}, incr_hiacc[] = {2,2,3,2,2,3};
  #elif defined(USE_AVX)
	const int incr_long[] = {14,14}, incr_med[] = {7,7,7,7}, incr_short[] = {3,4,3,4,3,4,3,4}, incr_hiacc[] = {0};
  #elif defined(USE_ARM_V8_SIMD)
	const int incr_long[] = {14,14,14,14}, incr_med[] = {7,7,7,7,7,7,7,7}, incr_short[] = {5,4,5,5,4,5,5,4,5,5,4,5}, incr_hiacc[] = {2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3};
  #else
	const int incr_long[] = {14,14,14,14}, incr_med[] = {7,7,7,7,7,7,7,7}, incr_short[] = {3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4}, incr_hiacc[] = {0};
  #endif
	// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
  #ifdef USE_AVX512
	// For nfold > 2,  RADIX/8 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
	int nfold = USE_SHORT_CY_CHAIN;
	const int nexec_long[] = {28}, nexec_med[] = {14,14}, nexec_short[] = {7,7,7,7}, nexec_hiacc[] = {4,3,4,3,4,3,4,3};	// nfold = 0-3
  #elif defined(USE_AVX)
	// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
	int nfold = USE_SHORT_CY_CHAIN + 1;
	const int nexec_long[] = {28,28}, nexec_med[] = {14,14,14,14}, nexec_short[] = {7,7,7,7,7,7,7,7}, nexec_hiacc[] = {4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3};	// nfold = 1-4
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		if(USE_SHORT_CY_CHAIN == 0)
			inc_arr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			inc_arr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			inc_arr = incr_short;
		else
			inc_arr = incr_hiacc;
	} else {	// MODULUS_TYPE_FERMAT:
	#ifdef USE_AVX
		if(USE_SHORT_CY_CHAIN == 0)
			inc_arr = nexec_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			inc_arr = nexec_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			inc_arr = nexec_short;
		else
			inc_arr = nexec_hiacc;
	#endif
	}
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p04 offset for loop control
#ifndef MULTITHREAD
	double re,im;
	static int t_offsets[32], dif_offsets[RADIX], dit_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13 elts:
	static int dif_p20_cperms[26], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
	static int dit_p20_cperms[26], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
#endif
	static double radix_inv, n2inv;
	uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;

#if !defined(USE_SSE2) && !defined(MULTITHREAD)
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;	// tmps for scalar-double DFT macros
#endif
// FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) means non-Nussbaumer radix-7, uses these sincos constants:
#if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/

#elif defined(USE_SSE2) || !LO_ADD	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:
	#warning LO_ADD = 0 defined at compile time ... using low-mul Nussbaumer-style 7-DFT. [NOTE: SSE2-SIMD build *requires* this]
	// SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation:
	const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
				 	cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
				 	cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
				 	cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
				/* Switch the sign of ss3 in these: */
				 	sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
				 	sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
				 	sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
				 	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#else
	#error Unhandled combination of preprocessor flags!	// Just in case I ever 'decide' to leave some holes in the above PP logic tree
#endif
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *addr, *addi;
	struct complex t[RADIX], *tptr;
	int *itmp,*itm2;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,k1,k2,m,m2;
  #ifdef USE_AVX512
	double t0,t1,t2,t3;
   #ifdef CARRY_16_WAY
	static struct uint32x16 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	static struct uint32x8  *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #endif
  #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it;
	static int ii[ODD_RADIX] = {-1,-1,-1,-1,-1,-1,-1};	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		,icycle[ODD_RADIX],ic_idx;
#ifdef USE_SSE2
	static int jcycle[ODD_RADIX],jc_idx;
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX],kc_idx;
	static int lcycle[ODD_RADIX],lc_idx;
  #endif
  #ifdef USE_AVX512
	static int mcycle[ODD_RADIX],mc_idx;
	static int ncycle[ODD_RADIX],nc_idx;
	static int ocycle[ODD_RADIX],oc_idx;
	static int pcycle[ODD_RADIX],pc_idx;
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
  #endif

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *two,*one,*sqrt2,*isrt2,*cc1,*ss1,*cc2,*ss2,*cc3,*ss3,	// radix-32 DFT trig consts
		*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,			// radix-7 DFT trig consts
		*max_err, *sse2_rnd, *half_arr,
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
		*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif
	vec_dbl
	#ifndef MULTITHREAD
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,
	#endif
		*tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
#else
	static int p0123[4];
#endif	// USE_SSE2

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy224_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
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
	static double *_cy_r[RADIX],*_cy_i[RADIX];
	if(!_jhi) {
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
	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	// v18: If use residue shift in context of Pépin test, need prp_mult = 2 whenever the 'shift = 2*shift + random[0,1]' update gets a 1-bit in the random slot
	if((TEST_TYPE == TEST_TYPE_PRIMALITY && MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	|| (TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;	ndivr_inv = (double)RADIX/n;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		fprintf(stderr,"ERROR: iter = %10d; NWT_BITS does not divide N/RADIX in %s.\n",iter,func);
		err = ERR_SKIP_RADIX_SET;
		return(err);
	}

	if(p != psave || n != nsave
	#ifdef USE_PTHREAD	// Oct 2021: cf. radix176_ditN_cy_dif1.c for why I added this
		|| (tdat != 0x0 && tdat[0].wt1 != wt1)
	#endif
	) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
		/* To-do: Support #thread change here! */
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		ASSERT(LO_ADD,"LO_ADD");
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		i = 16;
	   #else
		i = 8;
	   #endif
	  #elif defined(USE_AVX)	// AVX LOACC: Make CARRY_8_WAY default here:
		i = 8;
	  #elif defined(USE_SSE2)	// AVX and SSE2 modes use 4-way carry macros
		i = 4;
	  #else	// Scalar-double mode:
		i = 1;
	  #endif

		// For n a power of 2 don't need to worry about 32-bit integer overflow in the sw*NDIVR term,
		// but for non-power-of-2 n we must cast-to-uint64 to avoid such overflows fubaring the result:
		struct qfloat qt,qn;
		qt = i64_to_q(i*(uint64)sw*NDIVR % n);
		qn = i64_to_q((int64) n);
		qt = qfdiv(qt, qn);		// x = (sw*NDIVR (mod n))/n
		qt = qfmul(qt, QLN2);	// x*ln(2)...
		qt = qfexp(qt);			// ...and get 2^x via exp[x*ln(2)].
		wts_mult[0] = qfdbl(qt);		// a = 2^(x/n), with x = sw
		inv_mult[0] = qfdbl(qfinv(qt));	// Double-based inversion (1.0 / wts_mult_a[0]) often gets LSB wrong
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		//curr have w, 2/w, separate-mul-by-1-or-0.5 gives [w,w/2] and [1/w,2/w] for i = 0,1, resp:
		wts_mult[1] = 0.5*wts_mult[0];
		inv_mult[1] = 2.0*inv_mult[0];
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the largest
		power of 2 that is <= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else {
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> i;
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0 ... likely more threads than this leading radix can handle.", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD
		if(tdat == 0x0) {
			j = (uint32)sizeof(struct cy_thread_data_t);
			tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

			// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
			// so on that platform try to be clever and interleave main-thread and threadpool-work processing
			#if 0//def OS_TYPE_MACOSX

				if(CY_THREADS > 1) {
					main_work_units = CY_THREADS/2;
					pool_work_units = CY_THREADS - main_work_units;
					ASSERT(0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
					printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
				} else {
					main_work_units = 1;
					printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
				}

			#else

				pool_work_units = CY_THREADS;
				ASSERT(0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

			#endif

			fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
		}
	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].iter = iter;
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
		// Jan 2018: To support PRP-cofactor testing, must allow user to feed different input vectors
		// to this carry routine, so set at runtime rather than in init-only block:
		//	tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].wts_mult = wts_mult;
			tdat[ithread].inv_mult = inv_mult;
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < RE_IM_STRIDE; l++) {
				if( ((intptr_t)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy_r = &tdat[ithread].cy_dat[l];
					tdat[ithread].cy_i = tdat[ithread].cy_r + RADIX;
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",SZ_VD,l);
					break;
				}
			}
			ASSERT(l < RE_IM_STRIDE, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((intptr_t)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((intptr_t)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of radix224_creals_in_local_store vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix224_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix224_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1c0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x1c0;
		two  	= tmp + 0x0;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one 	= tmp + 0x1;
		sqrt2	= tmp + 0x2;
		isrt2	= tmp + 0x3;
		cc2		= tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		ss2		= tmp + 0x5;
		cc1		= tmp + 0x6;
		ss1		= tmp + 0x7;
		cc3		= tmp + 0x8;
		ss3		= tmp + 0x9;
		dc0		= tmp + 0x0a;	// radix-7 DFT trig consts - Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here
		ds0		= tmp + 0x0b;
		dc1		= tmp + 0x0c;
		ds1		= tmp + 0x0d;
		dc2		= tmp + 0x0e;
		ds2		= tmp + 0x0f;
		dc3  	= tmp + 0x10;
		ds3		= tmp + 0x11;	// Add extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro:
		tmp += 0x16;	// sc_ptr += 0x396
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x1c;	tmp += 2*0x1c;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x38;	tmp += 2*0x38;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;  +30 = 338
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(396 + 70 + 2) = 0x408; This is where the value of half_arr_offset224 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x70;	tmp += 2*0x70;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays; +60 = 368 complex
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(396 + e0 + 2) = 0x478; This is where the value of half_arr_offset224 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(half_arr_offset224 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix224_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix224_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(cc2, c16  );	VEC_DBL_INIT(ss2, s16  );
		VEC_DBL_INIT(cc1, c32_1);	VEC_DBL_INIT(ss1, s32_1);
		VEC_DBL_INIT(cc3, c32_3);	VEC_DBL_INIT(ss3, s32_3);
	  #if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD)
		// AVX2 (i.e. FMA) means non-Nussbaumer radix-7, uses these sincos constants:
		VEC_DBL_INIT(dc0, uc1  );	VEC_DBL_INIT(ds0, us1);
		VEC_DBL_INIT(dc1, uc2  );	VEC_DBL_INIT(ds1, us2);
		VEC_DBL_INIT(dc2, uc3  );	VEC_DBL_INIT(ds2, us3);
		VEC_DBL_INIT(dc3, 0.0  );	VEC_DBL_INIT(ds3, 0.0);	// Unused in non-Nussbaumer mode
	  #else
		/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(dc0, cx0-1);	VEC_DBL_INIT(ds0, sx0);
		VEC_DBL_INIT(dc1, cx1  );	VEC_DBL_INIT(ds1, sx1);
		VEC_DBL_INIT(dc2, cx2  );	VEC_DBL_INIT(ds2, sx2);
		VEC_DBL_INIT(dc3, cx3  );	VEC_DBL_INIT(ds3, sx3);
	  #endif
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)two;	// #bytes in 1st of above block of consts
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = SZ_VD;	// sse2_rnd is a solo (in the SIMD-vector) datum
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

	  #ifdef HIACC
		#ifdef USE_AVX512
		  #error Fermat-mod HIACC mode only available in AVX/AVX2 builds!
		#endif
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		// Simple qfloat-based loop to crank out the roots which populate the radix224_avx_negadwt_consts table:
			struct qfloat qc,qs,qx,qy,qt,qn,qmul;
			qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
			qc = qfcos(qt);	qs = qfsin(qt);
			qx = QONE;		qy = QZRO;
			for(j = 0; j < RADIX; j++) {
				printf("j = %3u: cos = %#16" PRIX64 "\n",j,qfdbl_as_uint64(qx));
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
		tmp64 = radix224_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix224_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix224_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix224_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix224_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix224_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix224_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
	   #ifdef USE_AVX512

		// Init exp(j*I*Pi/2/RADIX), for j = 0-7:
		tmp = base_negacyclic_root + 16;	// First 16 slots reserved for Re/Im parts of the 8 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix224_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      4];	tmp->d4 = *(double *)&tmp64;	// cos(04*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      5];	tmp->d5 = *(double *)&tmp64;	// cos(05*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      6];	tmp->d6 = *(double *)&tmp64;	// cos(06*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      7];	tmp->d7 = *(double *)&tmp64;	// cos(07*I*Pi/(2*RADIX))
		++tmp;
		tmp->d0 = 0.0;
		tmp64 = radix224_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-4];	tmp->d4 = *(double *)&tmp64;	// sin(04*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-5];	tmp->d5 = *(double *)&tmp64;	// sin(05*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-6];	tmp->d6 = *(double *)&tmp64;	// sin(06*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-7];	tmp->d7 = *(double *)&tmp64;	// sin(07*I*Pi/(2*RADIX))
		++tmp;	// 0x480(base_negacyclic_root)
		tmp64 = radix224_avx_negadwt_consts[      8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(08*I*Pi/(2*RADIX))
		++tmp;	// 0x4c0(base_negacyclic_root)
		tmp64 = radix224_avx_negadwt_consts[RADIX-8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(08*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 16;	// reset to point to start of above block

	   #elif defined(USE_AVX)

		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix224_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix224_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix224_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;	// 0x140(base_negacyclic_root)
		tmp64 = radix224_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;	// 0x160(base_negacyclic_root)
		tmp64 = radix224_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 8;	// reset to point to start of above block

	   #endif

		nbytes = 4*SZ_VD;	// 2 SIMD-complex data

	  #endif	// HIACC toggle

		// Propagate the above consts to the remaining threads:
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// AVX?
	}		/************************************************************************/
	else	/*                MODULUS_TYPE_MERSENNE:                                */
	{		/************************************************************************/
	#ifdef USE_AVX512
		// Each lookup-category in the 'mini-tables' used in AVX mode balloons from 16x32-bytes to 64x64-bytes,
		// so switch to an opmask-based scheme which starts with e.g. a broadcast constant and onditional doubling.
		// Here are the needed consts and opmasks:
		// [1] Fwd-wt multipliers: Init = 0.50 x 8, anytime AVX-style lookup into 1st table below would have bit = 0, double the corr. datum
		// [2] Inv-wt multipliers: Init = 0.25 x 8, anytime AVX-style lookup into 2nd table below would have bit = 0, double the corr. datum
		// [3] Fwd-base mults: Init = base[0] x 8, anytime AVX-style lookup into 3rd table below would have bit = 1, double the corr. datum
		// [4] Inv-base mults: Init = binv[1] x 8, anytime AVX-style lookup into 4th table below would have bit = 0, double the corr. datum
		// [5] [LOACC] Init = wts_mult[1] x 8, anytime AVX-style lookup into 5th table below would have bit = 0, double the corr. datum
		// [6] [LOACC] Init = inv_mult[0] x 8, anytime AVX-style lookup into 6th table below would have bit = 1, double the corr. datum
		nbytes = 0;
	#elif defined(USE_AVX)
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
		// In LOACC mode, put wts_mult and their inverses in the first 32 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[0];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[0];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[0];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[0];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		tmp->d0 = wts_mult[1];	tmp->d1 = wts_mult[1];	tmp->d2 = wts_mult[1];	tmp->d3 = wts_mult[1];	++tmp;
		/* inv_mult: */
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[0];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[0];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[0];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[0];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		tmp->d0 = inv_mult[1];	tmp->d1 = inv_mult[1];	tmp->d2 = inv_mult[1];	tmp->d3 = inv_mult[1];	++tmp;
		nbytes = 96 << L2_SZ_VD;

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
		// In LOACC mode, put wts_mult and their inverses in the first 8 slots below in place of the 1/2-stuff:
		/* wts_mult:*/
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[0];	++ctmp;
		ctmp->re = wts_mult[0];	ctmp->im = wts_mult[1];	++ctmp;
		ctmp->re = wts_mult[1];	ctmp->im = wts_mult[1];	++ctmp;
		/* inv_mult:*/
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[0];	++ctmp;
		ctmp->re = inv_mult[0];	ctmp->im = inv_mult[1];	++ctmp;
		ctmp->re = inv_mult[1];	ctmp->im = inv_mult[1];	++ctmp;
		nbytes = 24 << L2_SZ_VD;

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

	nbytes = 4 << L2_SZ_VD;

  #ifdef USE_AVX512
   #ifdef CARRY_16_WAY
	n_minus_sil   = (struct uint32x16*)sse_n + 1;
	n_minus_silp1 = (struct uint32x16*)sse_n + 2;
	sinwt         = (struct uint32x16*)sse_n + 3;
	sinwtm1       = (struct uint32x16*)sse_n + 4;
	nbytes += 256;
   #else
	n_minus_sil   = (struct uint32x8 *)sse_n + 1;
	n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
	sinwt         = (struct uint32x8 *)sse_n + 3;
	sinwtm1       = (struct uint32x8 *)sse_n + 4;
	nbytes += 128;
   #endif
  #elif defined(USE_AVX)
	n_minus_sil   = (struct uint32x4 *)sse_n + 1;
	n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
	sinwt         = (struct uint32x4 *)sse_n + 3;
	sinwtm1       = (struct uint32x4 *)sse_n + 4;
	nbytes += 64;
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
		NDIVR >>= 4;			pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
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

	#ifndef MULTITHREAD
		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 7 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

	/*** DIF indexing stuff: ***/

		dif_phi[0] =   0;		dit_phi[0] =   0;
		dif_phi[1] = pc0;		dit_phi[1] = p60;
		dif_phi[2] = pa0;		dit_phi[2] = pc0;
		dif_phi[3] = p80;		dit_phi[3] = p40;
		dif_phi[4] = p60;		dit_phi[4] = pa0;
		dif_phi[5] = p40;		dit_phi[5] = p20;
		dif_phi[6] = p20;		dit_phi[6] = p80;

		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0xc0<<1; dif_p20_cperms[l++] = 0xa0<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1; dif_p20_cperms[l++] = 0x20<<1; dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0xc0<<1; dif_p20_cperms[l++] = 0xa0<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0xd0<<1; dif_p20_cperms[l++] = 0xb0<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1; dif_p20_cperms[l++] = 0x10<<1; dif_p20_cperms[l++] = 0xd0<<1; dif_p20_cperms[l++] = 0xb0<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 5);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 5);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 6);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 6);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 0);

	   #else

		l = 0;
		dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40; dif_p20_cperms[l++] = p20; dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30; dif_p20_cperms[l++] = p10; dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 3) + 1);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 1);
		dif_p20_lo_offset[l++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 3) + 2);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 2);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 2);
		dif_p20_lo_offset[l++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 3) + 3);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 3);
		dif_p20_lo_offset[l++] = ((pe << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 4);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 4);
		dif_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 3) + 5);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 5);
		dif_p20_lo_offset[l++] = ((pf << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 3) + 6);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 6);
		dif_p20_lo_offset[l++] = ((pc << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 3) + 0);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 0);

	   #endif

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each radix-32 DFT.
		// Set 0: [0,1,2,3,4,5,6,7,9,8,b,a,d,c,f,e + p00],[2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p2+p10;	// <*** same rcol offsets as Set 6 lcol
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p3+p10;	// <*** Unique lcol offsets
		dif_offsets[0x02] = p2;		dif_offsets[0x12] = p1+p10;
		dif_offsets[0x03] = p3;		dif_offsets[0x13] =    p10;
		dif_offsets[0x04] = p4;		dif_offsets[0x14] = p6+p10;
		dif_offsets[0x05] = p5;		dif_offsets[0x15] = p7+p10;
		dif_offsets[0x06] = p6;		dif_offsets[0x16] = p5+p10;
		dif_offsets[0x07] = p7;		dif_offsets[0x17] = p4+p10;
		dif_offsets[0x08] = p9;		dif_offsets[0x18] = pb+p10;
		dif_offsets[0x09] = p8;		dif_offsets[0x19] = pa+p10;
		dif_offsets[0x0a] = pb;		dif_offsets[0x1a] = p8+p10;
		dif_offsets[0x0b] = pa;		dif_offsets[0x1b] = p9+p10;
		dif_offsets[0x0c] = pd;		dif_offsets[0x1c] = pf+p10;
		dif_offsets[0x0d] = pc;		dif_offsets[0x1d] = pe+p10;
		dif_offsets[0x0e] = pf;		dif_offsets[0x1e] = pc+p10;
		dif_offsets[0x0f] = pe;		dif_offsets[0x1f] = pd+p10;
		// Set 1: [d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pc0],[f,e,c,d,8,9,a,b,1,0,3,2,5,4,7,6 + pd0] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = pd;		dif_offsets[l+0x10] = pf+p10;	// <*** Unique rcol offsets
		dif_offsets[l+0x01] = pc;		dif_offsets[l+0x11] = pe+p10;	// <*** same lcol offsets as Set 2 rcol
		dif_offsets[l+0x02] = pf;		dif_offsets[l+0x12] = pc+p10;
		dif_offsets[l+0x03] = pe;		dif_offsets[l+0x13] = pd+p10;
		dif_offsets[l+0x04] = pb;		dif_offsets[l+0x14] = p8+p10;
		dif_offsets[l+0x05] = pa;		dif_offsets[l+0x15] = p9+p10;
		dif_offsets[l+0x06] = p8;		dif_offsets[l+0x16] = pa+p10;
		dif_offsets[l+0x07] = p9;		dif_offsets[l+0x17] = pb+p10;
		dif_offsets[l+0x08] = p2;		dif_offsets[l+0x18] = p1+p10;
		dif_offsets[l+0x09] = p3;		dif_offsets[l+0x19] =    p10;
		dif_offsets[l+0x0a] = p1;		dif_offsets[l+0x1a] = p3+p10;
		dif_offsets[l+0x0b] =  0;		dif_offsets[l+0x1b] = p2+p10;
		dif_offsets[l+0x0c] = p6;		dif_offsets[l+0x1c] = p5+p10;
		dif_offsets[l+0x0d] = p7;		dif_offsets[l+0x1d] = p4+p10;
		dif_offsets[l+0x0e] = p5;		dif_offsets[l+0x1e] = p7+p10;
		dif_offsets[l+0x0f] = p4;		dif_offsets[l+0x1f] = p6+p10;
		// Set 2: [6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + pb0],[d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pa0] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p6+p10;		dif_offsets[l+0x10] = pd;	// <*** same rcol offsets as Set 1 lcol
		dif_offsets[l+0x01] = p7+p10;		dif_offsets[l+0x11] = pc;	// <*** same lcol offsets as Set 3 rcol
		dif_offsets[l+0x02] = p5+p10;		dif_offsets[l+0x12] = pf;
		dif_offsets[l+0x03] = p4+p10;		dif_offsets[l+0x13] = pe;
		dif_offsets[l+0x04] = p1+p10;		dif_offsets[l+0x14] = pb;
		dif_offsets[l+0x05] =    p10;		dif_offsets[l+0x15] = pa;
		dif_offsets[l+0x06] = p3+p10;		dif_offsets[l+0x16] = p8;
		dif_offsets[l+0x07] = p2+p10;		dif_offsets[l+0x17] = p9;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = p2;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = p3;
		dif_offsets[l+0x0a] = pc+p10;		dif_offsets[l+0x1a] = p1;
		dif_offsets[l+0x0b] = pd+p10;		dif_offsets[l+0x1b] =  0;
		dif_offsets[l+0x0c] = p8+p10;		dif_offsets[l+0x1c] = p6;
		dif_offsets[l+0x0d] = p9+p10;		dif_offsets[l+0x1d] = p7;
		dif_offsets[l+0x0e] = pa+p10;		dif_offsets[l+0x1e] = p5;
		dif_offsets[l+0x0f] = pb+p10;		dif_offsets[l+0x1f] = p4;
		// Set 3: [4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p80],[6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + p90] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p4;		dif_offsets[l+0x10] = p6+p10;	// <*** same rcol offsets as Set 2 lcol
		dif_offsets[l+0x01] = p5;		dif_offsets[l+0x11] = p7+p10;	// <*** same lcol offsets as Set 4 rcol
		dif_offsets[l+0x02] = p6;		dif_offsets[l+0x12] = p5+p10;
		dif_offsets[l+0x03] = p7;		dif_offsets[l+0x13] = p4+p10;
		dif_offsets[l+0x04] = p2;		dif_offsets[l+0x14] = p1+p10;
		dif_offsets[l+0x05] = p3;		dif_offsets[l+0x15] =    p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = p3+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = p2+p10;
		dif_offsets[l+0x08] = pd;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = pc;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = pf;		dif_offsets[l+0x1a] = pc+p10;
		dif_offsets[l+0x0b] = pe;		dif_offsets[l+0x1b] = pd+p10;
		dif_offsets[l+0x0c] = pb;		dif_offsets[l+0x1c] = p8+p10;
		dif_offsets[l+0x0d] = pa;		dif_offsets[l+0x1d] = p9+p10;
		dif_offsets[l+0x0e] = p8;		dif_offsets[l+0x1e] = pa+p10;
		dif_offsets[l+0x0f] = p9;		dif_offsets[l+0x1f] = pb+p10;
		// Set 4: [b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p70],[4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p60] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = pb+p10;		dif_offsets[l+0x10] = p4;	// <*** same rcol offsets as Set 3 lcol
		dif_offsets[l+0x01] = pa+p10;		dif_offsets[l+0x11] = p5;	// <*** same lcol offsets as Set 5 rcol
		dif_offsets[l+0x02] = p8+p10;		dif_offsets[l+0x12] = p6;
		dif_offsets[l+0x03] = p9+p10;		dif_offsets[l+0x13] = p7;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = p2;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = p3;
		dif_offsets[l+0x06] = pc+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = pd+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = p6+p10;		dif_offsets[l+0x18] = pd;
		dif_offsets[l+0x09] = p7+p10;		dif_offsets[l+0x19] = pc;
		dif_offsets[l+0x0a] = p5+p10;		dif_offsets[l+0x1a] = pf;
		dif_offsets[l+0x0b] = p4+p10;		dif_offsets[l+0x1b] = pe;
		dif_offsets[l+0x0c] = p1+p10;		dif_offsets[l+0x1c] = pb;
		dif_offsets[l+0x0d] =    p10;		dif_offsets[l+0x1d] = pa;
		dif_offsets[l+0x0e] = p3+p10;		dif_offsets[l+0x1e] = p8;
		dif_offsets[l+0x0f] = p2+p10;		dif_offsets[l+0x1f] = p9;
		// Set 5: [9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p40],[b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p50] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p9;		dif_offsets[l+0x10] = pb+p10;	// <*** same rcol offsets as Set 4 lcol
		dif_offsets[l+0x01] = p8;		dif_offsets[l+0x11] = pa+p10;	// <*** same lcol offsets as Set 6 rcol
		dif_offsets[l+0x02] = pb;		dif_offsets[l+0x12] = p8+p10;
		dif_offsets[l+0x03] = pa;		dif_offsets[l+0x13] = p9+p10;
		dif_offsets[l+0x04] = pd;		dif_offsets[l+0x14] = pf+p10;
		dif_offsets[l+0x05] = pc;		dif_offsets[l+0x15] = pe+p10;
		dif_offsets[l+0x06] = pf;		dif_offsets[l+0x16] = pc+p10;
		dif_offsets[l+0x07] = pe;		dif_offsets[l+0x17] = pd+p10;
		dif_offsets[l+0x08] = p4;		dif_offsets[l+0x18] = p6+p10;
		dif_offsets[l+0x09] = p5;		dif_offsets[l+0x19] = p7+p10;
		dif_offsets[l+0x0a] = p6;		dif_offsets[l+0x1a] = p5+p10;
		dif_offsets[l+0x0b] = p7;		dif_offsets[l+0x1b] = p4+p10;
		dif_offsets[l+0x0c] = p2;		dif_offsets[l+0x1c] = p1+p10;
		dif_offsets[l+0x0d] = p3;		dif_offsets[l+0x1d] =    p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = p3+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = p2+p10;
		// Set 6: [2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p30],[9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p20] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p2+p10;		dif_offsets[l+0x10] = p9;	// <*** same rcol offsets as Set 5 lcol
		dif_offsets[l+0x01] = p3+p10;		dif_offsets[l+0x11] = p8;	// <*** same lcol offsets as Set 0 rcol
		dif_offsets[l+0x02] = p1+p10;		dif_offsets[l+0x12] = pb;
		dif_offsets[l+0x03] =    p10;		dif_offsets[l+0x13] = pa;
		dif_offsets[l+0x04] = p6+p10;		dif_offsets[l+0x14] = pd;
		dif_offsets[l+0x05] = p7+p10;		dif_offsets[l+0x15] = pc;
		dif_offsets[l+0x06] = p5+p10;		dif_offsets[l+0x16] = pf;
		dif_offsets[l+0x07] = p4+p10;		dif_offsets[l+0x17] = pe;
		dif_offsets[l+0x08] = pb+p10;		dif_offsets[l+0x18] = p4;
		dif_offsets[l+0x09] = pa+p10;		dif_offsets[l+0x19] = p5;
		dif_offsets[l+0x0a] = p8+p10;		dif_offsets[l+0x1a] = p6;
		dif_offsets[l+0x0b] = p9+p10;		dif_offsets[l+0x1b] = p7;
		dif_offsets[l+0x0c] = pf+p10;		dif_offsets[l+0x1c] = p2;
		dif_offsets[l+0x0d] = pe+p10;		dif_offsets[l+0x1d] = p3;
		dif_offsets[l+0x0e] = pc+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = pd+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/
		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0xc0<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0xa0<<1; dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0xc0<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x60<<1;
		dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0xb0<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x90<<1; dit_p20_cperms[l++] = 0xd0<<1; dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0xb0<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x90<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 6);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 1);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 2);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 3);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 5);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 6);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 5);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 6);

	   #else

		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 3) + 6);
		dit_p20_lo_offset[l++] = ((pe << 3) + 0);
		dit_p20_lo_offset[l++] = ((pd << 3) + 1);
		dit_p20_lo_offset[l++] = ((pc << 3) + 2);
		dit_p20_lo_offset[l++] = ((pb << 3) + 3);
		dit_p20_lo_offset[l++] = ((pa << 3) + 4);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 6);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 0);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 4);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 6);

	   #endif

	// dit_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		l = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p70],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p60] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pb+p10;		dit_offsets[l+0x10] = pb;
		dit_offsets[l+0x01] = pa+p10;		dit_offsets[l+0x11] = pa;
		dit_offsets[l+0x02] = p9+p10;		dit_offsets[l+0x12] = p9;
		dit_offsets[l+0x03] = p8+p10;		dit_offsets[l+0x13] = p8;
		dit_offsets[l+0x04] = pd+p10;		dit_offsets[l+0x14] = pd;
		dit_offsets[l+0x05] = pc+p10;		dit_offsets[l+0x15] = pc;
		dit_offsets[l+0x06] = pe+p10;		dit_offsets[l+0x16] = pe;
		dit_offsets[l+0x07] = pf+p10;		dit_offsets[l+0x17] = pf;
		dit_offsets[l+0x08] = p3+p10;		dit_offsets[l+0x18] = p3;
		dit_offsets[l+0x09] = p2+p10;		dit_offsets[l+0x19] = p2;
		dit_offsets[l+0x0a] = p1+p10;		dit_offsets[l+0x1a] = p1;
		dit_offsets[l+0x0b] =    p10;		dit_offsets[l+0x1b] =  0;
		dit_offsets[l+0x0c] = p5+p10;		dit_offsets[l+0x1c] = p5;
		dit_offsets[l+0x0d] = p4+p10;		dit_offsets[l+0x1d] = p4;
		dit_offsets[l+0x0e] = p6+p10;		dit_offsets[l+0x1e] = p6;
		dit_offsets[l+0x0f] = p7+p10;		dit_offsets[l+0x1f] = p7;
		// Set 2: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + pc0],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + pd0] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pd;		dit_offsets[l+0x10] = p5+p10;
		dit_offsets[l+0x01] = pc;		dit_offsets[l+0x11] = p4+p10;
		dit_offsets[l+0x02] = pe;		dit_offsets[l+0x12] = p6+p10;
		dit_offsets[l+0x03] = pf;		dit_offsets[l+0x13] = p7+p10;
		dit_offsets[l+0x04] = p9;		dit_offsets[l+0x14] = p1+p10;
		dit_offsets[l+0x05] = p8;		dit_offsets[l+0x15] =    p10;
		dit_offsets[l+0x06] = pa;		dit_offsets[l+0x16] = p2+p10;
		dit_offsets[l+0x07] = pb;		dit_offsets[l+0x17] = p3+p10;
		dit_offsets[l+0x08] = p5;		dit_offsets[l+0x18] = p9+p10;
		dit_offsets[l+0x09] = p4;		dit_offsets[l+0x19] = p8+p10;
		dit_offsets[l+0x0a] = p6;		dit_offsets[l+0x1a] = pa+p10;
		dit_offsets[l+0x0b] = p7;		dit_offsets[l+0x1b] = pb+p10;
		dit_offsets[l+0x0c] = p1;		dit_offsets[l+0x1c] = pe+p10;
		dit_offsets[l+0x0d] =  0;		dit_offsets[l+0x1d] = pf+p10;
		dit_offsets[l+0x0e] = p2;		dit_offsets[l+0x1e] = pc+p10;
		dit_offsets[l+0x0f] = p3;		dit_offsets[l+0x1f] = pd+p10;
		// Set 3: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p40],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p50] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p9;		dit_offsets[l+0x10] = p1+p10;
		dit_offsets[l+0x01] = p8;		dit_offsets[l+0x11] =    p10;
		dit_offsets[l+0x02] = pa;		dit_offsets[l+0x12] = p2+p10;
		dit_offsets[l+0x03] = pb;		dit_offsets[l+0x13] = p3+p10;
		dit_offsets[l+0x04] = pe;		dit_offsets[l+0x14] = p6+p10;
		dit_offsets[l+0x05] = pf;		dit_offsets[l+0x15] = p7+p10;
		dit_offsets[l+0x06] = pc;		dit_offsets[l+0x16] = p4+p10;
		dit_offsets[l+0x07] = pd;		dit_offsets[l+0x17] = p5+p10;
		dit_offsets[l+0x08] = p1;		dit_offsets[l+0x18] = pe+p10;
		dit_offsets[l+0x09] =  0;		dit_offsets[l+0x19] = pf+p10;
		dit_offsets[l+0x0a] = p2;		dit_offsets[l+0x1a] = pc+p10;
		dit_offsets[l+0x0b] = p3;		dit_offsets[l+0x1b] = pd+p10;
		dit_offsets[l+0x0c] = p6;		dit_offsets[l+0x1c] = pa+p10;
		dit_offsets[l+0x0d] = p7;		dit_offsets[l+0x1d] = pb+p10;
		dit_offsets[l+0x0e] = p4;		dit_offsets[l+0x1e] = p8+p10;
		dit_offsets[l+0x0f] = p5;		dit_offsets[l+0x1f] = p9+p10;
		// Set 4: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pb0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pa0] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p6+p10;		dit_offsets[l+0x10] = p6;
		dit_offsets[l+0x01] = p7+p10;		dit_offsets[l+0x11] = p7;
		dit_offsets[l+0x02] = p4+p10;		dit_offsets[l+0x12] = p4;
		dit_offsets[l+0x03] = p5+p10;		dit_offsets[l+0x13] = p5;
		dit_offsets[l+0x04] = p2+p10;		dit_offsets[l+0x14] = p2;
		dit_offsets[l+0x05] = p3+p10;		dit_offsets[l+0x15] = p3;
		dit_offsets[l+0x06] =    p10;		dit_offsets[l+0x16] =  0;
		dit_offsets[l+0x07] = p1+p10;		dit_offsets[l+0x17] = p1;
		dit_offsets[l+0x08] = pa+p10;		dit_offsets[l+0x18] = pa;
		dit_offsets[l+0x09] = pb+p10;		dit_offsets[l+0x19] = pb;
		dit_offsets[l+0x0a] = p8+p10;		dit_offsets[l+0x1a] = p8;
		dit_offsets[l+0x0b] = p9+p10;		dit_offsets[l+0x1b] = p9;
		dit_offsets[l+0x0c] = pc+p10;		dit_offsets[l+0x1c] = pc;
		dit_offsets[l+0x0d] = pd+p10;		dit_offsets[l+0x1d] = pd;
		dit_offsets[l+0x0e] = pf+p10;		dit_offsets[l+0x1e] = pf;
		dit_offsets[l+0x0f] = pe+p10;		dit_offsets[l+0x1f] = pe;
		// Set 5: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p20] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p2+p10;		dit_offsets[l+0x10] = p2;
		dit_offsets[l+0x01] = p3+p10;		dit_offsets[l+0x11] = p3;
		dit_offsets[l+0x02] =    p10;		dit_offsets[l+0x12] =  0;
		dit_offsets[l+0x03] = p1+p10;		dit_offsets[l+0x13] = p1;
		dit_offsets[l+0x04] = p4+p10;		dit_offsets[l+0x14] = p4;
		dit_offsets[l+0x05] = p5+p10;		dit_offsets[l+0x15] = p5;
		dit_offsets[l+0x06] = p7+p10;		dit_offsets[l+0x16] = p7;
		dit_offsets[l+0x07] = p6+p10;		dit_offsets[l+0x17] = p6;
		dit_offsets[l+0x08] = pc+p10;		dit_offsets[l+0x18] = pc;
		dit_offsets[l+0x09] = pd+p10;		dit_offsets[l+0x19] = pd;
		dit_offsets[l+0x0a] = pf+p10;		dit_offsets[l+0x1a] = pf;
		dit_offsets[l+0x0b] = pe+p10;		dit_offsets[l+0x1b] = pe;
		dit_offsets[l+0x0c] = p8+p10;		dit_offsets[l+0x1c] = p8;
		dit_offsets[l+0x0d] = p9+p10;		dit_offsets[l+0x1d] = p9;
		dit_offsets[l+0x0e] = pb+p10;		dit_offsets[l+0x1e] = pb;
		dit_offsets[l+0x0f] = pa+p10;		dit_offsets[l+0x1f] = pa;
		// Set 6: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p80],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p90] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p4;		dit_offsets[l+0x10] = p8+p10;
		dit_offsets[l+0x01] = p5;		dit_offsets[l+0x11] = p9+p10;
		dit_offsets[l+0x02] = p7;		dit_offsets[l+0x12] = pb+p10;
		dit_offsets[l+0x03] = p6;		dit_offsets[l+0x13] = pa+p10;
		dit_offsets[l+0x04] =  0;		dit_offsets[l+0x14] = pf+p10;
		dit_offsets[l+0x05] = p1;		dit_offsets[l+0x15] = pe+p10;
		dit_offsets[l+0x06] = p3;		dit_offsets[l+0x16] = pd+p10;
		dit_offsets[l+0x07] = p2;		dit_offsets[l+0x17] = pc+p10;
		dit_offsets[l+0x08] = p8;		dit_offsets[l+0x18] =    p10;
		dit_offsets[l+0x09] = p9;		dit_offsets[l+0x19] = p1+p10;
		dit_offsets[l+0x0a] = pb;		dit_offsets[l+0x1a] = p3+p10;
		dit_offsets[l+0x0b] = pa;		dit_offsets[l+0x1b] = p2+p10;
		dit_offsets[l+0x0c] = pf;		dit_offsets[l+0x1c] = p7+p10;
		dit_offsets[l+0x0d] = pe;		dit_offsets[l+0x1d] = p6+p10;
		dit_offsets[l+0x0e] = pd;		dit_offsets[l+0x1e] = p5+p10;
		dit_offsets[l+0x0f] = pc;		dit_offsets[l+0x1f] = p4+p10;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#endif	// MULTITHREAD?

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

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays!");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in %s.\n", func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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
		ASSERT(_bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

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
			ASSERT(wts_idx_incr != 0, "wts_idx_incr init failed!");

		#ifdef USE_SSE2
			wts_idx_inc2 = wts_idx_incr << (2*L2_SZ_VD - 3);	/* In the SIMD version, use icycle0-6 as actual address
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
				jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);	tmp->d1 = wt_arr[jcycle[i]];
			#ifdef USE_AVX
				kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);	tmp->d2 = wt_arr[kcycle[i]];
				lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);	tmp->d3 = wt_arr[lcycle[i]];
			#endif
			#ifdef USE_AVX512
				mcycle[i] = lcycle[i] + wts_idx_incr;	mcycle[i] += ( (-(mcycle[i] < 0)) & nwt);	tmp->d4 = wt_arr[mcycle[i]];
				ncycle[i] = mcycle[i] + wts_idx_incr;	ncycle[i] += ( (-(ncycle[i] < 0)) & nwt);	tmp->d5 = wt_arr[ncycle[i]];
				ocycle[i] = ncycle[i] + wts_idx_incr;	ocycle[i] += ( (-(ocycle[i] < 0)) & nwt);	tmp->d6 = wt_arr[ocycle[i]];
				pcycle[i] = ocycle[i] + wts_idx_incr;	pcycle[i] += ( (-(pcycle[i] < 0)) & nwt);	tmp->d7 = wt_arr[pcycle[i]];
			#endif
			}

			// Propagate the above wts-consts to the remaining threads:
			nbytes = ODD_RADIX*SZ_VD;
			tmp = half_arr;
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

			tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */

		  #ifdef USE_AVX512

			// Each transposed-data octet in the AVX-512 carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-8 index subsequences, which will be accessed via their head element
			// by the [ijklmnop]cycle* index octets in the respective carry-macro call:
			tm2 = tmp + ODD_RADIX;
			for(i = 0; i < ODD_RADIX; i++, tmp++, tm2++) {
				tmp->d0 = bs_arr   [i];	tmp->d1 = bs_arr   [(i+1)%ODD_RADIX];	tmp->d2 = bs_arr   [(i+2)%ODD_RADIX];	tmp->d3 = bs_arr   [(i+3)%ODD_RADIX];	tmp->d4 = bs_arr   [(i+4)%ODD_RADIX];	tmp->d5 = bs_arr   [(i+5)%ODD_RADIX];	tmp->d6 = bs_arr   [(i+6)%ODD_RADIX];	tmp->d7 = bs_arr   [(i+7)%ODD_RADIX];
				tm2->d0 = bsinv_arr[i];	tm2->d1 = bsinv_arr[(i+1)%ODD_RADIX];	tm2->d2 = bsinv_arr[(i+2)%ODD_RADIX];	tm2->d3 = bsinv_arr[(i+3)%ODD_RADIX];	tm2->d4 = bsinv_arr[(i+4)%ODD_RADIX];	tm2->d5 = bsinv_arr[(i+5)%ODD_RADIX];	tm2->d6 = bsinv_arr[(i+6)%ODD_RADIX];	tm2->d7 = bsinv_arr[(i+7)%ODD_RADIX];
			}

		  #elif defined(USE_AVX)

			// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-4 index subsequences, which will be accessed via their head element
			// by the [ijkl]cycle* index quartets in the respective carry-macro call:
			tm2 = tmp + ODD_RADIX;
			for(i = 0; i < ODD_RADIX; i++, tmp++, tm2++) {
				tmp->d0 = bs_arr   [i];	tmp->d1 = bs_arr   [(i+1)%ODD_RADIX];	tmp->d2 = bs_arr   [(i+2)%ODD_RADIX];	tmp->d3 = bs_arr   [(i+3)%ODD_RADIX];
				tm2->d0 = bsinv_arr[i];	tm2->d1 = bsinv_arr[(i+1)%ODD_RADIX];	tm2->d2 = bsinv_arr[(i+2)%ODD_RADIX];	tm2->d3 = bsinv_arr[(i+3)%ODD_RADIX];
			}

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
				icycle[i] <<= L2_SZ_VD;	jcycle[i] <<= L2_SZ_VD;
			#ifdef USE_AVX
				kcycle[i] <<= L2_SZ_VD;	lcycle[i] <<= L2_SZ_VD;
			#endif
			#ifdef USE_AVX512
				mcycle[i] <<= L2_SZ_VD;	ncycle[i] <<= L2_SZ_VD;	ocycle[i] <<= L2_SZ_VD;	pcycle[i] <<= L2_SZ_VD;
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
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
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
			  #endif
			  #ifdef USE_AVX
				tdat[0].kcycle[i] = kcycle[i];
				tdat[0].lcycle[i] = lcycle[i];
			  #endif
			  #ifdef USE_AVX512
				tdat[0].mcycle[i] = mcycle[i];
				tdat[0].ncycle[i] = ncycle[i];
				tdat[0].ocycle[i] = ocycle[i];
				tdat[0].pcycle[i] = pcycle[i];
			  #endif
			}
			// For remaining threads, simulate the loop-evolution of the above indices.
			// Note that the non-thread-associated *cycle[] arry data will get changed fom their above-inited
			// values in the loop here, but that's OK because in || mode only the thread-associated values matter:
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				jstart = 0;
				jhi = NDIVR/CY_THREADS;	// Earlier setting = NDIVR/CY_THREADS/2 was for simulating bjmodn evolution, need 2x that here
				// Get value of (negative) increment resulting from (jhi-jstart)/stride execs of *cycle[] += wts_idx_inc* (mod nwt*):
			#ifndef USE_SSE2
				j = ((int64)wts_idx_incr * ( (jhi-jstart)>>(L2_SZ_VD-2) ) % nwt  );	// []>>(L2_SZ_VD-2) is fast subst. for []/stride
			#else
				j = ((int64)wts_idx_inc2 * ( (jhi-jstart)>>(L2_SZ_VD-2) ) % nwt16);	// Cast wts_idx_inc* to signed 64-bit to avoid
						// overflow of product; further compute (jhi-jstart)/stride prior to multiply to gain more bits-to-spare.
			#endif
				// khi = 1 for Fermat-mod, thus no outer loop needed here
				for(i = 0; i < ODD_RADIX; i++) {
				#ifndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:
					icycle[i] += j;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt);
				#else
					icycle[i] += j;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt16);
					jcycle[i] += j;		jcycle[i] += ( (-(int)((uint32)jcycle[i] >> 31)) & nwt16);
				  #ifdef USE_AVX
					kcycle[i] += j;		kcycle[i] += ( (-(int)((uint32)kcycle[i] >> 31)) & nwt16);
					lcycle[i] += j;		lcycle[i] += ( (-(int)((uint32)lcycle[i] >> 31)) & nwt16);
				  #endif
				  #ifdef USE_AVX512
					mcycle[i] += j;		mcycle[i] += ( (-(int)((uint32)mcycle[i] >> 31)) & nwt16);
					ncycle[i] += j;		ncycle[i] += ( (-(int)((uint32)ncycle[i] >> 31)) & nwt16);
					ocycle[i] += j;		ocycle[i] += ( (-(int)((uint32)ocycle[i] >> 31)) & nwt16);
					pcycle[i] += j;		pcycle[i] += ( (-(int)((uint32)pcycle[i] >> 31)) & nwt16);
				  #endif
				#endif
				}
				for(i = 0; i < ODD_RADIX; i++) {
					tdat[ithread].icycle[i] = icycle[i];
				  #ifdef USE_SSE2
					tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
					tdat[ithread].jcycle[i] = jcycle[i];
				  #endif
				  #ifdef USE_AVX
					tdat[ithread].kcycle[i] = kcycle[i];
					tdat[ithread].lcycle[i] = lcycle[i];
				  #endif
				  #ifdef USE_AVX512
					tdat[ithread].mcycle[i] = mcycle[i];
					tdat[ithread].ncycle[i] = ncycle[i];
					tdat[ithread].ocycle[i] = ocycle[i];
					tdat[ithread].pcycle[i] = pcycle[i];
				  #endif
				}
			}
			// Restore the original loop-start values of the cycle arrays, since we use these for init of inv-wts below:
			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] = tdat[0].icycle[i];
			  #ifdef USE_SSE2
				jcycle[i] = tdat[0].jcycle[i];
			  #endif
			  #ifdef USE_AVX
				kcycle[i] = tdat[0].kcycle[i];
				lcycle[i] = tdat[0].lcycle[i];
			  #endif
			  #ifdef USE_AVX512
				mcycle[i] = tdat[0].mcycle[i];
				ncycle[i] = tdat[0].ncycle[i];
				ocycle[i] = tdat[0].ocycle[i];
				pcycle[i] = tdat[0].pcycle[i];
			  #endif
			}
		}

	#endif	// USE_PTHREAD

		first_entry=FALSE;
	}	/* endif(first_entry) */

	// Jun 2018: If LL test and shift applied, compute target index for data-processing loop.
	// Note that only 1 thread of the carry-processing set will hit the target, but all need the same logic to check for a hit:
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY) {
		if(RES_SHIFT) {
			itmp64 = shift_word(a, n, p, RES_SHIFT, 0.0);	// Note return value (specifically high 7 bytes thereof) is an unpadded index
			target_idx = (int)(itmp64 >>  8);	// This still needs to be (mod NDIVR)'ed, but first use unmodded form to compute needed DWT weights
			// Compute wt = 2^(target_idx*sw % n)/n and its reciprocal:
			uint32 sw_idx_modn = ((uint64)target_idx*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
			double target_wtfwd = pow(2.0, sw_idx_modn*0.5*n2inv);	// 0.5*n2inv = 0.5/(n/2) = 1.0/n
			target_set = target_idx*ndivr_inv;	// Which of the [RADIX] independent sub-carry-chains contains the target index?
			target_idx -= target_set*NDIVR;		// Fast computation of target_idx = (target_idx % NDIVR)
			// Now compute the doubles-pointer offset of the target double w.r.to the SIMD s1p00-... data layout:
			tidx_mod_stride = target_idx & (stride-1);	// Stride a power of 2, so can use AND-minus-1 for mod
			target_idx -= tidx_mod_stride;
		//	printf("Iter %d: cy_shift = %d, target_idx,tidx_mod_stride,target_set = %d,%d,%d\n",iter,(itmp64 & 255),target_idx,tidx_mod_stride,target_set);
		#ifdef USE_AVX512
			tidx_mod_stride = br16[tidx_mod_stride];
		#elif defined(USE_AVX)
			tidx_mod_stride = br8[tidx_mod_stride];
		#elif defined(USE_SSE2)
			tidx_mod_stride = br4[tidx_mod_stride];
		#endif
			target_set = (target_set<<(L2_SZ_VD-2)) + tidx_mod_stride;
			target_cy  = target_wtfwd * (-(int)(2u << (itmp64 & 255)));
		} else {
			target_idx = target_set = 0;
			target_cy = -2.0;
		}
	}

/*...The radix-224 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy_r[i][ithread] = 0;
			_cy_i[i][ithread] = 0;
		}
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r[0][0] = -2;
	}
  #endif
	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

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
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_mers;	/* Cleanup loop assumes carryins propagate at most 4 words up. */
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
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_ferm;	/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}

		// Now that full_pass-dependent scale factor known, init inverse weights tiny-table used for Fermat-mod
		for(i = 0; i < ODD_RADIX; i++) {
			wtinv_arr[i] = scale*wt1[ii[i]];
		}

	// In threaded mode, the master *cycle[] values are unmodified during main loop exec; only the thread-associated
	// copies of these index arrays get modified. In non-threaded mode we must separately store copies of the masters
	// in order to solve the save/restore issue. We start from the (static, unmodified during loop) ii[]-index values:
	#ifndef MULTITHREAD
		for(i = 0; i < ODD_RADIX; i++) {
			/* Reinit *cycle indices their proper starting values - recall in SIMD mode these all are ( << L2_SZ_VD): */
			icycle[i] = i;
		#ifdef USE_SSE2
			jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
		  #ifdef USE_AVX
			kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
			lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
		   #ifdef USE_AVX512
			mcycle[i] = lcycle[i] + wts_idx_incr;	mcycle[i] += ( (-(mcycle[i] < 0)) & nwt);
			ncycle[i] = mcycle[i] + wts_idx_incr;	ncycle[i] += ( (-(ncycle[i] < 0)) & nwt);
			ocycle[i] = ncycle[i] + wts_idx_incr;	ocycle[i] += ( (-(ocycle[i] < 0)) & nwt);
			pcycle[i] = ocycle[i] + wts_idx_incr;	pcycle[i] += ( (-(pcycle[i] < 0)) & nwt);
			mcycle[i] <<= L2_SZ_VD;	ncycle[i] <<= L2_SZ_VD;	ocycle[i] <<= L2_SZ_VD;	pcycle[i] <<= L2_SZ_VD;
		   #endif
			kcycle[i] <<= L2_SZ_VD;	lcycle[i] <<= L2_SZ_VD;
		  #endif
			icycle[i] <<= L2_SZ_VD;	jcycle[i] <<= L2_SZ_VD;
		#endif
		}
	#endif

	#ifdef USE_SSE2
		// Remember: *cycle[] entries all << L2_SZ_VD here - must left-shift-on-the-fly before using:
		tm2 = half_arr + ODD_RADIX;
		for(i = 0; i < ODD_RADIX; i++, tm2++) {
			tm2->d0 = wtinv_arr[icycle[i] >> L2_SZ_VD];
			tm2->d1 = wtinv_arr[jcycle[i] >> L2_SZ_VD];
		#ifdef USE_AVX
			tm2->d2 = wtinv_arr[kcycle[i] >> L2_SZ_VD];
			tm2->d3 = wtinv_arr[lcycle[i] >> L2_SZ_VD];
		#endif
		#ifdef USE_AVX512
			tm2->d4 = wtinv_arr[mcycle[i] >> L2_SZ_VD];
			tm2->d5 = wtinv_arr[ncycle[i] >> L2_SZ_VD];
			tm2->d6 = wtinv_arr[ocycle[i] >> L2_SZ_VD];
			tm2->d7 = wtinv_arr[pcycle[i] >> L2_SZ_VD];
		#endif
		}

		// Propagate the above inv-wts to the remaining threads - surrounding consts are unchanged:
		nbytes = ODD_RADIX*SZ_VD;
		tmp = half_arr + ODD_RADIX;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// USE_SSE2 ?
	}	// MODULUS_TYPE_FERMAT ?

#ifdef USE_PTHREAD
	for(ithread = 0; ithread < CY_THREADS; ++ithread) { tdat[ithread].iter = iter; }
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Carry-injection location for the shifted-residue -2 addend is only needed for full pass:
		if(full_pass) {
			tdat[0].target_idx = target_idx;
			tdat[0].target_set = target_set;
			tdat[0].target_cy  = target_cy;
		} else {
			tdat[0].target_idx = -1;
			tdat[0].target_set = 0;
			tdat[0].target_cy  = 0;
		}
		// Copy to the remaining threads:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			tdat[ithread].target_idx = tdat[0].target_idx;
			tdat[ithread].target_set = tdat[0].target_set;
			tdat[ithread].target_cy  = tdat[0].target_cy;
		}
	}
#endif
#ifdef USE_SSE2

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, SZ_VD);
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
		ASSERT(tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = 0.0;
		tdat[ithread].scale = scale;
		tdat[ithread].prp_mult = prp_mult;

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		ASSERT(tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		ASSERT(((tmp-1)->d0 == base[0] && (tmp-1)->d1 == baseinv[1] && (tmp-1)->d2 == wts_mult[1] && (tmp-1)->d3 == inv_mult[0]), "thread-local memcheck failed!");
	  #else
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #endif
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			/* No-Op */
		#elif defined(USE_AVX)
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
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
		if(full_pass) maxerr = 0.0;
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
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];
				tmp->d4 = _cy_r[l+4][ithread];
				tmp->d5 = _cy_r[l+5][ithread];
				tmp->d6 = _cy_r[l+6][ithread];
				tmp->d7 = _cy_r[l+7][ithread];
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
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
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				tmp->d0 = _cy_r[l  ][ithread];		tm2->d0 = _cy_i[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];		tm2->d1 = _cy_i[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];		tm2->d2 = _cy_i[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];		tm2->d3 = _cy_i[l+3][ithread];
				tmp->d4 = _cy_r[l+4][ithread];		tm2->d4 = _cy_i[l+4][ithread];
				tmp->d5 = _cy_r[l+5][ithread];		tm2->d5 = _cy_i[l+5][ithread];
				tmp->d6 = _cy_r[l+6][ithread];		tm2->d6 = _cy_i[l+6][ithread];
				tmp->d7 = _cy_r[l+7][ithread];		tm2->d7 = _cy_i[l+7][ithread];
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
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
		#include "radix224_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
				_cy_r[l+2][ithread] = tmp->d2;
				_cy_r[l+3][ithread] = tmp->d3;
				_cy_r[l+4][ithread] = tmp->d4;
				_cy_r[l+5][ithread] = tmp->d5;
				_cy_r[l+6][ithread] = tmp->d6;
				_cy_r[l+7][ithread] = tmp->d7;
			}
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
				_cy_r[l+2][ithread] = tmp->d2;
				_cy_r[l+3][ithread] = tmp->d3;
			}
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
			}
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];
			}
		#endif
		}
		else
		{
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				_cy_r[l  ][ithread] = tmp->d0;		_cy_i[l  ][ithread] = tm2->d0;
				_cy_r[l+1][ithread] = tmp->d1;		_cy_i[l+1][ithread] = tm2->d1;
				_cy_r[l+2][ithread] = tmp->d2;		_cy_i[l+2][ithread] = tm2->d2;
				_cy_r[l+3][ithread] = tmp->d3;		_cy_i[l+3][ithread] = tm2->d3;
				_cy_r[l+4][ithread] = tmp->d4;		_cy_i[l+4][ithread] = tm2->d4;
				_cy_r[l+5][ithread] = tmp->d5;		_cy_i[l+5][ithread] = tm2->d5;
				_cy_r[l+6][ithread] = tmp->d6;		_cy_i[l+6][ithread] = tm2->d6;
				_cy_r[l+7][ithread] = tmp->d7;		_cy_i[l+7][ithread] = tm2->d7;
			}
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				_cy_r[l  ][ithread] = tmp->d0;		_cy_i[l  ][ithread] = tm2->d0;
				_cy_r[l+1][ithread] = tmp->d1;		_cy_i[l+1][ithread] = tm2->d1;
				_cy_r[l+2][ithread] = tmp->d2;		_cy_i[l+2][ithread] = tm2->d2;
				_cy_r[l+3][ithread] = tmp->d3;		_cy_i[l+3][ithread] = tm2->d3;
			}
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				_cy_r[l][ithread] = tmp->d0;		_cy_i[l][ithread] = tmp->d1;
			}
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];		_cy_i[l][ithread] = cy_i[l];
			}
		#endif
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(0x0 == cy224_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
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
	scale = prp_mult = 1;	// Don't scale by prp_mult on cleanup-pass, since that multiplier has already been applied to every residue word.

	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large:
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
		j_jhi = jhi_wrap_ferm;
	else
		j_jhi = jhi_wrap_mers;

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
		*fracmax = maxerr;
	}
	if(dtmp != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		mlucas_fprint(cbuf,INTERACT);
		err = ERR_CARRY;
		return(err);
	}
	return(0);
}

/****************/

void radix224_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-224 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int l,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0, first_entry=TRUE;
	static int t_offsets[32], dif_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13 elts:
	static int dif_p20_cperms[26], dif_p20_lo_offset[32], dif_phi[7];
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;

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
		NDIVR >>= 4;			pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 7 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

		dif_phi[0] =   0;
		dif_phi[1] = pc0;
		dif_phi[2] = pa0;
		dif_phi[3] = p80;
		dif_phi[4] = p60;
		dif_phi[5] = p40;
		dif_phi[6] = p20;

		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13
		l = 0;
		dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40; dif_p20_cperms[l++] = p20; dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30; dif_p20_cperms[l++] = p10; dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 3) + 1);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 1);
		dif_p20_lo_offset[l++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 3) + 2);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 2);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 2);
		dif_p20_lo_offset[l++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 3) + 3);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 3);
		dif_p20_lo_offset[l++] = ((pe << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 4);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 4);
		dif_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 3) + 5);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 5);
		dif_p20_lo_offset[l++] = ((pf << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 3) + 6);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 6);
		dif_p20_lo_offset[l++] = ((pc << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 3) + 0);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 0);

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,2,3,4,5,6,7,9,8,b,a,d,c,f,e + p00],[2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p2+p10;	// <*** same rcol offsets as Set 6 lcol
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p3+p10;	// <*** Unique lcol offsets
		dif_offsets[0x02] = p2;		dif_offsets[0x12] = p1+p10;
		dif_offsets[0x03] = p3;		dif_offsets[0x13] =    p10;
		dif_offsets[0x04] = p4;		dif_offsets[0x14] = p6+p10;
		dif_offsets[0x05] = p5;		dif_offsets[0x15] = p7+p10;
		dif_offsets[0x06] = p6;		dif_offsets[0x16] = p5+p10;
		dif_offsets[0x07] = p7;		dif_offsets[0x17] = p4+p10;
		dif_offsets[0x08] = p9;		dif_offsets[0x18] = pb+p10;
		dif_offsets[0x09] = p8;		dif_offsets[0x19] = pa+p10;
		dif_offsets[0x0a] = pb;		dif_offsets[0x1a] = p8+p10;
		dif_offsets[0x0b] = pa;		dif_offsets[0x1b] = p9+p10;
		dif_offsets[0x0c] = pd;		dif_offsets[0x1c] = pf+p10;
		dif_offsets[0x0d] = pc;		dif_offsets[0x1d] = pe+p10;
		dif_offsets[0x0e] = pf;		dif_offsets[0x1e] = pc+p10;
		dif_offsets[0x0f] = pe;		dif_offsets[0x1f] = pd+p10;
		// Set 1: [d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pc0],[f,e,c,d,8,9,a,b,1,0,3,2,5,4,7,6 + pd0] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = pd;		dif_offsets[l+0x10] = pf+p10;	// <*** Unique rcol offsets
		dif_offsets[l+0x01] = pc;		dif_offsets[l+0x11] = pe+p10;	// <*** same lcol offsets as Set 2 rcol
		dif_offsets[l+0x02] = pf;		dif_offsets[l+0x12] = pc+p10;
		dif_offsets[l+0x03] = pe;		dif_offsets[l+0x13] = pd+p10;
		dif_offsets[l+0x04] = pb;		dif_offsets[l+0x14] = p8+p10;
		dif_offsets[l+0x05] = pa;		dif_offsets[l+0x15] = p9+p10;
		dif_offsets[l+0x06] = p8;		dif_offsets[l+0x16] = pa+p10;
		dif_offsets[l+0x07] = p9;		dif_offsets[l+0x17] = pb+p10;
		dif_offsets[l+0x08] = p2;		dif_offsets[l+0x18] = p1+p10;
		dif_offsets[l+0x09] = p3;		dif_offsets[l+0x19] =    p10;
		dif_offsets[l+0x0a] = p1;		dif_offsets[l+0x1a] = p3+p10;
		dif_offsets[l+0x0b] =  0;		dif_offsets[l+0x1b] = p2+p10;
		dif_offsets[l+0x0c] = p6;		dif_offsets[l+0x1c] = p5+p10;
		dif_offsets[l+0x0d] = p7;		dif_offsets[l+0x1d] = p4+p10;
		dif_offsets[l+0x0e] = p5;		dif_offsets[l+0x1e] = p7+p10;
		dif_offsets[l+0x0f] = p4;		dif_offsets[l+0x1f] = p6+p10;
		// Set 2: [6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + pb0],[d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pa0] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = p6+p10;		dif_offsets[l+0x10] = pd;	// <*** same rcol offsets as Set 1 lcol
		dif_offsets[l+0x01] = p7+p10;		dif_offsets[l+0x11] = pc;	// <*** same lcol offsets as Set 3 rcol
		dif_offsets[l+0x02] = p5+p10;		dif_offsets[l+0x12] = pf;
		dif_offsets[l+0x03] = p4+p10;		dif_offsets[l+0x13] = pe;
		dif_offsets[l+0x04] = p1+p10;		dif_offsets[l+0x14] = pb;
		dif_offsets[l+0x05] =    p10;		dif_offsets[l+0x15] = pa;
		dif_offsets[l+0x06] = p3+p10;		dif_offsets[l+0x16] = p8;
		dif_offsets[l+0x07] = p2+p10;		dif_offsets[l+0x17] = p9;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = p2;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = p3;
		dif_offsets[l+0x0a] = pc+p10;		dif_offsets[l+0x1a] = p1;
		dif_offsets[l+0x0b] = pd+p10;		dif_offsets[l+0x1b] =  0;
		dif_offsets[l+0x0c] = p8+p10;		dif_offsets[l+0x1c] = p6;
		dif_offsets[l+0x0d] = p9+p10;		dif_offsets[l+0x1d] = p7;
		dif_offsets[l+0x0e] = pa+p10;		dif_offsets[l+0x1e] = p5;
		dif_offsets[l+0x0f] = pb+p10;		dif_offsets[l+0x1f] = p4;
		// Set 3: [4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p80],[6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + p90] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = p4;		dif_offsets[l+0x10] = p6+p10;	// <*** same rcol offsets as Set 2 lcol
		dif_offsets[l+0x01] = p5;		dif_offsets[l+0x11] = p7+p10;	// <*** same lcol offsets as Set 4 rcol
		dif_offsets[l+0x02] = p6;		dif_offsets[l+0x12] = p5+p10;
		dif_offsets[l+0x03] = p7;		dif_offsets[l+0x13] = p4+p10;
		dif_offsets[l+0x04] = p2;		dif_offsets[l+0x14] = p1+p10;
		dif_offsets[l+0x05] = p3;		dif_offsets[l+0x15] =    p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = p3+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = p2+p10;
		dif_offsets[l+0x08] = pd;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = pc;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = pf;		dif_offsets[l+0x1a] = pc+p10;
		dif_offsets[l+0x0b] = pe;		dif_offsets[l+0x1b] = pd+p10;
		dif_offsets[l+0x0c] = pb;		dif_offsets[l+0x1c] = p8+p10;
		dif_offsets[l+0x0d] = pa;		dif_offsets[l+0x1d] = p9+p10;
		dif_offsets[l+0x0e] = p8;		dif_offsets[l+0x1e] = pa+p10;
		dif_offsets[l+0x0f] = p9;		dif_offsets[l+0x1f] = pb+p10;
		// Set 4: [b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p70],[4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p60] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = pb+p10;		dif_offsets[l+0x10] = p4;	// <*** same rcol offsets as Set 3 lcol
		dif_offsets[l+0x01] = pa+p10;		dif_offsets[l+0x11] = p5;	// <*** same lcol offsets as Set 5 rcol
		dif_offsets[l+0x02] = p8+p10;		dif_offsets[l+0x12] = p6;
		dif_offsets[l+0x03] = p9+p10;		dif_offsets[l+0x13] = p7;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = p2;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = p3;
		dif_offsets[l+0x06] = pc+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = pd+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = p6+p10;		dif_offsets[l+0x18] = pd;
		dif_offsets[l+0x09] = p7+p10;		dif_offsets[l+0x19] = pc;
		dif_offsets[l+0x0a] = p5+p10;		dif_offsets[l+0x1a] = pf;
		dif_offsets[l+0x0b] = p4+p10;		dif_offsets[l+0x1b] = pe;
		dif_offsets[l+0x0c] = p1+p10;		dif_offsets[l+0x1c] = pb;
		dif_offsets[l+0x0d] =    p10;		dif_offsets[l+0x1d] = pa;
		dif_offsets[l+0x0e] = p3+p10;		dif_offsets[l+0x1e] = p8;
		dif_offsets[l+0x0f] = p2+p10;		dif_offsets[l+0x1f] = p9;
		// Set 5: [9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p40],[b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p50] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = p9;		dif_offsets[l+0x10] = pb+p10;	// <*** same rcol offsets as Set 4 lcol
		dif_offsets[l+0x01] = p8;		dif_offsets[l+0x11] = pa+p10;	// <*** same lcol offsets as Set 6 rcol
		dif_offsets[l+0x02] = pb;		dif_offsets[l+0x12] = p8+p10;
		dif_offsets[l+0x03] = pa;		dif_offsets[l+0x13] = p9+p10;
		dif_offsets[l+0x04] = pd;		dif_offsets[l+0x14] = pf+p10;
		dif_offsets[l+0x05] = pc;		dif_offsets[l+0x15] = pe+p10;
		dif_offsets[l+0x06] = pf;		dif_offsets[l+0x16] = pc+p10;
		dif_offsets[l+0x07] = pe;		dif_offsets[l+0x17] = pd+p10;
		dif_offsets[l+0x08] = p4;		dif_offsets[l+0x18] = p6+p10;
		dif_offsets[l+0x09] = p5;		dif_offsets[l+0x19] = p7+p10;
		dif_offsets[l+0x0a] = p6;		dif_offsets[l+0x1a] = p5+p10;
		dif_offsets[l+0x0b] = p7;		dif_offsets[l+0x1b] = p4+p10;
		dif_offsets[l+0x0c] = p2;		dif_offsets[l+0x1c] = p1+p10;
		dif_offsets[l+0x0d] = p3;		dif_offsets[l+0x1d] =    p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = p3+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = p2+p10;
		// Set 6: [2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p30],[9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p20] (mod p20):
		l +=32;
		dif_offsets[l+0x00] = p2+p10;		dif_offsets[l+0x10] = p9;	// <*** same rcol offsets as Set 5 lcol
		dif_offsets[l+0x01] = p3+p10;		dif_offsets[l+0x11] = p8;	// <*** same lcol offsets as Set 0 rcol
		dif_offsets[l+0x02] = p1+p10;		dif_offsets[l+0x12] = pb;
		dif_offsets[l+0x03] =    p10;		dif_offsets[l+0x13] = pa;
		dif_offsets[l+0x04] = p6+p10;		dif_offsets[l+0x14] = pd;
		dif_offsets[l+0x05] = p7+p10;		dif_offsets[l+0x15] = pc;
		dif_offsets[l+0x06] = p5+p10;		dif_offsets[l+0x16] = pf;
		dif_offsets[l+0x07] = p4+p10;		dif_offsets[l+0x17] = pe;
		dif_offsets[l+0x08] = pb+p10;		dif_offsets[l+0x18] = p4;
		dif_offsets[l+0x09] = pa+p10;		dif_offsets[l+0x19] = p5;
		dif_offsets[l+0x0a] = p8+p10;		dif_offsets[l+0x1a] = p6;
		dif_offsets[l+0x0b] = p9+p10;		dif_offsets[l+0x1b] = p7;
		dif_offsets[l+0x0c] = pf+p10;		dif_offsets[l+0x1c] = p2;
		dif_offsets[l+0x0d] = pe+p10;		dif_offsets[l+0x1d] = p3;
		dif_offsets[l+0x0e] = pc+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = pd+p10;		dif_offsets[l+0x1f] =  0;
	}

/*...The radix-224 pass is here.	*/
	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 += ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 32 radix-7 transforms. */
	/*
	Twiddleless version arranges 32 sets of radix-7 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 224 = 0xe0) 32 (= 0x20) horizontally and 7 vertically. Indexing in hex for clarity and using
	[evn|odd]0-6 notation in the rightmost column to flag reusable 7-perms [in fact simple circular (0-6)-element
	shifts of a basic pattern] of (0-6)*p20 and p10 + (0-6)*p20:

		00,c0,a0,80,60,40,20		00,c0,a0,80,60,40,20 + p0		[evn0] + p0
		d9,b9,99,79,59,39,19		d0,b0,90,70,50,30,10 + p9		[odd0] + p9
		d2,b2,92,72,52,32,12		d0,b0,90,70,50,30,10 + p2		[odd0] + p2
		cb,ab,8b,6b,4b,2b,0b		c0,a0,80,60,40,20,00 + pb		[evn1] + pb
		c4,a4,84,64,44,24,04		c0,a0,80,60,40,20,00 + p4		[evn1] + p4
		bd,9d,7d,5d,3d,1d,dd		b0,90,70,50,30,10,d0 + pd		[odd1] + pd
		b6,96,76,56,36,16,d6		b0,90,70,50,30,10,d0 + p6		[odd1] + p6
		af,8f,6f,4f,2f,0f,cf		a0,80,60,40,20,00,c0 + pf		[evn2] + pf
		a8,88,68,48,28,08,c8		a0,80,60,40,20,00,c0 + p8		[evn2] + p8
		a1,81,61,41,21,01,c1		a0,80,60,40,20,00,c0 + p1		[evn2] + p1
		9a,7a,5a,3a,1a,da,ba		90,70,50,30,10,d0,b0 + pa		[odd2] + pa
		93,73,53,33,13,d3,b3		90,70,50,30,10,d0,b0 + p3		[odd2] + p3
		8c,6c,4c,2c,0c,cc,ac		80,60,40,20,00,c0,a0 + pc		[evn3] + pc
		85,65,45,25,05,c5,a5		80,60,40,20,00,c0,a0 + p5		[evn3] + p5
		7e,5e,3e,1e,de,be,9e		70,50,30,10,d0,b0,90 + pe		[odd3] + pe
		77,57,37,17,d7,b7,97		70,50,30,10,d0,b0,90 + p7		[odd3] + p7
		70,50,30,10,d0,b0,90	=	70,50,30,10,d0,b0,90 + p0	=	[odd3] + p0
		69,49,29,09,c9,a9,89		60,40,20,00,c0,a0,80 + p9		[evn4] + p9
		62,42,22,02,c2,a2,82		60,40,20,00,c0,a0,80 + p2		[evn4] + p2
		5b,3b,1b,db,bb,9b,7b		50,30,10,d0,b0,90,70 + pb		[odd4] + pb
		54,34,14,d4,b4,94,74		50,30,10,d0,b0,90,70 + p4		[odd4] + p4
		4d,2d,0d,cd,ad,8d,6d		40,20,00,c0,a0,80,60 + pd		[evn5] + pd
		46,26,06,c6,a6,86,66		40,20,00,c0,a0,80,60 + p6		[evn5] + p6
		3f,1f,df,bf,9f,7f,5f		30,10,d0,b0,90,70,50 + pf		[odd5] + pf
		38,18,d8,b8,98,78,58		30,10,d0,b0,90,70,50 + p8		[odd5] + p8
		31,11,d1,b1,91,71,51		30,10,d0,b0,90,70,50 + p1		[odd5] + p1
		2a,0a,ca,aa,8a,6a,4a		20,00,c0,a0,80,60,40 + pa		[evn6] + pa
		23,03,c3,a3,83,63,43		20,00,c0,a0,80,60,40 + p3		[evn6] + p3
		1c,dc,bc,9c,7c,5c,3c		10,d0,b0,90,70,50,30 + pc		[odd6] + pc
		15,d5,b5,95,75,55,35		10,d0,b0,90,70,50,30 + p5		[odd6] + p5
		0e,ce,ae,8e,6e,4e,2e		00,c0,a0,80,60,40,20 + pe		[evn0] + pe
		07,c7,a7,87,67,47,27		00,c0,a0,80,60,40,20 + p7		[evn0] + p7
	*/
		tptr = t;
		for(l = 0; l < 32; l++) {
			int k = dif_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of k; the "which length-13 half of the dif_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 13)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/13)
						+ (k & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4], k5 = dif_p20_cperms[ic+5], k6 = dif_p20_cperms[ic+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			k = (k & 0x7fffffff) >> 3;
			jt = j1+k; jp = j2+k;
			RADIX_07_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im,
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++;
		}

	/*...and now do 7 radix-32 transforms;
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the radix-32 DFTs to the required ordering, which in terms of our p-offsets is

		00,01,02,03,04,05,06,07,09,08,0b,0a,0d,0c,0f,0e,12,13,11,10,16,17,15,14,1b,1a,18,19,1f,1e,1c,1d
		cd,cc,cf,ce,cb,ca,c8,c9,c2,c3,c1,c0,c6,c7,c5,c4,df,de,dc,dd,d8,d9,da,db,d1,d0,d3,d2,d5,d4,d7,d6
		b6,b7,b5,b4,b1,b0,b3,b2,bf,be,bc,bd,b8,b9,ba,bb,ad,ac,af,ae,ab,aa,a8,a9,a2,a3,a1,a0,a6,a7,a5,a4
		84,85,86,87,82,83,81,80,8d,8c,8f,8e,8b,8a,88,89,96,97,95,94,91,90,93,92,9f,9e,9c,9d,98,99,9a,9b
		7b,7a,78,79,7f,7e,7c,7d,76,77,75,74,71,70,73,72,64,65,66,67,62,63,61,60,6d,6c,6f,6e,6b,6a,68,69
		49,48,4b,4a,4d,4c,4f,4e,44,45,46,47,42,43,41,40,5b,5a,58,59,5f,5e,5c,5d,56,57,55,54,51,50,53,52
		32,33,31,30,36,37,35,34,3b,3a,38,39,3f,3e,3c,3d,29,28,2b,2a,2d,2c,2f,2e,24,25,26,27,22,23,21,20

		[0,1,2,3,4,5,6,7,9,8,b,a,d,c,f,e + p00],[2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p10]
		[d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pc0],[f,e,c,d,8,9,a,b,1,0,3,2,5,4,7,6 + pd0]
		[6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + pb0],[d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pa0]
	=	[4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p80],[6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + p90]
		[b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p70],[4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p60]
		[9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p40],[b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p50]
		[2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p30],[9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p20]
	*/
		tptr = t;
	// In scalar-double mode this makes no objcode size difference because the radix-32 DFTs
	// are separate functions defined in dft_macro.c,but use same loop template in SIMD mode
		for(l = 0; l < 7; l++) {
			jt = j1+dif_phi[l]; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+(l<<5),RE_IM_STRIDE);	tptr += 32;
		}
	}
}

/***************/

void radix224_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-224 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int l,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0, first_entry=TRUE;
	static int t_offsets[32], dit_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13 elts:
	static int dit_p20_cperms[26], dit_p20_lo_offset[32], dit_phi[7];
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;

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
		NDIVR >>= 4;			pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-32 outputs.
		// t_offsets w.r.to: t-array, same for all 7 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

		dit_phi[0] =   0;
		dit_phi[1] = p60;
		dit_phi[2] = pc0;
		dit_phi[3] = p40;
		dit_phi[4] = pa0;
		dit_phi[5] = p20;
		dit_phi[6] = p80;

		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13
		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 3) + 6);
		dit_p20_lo_offset[l++] = ((pe << 3) + 0);
		dit_p20_lo_offset[l++] = ((pd << 3) + 1);
		dit_p20_lo_offset[l++] = ((pc << 3) + 2);
		dit_p20_lo_offset[l++] = ((pb << 3) + 3);
		dit_p20_lo_offset[l++] = ((pa << 3) + 4);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 6);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 0);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 4);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 6);

	// dit_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		l = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p70],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p60] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = pb+p10;		dit_offsets[l+0x10] = pb;
		dit_offsets[l+0x01] = pa+p10;		dit_offsets[l+0x11] = pa;
		dit_offsets[l+0x02] = p9+p10;		dit_offsets[l+0x12] = p9;
		dit_offsets[l+0x03] = p8+p10;		dit_offsets[l+0x13] = p8;
		dit_offsets[l+0x04] = pd+p10;		dit_offsets[l+0x14] = pd;
		dit_offsets[l+0x05] = pc+p10;		dit_offsets[l+0x15] = pc;
		dit_offsets[l+0x06] = pe+p10;		dit_offsets[l+0x16] = pe;
		dit_offsets[l+0x07] = pf+p10;		dit_offsets[l+0x17] = pf;
		dit_offsets[l+0x08] = p3+p10;		dit_offsets[l+0x18] = p3;
		dit_offsets[l+0x09] = p2+p10;		dit_offsets[l+0x19] = p2;
		dit_offsets[l+0x0a] = p1+p10;		dit_offsets[l+0x1a] = p1;
		dit_offsets[l+0x0b] =    p10;		dit_offsets[l+0x1b] =  0;
		dit_offsets[l+0x0c] = p5+p10;		dit_offsets[l+0x1c] = p5;
		dit_offsets[l+0x0d] = p4+p10;		dit_offsets[l+0x1d] = p4;
		dit_offsets[l+0x0e] = p6+p10;		dit_offsets[l+0x1e] = p6;
		dit_offsets[l+0x0f] = p7+p10;		dit_offsets[l+0x1f] = p7;
		// Set 2: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + pc0],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + pd0] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = pd;		dit_offsets[l+0x10] = p5+p10;
		dit_offsets[l+0x01] = pc;		dit_offsets[l+0x11] = p4+p10;
		dit_offsets[l+0x02] = pe;		dit_offsets[l+0x12] = p6+p10;
		dit_offsets[l+0x03] = pf;		dit_offsets[l+0x13] = p7+p10;
		dit_offsets[l+0x04] = p9;		dit_offsets[l+0x14] = p1+p10;
		dit_offsets[l+0x05] = p8;		dit_offsets[l+0x15] =    p10;
		dit_offsets[l+0x06] = pa;		dit_offsets[l+0x16] = p2+p10;
		dit_offsets[l+0x07] = pb;		dit_offsets[l+0x17] = p3+p10;
		dit_offsets[l+0x08] = p5;		dit_offsets[l+0x18] = p9+p10;
		dit_offsets[l+0x09] = p4;		dit_offsets[l+0x19] = p8+p10;
		dit_offsets[l+0x0a] = p6;		dit_offsets[l+0x1a] = pa+p10;
		dit_offsets[l+0x0b] = p7;		dit_offsets[l+0x1b] = pb+p10;
		dit_offsets[l+0x0c] = p1;		dit_offsets[l+0x1c] = pe+p10;
		dit_offsets[l+0x0d] =  0;		dit_offsets[l+0x1d] = pf+p10;
		dit_offsets[l+0x0e] = p2;		dit_offsets[l+0x1e] = pc+p10;
		dit_offsets[l+0x0f] = p3;		dit_offsets[l+0x1f] = pd+p10;
		// Set 3: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p40],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p50] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = p9;		dit_offsets[l+0x10] = p1+p10;
		dit_offsets[l+0x01] = p8;		dit_offsets[l+0x11] =    p10;
		dit_offsets[l+0x02] = pa;		dit_offsets[l+0x12] = p2+p10;
		dit_offsets[l+0x03] = pb;		dit_offsets[l+0x13] = p3+p10;
		dit_offsets[l+0x04] = pe;		dit_offsets[l+0x14] = p6+p10;
		dit_offsets[l+0x05] = pf;		dit_offsets[l+0x15] = p7+p10;
		dit_offsets[l+0x06] = pc;		dit_offsets[l+0x16] = p4+p10;
		dit_offsets[l+0x07] = pd;		dit_offsets[l+0x17] = p5+p10;
		dit_offsets[l+0x08] = p1;		dit_offsets[l+0x18] = pe+p10;
		dit_offsets[l+0x09] =  0;		dit_offsets[l+0x19] = pf+p10;
		dit_offsets[l+0x0a] = p2;		dit_offsets[l+0x1a] = pc+p10;
		dit_offsets[l+0x0b] = p3;		dit_offsets[l+0x1b] = pd+p10;
		dit_offsets[l+0x0c] = p6;		dit_offsets[l+0x1c] = pa+p10;
		dit_offsets[l+0x0d] = p7;		dit_offsets[l+0x1d] = pb+p10;
		dit_offsets[l+0x0e] = p4;		dit_offsets[l+0x1e] = p8+p10;
		dit_offsets[l+0x0f] = p5;		dit_offsets[l+0x1f] = p9+p10;
		// Set 4: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pb0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pa0] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = p6+p10;		dit_offsets[l+0x10] = p6;
		dit_offsets[l+0x01] = p7+p10;		dit_offsets[l+0x11] = p7;
		dit_offsets[l+0x02] = p4+p10;		dit_offsets[l+0x12] = p4;
		dit_offsets[l+0x03] = p5+p10;		dit_offsets[l+0x13] = p5;
		dit_offsets[l+0x04] = p2+p10;		dit_offsets[l+0x14] = p2;
		dit_offsets[l+0x05] = p3+p10;		dit_offsets[l+0x15] = p3;
		dit_offsets[l+0x06] =    p10;		dit_offsets[l+0x16] =  0;
		dit_offsets[l+0x07] = p1+p10;		dit_offsets[l+0x17] = p1;
		dit_offsets[l+0x08] = pa+p10;		dit_offsets[l+0x18] = pa;
		dit_offsets[l+0x09] = pb+p10;		dit_offsets[l+0x19] = pb;
		dit_offsets[l+0x0a] = p8+p10;		dit_offsets[l+0x1a] = p8;
		dit_offsets[l+0x0b] = p9+p10;		dit_offsets[l+0x1b] = p9;
		dit_offsets[l+0x0c] = pc+p10;		dit_offsets[l+0x1c] = pc;
		dit_offsets[l+0x0d] = pd+p10;		dit_offsets[l+0x1d] = pd;
		dit_offsets[l+0x0e] = pf+p10;		dit_offsets[l+0x1e] = pf;
		dit_offsets[l+0x0f] = pe+p10;		dit_offsets[l+0x1f] = pe;
		// Set 5: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p20] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = p2+p10;		dit_offsets[l+0x10] = p2;
		dit_offsets[l+0x01] = p3+p10;		dit_offsets[l+0x11] = p3;
		dit_offsets[l+0x02] =    p10;		dit_offsets[l+0x12] =  0;
		dit_offsets[l+0x03] = p1+p10;		dit_offsets[l+0x13] = p1;
		dit_offsets[l+0x04] = p4+p10;		dit_offsets[l+0x14] = p4;
		dit_offsets[l+0x05] = p5+p10;		dit_offsets[l+0x15] = p5;
		dit_offsets[l+0x06] = p7+p10;		dit_offsets[l+0x16] = p7;
		dit_offsets[l+0x07] = p6+p10;		dit_offsets[l+0x17] = p6;
		dit_offsets[l+0x08] = pc+p10;		dit_offsets[l+0x18] = pc;
		dit_offsets[l+0x09] = pd+p10;		dit_offsets[l+0x19] = pd;
		dit_offsets[l+0x0a] = pf+p10;		dit_offsets[l+0x1a] = pf;
		dit_offsets[l+0x0b] = pe+p10;		dit_offsets[l+0x1b] = pe;
		dit_offsets[l+0x0c] = p8+p10;		dit_offsets[l+0x1c] = p8;
		dit_offsets[l+0x0d] = p9+p10;		dit_offsets[l+0x1d] = p9;
		dit_offsets[l+0x0e] = pb+p10;		dit_offsets[l+0x1e] = pb;
		dit_offsets[l+0x0f] = pa+p10;		dit_offsets[l+0x1f] = pa;
		// Set 6: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p80],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p90] (mod p20):
		l +=32;
		dit_offsets[l+0x00] = p4;		dit_offsets[l+0x10] = p8+p10;
		dit_offsets[l+0x01] = p5;		dit_offsets[l+0x11] = p9+p10;
		dit_offsets[l+0x02] = p7;		dit_offsets[l+0x12] = pb+p10;
		dit_offsets[l+0x03] = p6;		dit_offsets[l+0x13] = pa+p10;
		dit_offsets[l+0x04] =  0;		dit_offsets[l+0x14] = pf+p10;
		dit_offsets[l+0x05] = p1;		dit_offsets[l+0x15] = pe+p10;
		dit_offsets[l+0x06] = p3;		dit_offsets[l+0x16] = pd+p10;
		dit_offsets[l+0x07] = p2;		dit_offsets[l+0x17] = pc+p10;
		dit_offsets[l+0x08] = p8;		dit_offsets[l+0x18] =    p10;
		dit_offsets[l+0x09] = p9;		dit_offsets[l+0x19] = p1+p10;
		dit_offsets[l+0x0a] = pb;		dit_offsets[l+0x1a] = p3+p10;
		dit_offsets[l+0x0b] = pa;		dit_offsets[l+0x1b] = p2+p10;
		dit_offsets[l+0x0c] = pf;		dit_offsets[l+0x1c] = p7+p10;
		dit_offsets[l+0x0d] = pe;		dit_offsets[l+0x1d] = p6+p10;
		dit_offsets[l+0x0e] = pd;		dit_offsets[l+0x1e] = p5+p10;
		dit_offsets[l+0x0f] = pc;		dit_offsets[l+0x1f] = p4+p10;
	}

/*...The radix-224 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX512
		j1 = (j & mask03) + br16[j&15];
	#elif defined(USE_AVX)
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 += ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so using output of test_fft_radix(),
	store the 7 index-offset 32-tets going into the radix-32 DFTs in the dit_offset array.

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially):
	Combined DIT input-scramble array:

		[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p70],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p60]
		[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + pc0],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + pd0]
		[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p40],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p50]
		[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pb0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pa0]
		[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p20]
		[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p80],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p90]
	*/
	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 7 radix-32 transforms...*/
		tptr = t;
	// In scalar-double mode this makes no objcode size difference because the radix-32 DFTs
	// are separate functions defined in dft_macro.c,but use same loop template in SIMD mode
		for(l = 0; l < 7; l++) {
			jt = j1+dit_phi[l]; RADIX_32_DIT((a+jt),dit_offsets+(l<<5),RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		}

	/*...and now do 32 radix-7 transforms, with the columns of t*[r,i] output pairs in the above 7x radix-32 set now acting as input rows.
	Since our first-look oindex ordering was +p0x[0,20,40,60,80,a0,c0] for each radix-7 and incrementing += p1 between those DFTs,
	arrange resulting mismatched-data-sorted index permutation into 7 vertical 32-entry columns to get needed oindex patterns.
	Indexing in hex for clarity and using [evn|odd]0-6 notation in the rightmost column to flag reusable 7-perms
	[in fact simple circular (0-6)-element shifts of a basic pattern] of (0-6)*p20 and p10 + (0-6)*p20:

		00,40,80,c0,20,60,a0		00,40,80,c0,20,60,a0 + p0		[evn0] + p0
		3f,7f,bf,1f,5f,9f,df		30,70,b0,10,50,90,d0 + pf		[odd0] + pf
		7e,be,1e,5e,9e,de,3e		70,b0,10,50,90,d0,30 + pe		[odd1] + pe
		bd,1d,5d,9d,dd,3d,7d		b0,10,50,90,d0,30,70 + pd		[odd2] + pd
		1c,5c,9c,dc,3c,7c,bc		10,50,90,d0,30,70,b0 + pc		[odd3] + pc
		5b,9b,db,3b,7b,bb,1b		50,90,d0,30,70,b0,10 + pb		[odd4] + pb
		9a,da,3a,7a,ba,1a,5a		90,d0,30,70,b0,10,50 + pa		[odd5] + pa
		d9,39,79,b9,19,59,99		d0,30,70,b0,10,50,90 + p9		[odd6] + p9
		38,78,b8,18,58,98,d8		30,70,b0,10,50,90,d0 + p8		[odd0] + p8
		77,b7,17,57,97,d7,37		70,b0,10,50,90,d0,30 + p7		[odd1] + p7
		b6,16,56,96,d6,36,76		b0,10,50,90,d0,30,70 + p6		[odd2] + p6
		15,55,95,d5,35,75,b5		10,50,90,d0,30,70,b0 + p5		[odd3] + p5
		54,94,d4,34,74,b4,14		50,90,d0,30,70,b0,10 + p4		[odd4] + p4
		93,d3,33,73,b3,13,53		90,d0,30,70,b0,10,50 + p3		[odd5] + p3
		d2,32,72,b2,12,52,92		d0,30,70,b0,10,50,90 + p2		[odd6] + p2
		31,71,b1,11,51,91,d1		30,70,b0,10,50,90,d0 + p1		[odd0] + p1
		70,b0,10,50,90,d0,30	=	70,b0,10,50,90,d0,30 + p0	=	[odd1] + p0
		af,0f,4f,8f,cf,2f,6f		a0,00,40,80,c0,20,60 + pf		[evn6] + pf
		0e,4e,8e,ce,2e,6e,ae		00,40,80,c0,20,60,a0 + pe		[evn0] + pe
		4d,8d,cd,2d,6d,ad,0d		40,80,c0,20,60,a0,00 + pd		[evn1] + pd
		8c,cc,2c,6c,ac,0c,4c		80,c0,20,60,a0,00,40 + pc		[evn2] + pc
		cb,2b,6b,ab,0b,4b,8b		c0,20,60,a0,00,40,80 + pb		[evn3] + pb
		2a,6a,aa,0a,4a,8a,ca		20,60,a0,00,40,80,c0 + pa		[evn4] + pa
		69,a9,09,49,89,c9,29		60,a0,00,40,80,c0,20 + p9		[evn5] + p9
		a8,08,48,88,c8,28,68		a0,00,40,80,c0,20,60 + p8		[evn6] + p8
		07,47,87,c7,27,67,a7		00,40,80,c0,20,60,a0 + p7		[evn0] + p7
		46,86,c6,26,66,a6,06		40,80,c0,20,60,a0,00 + p6		[evn1] + p6
		85,c5,25,65,a5,05,45		80,c0,20,60,a0,00,40 + p5		[evn2] + p5
		c4,24,64,a4,04,44,84		c0,20,60,a0,00,40,80 + p4		[evn3] + p4
		23,63,a3,03,43,83,c3		20,60,a0,00,40,80,c0 + p3		[evn4] + p3
		62,a2,02,42,82,c2,22		60,a0,00,40,80,c0,20 + p2		[evn5] + p2
		a1,01,41,81,c1,21,61		a0,00,40,80,c0,20,60 + p1		[evn6] + p1
	*/
		tptr = t;
		for(l = 0; l < 32; l++) {
			int k = dit_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of k; the "which length-13 half of the dit_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 13)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/13)
						+ (k & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4], k5 = dit_p20_cperms[ic+5], k6 = dit_p20_cperms[ic+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			k = (k & 0x7fffffff) >> 3;
			jt = j1+k; jp = j2+k;
			RADIX_07_DFT(
				tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im,
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy224_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		const char func[] = "cy224_process_chunk";
	  // LO_ADD = 1 in masterdefs.h means non-Nussbaumer radix-7, uses these sincos constants:
	  #if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/

	  #elif defined(USE_SSE2) || !LO_ADD	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:
		const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
						cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
						cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
						cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
					/* Switch the sign of ss3 in these: */
						sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
						sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
						sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
						sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
	  #endif
		double *addr,*addi;
		struct complex *tptr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
			,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0;
		int poff[RADIX>>2], t_offsets[32], dif_offsets[RADIX], dit_offsets[RADIX];
		// Need storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13 elts:
		int dif_p20_cperms[26], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
		int dit_p20_cperms[26], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
		int j,j1,j2,k,l,ntmp;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 14|28|56 for avx512,avx,sse, respectively.
		// But fixed-incr too restrictive here, so 'divide 14|28|56 into pieces' via increment-array whose elts sum to 14|28|56:
		const int *incr,*inc_arr = 0x0;	// Fermat-mod only uses inc_arr in AVX+ mode, init = 0 to quiet 'uninit' warnings.
	  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
		const int incr_long[] = {14}, incr_med[] = {7,7}, incr_short[] = {3,4,3,4}, incr_hiacc[] = {2,2,3,2,2,3};
	  #elif defined(USE_AVX)
		const int incr_long[] = {14,14}, incr_med[] = {7,7,7,7}, incr_short[] = {3,4,3,4,3,4,3,4}, incr_hiacc[] = {0};
	  #elif defined(USE_ARM_V8_SIMD)
		const int incr_long[] = {14,14,14,14}, incr_med[] = {7,7,7,7,7,7,7,7}, incr_short[] = {5,4,5,5,4,5,5,4,5,5,4,5}, incr_hiacc[] = {2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3,2,2,3};
	  #else
		const int incr_long[] = {14,14,14,14}, incr_med[] = {7,7,7,7,7,7,7,7}, incr_short[] = {3,4,3,4,3,4,3,4,3,4,3,4,3,4,3,4}, incr_hiacc[] = {0};
	  #endif
		// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
	  #ifdef USE_AVX512
		// For nfold > 2,  RADIX/8 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
//		int nfold = USE_SHORT_CY_CHAIN;
//		const int nexec_long[] = {28}, nexec_med[] = {14,14}, nexec_short[] = {7,7,7,7}, nexec_hiacc[] = {4,3,4,3,4,3,4,3};	// nfold = 0-3
		int nfold = USE_SHORT_CY_CHAIN + 1;
		const int nexec_long[] = {14,14}, nexec_med[] = {7,7,7,7}, nexec_short[] = {4,3,4,3,4,3,4,3}, nexec_hiacc[] = {2,2,2,1,2,2,2,1,2,2,2,1,2,2,2,1};	// nfold = 1-4
	  #elif defined(USE_AVX)
		// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
		int nfold = USE_SHORT_CY_CHAIN + 1;
		const int nexec_long[] = {28,28}, nexec_med[] = {14,14,14,14}, nexec_short[] = {7,7,7,7,7,7,7,7}, nexec_hiacc[] = {4,3,4,3,4,3,4,3,4,3,4,3,4,3,4,3};	// nfold = 1-4
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			if(USE_SHORT_CY_CHAIN == 0)
				inc_arr = incr_long;
			else if(USE_SHORT_CY_CHAIN == 1)
				inc_arr = incr_med;
			else if(USE_SHORT_CY_CHAIN == 2)
				inc_arr = incr_short;
			else
				inc_arr = incr_hiacc;
		} else {	// MODULUS_TYPE_FERMAT:
		#ifdef USE_AVX
			if(USE_SHORT_CY_CHAIN == 0)
				inc_arr = nexec_long;
			else if(USE_SHORT_CY_CHAIN == 1)
				inc_arr = nexec_med;
			else if(USE_SHORT_CY_CHAIN == 2)
				inc_arr = nexec_short;
			else
				inc_arr = nexec_hiacc;
		#endif
		}

	#ifdef USE_AVX512
		double t0,t1,t2,t3;
	  #ifdef CARRY_16_WAY
		struct uint32x16 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	  #else
		struct uint32x8  *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	  #endif
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#endif
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
		double rt,it;
		int k1,k2;

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2,	// utility ptrs
			*va0,*va1,*va2,*va3,*va4,*va5,*va6,
			*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6;
		int *itmp,*itm2;			// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *two,*one,*sqrt2,*isrt2,*cc1,*ss1,*cc2,*ss2,*cc3,*ss3,	// radix-32 DFT trig consts
			*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,			// radix-7 DFT trig consts
			*max_err, *sse2_rnd, *half_arr,
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
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

		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		int p0123[4];
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy_r = thread_arg->cy_r,*cy_i = thread_arg->cy_i, temp,frac;
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		double re,im,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;	// tmps for scalar-double DFT macros
		struct complex t[RADIX];
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
		int thr_id = thread_arg->tid;
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int target_idx = thread_arg->target_idx;
		int target_set = thread_arg->target_set;
		double target_cy  = thread_arg->target_cy;
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
		double prp_mult = thread_arg->prp_mult;
	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		double *wts_mult = thread_arg->wts_mult;	// Const Intra-block wts-multiplier...
		double *inv_mult = thread_arg->inv_mult;	// ...and 2*(its multiplicative inverse).
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
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
		NDIVR >>= 4;			pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
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

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 7 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

	/*** DIF indexing stuff: ***/

		dif_phi[0] =   0;		dit_phi[0] =   0;
		dif_phi[1] = pc0;		dit_phi[1] = p60;
		dif_phi[2] = pa0;		dit_phi[2] = pc0;
		dif_phi[3] = p80;		dit_phi[3] = p40;
		dif_phi[4] = p60;		dit_phi[4] = pa0;
		dif_phi[5] = p40;		dit_phi[5] = p20;
		dif_phi[6] = p20;		dit_phi[6] = p80;

		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0xc0<<1; dif_p20_cperms[l++] = 0xa0<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1; dif_p20_cperms[l++] = 0x20<<1; dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0xc0<<1; dif_p20_cperms[l++] = 0xa0<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0xd0<<1; dif_p20_cperms[l++] = 0xb0<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1; dif_p20_cperms[l++] = 0x10<<1; dif_p20_cperms[l++] = 0xd0<<1; dif_p20_cperms[l++] = 0xb0<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:

		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 5);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 5);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 6);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 6);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 0);

	   #else

		l = 0;
		dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40; dif_p20_cperms[l++] = p20; dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = pc0; dif_p20_cperms[l++] = pa0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30; dif_p20_cperms[l++] = p10; dif_p20_cperms[l++] = pd0; dif_p20_cperms[l++] = pb0; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 3) + 1);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 1);
		dif_p20_lo_offset[l++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 3) + 2);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 2);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 2);
		dif_p20_lo_offset[l++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 3) + 3);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 3);
		dif_p20_lo_offset[l++] = ((pe << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 4);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 4);
		dif_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 3) + 5);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 5);
		dif_p20_lo_offset[l++] = ((pf << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 5) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 3) + 6);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 6);
		dif_p20_lo_offset[l++] = ((pc << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 6) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 3) + 0);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 0);

	   #endif

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each radix-32 DFT.
		// Set 0: [0,1,2,3,4,5,6,7,9,8,b,a,d,c,f,e + p00],[2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p2+p10;	// <*** same rcol offsets as Set 6 lcol
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p3+p10;	// <*** Unique lcol offsets
		dif_offsets[0x02] = p2;		dif_offsets[0x12] = p1+p10;
		dif_offsets[0x03] = p3;		dif_offsets[0x13] =    p10;
		dif_offsets[0x04] = p4;		dif_offsets[0x14] = p6+p10;
		dif_offsets[0x05] = p5;		dif_offsets[0x15] = p7+p10;
		dif_offsets[0x06] = p6;		dif_offsets[0x16] = p5+p10;
		dif_offsets[0x07] = p7;		dif_offsets[0x17] = p4+p10;
		dif_offsets[0x08] = p9;		dif_offsets[0x18] = pb+p10;
		dif_offsets[0x09] = p8;		dif_offsets[0x19] = pa+p10;
		dif_offsets[0x0a] = pb;		dif_offsets[0x1a] = p8+p10;
		dif_offsets[0x0b] = pa;		dif_offsets[0x1b] = p9+p10;
		dif_offsets[0x0c] = pd;		dif_offsets[0x1c] = pf+p10;
		dif_offsets[0x0d] = pc;		dif_offsets[0x1d] = pe+p10;
		dif_offsets[0x0e] = pf;		dif_offsets[0x1e] = pc+p10;
		dif_offsets[0x0f] = pe;		dif_offsets[0x1f] = pd+p10;
		// Set 1: [d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pc0],[f,e,c,d,8,9,a,b,1,0,3,2,5,4,7,6 + pd0] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = pd;		dif_offsets[l+0x10] = pf+p10;	// <*** Unique rcol offsets
		dif_offsets[l+0x01] = pc;		dif_offsets[l+0x11] = pe+p10;	// <*** same lcol offsets as Set 2 rcol
		dif_offsets[l+0x02] = pf;		dif_offsets[l+0x12] = pc+p10;
		dif_offsets[l+0x03] = pe;		dif_offsets[l+0x13] = pd+p10;
		dif_offsets[l+0x04] = pb;		dif_offsets[l+0x14] = p8+p10;
		dif_offsets[l+0x05] = pa;		dif_offsets[l+0x15] = p9+p10;
		dif_offsets[l+0x06] = p8;		dif_offsets[l+0x16] = pa+p10;
		dif_offsets[l+0x07] = p9;		dif_offsets[l+0x17] = pb+p10;
		dif_offsets[l+0x08] = p2;		dif_offsets[l+0x18] = p1+p10;
		dif_offsets[l+0x09] = p3;		dif_offsets[l+0x19] =    p10;
		dif_offsets[l+0x0a] = p1;		dif_offsets[l+0x1a] = p3+p10;
		dif_offsets[l+0x0b] =  0;		dif_offsets[l+0x1b] = p2+p10;
		dif_offsets[l+0x0c] = p6;		dif_offsets[l+0x1c] = p5+p10;
		dif_offsets[l+0x0d] = p7;		dif_offsets[l+0x1d] = p4+p10;
		dif_offsets[l+0x0e] = p5;		dif_offsets[l+0x1e] = p7+p10;
		dif_offsets[l+0x0f] = p4;		dif_offsets[l+0x1f] = p6+p10;
		// Set 2: [6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + pb0],[d,c,f,e,b,a,8,9,2,3,1,0,6,7,5,4 + pa0] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p6+p10;		dif_offsets[l+0x10] = pd;	// <*** same rcol offsets as Set 1 lcol
		dif_offsets[l+0x01] = p7+p10;		dif_offsets[l+0x11] = pc;	// <*** same lcol offsets as Set 3 rcol
		dif_offsets[l+0x02] = p5+p10;		dif_offsets[l+0x12] = pf;
		dif_offsets[l+0x03] = p4+p10;		dif_offsets[l+0x13] = pe;
		dif_offsets[l+0x04] = p1+p10;		dif_offsets[l+0x14] = pb;
		dif_offsets[l+0x05] =    p10;		dif_offsets[l+0x15] = pa;
		dif_offsets[l+0x06] = p3+p10;		dif_offsets[l+0x16] = p8;
		dif_offsets[l+0x07] = p2+p10;		dif_offsets[l+0x17] = p9;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = p2;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = p3;
		dif_offsets[l+0x0a] = pc+p10;		dif_offsets[l+0x1a] = p1;
		dif_offsets[l+0x0b] = pd+p10;		dif_offsets[l+0x1b] =  0;
		dif_offsets[l+0x0c] = p8+p10;		dif_offsets[l+0x1c] = p6;
		dif_offsets[l+0x0d] = p9+p10;		dif_offsets[l+0x1d] = p7;
		dif_offsets[l+0x0e] = pa+p10;		dif_offsets[l+0x1e] = p5;
		dif_offsets[l+0x0f] = pb+p10;		dif_offsets[l+0x1f] = p4;
		// Set 3: [4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p80],[6,7,5,4,1,0,3,2,f,e,c,d,8,9,a,b + p90] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p4;		dif_offsets[l+0x10] = p6+p10;	// <*** same rcol offsets as Set 2 lcol
		dif_offsets[l+0x01] = p5;		dif_offsets[l+0x11] = p7+p10;	// <*** same lcol offsets as Set 4 rcol
		dif_offsets[l+0x02] = p6;		dif_offsets[l+0x12] = p5+p10;
		dif_offsets[l+0x03] = p7;		dif_offsets[l+0x13] = p4+p10;
		dif_offsets[l+0x04] = p2;		dif_offsets[l+0x14] = p1+p10;
		dif_offsets[l+0x05] = p3;		dif_offsets[l+0x15] =    p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = p3+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = p2+p10;
		dif_offsets[l+0x08] = pd;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = pc;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = pf;		dif_offsets[l+0x1a] = pc+p10;
		dif_offsets[l+0x0b] = pe;		dif_offsets[l+0x1b] = pd+p10;
		dif_offsets[l+0x0c] = pb;		dif_offsets[l+0x1c] = p8+p10;
		dif_offsets[l+0x0d] = pa;		dif_offsets[l+0x1d] = p9+p10;
		dif_offsets[l+0x0e] = p8;		dif_offsets[l+0x1e] = pa+p10;
		dif_offsets[l+0x0f] = p9;		dif_offsets[l+0x1f] = pb+p10;
		// Set 4: [b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p70],[4,5,6,7,2,3,1,0,d,c,f,e,b,a,8,9 + p60] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = pb+p10;		dif_offsets[l+0x10] = p4;	// <*** same rcol offsets as Set 3 lcol
		dif_offsets[l+0x01] = pa+p10;		dif_offsets[l+0x11] = p5;	// <*** same lcol offsets as Set 5 rcol
		dif_offsets[l+0x02] = p8+p10;		dif_offsets[l+0x12] = p6;
		dif_offsets[l+0x03] = p9+p10;		dif_offsets[l+0x13] = p7;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = p2;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = p3;
		dif_offsets[l+0x06] = pc+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = pd+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = p6+p10;		dif_offsets[l+0x18] = pd;
		dif_offsets[l+0x09] = p7+p10;		dif_offsets[l+0x19] = pc;
		dif_offsets[l+0x0a] = p5+p10;		dif_offsets[l+0x1a] = pf;
		dif_offsets[l+0x0b] = p4+p10;		dif_offsets[l+0x1b] = pe;
		dif_offsets[l+0x0c] = p1+p10;		dif_offsets[l+0x1c] = pb;
		dif_offsets[l+0x0d] =    p10;		dif_offsets[l+0x1d] = pa;
		dif_offsets[l+0x0e] = p3+p10;		dif_offsets[l+0x1e] = p8;
		dif_offsets[l+0x0f] = p2+p10;		dif_offsets[l+0x1f] = p9;
		// Set 5: [9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p40],[b,a,8,9,f,e,c,d,6,7,5,4,1,0,3,2 + p50] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p9;		dif_offsets[l+0x10] = pb+p10;	// <*** same rcol offsets as Set 4 lcol
		dif_offsets[l+0x01] = p8;		dif_offsets[l+0x11] = pa+p10;	// <*** same lcol offsets as Set 6 rcol
		dif_offsets[l+0x02] = pb;		dif_offsets[l+0x12] = p8+p10;
		dif_offsets[l+0x03] = pa;		dif_offsets[l+0x13] = p9+p10;
		dif_offsets[l+0x04] = pd;		dif_offsets[l+0x14] = pf+p10;
		dif_offsets[l+0x05] = pc;		dif_offsets[l+0x15] = pe+p10;
		dif_offsets[l+0x06] = pf;		dif_offsets[l+0x16] = pc+p10;
		dif_offsets[l+0x07] = pe;		dif_offsets[l+0x17] = pd+p10;
		dif_offsets[l+0x08] = p4;		dif_offsets[l+0x18] = p6+p10;
		dif_offsets[l+0x09] = p5;		dif_offsets[l+0x19] = p7+p10;
		dif_offsets[l+0x0a] = p6;		dif_offsets[l+0x1a] = p5+p10;
		dif_offsets[l+0x0b] = p7;		dif_offsets[l+0x1b] = p4+p10;
		dif_offsets[l+0x0c] = p2;		dif_offsets[l+0x1c] = p1+p10;
		dif_offsets[l+0x0d] = p3;		dif_offsets[l+0x1d] =    p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = p3+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = p2+p10;
		// Set 6: [2,3,1,0,6,7,5,4,b,a,8,9,f,e,c,d + p30],[9,8,b,a,d,c,f,e,4,5,6,7,2,3,1,0 + p20] (mod p20):
		l += 32;
		dif_offsets[l+0x00] = p2+p10;		dif_offsets[l+0x10] = p9;	// <*** same rcol offsets as Set 5 lcol
		dif_offsets[l+0x01] = p3+p10;		dif_offsets[l+0x11] = p8;	// <*** same lcol offsets as Set 0 rcol
		dif_offsets[l+0x02] = p1+p10;		dif_offsets[l+0x12] = pb;
		dif_offsets[l+0x03] =    p10;		dif_offsets[l+0x13] = pa;
		dif_offsets[l+0x04] = p6+p10;		dif_offsets[l+0x14] = pd;
		dif_offsets[l+0x05] = p7+p10;		dif_offsets[l+0x15] = pc;
		dif_offsets[l+0x06] = p5+p10;		dif_offsets[l+0x16] = pf;
		dif_offsets[l+0x07] = p4+p10;		dif_offsets[l+0x17] = pe;
		dif_offsets[l+0x08] = pb+p10;		dif_offsets[l+0x18] = p4;
		dif_offsets[l+0x09] = pa+p10;		dif_offsets[l+0x19] = p5;
		dif_offsets[l+0x0a] = p8+p10;		dif_offsets[l+0x1a] = p6;
		dif_offsets[l+0x0b] = p9+p10;		dif_offsets[l+0x1b] = p7;
		dif_offsets[l+0x0c] = pf+p10;		dif_offsets[l+0x1c] = p2;
		dif_offsets[l+0x0d] = pe+p10;		dif_offsets[l+0x1d] = p3;
		dif_offsets[l+0x0e] = pc+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = pd+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/

		// Init storage for 2 circular-shifts perms of a basic 7-vector, with shift count in [0,6] that means 2*13

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0xc0<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0xa0<<1; dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0xc0<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x60<<1;
		dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0xb0<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x90<<1; dit_p20_cperms[l++] = 0xd0<<1; dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0xb0<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x90<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 6);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 1);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 2);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 3);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 5);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 6);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 5);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 6);

	   #else

		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p90;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-6; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-6 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 3) + 6);
		dit_p20_lo_offset[l++] = ((pe << 3) + 0);
		dit_p20_lo_offset[l++] = ((pd << 3) + 1);
		dit_p20_lo_offset[l++] = ((pc << 3) + 2);
		dit_p20_lo_offset[l++] = ((pb << 3) + 3);
		dit_p20_lo_offset[l++] = ((pa << 3) + 4);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 6);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 0);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 4);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 5);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 6);

	   #endif

	// dit_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		l = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p70],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p60] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pb+p10;		dit_offsets[l+0x10] = pb;
		dit_offsets[l+0x01] = pa+p10;		dit_offsets[l+0x11] = pa;
		dit_offsets[l+0x02] = p9+p10;		dit_offsets[l+0x12] = p9;
		dit_offsets[l+0x03] = p8+p10;		dit_offsets[l+0x13] = p8;
		dit_offsets[l+0x04] = pd+p10;		dit_offsets[l+0x14] = pd;
		dit_offsets[l+0x05] = pc+p10;		dit_offsets[l+0x15] = pc;
		dit_offsets[l+0x06] = pe+p10;		dit_offsets[l+0x16] = pe;
		dit_offsets[l+0x07] = pf+p10;		dit_offsets[l+0x17] = pf;
		dit_offsets[l+0x08] = p3+p10;		dit_offsets[l+0x18] = p3;
		dit_offsets[l+0x09] = p2+p10;		dit_offsets[l+0x19] = p2;
		dit_offsets[l+0x0a] = p1+p10;		dit_offsets[l+0x1a] = p1;
		dit_offsets[l+0x0b] =    p10;		dit_offsets[l+0x1b] =  0;
		dit_offsets[l+0x0c] = p5+p10;		dit_offsets[l+0x1c] = p5;
		dit_offsets[l+0x0d] = p4+p10;		dit_offsets[l+0x1d] = p4;
		dit_offsets[l+0x0e] = p6+p10;		dit_offsets[l+0x1e] = p6;
		dit_offsets[l+0x0f] = p7+p10;		dit_offsets[l+0x1f] = p7;
		// Set 2: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + pc0],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + pd0] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pd;		dit_offsets[l+0x10] = p5+p10;
		dit_offsets[l+0x01] = pc;		dit_offsets[l+0x11] = p4+p10;
		dit_offsets[l+0x02] = pe;		dit_offsets[l+0x12] = p6+p10;
		dit_offsets[l+0x03] = pf;		dit_offsets[l+0x13] = p7+p10;
		dit_offsets[l+0x04] = p9;		dit_offsets[l+0x14] = p1+p10;
		dit_offsets[l+0x05] = p8;		dit_offsets[l+0x15] =    p10;
		dit_offsets[l+0x06] = pa;		dit_offsets[l+0x16] = p2+p10;
		dit_offsets[l+0x07] = pb;		dit_offsets[l+0x17] = p3+p10;
		dit_offsets[l+0x08] = p5;		dit_offsets[l+0x18] = p9+p10;
		dit_offsets[l+0x09] = p4;		dit_offsets[l+0x19] = p8+p10;
		dit_offsets[l+0x0a] = p6;		dit_offsets[l+0x1a] = pa+p10;
		dit_offsets[l+0x0b] = p7;		dit_offsets[l+0x1b] = pb+p10;
		dit_offsets[l+0x0c] = p1;		dit_offsets[l+0x1c] = pe+p10;
		dit_offsets[l+0x0d] =  0;		dit_offsets[l+0x1d] = pf+p10;
		dit_offsets[l+0x0e] = p2;		dit_offsets[l+0x1e] = pc+p10;
		dit_offsets[l+0x0f] = p3;		dit_offsets[l+0x1f] = pd+p10;
		// Set 3: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p40],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p50] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p9;		dit_offsets[l+0x10] = p1+p10;
		dit_offsets[l+0x01] = p8;		dit_offsets[l+0x11] =    p10;
		dit_offsets[l+0x02] = pa;		dit_offsets[l+0x12] = p2+p10;
		dit_offsets[l+0x03] = pb;		dit_offsets[l+0x13] = p3+p10;
		dit_offsets[l+0x04] = pe;		dit_offsets[l+0x14] = p6+p10;
		dit_offsets[l+0x05] = pf;		dit_offsets[l+0x15] = p7+p10;
		dit_offsets[l+0x06] = pc;		dit_offsets[l+0x16] = p4+p10;
		dit_offsets[l+0x07] = pd;		dit_offsets[l+0x17] = p5+p10;
		dit_offsets[l+0x08] = p1;		dit_offsets[l+0x18] = pe+p10;
		dit_offsets[l+0x09] =  0;		dit_offsets[l+0x19] = pf+p10;
		dit_offsets[l+0x0a] = p2;		dit_offsets[l+0x1a] = pc+p10;
		dit_offsets[l+0x0b] = p3;		dit_offsets[l+0x1b] = pd+p10;
		dit_offsets[l+0x0c] = p6;		dit_offsets[l+0x1c] = pa+p10;
		dit_offsets[l+0x0d] = p7;		dit_offsets[l+0x1d] = pb+p10;
		dit_offsets[l+0x0e] = p4;		dit_offsets[l+0x1e] = p8+p10;
		dit_offsets[l+0x0f] = p5;		dit_offsets[l+0x1f] = p9+p10;
		// Set 4: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pb0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pa0] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p6+p10;		dit_offsets[l+0x10] = p6;
		dit_offsets[l+0x01] = p7+p10;		dit_offsets[l+0x11] = p7;
		dit_offsets[l+0x02] = p4+p10;		dit_offsets[l+0x12] = p4;
		dit_offsets[l+0x03] = p5+p10;		dit_offsets[l+0x13] = p5;
		dit_offsets[l+0x04] = p2+p10;		dit_offsets[l+0x14] = p2;
		dit_offsets[l+0x05] = p3+p10;		dit_offsets[l+0x15] = p3;
		dit_offsets[l+0x06] =    p10;		dit_offsets[l+0x16] =  0;
		dit_offsets[l+0x07] = p1+p10;		dit_offsets[l+0x17] = p1;
		dit_offsets[l+0x08] = pa+p10;		dit_offsets[l+0x18] = pa;
		dit_offsets[l+0x09] = pb+p10;		dit_offsets[l+0x19] = pb;
		dit_offsets[l+0x0a] = p8+p10;		dit_offsets[l+0x1a] = p8;
		dit_offsets[l+0x0b] = p9+p10;		dit_offsets[l+0x1b] = p9;
		dit_offsets[l+0x0c] = pc+p10;		dit_offsets[l+0x1c] = pc;
		dit_offsets[l+0x0d] = pd+p10;		dit_offsets[l+0x1d] = pd;
		dit_offsets[l+0x0e] = pf+p10;		dit_offsets[l+0x1e] = pf;
		dit_offsets[l+0x0f] = pe+p10;		dit_offsets[l+0x1f] = pe;
		// Set 5: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p20] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p2+p10;		dit_offsets[l+0x10] = p2;
		dit_offsets[l+0x01] = p3+p10;		dit_offsets[l+0x11] = p3;
		dit_offsets[l+0x02] =    p10;		dit_offsets[l+0x12] =  0;
		dit_offsets[l+0x03] = p1+p10;		dit_offsets[l+0x13] = p1;
		dit_offsets[l+0x04] = p4+p10;		dit_offsets[l+0x14] = p4;
		dit_offsets[l+0x05] = p5+p10;		dit_offsets[l+0x15] = p5;
		dit_offsets[l+0x06] = p7+p10;		dit_offsets[l+0x16] = p7;
		dit_offsets[l+0x07] = p6+p10;		dit_offsets[l+0x17] = p6;
		dit_offsets[l+0x08] = pc+p10;		dit_offsets[l+0x18] = pc;
		dit_offsets[l+0x09] = pd+p10;		dit_offsets[l+0x19] = pd;
		dit_offsets[l+0x0a] = pf+p10;		dit_offsets[l+0x1a] = pf;
		dit_offsets[l+0x0b] = pe+p10;		dit_offsets[l+0x1b] = pe;
		dit_offsets[l+0x0c] = p8+p10;		dit_offsets[l+0x1c] = p8;
		dit_offsets[l+0x0d] = p9+p10;		dit_offsets[l+0x1d] = p9;
		dit_offsets[l+0x0e] = pb+p10;		dit_offsets[l+0x1e] = pb;
		dit_offsets[l+0x0f] = pa+p10;		dit_offsets[l+0x1f] = pa;
		// Set 6: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p80],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p90] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p4;		dit_offsets[l+0x10] = p8+p10;
		dit_offsets[l+0x01] = p5;		dit_offsets[l+0x11] = p9+p10;
		dit_offsets[l+0x02] = p7;		dit_offsets[l+0x12] = pb+p10;
		dit_offsets[l+0x03] = p6;		dit_offsets[l+0x13] = pa+p10;
		dit_offsets[l+0x04] =  0;		dit_offsets[l+0x14] = pf+p10;
		dit_offsets[l+0x05] = p1;		dit_offsets[l+0x15] = pe+p10;
		dit_offsets[l+0x06] = p3;		dit_offsets[l+0x16] = pd+p10;
		dit_offsets[l+0x07] = p2;		dit_offsets[l+0x17] = pc+p10;
		dit_offsets[l+0x08] = p8;		dit_offsets[l+0x18] =    p10;
		dit_offsets[l+0x09] = p9;		dit_offsets[l+0x19] = p1+p10;
		dit_offsets[l+0x0a] = pb;		dit_offsets[l+0x1a] = p3+p10;
		dit_offsets[l+0x0b] = pa;		dit_offsets[l+0x1b] = p2+p10;
		dit_offsets[l+0x0c] = pf;		dit_offsets[l+0x1c] = p7+p10;
		dit_offsets[l+0x0d] = pe;		dit_offsets[l+0x1d] = p6+p10;
		dit_offsets[l+0x0e] = pd;		dit_offsets[l+0x1e] = p5+p10;
		dit_offsets[l+0x0f] = pc;		dit_offsets[l+0x1f] = p4+p10;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		tmp = r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1c0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x1c0;
		two  	= tmp + 0x0;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one 	= tmp + 0x1;
		sqrt2	= tmp + 0x2;
		isrt2	= tmp + 0x3;
		cc2		= tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		ss2		= tmp + 0x5;
		cc1		= tmp + 0x6;
		ss1		= tmp + 0x7;
		cc3		= tmp + 0x8;
		ss3		= tmp + 0x9;
		dc0		= tmp + 0x0a;	// radix-7 DFT trig consts - Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here
		ds0		= tmp + 0x0b;
		dc1		= tmp + 0x0c;
		ds1		= tmp + 0x0d;
		dc2		= tmp + 0x0e;
		ds2		= tmp + 0x0f;
		dc3  	= tmp + 0x10;
		ds3		= tmp + 0x11;	// Add extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro:
		tmp += 0x16;	// sc_ptr += 0x396
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x1c;	tmp += 2*0x1c;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x38;	tmp += 0x70;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;  +30 = 338
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(396 + 70 + 2) = 0x408; This is where the value of half_arr_offset224 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #else
		cy_r = tmp;	cy_i = tmp+0x70;	tmp += 0xe0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays; +60 = 368 complex
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(396 + e0 + 2) = 0x476; This is where the value of half_arr_offset224 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif

		ASSERT((two->d0 == 2.0 && two->d1 == 2.0), "thread-local memcheck failed!");
	  #if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD)
		// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
		ASSERT((ds3->d0 == 0.0 && ds3->d1 == 0.0), "thread-local memcheck failed!");
	  #else
		/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
		ASSERT((ds3->d0 == sx3 && ds3->d1 == sx3), "thread-local memcheck failed!");
	  #endif
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
	  #ifndef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts:
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	  #endif
		tmp = half_arr;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	  #ifdef USE_AVX512
		/* No-Op */
	  #elif defined(USE_AVX)
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix224_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		n_minus_sil   = (struct uint32x16*)sse_n + 1;
		n_minus_silp1 = (struct uint32x16*)sse_n + 2;
		sinwt         = (struct uint32x16*)sse_n + 3;
		sinwtm1       = (struct uint32x16*)sse_n + 4;
	   #else
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
	   #endif
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
	  #endif
	  #ifdef USE_AVX
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
		// is otherwise unused in Fermat-Mod mode - for local storage of these cycle tables.
		/*** Pointer-inits: ***/
		int *icycle = bjmodn,ic_idx;
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int *jcycle = icycle + ODD_RADIX,jc_idx;
	  #ifdef USE_AVX
		int *kcycle = jcycle + ODD_RADIX,kc_idx;
		int *lcycle = kcycle + ODD_RADIX,lc_idx;
	  #endif
	  #ifdef USE_AVX512
		int *mcycle = lcycle + ODD_RADIX,mc_idx;
		int *ncycle = mcycle + ODD_RADIX,nc_idx;
		int *ocycle = ncycle + ODD_RADIX,oc_idx;
		int *pcycle = ocycle + ODD_RADIX,pc_idx;
	  #endif
	#endif
		/*** Value-inits: ***/
		for(j = 0; j < ODD_RADIX; j++) {
			icycle[j] = thread_arg->icycle[j];
		#ifdef USE_SSE2
			jcycle[j] = thread_arg->jcycle[j];
		  #ifdef USE_AVX
			kcycle[j] = thread_arg->kcycle[j];
			lcycle[j] = thread_arg->lcycle[j];
		  #endif
		  #ifdef USE_AVX512
			mcycle[j] = thread_arg->mcycle[j];
			ncycle[j] = thread_arg->ncycle[j];
			ocycle[j] = thread_arg->ocycle[j];
			pcycle[j] = thread_arg->pcycle[j];
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
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				tmp->d0 = *(addr+l  );
				tmp->d1 = *(addr+l+1);
				tmp->d2 = *(addr+l+2);
				tmp->d3 = *(addr+l+3);
				tmp->d4 = *(addr+l+4);
				tmp->d5 = *(addr+l+5);
				tmp->d6 = *(addr+l+6);
				tmp->d7 = *(addr+l+7);
			}
		#elif defined(USE_AVX)
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
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				tmp->d0 = *(addr+l  );		tm2->d0 = *(addi+l  );
				tmp->d1 = *(addr+l+1);		tm2->d1 = *(addi+l+1);
				tmp->d2 = *(addr+l+2);		tm2->d2 = *(addi+l+2);
				tmp->d3 = *(addr+l+3);		tm2->d3 = *(addi+l+3);
				tmp->d4 = *(addr+l+4);		tm2->d4 = *(addi+l+4);
				tmp->d5 = *(addr+l+5);		tm2->d5 = *(addi+l+5);
				tmp->d6 = *(addr+l+6);		tm2->d6 = *(addi+l+6);
				tmp->d7 = *(addr+l+7);		tm2->d7 = *(addi+l+7);
			}
		#elif defined(USE_AVX)
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
		#include "radix224_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			addr = thread_arg->cy_r;
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				*(addr+l  ) = tmp->d0;
				*(addr+l+1) = tmp->d1;
				*(addr+l+2) = tmp->d2;
				*(addr+l+3) = tmp->d3;
				*(addr+l+4) = tmp->d4;
				*(addr+l+5) = tmp->d5;
				*(addr+l+6) = tmp->d6;
				*(addr+l+7) = tmp->d7;
			}
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		#elif defined(USE_AVX)
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
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tm2->d0;
				*(addr+l+1) = tmp->d1;		*(addi+l+1) = tm2->d1;
				*(addr+l+2) = tmp->d2;		*(addi+l+2) = tm2->d2;
				*(addr+l+3) = tmp->d3;		*(addi+l+3) = tm2->d3;
				*(addr+l+4) = tmp->d4;		*(addi+l+4) = tm2->d4;
				*(addr+l+5) = tmp->d5;		*(addi+l+5) = tm2->d5;
				*(addr+l+6) = tmp->d6;		*(addi+l+6) = tm2->d6;
				*(addr+l+7) = tmp->d7;		*(addi+l+7) = tm2->d7;
			}
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		#elif defined(USE_AVX)
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

#undef RADIX
#undef ODD_RADIX
#undef PFETCH_DIST
