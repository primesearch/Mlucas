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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef Mlucas_h_included
#define Mlucas_h_included

#include "align.h"
#include "carry.h"
#include "dft_macro.h"
#include "factor.h"
#include "prefetch.h"
#include "util.h"

#ifdef INCLUDE_PM1
	#include "gcd_lehmer.h"
#endif

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* Mlucas.c: */
void	Mlucas_init(void);
uint32	ernstMain
(
	int		mod_type,
	int		test_type,
	uint64	exponent,
	uint32	fft_length,
	int		radix_set,
	uint32	maxFFT,
	uint32	iterations,	/* Use to store log2[max factor depth] in TF mode */
	int		error_checking,
	uint64	*sh0,	/* Reference/Computed mod-(2^64) residue */
	uint64	*sh1,	/* Reference/Computed Selfridge/Hurwitz Residue mod 2^35-1 */
	uint64	*sh2,	/* Reference/Computed Selfridge/Hurwitz Residue mod 2^36-1 */
	int		scrnFlag,
	double	*runtime
);

int		cfgNeedsUpdating(char*in_line);
void	printMlucasErrCode(int retVal);
uint64	res64	(double a[], int n, const uint64 p, int *nbits, char *hex_res);
void	resSH	(double a[], int n, const uint64 p, uint64*Res35m1, uint64*Res36m1);
void	hex_res_printtofile(double a[], int n, const uint64 p, int timing_test_iters, FILE *fp);
int		convert_LL_savefiles(uint64 psave, FILE*fp, uint32*ilo, uint32 ndim, int32 arr_scratch[], double a[]);
int		read_ppm1_savefiles	(uint64 p, FILE*fp, uint32*ilo, uint8 arr_tmp[], uint64*Res64, uint64*Res35m1, uint64*Res36m1);
/* Moved this prot to Mlucas.c:
void	write_ppm1_savefiles(uint64 p, FILE*fp, uint32 ihi, uint8 arr_tmp[], uint64 Res64, uint64 Res35m1, uint64 Res36m1);
*/
int		convert_res_bytewise_FP(const uint8 arr_tmp[], double a[], int n, const uint64 p, const uint64 Res64, const uint64 Res35m1, const uint64 Res36m1);
void	convert_res_FP_bytewise(const double a[], uint8 arr_tmp[], int n, const uint64 p,       uint64*Res64,       uint64*Res35m1,       uint64*Res36m1);
uint32	get_default_factoring_depth(uint64 p);
void	write_fft_debug_data(double a[], int jlo, int jhi);

/* br.c: */
int		reverse(uint32 i, uint32 n);
void	bit_reverse_int(int vec[], int n, int nradices, int radix[], int incr, int*scratch);

/* get_fft_radices.c: */
int		get_fft_radices			(uint32 kblocks, int radix_set, int *nradices, int radix_vec[], int radix_vec_dim);
void	test_fft_radixtables	(void);
uint32	get_default_fft_length	(uint32 p);
uint32	get_nextlarger_fft_length	(uint32 n);
uint32	given_N_get_maxP		(uint32 N);

/* get_preferred_fft_radix.c: */
uint32	get_preferred_fft_radix(uint32 kblocks);
uint32	extractFFTlengthFrom32Bit (uint32 n);
void	extractFFTradicesFrom32Bit(uint32 n);

/* radix*_ditN_cy_dif1.c: */
void	radix5_dif_pass1	(double a[], int n);
void	radix6_dif_pass1	(double a[], int n);
void	radix6_dif_pass1A	(double a[], int n);
void	radix6_dif_pass1B	(double a[], int n);
void	radix7_dif_pass1	(double a[], int n);
void	radix8_dif_pass1	(double a[], int n);
void	radix9_dif_pass1	(double a[], int n);
void	radix10_dif_pass1	(double a[], int n);
void	radix11_dif_pass1	(double a[], int n);
void	radix12_dif_pass1	(double a[], int n);
void	radix13_dif_pass1	(double a[], int n);
void	radix14_dif_pass1	(double a[], int n);
void	radix15_dif_pass1	(double a[], int n);
void	radix16_dif_pass1	(double a[], int n);
void	radix18_dif_pass1	(double a[], int n);
void	radix20_dif_pass1	(double a[], int n);
void	radix22_dif_pass1	(double a[], int n);
void	radix24_dif_pass1	(double a[], int n);
void	radix26_dif_pass1	(double a[], int n);
void	radix28_dif_pass1	(double a[], int n);
void	radix30_dif_pass1	(double a[], int n);
void	radix32_dif_pass1	(double a[], int n);
void	radix36_dif_pass1	(double a[], int n);
void	radix40_dif_pass1	(double a[], int n);
void	radix44_dif_pass1	(double a[], int n);
void	radix48_dif_pass1	(double a[], int n);
void	radix52_dif_pass1	(double a[], int n);
void	radix56_dif_pass1	(double a[], int n);
void	radix60_dif_pass1	(double a[], int n);
void	radix64_dif_pass1	(double a[], int n);
void	radix72_dif_pass1	(double a[], int n);
void	radix80_dif_pass1	(double a[], int n);

void	radix5_dit_pass1	(double a[], int n);
void	radix6_dit_pass1	(double a[], int n);
void	radix7_dit_pass1	(double a[], int n);
void	radix8_dit_pass1	(double a[], int n);
void	radix9_dit_pass1	(double a[], int n);
void	radix10_dit_pass1	(double a[], int n);
void	radix11_dit_pass1	(double a[], int n);
void	radix12_dit_pass1	(double a[], int n);
void	radix13_dit_pass1	(double a[], int n);
void	radix13_dit_pass1A	(double a[], int n);
void	radix13_dit_pass1B	(double a[], int n);
void	radix14_dit_pass1	(double a[], int n);
void	radix15_dit_pass1	(double a[], int n);
void	radix16_dit_pass1	(double a[], int n);
void	radix18_dit_pass1	(double a[], int n);
void	radix20_dit_pass1	(double a[], int n);
void	radix22_dit_pass1	(double a[], int n);
void	radix24_dit_pass1	(double a[], int n);
void	radix26_dit_pass1	(double a[], int n);
void	radix28_dit_pass1	(double a[], int n);
void	radix30_dit_pass1	(double a[], int n);
void	radix32_dit_pass1	(double a[], int n);
void	radix36_dit_pass1	(double a[], int n);
void	radix40_dit_pass1	(double a[], int n);
void	radix44_dit_pass1	(double a[], int n);
void	radix48_dit_pass1	(double a[], int n);
void	radix52_dit_pass1	(double a[], int n);
void	radix56_dit_pass1	(double a[], int n);
void	radix60_dit_pass1	(double a[], int n);
void	radix64_dit_pass1	(double a[], int n);
void	radix72_dit_pass1	(double a[], int n);
void	radix80_dit_pass1	(double a[], int n);

void	radix8_dif_pass		(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix16_dif_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix32_dif_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix64_dif_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);

void	radix8_dit_pass		(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix16_dit_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix32_dit_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);
void	radix64_dit_pass	(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr);

/* mers_mod_square.c: */
int		mers_mod_square		(double a[], int arr_scratch[], int n, int ilo, int ihi, uint64 p, uint32 *err_iter, int scrnFlag, double *tdiff);
void	mers_process_chunk	(double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[], int index[], int block_index[], int ii, int nradices_prim, int radix_prim[], int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[]);

	/* radix{16|32|64}_wrapper_square.c: */
	void	pair_square  (double *x1, double *y1, double *x2, double *y2, double c, double s);
	void	pair_square2A(double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, double *x4, double *y4, double c, double s);
	void	pair_square2B(double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, double *x4, double *y4, double c, double s);

	void	radix16_wrapper_ini		(int n, int radix0, int iblock, int nradices_prim, int radix_prim[], int i[], int j1[], int j2[], int j2_start[], int k[], int m[], int blocklen[], int blocklen_sum[]);
	void	radix16_wrapper_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum, int INIT_MODE);

	void	radix32_wrapper_ini		(int n, int radix0, int iblock, int nradices_prim, int radix_prim[], int i[], int j1[], int j2[], int j2_start[], int k[], int m[], int blocklen[], int blocklen_sum[]);
	void	radix32_wrapper_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum, int INIT_MODE);

	void	radix64_wrapper_ini		(int n, int radix0, int iblock, int nradices_prim, int radix_prim[], int i[], int j1[], int j2[], int j2_start[], int k[], int m[], int blocklen[], int blocklen_sum[]);
	void	radix64_wrapper_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum, int INIT_MODE);

/* fermat_mod_square.c: */
int 	fermat_mod_square	(double a[], int arr_scratch[], int n, int ilo, int ihi, uint64 p, uint32 *err_iter, int scrnFlag, double *tdiff);
void	fermat_process_chunk(double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[], int index[],                    int ii, int nradices_prim, int radix_prim[], int skip_square);

	/* radix{16|32|64}_dyadic_square.c: */
	void	radix16_dyadic_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr);
	void	radix32_dyadic_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr);
	void	radix64_dyadic_square	(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr);

#ifdef INCLUDE_PM1
	/* pairFFT_mul.c: */
	void	pairFFT_mul(double x[], double y[], int n, int INIT_ARRAYS, int FORWARD_FFT_ONLY);
	void	pairFFT_mul_process_chunk(double a[], double ab_mul[], double cd_mul[], int n, int nradices, int radix_vec[], struct complex rt0[], struct complex rt1[], int index[], int ii, int nradices_prim, int radix_prim[], int FORWARD_FFT_ONLY, int skip_square);

	void	radix16_pairFFT_mul		(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int INIT_ARRAYS, int FORWARD_FFT_ONLY, int skip_square);
	void	radix32_pairFFT_mul		(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int INIT_ARRAYS, int FORWARD_FFT_ONLY, int skip_square);
#endif

/* Each of the final-DIT-pass/propagate-carries/initial-DIF-pass routines comes in both an error-checking and non-error-checking version: */

	int	radix5_ditN_cy_dif1 		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix6_ditN_cy_dif1 		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int radix7_ditN_cy_dif1 		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int radix8_ditN_cy_dif1 		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix9_ditN_cy_dif1 		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix10_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix11_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix12_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix13_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix14_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix15_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int radix16_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix18_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix20_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix22_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix24_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix26_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix28_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix30_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix32_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix36_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix40_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix44_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix48_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix52_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix56_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix60_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix64_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix72_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix80_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);
	int	radix88_ditN_cy_dif1		(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter, double *fracmax, uint64 p);

	/* Only define 2nd version of carry routine[s] with ROE checking disabled in non-SSE2 mode, as SSE2 ROE checking is cheap.
	For large radices (> 36) also define just a single (ROE-enabled) version of the carry routine and use it irrespective of the command-line ROE-check flag value: */
  #ifndef USE_SSE2
	int	radix5_ditN_cy_dif1_nochk 	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix6_ditN_cy_dif1_nochk 	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix7_ditN_cy_dif1_nochk 	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int radix8_ditN_cy_dif1_nochk 	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int	radix9_ditN_cy_dif1_nochk 	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix10_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix11_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix12_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix13_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix14_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int	radix15_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int radix16_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int	radix18_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix20_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix22_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix24_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix26_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[],                                             double base[], double baseinv[], int iter,                  uint64 p);
	int	radix28_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int	radix30_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int radix32_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
	int	radix36_ditN_cy_dif1_nochk	(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter,                  uint64 p);
  #endif

#endif	/* Mlucas_h_included */
