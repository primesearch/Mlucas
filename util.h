/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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
#ifndef util_h_included
#define util_h_included

#include "masterdefs.h"
#include "types.h"
#include "float_intrin.h"
#include "mi64.h"
#include "prefetch.h"	/* Since this file doesn't directly include Mlucas.h, need this here */
#include "qfloat.h"
#include "rng_isaac.h"

#include "Mdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* get_fp_rnd_const.c: */
void	get_fp_rnd_const(double*RND_A, double*RND_B);

/* test_fft_radix.c: */
void	test_fft_radix(void);

/* util.c: */
void	host_init(void);	/* This one is a wrapper for calls to the next few: */
char*	get_time_str(double tdiff);

/* get_cpuid.c: x86-style CPUs */
void	get_cpu(void);
uint32	has_sse  (void);
uint32	has_sse2 (void);
uint32	has_sse3 (void);
void	cpu_details(void);

void	print_host_info(void);
void	set_x87_fpu_params(unsigned short FPU_MODE);
void	info_x87_fpu_ctrl(uint64 FPUCTRL);
void	check_nbits_in_types(void);
int		test_mul(void);	/* This one is actually defined in imul_macro.c */
double	errprint_sincos	(double *x, double *y, double theta);

/* Multithreading-related: */
int		get_num_cores(void);
int		test_pthreads(int ncpu, int verbose);
void* 	ex_loop(void* data);
void*	PrintHello(void *threadid);
void*	do_loop(void*targ);

int		file_valid(FILE*fp);
/* Need a portable way to implement these:
int		file_valid_for_read	(FILE*fp);
int		file_valid_for_write(FILE*fp);
*/
void	INFO	(long line, char*file, char*info_string, char*info_file, int copy2stderr);
void	WARN	(long line, char*file, char*warn_string, char*warn_file, int copy2stderr);
void	ASSERT	(long line, char*file, int expr, char*assert_string);
void	VAR_WARN(char *typelist, ...);

uint32	trailz32(uint32 x);
uint32	trailz64(uint64 x);

uint32	leadz32	(uint32  i);
uint32	leadz64	(uint64  i);
uint32	leadz128(uint128 i);
uint32	leadz192(uint192 i);
uint32	leadz256(uint256 i);

uint32	isPow2(uint32 i32);
uint32	isPow4(uint32 i32);
uint32	popcount32(uint32 num);
uint64	popcount64(uint64 num);

uint32	ibits32	(uint32 i, uint32 beg, uint32 nbits);
uint64	ibits64	(uint64 i, uint32 beg, uint32 nbits);
uint64	getbits64(uint64 x, uint32 src_bit_start, uint32 nbits,           uint32 tgt_bit_start);
void	mvbits64 (uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start);
int		pprimeF	(uint32 p, uint32 z);
int		isPRP	(uint32 p);
int		pprimeF_64	(uint64 p, uint64 z);
int		isPRP_64	(uint64 p);
uint32	twompmodq32(uint32 p, uint32 q);	// 2^-p % q
int		twopmodq32(uint32 p, uint32 q);		// (2^-p % q) == 0
int		twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7);
uint32	gcd32(uint32 x, uint32 y);
uint64	gcd64(uint64 x, uint64 y);
uint32	egcd32	(uint32 *x, uint32 *y);
int		modinv32(uint32 z, uint32 n);

struct	complex cmul(struct complex *, struct complex *);

uint128 xmody128(uint128 x, uint128 y);
uint192 xmody192(const uint192 x, const uint192 y, uint192*quot);
uint256 xmody256(const uint256 x, const uint256 y, uint256*quot);

uint32	x128_div_y32(uint128 *x, uint32 y);
uint32	x192_div_y32(uint192 *x, uint32 y);
uint32	x256_div_y32(uint256 *x, uint32 y);

uint32	x128_mod_y32(uint128  x, uint32 y);
uint32	x192_mod_y32(uint192  x, uint32 y);
uint32	x256_mod_y32(uint256  x, uint32 y);

/* Need the first of these because some compilers (e.g. .NET) don't properly print 64-bit ints: */
int		convert_uint64_base10_char (char*char_buf, uint64  q   );
int		convert_uint64_base16_char (char*char_buf, uint64  q   );
int		convert_uint96_base10_char (char*char_buf, uint96  q96 );
int		convert_uint96ptr_base10_char(char*char_buf, uint96*q96);
int		convert_uint128_base10_char(char*char_buf, uint128 q128);
int		convert_uint192_base10_char(char*char_buf, uint192 q192);
int		convert_uint256_base10_char(char*char_buf, uint256 q256);

double	convert_base10_char_double (const char*char_buf);
uint64	convert_base10_char_uint64 (const char*char_buf);
uint96	convert_base10_char_uint96 (const char*char_buf);
uint128	convert_base10_char_uint128(const char*char_buf);
uint192	convert_base10_char_uint192(const char*char_buf);
uint256	convert_base10_char_uint256(const char*char_buf);

uint64	TEST_BIT96 (uint96  __x, uint32 __bit);
uint64	TEST_BIT128(uint128 __x, uint32 __bit);
uint64	TEST_BIT160(uint160 __x, uint32 __bit);
uint64	TEST_BIT192(uint192 __x, uint32 __bit);
uint64	TEST_BIT256(uint256 __x, uint32 __bit);

double	finvest  (double x, uint32 numbits);
double	fisqrtest(double x, uint32 numbits);

/* imul50_macro.c: */
/*********** 50x50==>100-bit fmadd-based product algorithms ************/
/*********** Should only be called if USE_FMADD defined! *******/
void	test_mul50x50();
void	mul50x50_debug(double a, double b, double *lo, double *hi);
void	mul52x52_debug(double a, double b, double *lo, double *hi);
void	mul53x53_debug(double a, double b, double *lo, double *hi);

/*********************** 50+ BIT MULTIPLY MACROS **********************/

#ifdef __cplusplus
}
#endif

#endif	/* util_h_included */

