/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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

// Control over floating type used in DNINT emulation:
#ifdef X64_ASM
	#if FP_MANTISSA_BITS_DOUBLE != 64
		#error x86_64 asm requires FP_MANTISSA_BITS_DOUBLE == 64!
	#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __CUDA_ARCH__
	#define DEV __device__
#else
	#define DEV /* */
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

/* getRealTime.c: */
double	 getRealTime();

/* get_cpuid.c: x86-style CPUs */
void	get_cpu(void);
uint32	has_sse  (void);
uint32	has_sse2 (void);
uint32	has_sse3 (void);
uint32	has_sse3e(void);
uint32	has_sse41(void);
uint32	has_sse42(void);
uint32	has_avx  (void);
uint32	has_avx2 (void);
uint32	has_avx512(void);
void	cpu_details(void);

/* util.c: */
void	host_init(void);	/* This one is a wrapper for calls to the next few: */
char*	get_time_str(double tdiff);
void	set_stacklimit_restart(char *argv[]);
void	print_host_info(void);
void	set_x87_fpu_params(unsigned short FPU_MODE);
void	info_x87_fpu_ctrl(uint64 FPUCTRL);
void	check_nbits_in_types(void);
int		test_mul(void);	/* This one is actually defined in imul_macro.c */
double	errprint_sincos	(double *x, double *y, double theta);

// EWM: Jun 2015: Path-related code due to Alex Vong,
// as part of his Debian-freeware-packaging-of-Mlucas project:
extern char *MLUCAS_PATH; /* MLUCAS_PATH is set by set_mlucas_path()  */
void	set_mlucas_path(void); /* Set MLUCAS_PATH  */
char	*quote_spaces(char *dest, char *src); /* Double-quote spaces in string  */
int		mkdir_p(char *path); /* Emulate `mkdir -p path'  */
char	*shell_quote(char *dest, char *src); /* Escape shell meta characters  */
FILE	*mlucas_fopen(const char *path, const char *mode); /* fopen() wrapper  */

#ifdef USE_GPU
//#if defined(USE_GPU) && defined(__CUDACC__)
	// GPU-related diagnostics and basic self-tests:
	void cudaVecModpowTest64();
	void cudaVecModpowTest78_0();
	void cudaVecModpowTest78();
	void cudaVecModpowTest96();
#endif

/* Multithreading-specific stuff: */
#ifdef MULTITHREAD

	int		get_num_cores(void);
	int		test_pthreads(int ncpu, int verbose);
	void* 	ex_loop(void* data);
	void*	PrintHello(void *threadid);
	void*	do_loop(void*targ);
	uint32	parseAffinityTriplet(char*istr);
	void	parseAffinityString(char*istr);

#endif	// MULTITHREAD ?

int		file_valid(FILE*fp);
/* Need a portable way to implement these:
int		file_valid_for_read	(FILE*fp);
int		file_valid_for_write(FILE*fp);
*/
void	INFO	(long line, char*file, char*info_string, char*info_file, int copy2stderr);
void	WARN	(long line, char*file, char*warn_string, char*warn_file, int copy2stderr);
#ifdef __CUDA_ARCH__
	__device__
	void	ASSERT(long line, char*file, int expr, char*assert_string);
#else
	void	ASSERT	(long line, char*file, int expr, char*assert_string);
#endif
void	VAR_WARN(char *typelist, ...);

void	byte_bitstr(const uint8  byte, char*ostr);
void	ui32_bitstr(const uint32 ui32, char*ostr);
void	ui64_bitstr(const uint64 ui64, char*ostr);

int		reverse(uint32 i, uint32 n);

// 32 and 64-bit analogs of the F90 intrinsic ISHFT function:
DEV uint32	ishft32(uint32 x, int shift);
DEV uint64	ishft64(uint64 x, int shift);

DEV uint32	trailz32(uint32 x);
DEV uint32	trailz64(uint64 x);

DEV uint32	leadz32	(uint32  i);
DEV uint32	leadz64	(uint64  i);
DEV uint32	leadz128(uint128 i);
DEV uint32	leadz192(uint192 i);
DEV uint32	leadz256(uint256 i);

DEV uint32	isPow2(uint32 i32);
DEV uint32	isPow4(uint32 i32);
DEV uint64	isPow2_64(uint64 i64);
DEV uint64	isPow4_64(uint64 i64);
DEV uint32	popcount32(uint32 num);
DEV uint32	popcount64(uint64 num);
DEV int		ith_set_bit32(uint32 num, uint32 i);
DEV int		ith_set_bit64(uint64 num, uint32 i);

DEV uint32	ibits32	(uint32 i, uint32 beg, uint32 nbits);
DEV uint64	ibits64	(uint64 i, uint32 beg, uint32 nbits);
DEV uint64	getbits64(uint64 x, uint32 src_bit_start, uint32 nbits,           uint32 tgt_bit_start);
DEV void	mvbits64 (uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start);
DEV int		pprimeF	(uint32 p, uint32 z);
DEV int		isPRP	(uint32 p);
DEV int		pprimeF64(uint64 p, uint64 z);
DEV int		isPRP64	(uint64 p);
DEV uint32	twompmodq32(uint32 p, uint32 q);	// 2^-p % q
DEV int		twopmodq32(uint32 p, uint32 q);		// (2^-p % q) == 0
DEV int		twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7);
DEV uint64	twompmodq64(uint64 p, uint64 q);	// 2^-p % q
DEV uint32	gcd32(uint32 x, uint32 y);
DEV uint64	gcd64(uint64 x, uint64 y);
DEV uint32	egcd32	(uint32 *x, uint32 *y);
DEV int		modinv32(uint32 z, uint32 n);

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
int		convert_uint64_base2_char  (char*char_buf, uint64  q   );
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

/* To_Do: Collect in imul50_macro.[h|c] */
/*********** 50x50==>100-bit fmadd-based product algorithms - Should only be called if USE_FMADD defined! *******/
uint32	test_mul50x50();
void	mul50x50_debug(double a, double b, double *lo, double *hi);
int		cmp_fma_lohi_vs_exact(double dx, double dy, double dhi, double dlo, uint64 ix, uint64 iy, uint64 ihi, uint64 ilo);

/*********************** 50+ BIT MULTIPLY MACROS **********************/
	/* to-do! */

/*********************** Bit-level unitilities: *************************************/
void bit_clr32(uint32*arr, uint32 bit);
void bit_clr64(uint64*arr, uint64 bit);
void bit_clr32_x4(uint32*arr, uint32 bit1, uint32 bit2, uint32 bit3, uint32 bit4);
void bit_clr64_x4(uint64*arr, uint64 bit1, uint64 bit2, uint64 bit3, uint64 bit4);

/**** Core-FFT-routines functionality tests wrapped in timed repeater loops in util.c: ****/
int	test_radix4_dft();
int	test_radix8_dft();
int	test_radix16_dft();
int	test_radix32_dft();

// In transpose-tests, NxN refers to linear memory treated as a NxN matrix-of-doubles:
#ifdef USE_AVX
int	test_simd_transpose_4x4();
#endif
#ifdef USE_AVX512
int	test_simd_transpose_8x8();
#endif
#ifdef USE_AVX1024
int	test_simd_transpose_16x16();
#endif

#ifdef USE_FGT61
	#include "fgt_m61.h"	// Mod-q utility macros
	/****** Prototypes for functions defined in fgt_m61.c are collected here: *******/
	uint64 prodq8(const uint64 x, const uint64 y);
	uint64 mul_by_3bit(const uint64 a, const uint64 x);
	uint64 rmul_modq(const uint64 x, const uint64 y);
	void cmul_modq (const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout, uint64*yout);
	void cmul_modq8(const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout, uint64*yout);
	void csqr_modq (const uint64 a0, const uint64 a1, uint64*xout, uint64*yout);
	void cmul_conj_modq(const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout1, uint64*yout1, uint64*xout2, uint64*yout2);
	void prim_root_q(const uint64 ord, uint64*root_re, uint64*root_im);
	void pow_modq(const uint64 power, const uint64 re_in, const uint64 im_in, uint64*re_out, uint64*im_out);
#endif

#ifdef __cplusplus
}
#endif

#endif	/* util_h_included */

