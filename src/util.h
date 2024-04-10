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

static const char bytestr[256][9] = {
	"00000000","00000001","00000010","00000011","00000100","00000101","00000110","00000111","00001000","00001001","00001010","00001011","00001100","00001101","00001110","00001111",
	"00010000","00010001","00010010","00010011","00010100","00010101","00010110","00010111","00011000","00011001","00011010","00011011","00011100","00011101","00011110","00011111",
	"00100000","00100001","00100010","00100011","00100100","00100101","00100110","00100111","00101000","00101001","00101010","00101011","00101100","00101101","00101110","00101111",
	"00110000","00110001","00110010","00110011","00110100","00110101","00110110","00110111","00111000","00111001","00111010","00111011","00111100","00111101","00111110","00111111",
	"01000000","01000001","01000010","01000011","01000100","01000101","01000110","01000111","01001000","01001001","01001010","01001011","01001100","01001101","01001110","01001111",
	"01010000","01010001","01010010","01010011","01010100","01010101","01010110","01010111","01011000","01011001","01011010","01011011","01011100","01011101","01011110","01011111",
	"01100000","01100001","01100010","01100011","01100100","01100101","01100110","01100111","01101000","01101001","01101010","01101011","01101100","01101101","01101110","01101111",
	"01110000","01110001","01110010","01110011","01110100","01110101","01110110","01110111","01111000","01111001","01111010","01111011","01111100","01111101","01111110","01111111",
	"10000000","10000001","10000010","10000011","10000100","10000101","10000110","10000111","10001000","10001001","10001010","10001011","10001100","10001101","10001110","10001111",
	"10010000","10010001","10010010","10010011","10010100","10010101","10010110","10010111","10011000","10011001","10011010","10011011","10011100","10011101","10011110","10011111",
	"10100000","10100001","10100010","10100011","10100100","10100101","10100110","10100111","10101000","10101001","10101010","10101011","10101100","10101101","10101110","10101111",
	"10110000","10110001","10110010","10110011","10110100","10110101","10110110","10110111","10111000","10111001","10111010","10111011","10111100","10111101","10111110","10111111",
	"11000000","11000001","11000010","11000011","11000100","11000101","11000110","11000111","11001000","11001001","11001010","11001011","11001100","11001101","11001110","11001111",
	"11010000","11010001","11010010","11010011","11010100","11010101","11010110","11010111","11011000","11011001","11011010","11011011","11011100","11011101","11011110","11011111",
	"11100000","11100001","11100010","11100011","11100100","11100101","11100110","11100111","11101000","11101001","11101010","11101011","11101100","11101101","11101110","11101111",
	"11110000","11110001","11110010","11110011","11110100","11110101","11110110","11110111","11111000","11111001","11111010","11111011","11111100","11111101","11111110","11111111"
};

/*** Const-arrays used in POPCNT and related bit-counting functions: ***/
// Bytewise POPCNT:
static const uint8 pop8[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
// Bytewise LEADZ:
static const uint8 lz8[256] = {
	8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
// Bytewise TRAILZ:
static const uint8 tz8[256] = {
	8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
	6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0
};

/* Bytewise ith_set_bit - the 8 4-bit subfields of each uint32 are numbered 1-8 and store the
location [0-7] of the 1st-8th set bit in each byte, with hex-F denoting 'no such set bit'.
Here is the simple C-snippet used to auto-generate the array data:
uint32 i,j,bit,x32;
for(i = 0; i < 256; i++) {
	x32 = bit = 0;
	for(j = 0; (j < 8) && (bit < pop8[i]); j++) {
		if((i>>j)&1) x32 += j<<(4*bit++);
	}
	for(j = bit ; j < 8; j++) {
		x32 += 0xf<<(4*j);
	}
	printf("0x%8X,",x32);
}
printf("\n");
*/
static const uint32 ith_set_bit8[256] = {
	0xFFFFFFFF,0xFFFFFFF0,0xFFFFFFF1,0xFFFFFF10,0xFFFFFFF2,0xFFFFFF20,0xFFFFFF21,0xFFFFF210,0xFFFFFFF3,0xFFFFFF30,0xFFFFFF31,0xFFFFF310,0xFFFFFF32,0xFFFFF320,0xFFFFF321,0xFFFF3210,
	0xFFFFFFF4,0xFFFFFF40,0xFFFFFF41,0xFFFFF410,0xFFFFFF42,0xFFFFF420,0xFFFFF421,0xFFFF4210,0xFFFFFF43,0xFFFFF430,0xFFFFF431,0xFFFF4310,0xFFFFF432,0xFFFF4320,0xFFFF4321,0xFFF43210,
	0xFFFFFFF5,0xFFFFFF50,0xFFFFFF51,0xFFFFF510,0xFFFFFF52,0xFFFFF520,0xFFFFF521,0xFFFF5210,0xFFFFFF53,0xFFFFF530,0xFFFFF531,0xFFFF5310,0xFFFFF532,0xFFFF5320,0xFFFF5321,0xFFF53210,
	0xFFFFFF54,0xFFFFF540,0xFFFFF541,0xFFFF5410,0xFFFFF542,0xFFFF5420,0xFFFF5421,0xFFF54210,0xFFFFF543,0xFFFF5430,0xFFFF5431,0xFFF54310,0xFFFF5432,0xFFF54320,0xFFF54321,0xFF543210,
	0xFFFFFFF6,0xFFFFFF60,0xFFFFFF61,0xFFFFF610,0xFFFFFF62,0xFFFFF620,0xFFFFF621,0xFFFF6210,0xFFFFFF63,0xFFFFF630,0xFFFFF631,0xFFFF6310,0xFFFFF632,0xFFFF6320,0xFFFF6321,0xFFF63210,
	0xFFFFFF64,0xFFFFF640,0xFFFFF641,0xFFFF6410,0xFFFFF642,0xFFFF6420,0xFFFF6421,0xFFF64210,0xFFFFF643,0xFFFF6430,0xFFFF6431,0xFFF64310,0xFFFF6432,0xFFF64320,0xFFF64321,0xFF643210,
	0xFFFFFF65,0xFFFFF650,0xFFFFF651,0xFFFF6510,0xFFFFF652,0xFFFF6520,0xFFFF6521,0xFFF65210,0xFFFFF653,0xFFFF6530,0xFFFF6531,0xFFF65310,0xFFFF6532,0xFFF65320,0xFFF65321,0xFF653210,
	0xFFFFF654,0xFFFF6540,0xFFFF6541,0xFFF65410,0xFFFF6542,0xFFF65420,0xFFF65421,0xFF654210,0xFFFF6543,0xFFF65430,0xFFF65431,0xFF654310,0xFFF65432,0xFF654320,0xFF654321,0xF6543210,
	0xFFFFFFF7,0xFFFFFF70,0xFFFFFF71,0xFFFFF710,0xFFFFFF72,0xFFFFF720,0xFFFFF721,0xFFFF7210,0xFFFFFF73,0xFFFFF730,0xFFFFF731,0xFFFF7310,0xFFFFF732,0xFFFF7320,0xFFFF7321,0xFFF73210,
	0xFFFFFF74,0xFFFFF740,0xFFFFF741,0xFFFF7410,0xFFFFF742,0xFFFF7420,0xFFFF7421,0xFFF74210,0xFFFFF743,0xFFFF7430,0xFFFF7431,0xFFF74310,0xFFFF7432,0xFFF74320,0xFFF74321,0xFF743210,
	0xFFFFFF75,0xFFFFF750,0xFFFFF751,0xFFFF7510,0xFFFFF752,0xFFFF7520,0xFFFF7521,0xFFF75210,0xFFFFF753,0xFFFF7530,0xFFFF7531,0xFFF75310,0xFFFF7532,0xFFF75320,0xFFF75321,0xFF753210,
	0xFFFFF754,0xFFFF7540,0xFFFF7541,0xFFF75410,0xFFFF7542,0xFFF75420,0xFFF75421,0xFF754210,0xFFFF7543,0xFFF75430,0xFFF75431,0xFF754310,0xFFF75432,0xFF754320,0xFF754321,0xF7543210,
	0xFFFFFF76,0xFFFFF760,0xFFFFF761,0xFFFF7610,0xFFFFF762,0xFFFF7620,0xFFFF7621,0xFFF76210,0xFFFFF763,0xFFFF7630,0xFFFF7631,0xFFF76310,0xFFFF7632,0xFFF76320,0xFFF76321,0xFF763210,
	0xFFFFF764,0xFFFF7640,0xFFFF7641,0xFFF76410,0xFFFF7642,0xFFF76420,0xFFF76421,0xFF764210,0xFFFF7643,0xFFF76430,0xFFF76431,0xFF764310,0xFFF76432,0xFF764320,0xFF764321,0xF7643210,
	0xFFFFF765,0xFFFF7650,0xFFFF7651,0xFFF76510,0xFFFF7652,0xFFF76520,0xFFF76521,0xFF765210,0xFFFF7653,0xFFF76530,0xFFF76531,0xFF765310,0xFFF76532,0xFF765320,0xFF765321,0xF7653210,
	0xFFFF7654,0xFFF76540,0xFFF76541,0xFF765410,0xFFF76542,0xFF765420,0xFF765421,0xF7654210,0xFFF76543,0xFF765430,0xFF765431,0xF7654310,0xFF765432,0xF7654320,0xF7654321,0x76543210
};

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

/* get_cpuid.c: x86 and other inline-ASM-targeted CPUs */
#ifdef CPU_IS_ARM_EABI
	int has_asimd(void);
#elif(defined(CPU_IS_X86) || defined(CPU_IS_IA64) || defined(CPU_IS_X86_64) || defined(CPU_IS_K1OM))
	void	get_cpu(void);
	uint32	has_sse  (void);
	uint32	has_sse2 (void);
	uint32	has_sse3 (void);
	uint32	has_sse3e(void);
	uint32	has_sse41(void);
	uint32	has_sse42(void);
	uint32	has_avx  (void);
	uint32	has_avx2 (void);
	uint32	has_imci512(void);	// 1st-gen Xeon Phi (Knights Ferry, Knights Corner)
	uint32	has_avx512(void);
	void	cpu_details(void);
#endif

/* util.c: */
int		mlucas_nanosleep(const struct timespec *req);
void	host_init(void);	/* This one is a wrapper for calls to the next few: */
double	get_time    (double tdiff);
char*	get_time_str(double tdiff);
void	set_stacklimit_restart(char *argv[]);
uint32	get_system_ram(void);
void	print_host_info(void);
uint32	x86_simd_mxcsr_getval(void);
uint32	x86_simd_mxcsr_setval(uint32 MXCSR_VALUE);
uint32	x86_simd_mxcsr_toggle(uint32 MXCSR_MASK);
void	set_x87_fpu_params(unsigned short FPU_MODE);	// Arg is a full FPU-control word, i.e. its value replaces the current one
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
// v20: Add simple utility to print the input string to the current-assignment logfile and/or to stderr:
void	mlucas_fprint(char*const cstr, uint32 echo_to_stderr);
double	mlucas_getOptVal(const char*fname, char*optname);

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
  #if INCLUDE_HWLOC
	#include <hwloc.h>	// Not clear to me why needed to add this redundant include here ... above include of Mdata.h should include it.
	int		num_sockets_of_core_set(hwloc_topology_t topology, int lidx_lo, int lidx_hi);
  #endif
	uint32	parseAffinityTriplet(char*istr, int hwloc_topo);
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
	// void	ASSERT	(long line, char*file, int expr, char*assert_string);
	void _ASSERT(const char*assertion, const char*file, long line, const char*func, bool expr, const char*assert_string);
#endif

#define ASSERT(expr, assert_string) _ASSERT(#expr, __FILE__, __LINE__, __func__, (expr), assert_string)

void	VAR_WARN(char *typelist, ...);

void	byte_bitstr(const uint8  byte, char*ostr);
void	ui32_bitstr(const uint32 ui32, char*ostr);
void	ui64_bitstr(const uint64 ui64, char*ostr);

uint32	reverse(uint32 i, uint32 nbits);
uint64	reverse64(uint64 i, uint32 nbits);

// 32 and 64-bit leftward circular shift, shift count n assumed unsigned < #bits-in-type:
DEV uint32 cshft32(uint32 x, uint32 n);
DEV uint64 cshft64(uint64 x, uint64 n);
// 32 and 64-bit analogs of the F90 intrinsic ISHFT function:
DEV uint32	ishft32(uint32 x, int shift);
DEV uint64	ishft64(uint64 x, int shift);

DEV uint32	trailz32(uint32 x);
DEV uint32	trailz64(uint64 x);
DEV uint32	trailz128(uint128 i);
DEV uint32	trailz192(uint192 i);
DEV uint32	trailz256(uint256 i);

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

DEV uint32	nbits32(uint32 i);
DEV uint64	nbits64(uint64 i);
DEV uint32	ibits32	(uint32 i, uint32 beg, uint32 nbits);
DEV uint64	ibits64	(uint64 i, uint32 beg, uint32 nbits);
DEV uint64	getbits64(uint64 x, uint32 src_bit_start, uint32 nbits,           uint32 tgt_bit_start);
DEV void	mvbits64 (uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start);
DEV int		pprimeF	(uint32 p, uint32 z);
DEV uint32	is_f2psp(uint32 n, uint32*idx_next_psp);
DEV uint32	is_prime(uint32 n);
DEV uint32	next_prime(uint32 n, int dir);
DEV uint32	nprimes_in_range(uint32 b1, uint32 b2);
DEV int		pprimeF64(uint64 p, uint64 z);
DEV int		isPRP64	(uint64 p);
DEV uint32	twompmodq32(uint32 p, uint32 q);	// 2^-p % q
DEV int		twopmodq32(uint32 p, uint32 q);		// (2^-p % q) == 0
DEV int		twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7);
DEV uint64	twompmodq64(uint64 p, uint64 q);	// 2^-p % q
DEV uint32	gcd32(uint32 x, uint32 y);
DEV uint64	gcd64(uint64 x, uint64 y);
DEV uint32	egcd32	(uint32 *x, uint32 *y);
DEV uint64	egcd64	(uint64 *x, uint64 *y);
DEV int		modinv32(uint32 z, uint32 n);
DEV int64	modinv64(uint64 z, uint64 n);

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

