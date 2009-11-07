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

#include "util.h"
#include "imul_macro.h"

#undef PLATFORM_SKIP_RND_CONST_ENFORCEMENT

/*********************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h:           */
/*********************************************************************************/

/* These externs yield 64-mantissa-bit register mode (bits <9:8> = 3),
with IEEE and truncating rounding mode set via bits <11:10> = 0 and 3,
respectively. The other 12 bits are identical to the MSVC defaults: */
unsigned short FPU_64RND = 0x037f, FPU_64CHOP = 0x0f7f;

const int CHAROFFSET = '0';

double RND_A, RND_B;	/* Used for fast NINT emulation; set in util.c. */

double             TWO13FLINV;	/* (double)2^13 inverse */
double TWO26FLOAT, TWO26FLINV;	/* (double)2^26 and inverse */
double TWO50FLOAT, TWO50FLINV;	/* (double)2^50 and inverse */
double TWO51FLOAT, TWO51FLINV;	/* (double)2^51 and inverse */
double TWO52FLOAT, TWO52FLINV;	/* (double)2^52 and inverse */
double TWO53FLOAT, TWO53FLINV;	/* (double)2^53 and inverse */
double TWO54FLOAT;	/* (double)2^54 */
double TWO64FLOAT, TWO64FLINV;	/* (double)2^64 and inverse */

int32 DAT_BITS, PAD_BITS;	/* Array padding parameters */

/* Fixed-size (but only necessarily constant during a given FFT-based MUL)
   base for generic FFT-based mul:

	FFT_MUL_BASE = 2^(FFT_MUL_BITS), where FFT_MUL_BITS #def'ed in Mdata.h
*/
double FFT_MUL_BASE, FFT_MUL_BASE_INV;

double ISRT2;	/* 1/sqrt(2) */

/***********************/

/*** 11/23/05: MSVC/.NET buggered things up with the second of these tables
     when each table was local to its respective calling function, so moved 'em here: ***/

	/* Table of approximate byte-inverses of 1.{byteval} is here. Since we know
	the input is in [1, 2), we know the multiplicative inverse is in (0.5, 1],
	i.e. we know that the MSB of the inverse (the one immediately right of the
	binary point) is 1. Thus we can use the hidden-bit-is-1 property of inputs
	to also gain a bit of precision in the bytewise approximate inverses, by
	neglecting the leading-order bit - since that one would get stored in the
	hidden-bit slot of the output anyway, this also makes our work easier. */
	/* Unix bc code:
	bc -l
	ibase=2
	obase=2
	d=0.00000001;
	x=1.000000001-d;
	x+=d;1/x
	{256 of these, round 10th bit into MS9, replace MSB by '0x', convert rest to hex}
	*/
	static const uint8 byte_lookup_finvest[256] = {
	0xff,0xfd,0xfb,0xf9,0xf7,0xf5,0xf3,0xf1,0xf0,0xee,0xec,0xea,0xe8,0xe6,0xe5,0xe3,
	0xe1,0xdf,0xdd,0xdc,0xda,0xd8,0xd7,0xd5,0xd3,0xd2,0xd0,0xce,0xcd,0xcb,0xc9,0xc8,
	0xc6,0xc5,0xc3,0xc2,0xc0,0xbf,0xbd,0xbc,0xba,0xb9,0xb7,0xb6,0xb4,0xb3,0xb1,0xb0,
	0xae,0xad,0xac,0xaa,0xa9,0xa7,0xa6,0xa5,0xa3,0xa2,0xa1,0x9f,0x9e,0x9d,0x9c,0x9a,
	0x99,0x98,0x96,0x95,0x94,0x93,0x91,0x90,0x8f,0x8e,0x8d,0x8b,0x8a,0x89,0x88,0x87,
	0x86,0x84,0x83,0x82,0x81,0x80,0x7f,0x7e,0x7c,0x7b,0x7a,0x79,0x78,0x77,0x76,0x75,
	0x74,0x73,0x72,0x71,0x70,0x6f,0x6e,0x6d,0x6c,0x6b,0x6a,0x69,0x68,0x67,0x66,0x65,
	0x64,0x63,0x62,0x61,0x60,0x5f,0x5e,0x5d,0x5c,0x5b,0x5a,0x59,0x58,0x58,0x57,0x56,
	0x55,0x54,0x53,0x52,0x51,0x51,0x50,0x4f,0x4e,0x4d,0x4c,0x4b,0x4b,0x4a,0x49,0x48,
	0x47,0x46,0x46,0x45,0x44,0x43,0x42,0x42,0x41,0x40,0x3f,0x3f,0x3e,0x3d,0x3c,0x3b,
	0x3b,0x3a,0x39,0x38,0x38,0x37,0x36,0x35,0x35,0x34,0x33,0x33,0x32,0x31,0x30,0x30,
	0x2f,0x2e,0x2e,0x2d,0x2c,0x2c,0x2b,0x2a,0x2a,0x29,0x28,0x28,0x27,0x26,0x26,0x25,
	0x24,0x24,0x23,0x22,0x22,0x21,0x20,0x20,0x1f,0x1e,0x1e,0x1d,0x1d,0x1c,0x1b,0x1b,
	0x1a,0x1a,0x19,0x18,0x18,0x17,0x17,0x16,0x15,0x15,0x14,0x14,0x13,0x12,0x12,0x11,
	0x11,0x10,0x10,0x0f,0x0f,0x0e,0x0d,0x0d,0x0c,0x0c,0x0b,0x0b,0x0a,0x0a,0x09,0x09,
	0x08,0x07,0x07,0x06,0x06,0x05,0x05,0x04,0x04,0x03,0x03,0x02,0x02,0x01,0x01,0x00
	};

/***********************/

	/* Table of approximate byte-inverses of 2 * 1.{byteval} is here. Since we know
	the input is in [1, 4), we know the inverse-square-rootis in (0.5, 1],
	i.e. we know that the MSB of the ISQRT (the one immediately right of the
	binary point) is 1.	We cheat a little on the 0 element of the byte table,
	since sqrt(1.000000001) really should give 0x100, not 0xff. But the
	alternative is using uint16s, which doubles the size of the table. */
	/* Unix bc code:
	bc -l
	ibase=2
	obase=2
	d=0.00000001;
	x=1.000000001-d;
	x+=d;1/sqrt(x)
	{768 of these, round 10th bit into MS9, replace MSB by '0x', convert rest to hex}
	*/
	/* Used to store MS 8 non-hidden mantissa bits. We'd need to use a 16-bit int
	to allow for the possibility of a carryout (i.e. result = 256) from rounding
	the 9th-most-significant NHB into the upper 8 (which would involve
	additional logic to handle), we instead deal with the issue of rounding
	by assuming the midpoint - e.g. if truncating to the MS 8 NHBs yields
	a certain integer in [0,255], we assume the resulting roundoff error
	is always 0.5, i.e. our precomputed 1/x values are approximations to
	the resulting midpoints. This also avoids our having to treat an input
	of 1.00000000 as a special case, since we munge that to 1.000000001,
	whose inverse is < 1.0: */
	static const uint8 byte_lookup_fisqrtest[768] = {
	0xff,0xff,0xfe,0xfd,0xfc,0xfb,0xfa,0xf9,0xf8,0xf7,0xf6,0xf5,0xf4,0xf3,0xf2,0xf1,
	0xf0,0xef,0xee,0xee,0xed,0xec,0xeb,0xea,0xe9,0xe8,0xe7,0xe7,0xe6,0xe5,0xe4,0xe3,
	0xe2,0xe1,0xe1,0xe0,0xdf,0xde,0xdd,0xdd,0xdc,0xdb,0xda,0xd9,0xd9,0xd8,0xd7,0xd6,
	0xd5,0xd5,0xd4,0xd3,0xd2,0xd2,0xd1,0xd0,0xcf,0xcf,0xce,0xcd,0xcc,0xcc,0xcb,0xca,
	0xca,0xc9,0xc8,0xc7,0xc7,0xc6,0xc5,0xc5,0xc4,0xc3,0xc3,0xc2,0xc1,0xc1,0xc0,0xbf,
	0xbf,0xbe,0xbd,0xbd,0xbc,0xbb,0xbb,0xba,0xb9,0xb9,0xb8,0xb7,0xb7,0xb6,0xb6,0xb5,
	0xb4,0xb4,0xb3,0xb2,0xb2,0xb1,0xb1,0xb0,0xaf,0xaf,0xae,0xae,0xad,0xac,0xac,0xab,
	0xab,0xaa,0xaa,0xa9,0xa8,0xa8,0xa7,0xa7,0xa6,0xa6,0xa5,0xa5,0xa4,0xa3,0xa3,0xa2,
	0xa2,0xa1,0xa1,0xa0,0xa0,0x9f,0x9f,0x9e,0x9d,0x9d,0x9c,0x9c,0x9b,0x9b,0x9a,0x9a,
	0x99,0x99,0x98,0x98,0x97,0x97,0x96,0x96,0x95,0x95,0x94,0x94,0x93,0x93,0x92,0x92,
	0x91,0x91,0x90,0x90,0x8f,0x8f,0x8f,0x8e,0x8e,0x8d,0x8d,0x8c,0x8c,0x8b,0x8b,0x8a,
	0x8a,0x89,0x89,0x89,0x88,0x88,0x87,0x87,0x86,0x86,0x85,0x85,0x85,0x84,0x84,0x83,
	0x83,0x82,0x82,0x82,0x81,0x81,0x80,0x80,0x7f,0x7f,0x7f,0x7e,0x7e,0x7d,0x7d,0x7d,
	0x7c,0x7c,0x7b,0x7b,0x7a,0x7a,0x7a,0x79,0x79,0x78,0x78,0x78,0x77,0x77,0x76,0x76,
	0x76,0x75,0x75,0x75,0x74,0x74,0x73,0x73,0x73,0x72,0x72,0x72,0x71,0x71,0x70,0x70,
	0x70,0x6f,0x6f,0x6f,0x6e,0x6e,0x6d,0x6d,0x6d,0x6c,0x6c,0x6c,0x6b,0x6b,0x6b,0x6a,
	0x6a,0x6a,0x69,0x69,0x68,0x68,0x68,0x67,0x67,0x67,0x66,0x66,0x66,0x65,0x65,0x65,
	0x64,0x64,0x64,0x63,0x63,0x63,0x62,0x62,0x62,0x61,0x61,0x61,0x60,0x60,0x60,0x5f,
	0x5f,0x5f,0x5e,0x5e,0x5e,0x5d,0x5d,0x5d,0x5d,0x5c,0x5c,0x5c,0x5b,0x5b,0x5b,0x5a,
	0x5a,0x5a,0x59,0x59,0x59,0x58,0x58,0x58,0x58,0x57,0x57,0x57,0x56,0x56,0x56,0x55,
	0x55,0x55,0x55,0x54,0x54,0x54,0x53,0x53,0x53,0x53,0x52,0x52,0x52,0x51,0x51,0x51,
	0x51,0x50,0x50,0x50,0x4f,0x4f,0x4f,0x4f,0x4e,0x4e,0x4e,0x4d,0x4d,0x4d,0x4d,0x4c,
	0x4c,0x4c,0x4c,0x4b,0x4b,0x4b,0x4a,0x4a,0x4a,0x4a,0x49,0x49,0x49,0x49,0x48,0x48,
	0x48,0x48,0x47,0x47,0x47,0x47,0x46,0x46,0x46,0x45,0x45,0x45,0x45,0x44,0x44,0x44,
	0x44,0x43,0x43,0x43,0x43,0x42,0x42,0x42,0x42,0x41,0x41,0x41,0x41,0x40,0x40,0x40,
	0x40,0x3f,0x3f,0x3f,0x3f,0x3f,0x3e,0x3e,0x3e,0x3e,0x3d,0x3d,0x3d,0x3d,0x3c,0x3c,
	0x3c,0x3c,0x3b,0x3b,0x3b,0x3b,0x3a,0x3a,0x3a,0x3a,0x3a,0x39,0x39,0x39,0x39,0x38,
	0x38,0x38,0x38,0x38,0x37,0x37,0x37,0x37,0x36,0x36,0x36,0x36,0x36,0x35,0x35,0x35,
	0x35,0x34,0x34,0x34,0x34,0x34,0x33,0x33,0x33,0x33,0x32,0x32,0x32,0x32,0x32,0x31,
	0x31,0x31,0x31,0x31,0x30,0x30,0x30,0x30,0x30,0x2f,0x2f,0x2f,0x2f,0x2e,0x2e,0x2e,
	0x2e,0x2e,0x2d,0x2d,0x2d,0x2d,0x2d,0x2c,0x2c,0x2c,0x2c,0x2c,0x2b,0x2b,0x2b,0x2b,
	0x2b,0x2a,0x2a,0x2a,0x2a,0x2a,0x29,0x29,0x29,0x29,0x29,0x28,0x28,0x28,0x28,0x28,
	0x28,0x27,0x27,0x27,0x27,0x27,0x26,0x26,0x26,0x26,0x26,0x25,0x25,0x25,0x25,0x25,
	0x24,0x24,0x24,0x24,0x24,0x24,0x23,0x23,0x23,0x23,0x23,0x22,0x22,0x22,0x22,0x22,
	0x22,0x21,0x21,0x21,0x21,0x21,0x20,0x20,0x20,0x20,0x20,0x20,0x1f,0x1f,0x1f,0x1f,
	0x1f,0x1f,0x1e,0x1e,0x1e,0x1e,0x1e,0x1d,0x1d,0x1d,0x1d,0x1d,0x1d,0x1c,0x1c,0x1c,
	0x1c,0x1c,0x1c,0x1b,0x1b,0x1b,0x1b,0x1b,0x1b,0x1a,0x1a,0x1a,0x1a,0x1a,0x1a,0x19,
	0x19,0x19,0x19,0x19,0x19,0x18,0x18,0x18,0x18,0x18,0x18,0x17,0x17,0x17,0x17,0x17,
	0x17,0x16,0x16,0x16,0x16,0x16,0x16,0x15,0x15,0x15,0x15,0x15,0x15,0x15,0x14,0x14,
	0x14,0x14,0x14,0x14,0x13,0x13,0x13,0x13,0x13,0x13,0x13,0x12,0x12,0x12,0x12,0x12,
	0x12,0x11,0x11,0x11,0x11,0x11,0x11,0x11,0x10,0x10,0x10,0x10,0x10,0x10,0x0f,0x0f,
	0x0f,0x0f,0x0f,0x0f,0x0f,0x0e,0x0e,0x0e,0x0e,0x0e,0x0e,0x0e,0x0d,0x0d,0x0d,0x0d,
	0x0d,0x0d,0x0d,0x0c,0x0c,0x0c,0x0c,0x0c,0x0c,0x0c,0x0b,0x0b,0x0b,0x0b,0x0b,0x0b,
	0x0b,0x0a,0x0a,0x0a,0x0a,0x0a,0x0a,0x0a,0x09,0x09,0x09,0x09,0x09,0x09,0x09,0x08,
	0x08,0x08,0x08,0x08,0x08,0x08,0x08,0x07,0x07,0x07,0x07,0x07,0x07,0x07,0x06,0x06,
	0x06,0x06,0x06,0x06,0x06,0x05,0x05,0x05,0x05,0x05,0x05,0x05,0x05,0x04,0x04,0x04,
	0x04,0x04,0x04,0x04,0x04,0x03,0x03,0x03,0x03,0x03,0x03,0x03,0x02,0x02,0x02,0x02,
	0x02,0x02,0x02,0x02,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x00,0x00,0x00,0x00
	};


/********FILE-ACCESS-RELATED:**********/

int		file_valid(FILE*fp)
{
	if(fp && !ferror(fp) && !feof(fp))
		return TRUE;
	else
		return FALSE;
}

/* Apparently each implementation uses its own names for these file-internals flags -
TODO: need a portable way of checking file-access mode here.
int		file_valid_for_read (FILE*fp)
{
	blah... (similar to below)
}

int		file_valid_for_write(FILE*fp)
{
	if(fp && (fp->flag & (_EOF | _ERR) == 0) && (fp->flag & _WRITE != 0) )
		return TRUE;
	else
		return FALSE;
}
*/
/******************/

/* Print key platform info, (on x86) set FPU mode, do some basic self-tests: */
void host_init(void)
{
	double dbl;
	print_host_info();
	set_x87_fpu_params(FPU_64RND);
	check_nbits_in_types();

	/* Test wide-mul routines: */
	printf("INFO: testing IMUL routines...\n");
	ASSERT(HERE, test_mul() == 0, "test_mul() returns nonzero!");

	/* Check qfloat routines (this call is also needed to init various qfloat global constants): */
	printf("INFO: testing qfloat routines...\n");
	qtest();

	/* Use qfloat routines to set the global floating-point constant 1/sqrt(2): */
	ISRT2 = qfdbl(qfmul_pow2(QSQRT2, -1));		/* 1/sqrt2	*/

	/* Various useful precomputed powers of 2 in floating-double form: */
	TWO26FLOAT = (double)0x04000000;
	TWO26FLINV = 1.0/TWO26FLOAT;
	dbl = qfdbl(qfmul_pow2(QONE, -26));
	ASSERT(HERE, TWO26FLINV == dbl, "TWO26FLINV!");
	TWO13FLINV = qfdbl(qfmul_pow2(QONE, -13));

	TWO50FLOAT = (double)1.0*0x01000000*0x04000000;
	TWO50FLINV = 1.0/TWO50FLOAT;

	TWO51FLOAT = (double)1.0*0x02000000*0x04000000;
	TWO51FLINV = 1.0/TWO51FLOAT;

	TWO52FLOAT = (double)1.0*0x04000000*0x04000000;
	TWO52FLINV = 1.0/TWO52FLOAT;

	TWO53FLOAT = (double)1.0*0x08000000*0x04000000;
	TWO53FLINV = 1.0/TWO53FLOAT;

	TWO54FLOAT = (double)1.0*0x08000000*0x08000000;	/* (double)2^54 */

	TWO64FLOAT = (double)4.0*0x80000000*0x80000000;
	TWO64FLINV = 1.0/TWO64FLOAT;

#ifndef USE_SSE2
//	test_fft_radix();
//	exit(0);
#endif
}

/***The following 3 routines MUST BE CALLED IN THE SAME ORDER AS IN host_init()!***/

void print_host_info(void)
{
#if EWM_DEBUG
	printf("INFO: Program compiled with debugging diagnostics ON.\n");
#endif

	printf("CPU Family = %s, OS = %s, %2d-bit Version, compiled with %s, Version %s.\n", CPU_NAME, OS_NAME, OS_BITS, COMPILER_NAME, COMPILER_VERSION);

#if(defined(CPU_TYPE_IA32) || defined(CPU_TYPE_AMD64) || defined(CPU_TYPE_X86_64))

//	get_cpu();

/* Enable this call to get gory details: */
	#if(1)
		/* if(1) --> if(0) Enables section below */
	#elif(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_ICC))
		cpu_details();
	#endif

  #if(defined(USE_SSE3))

	if(has_sse3())
	{
		printf("INFO: CPU supports SSE3 instruction set.\n");
		ASSERT(HERE, has_sse2(), "SSE3, but not SSE2 support detected on this CPU! Check get_cpuid functionality.");
	}
	else if(has_sse2())
	{
		ASSERT(HERE, 0, "#define USE_SSE3 invoked but only SSE2 support detected on this CPU! Try rebuilding with USE_SSE2 instead.\n");
	}
	else
	{
		ASSERT(HERE, 0, "#define USE_SSE3 invoked but no SSE2/3 support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #elif(defined(USE_SSE2))

	if(has_sse3())
		printf("INFO: CPU supports SSE3 instruction set, but only USE_SSE2 used in current Mlucas version.\n");
	else if(has_sse2())
		printf("INFO: CPU supports SSE2 instruction set.\n");
	else
	{
		ASSERT(HERE, 0, "#define USE_SSE2 invoked but no SSE2 support detected on this CPU! Check get_cpuid functionality and CPU type.\n");
	}

  #else

	if(has_sse3())
		printf("Info: CPU supports SSE2/3 instruction set, but USE_SSE2/3 not invoked for build.\n");
	else if(has_sse2())
		printf("Info: CPU supports SSE2 instruction set, but USE_SSE2 not invoked for build.\n");

  #endif

#endif

#if PFETCH
	printf("INFO: Using prefetch.\n");
#endif

#ifdef MUL_LOHI64_SUBROUTINE
	printf("INFO: Using subroutine form of MUL_LOHI64.\n");
#else
	printf("INFO: Using inline-macro form of MUL_LOHI64.\n");
#endif

#ifdef USE_FMADD
	printf("INFO: Using FMADD-based 100-bit modmul routines for factoring.\n");
#elif(defined(USE_FLOAT))
	printf("INFO: Using floating-double-based modmul routines for factoring.\n");
#endif

}

/*
80x87 FPU Control Register stuff: FPUCTRL is a 16-bit register. Counting from 0-15,
the 4 bits of interest in the present application are bits 8 and 9 which control the
FPU precision, and bits 10 and 11, which control the FPU rounding mode:

	Bits <9:8>:	Value	Precision		<11:10>:	Value	Rounding Mode
				-----	---------					-----	---------------
				00		24 bits						00		==> nearest (i.e. IEEE)
				01		Reserved					01		==> -oo
				10		53 bits						10		==> +oo
				11		64 bits						11		==> 0 (i.e. truncate)

For the purpose of completeness, the other FPU control bits are as follows
( Adapted from http://maven.smith.edu/~thiebaut/ArtOfAssembly/CH14/CH14-3.html ):

	Bits zero through five are the exception masks. These are similar to the
	interrupt enable bit in the 80x86's flags register. If these bits contain
	a one, the corresponding condition is ignored by the 80x87 FPU. However,
	if any bit contains zero, and the corresponding condition occurs, then
	the FPU immediately generates an interrupt so the program can handle the
	degenerate condition.

	Bit zero corresponds to an invalid operation error. This generally occurs
	as the result of a programming error. Problem which raise the invalid
	operation exception include pushing more than eight items onto the stack
	or attempting to pop an item off an empty stack, taking the square root
	of a negative number, or loading a non-empty register.

	Bit one masks the denormalized interrupt which occurs whenever you try to
	manipulate denormalized values. Denormalized values generally occur when
	you load arbitrary extended precision values into the FPU or work with
	very small numbers just beyond the range of the FPU's capabilities.
	Normally, you would probably not enable this exception.

	Bit two masks the zero divide exception. If this bit contains zero, the FPU
	will generate an interrupt if you attempt to divide a nonzero value by zero.
	If you do not enable the zero division exception, the FPU will produce NaN
	(not a number) whenever you perform a zero division.

	Bit three masks the overflow exception. The FPU will raise the overflow
	exception if a calculation overflows or if you attempt to store a value
	which is too large to fit into a destination operand (e.g., storing a large
	extended precision value into a single precision variable).

	Bit four, if set, masks the underflow exception. Underflow occurs when the
	result is too small to fit in the desintation operand. Like overflow, this
	exception can occur whenever you store a small extended precision value into
	a smaller variable (single or double precision) or when the result of
	computation is too small for extended precision.

	Bit five controls whether the precision exception can occur. A precision
	exception occurs whenever the FPU produces an imprecise result, generally
	the result of an internal rounding operation. Although many operations will
	produce an exact result, many more will not. For example, dividing one by
	ten will produce an inexact result. Therefore, this bit is usually one since
	inexact results are very common.

	Bits six and thirteen through fifteen in the control register are currently
	undefined and reserved for future use. Bit seven is the interrupt enable mask,
	but it is only active on the 8087 FPU; a zero in this bit enables 8087
	interrupts and a one disables FPU interrupts.

	The 80x87 provides two instructions, FLDCW (load control word) and FSTCW (store
	control word), that let you load and store the contents of the control register.
	The single operand to these instructions must be a 16 bit memory location. The
	FLDCW instruction loads the control register from the specified memory location,
	FSTCW stores the control register into the specified memory location.
*/

/* On the Itanium, the FPU Control Register is 64-bit, and divided into 5 subfields:

	- A single 6-bit subfield (bits <5:0> containing the 6 trap disable bits;

	- A quartet of 13-bit Floating-Point Status bitfields (numbered 0 through 3),
	all having the same format
	(A summary the precise meaning of the bits seems exorbitantly difficult to find -
	apparently Intel doesn't want users monkeying with these),
	but with only FPSB 0 accessible to the user, via the following intrinsics
	defined in the <xmmintrin.h> header file:

	extern unsigned __int64 _mm_getfpsr(void);
	extern void _mm_setfpsr(unsigned __int64);	(maps to the Maps to the fsetc.sf0 hardware instruction)

	Using ICC v9.0, the getfpsr function returns 0x8A70037F,
	The least-signficant 12 bits of which are identical to the
	default ones we use for the x87, so whatever the internal order
	of the bitfields, it appears that the _mm_getfpsr and _mm_setfpsr
	intrinsics munge things so the user-settable FPUCTRL bits are in
	the same order and locations as on x87.
*/
void set_x87_fpu_params(unsigned short FPU_MODE)
{
	#ifdef CPU_TYPE_IA64
		int64 FPUCTRL;
	#else
		unsigned short FPUCTRL;
	#endif

	ASSERT(HERE, (FPU_64RND == FPU_MODE) || (FPU_64CHOP == 0x0f7f), "Illegal value of FPU_MODE");

	/* Copy the FPU control word set by the compiler to a local variable
	(mainly so can see what the compiler sets), then overwrite with one of the above:
	*/
#if(defined(COMPILER_TYPE_ICC) && defined(CPU_TYPE_IA32))
	/**/
#elif(FP_MANTISSA_BITS_DOUBLE == 64)/* (defined(CPU_TYPE_IA32) || defined(CPU_TYPE_IA64))) */

	#ifdef CPU_TYPE_IA64

		#ifndef COMPILER_TYPE_ICC
			#error unsupported compiler type for ia64!
		#endif
		FPUCTRL = _mm_getfpsr();
		info_x87_fpu_ctrl(FPUCTRL);
		/* Just use the same full-16-bit constant on all platforms, to ensure that there
		are no compiler-based differences in the other 12 bits, either: */
		FPUCTRL = (FPUCTRL & 0x0000) + FPU_MODE;
		_mm_setfpsr(FPUCTRL);

	#else

		#if(defined(COMPILER_TYPE_MWERKS) || defined(COMPILER_TYPE_MSVC))

			__asm	fstcw	FPUCTRL
			/*_controlfp(_PC_64, _MCW_PC);*/

		#elif(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

			__asm__ volatile ("fstcw %0" : "=m" (FPUCTRL) );

		#else
			#error unsupported compiler type for x87!
		#endif
		info_x87_fpu_ctrl((uint64)FPUCTRL);

		/* Irrespective of what values the compiler set for bitfields <9:8> and <11:10>
		of FPUCTRL, set <9:8> = 11 and <11:10> = 00 to get the full 64-mantissa-bit
		precision available on the x87 and to ensure IEEE rounding mode, respectively:
		*/
	  #if 1
		/* Just use the same full-16-bit constant on all platforms, to ensure that there
		are no compiler-based differences in the other 12 bits, either: */
		FPUCTRL = (FPUCTRL & 0x0000) + FPU_MODE;
	  #else
		***obsolete:***
		FPUCTRL &= 0xf0ff;	/* Clear bits 8:11... */
		FPUCTRL |= 0x0300;	/* And set them to the desired value. */
	  #endif

		/* ...and then reload the FPU control word for the changes to take effect. */
		#if(defined(COMPILER_TYPE_MWERKS) || defined(COMPILER_TYPE_MSVC))
			__asm	fldcw	FPUCTRL
		#elif(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
			__asm__ volatile ("fldcw %0" :: "m" (FPUCTRL) );
		#endif

	#endif	/* endif(CPU_TYPE_IA32...) */

#endif	/* endif(FP_MANTISSA_BITS_DOUBLE) */
}

void info_x87_fpu_ctrl(uint64 FPUCTRL)
{
#if EWM_DEBUG
	printf("INFO: x87 FPU Control Word = %16X.\n", (uint64)FPUCTRL);
#endif

	/* Check bits <9:8>, and warn if the compiler isn't specifying 64-bit precision: */
	switch ((FPUCTRL >> 8) & 0x3) {
	case 0x3:
		break;
	case 0x2:
		printf("INFO: compiler sets x87 FPU to 53-bit mantissa mode. Overriding...Setting to 64-bit mode.\n");
		break;
	case 0x0:
		printf("INFO: compiler sets x87 FPU to 24-bit mantissa mode. Overriding...Setting to 64-bit mode.\n");
		break;
	default:
		printf("INFO: compiler sets x87 FPU to unknown precision. Overriding...Setting to 64-bit mode.\n");
	}

	/* Check bits <11:10>, and warn if the compiler isn't specifying 64-bit precision: */
	switch ((FPUCTRL >> 10) & 0x3) {
	case 0x0:
		break;
	case 0x1:
		printf("INFO: compiler sets x87 FPU to [round ==> -oo] rounding mode. Overriding...Setting to [round ==> nearest].\n");
		break;
	case 0x2:
		printf("INFO: compiler sets x87 FPU to [round ==> +oo] rounding mode. Overriding...Setting to [round ==> nearest].\n");
		break;
	case 0x3:
		printf("INFO: compiler sets x87 FPU to [round ==> 0] (truncate) rounding mode. Overriding...Setting to [round ==> nearest].\n");
		break;
	default:
		ASSERT(HERE, 0,"0");
	}
}

/******* DEFINE GLOBALS AND TEST RND-CONST FAST-NINT: *******/
void check_nbits_in_types(void)
{
	uint32 i;
	double ftmp, fran, ferr, finv, fsrt;
	double tpi = 3.1415926535897932384;
	double ln2 = LOG2;

	/* Make sure TRUE and FALSE behave as required: */
	ASSERT(HERE, !FALSE && TRUE, "TRUE and FALSE do not behave as required in check_nbits_in_types");

/* Check lengths of basic data types: */
    ASSERT(HERE, sizeof( int8 ) == 1          , "sizeof( int8 ) != 1          ");
    ASSERT(HERE, sizeof(uint8 ) == 1          , "sizeof(uint8 ) != 1          ");
    ASSERT(HERE, sizeof( int16) == 2          , "sizeof( int16) != 2          ");
    ASSERT(HERE, sizeof(uint16) == 2          , "sizeof(uint16) != 2          ");
    ASSERT(HERE, sizeof( int32) == 4          , "sizeof( int32) != 4          ");
    ASSERT(HERE, sizeof(uint32) == 4          , "sizeof(uint32) != 4          ");
    ASSERT(HERE, sizeof( int64) == 8          , "sizeof( int64) != 8          ");
    ASSERT(HERE, sizeof(uint64) == 8          , "sizeof(uint64) != 8          ");
    ASSERT(HERE, sizeof(long) == sizeof(void*), "sizeof(long) != sizeof(void*)");    /* ALIGN_DOUBLES assumes this. */

/* AltiVec vector types: */
#if(CPU_HAS_ALTIVEC || CPU_IS_CELL)
	ASSERT(HERE, sizeof(vec_uint8X16) == 16 , "sizeof(vec_uint8X16) != 16 ");
	ASSERT(HERE, sizeof(vec_uint16X8) == 16 , "sizeof(vec_uint16x8) != 16 ");
	ASSERT(HERE, sizeof(vec_uint32X4) == 16 , "sizeof(vec_uint32x4) != 16 ");
#endif

/* We move this into a function defined in a separate file in an
attempt to avoid too-clever compilers realizing that RND_A and RND_B
have the same value and optimizing the +RND_A-RND_B sequences below away:
*/
	get_fp_rnd_const(&RND_A, &RND_B);

	/* Attempted workaround for the ICC v9.0 weirdness here: */
	if(DNINT(tpi) != 3.0 || DNINT(ln2) != 1.0)
	{
		if(FP_MANTISSA_BITS_DOUBLE == 64)
		{
			sprintf(cbuf, "WARN: 64-bit rounding constant not behaving as expected - trying 53-bit version.\n");
			DBG_WARN(HERE, cbuf, "", 0);
			/*
			If there are any platforms which fail here not because they need
			the 53-bit form of the rounding constant but rather because they inevitably
			optimize away the rounding add/sub sequence, but which appear to do the
			right thing in the actual ditN_cy_dif1 routines, indicate in the #if here:
			*/
		#if(1)
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
		#else
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
			RND_A = 3.0*0x4000000*0x2000000;
			RND_B =12.0*0x2000000*0x1000000;
		#endif
		}
		else
		{
			sprintf(cbuf, "WARN: 53-bit rounding constant not behaving as expected - trying 64-bit version.\n");
			DBG_WARN(HERE, cbuf, "", 0);
			/*
			If there are any platforms which fail here not because they need
			the 53-bit form of the rounding constant but rather because they inevitably
			optimize away the rounding add/sub sequence, but which appear to do the
			right thing in the actual ditN_cy_dif1 routines, indicate in the #if here:
			*/
		#if(1)
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
		#else
		  #define PLATFORM_SKIP_RND_CONST_ENFORCEMENT
			RND_A = 3.0*0x4000000*0x2000000*0x800;
			RND_B =12.0*0x2000000*0x1000000*0x800;
		#endif
		}
	}

#ifdef PLATFORM_SKIP_RND_CONST_ENFORCEMENT

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, pi  = %20.3f,  DNINT(pi ) = %20.3f\n", RND_A, tpi, DNINT(tpi));
	if(DNINT(tpi) != 3.0)
		DBG_WARN(HERE, cbuf, "", TRUE);

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, ln2 = %20.3f,  DNINT(ln2) = %20.3f\n", RND_A, ln2, DNINT(ln2));
	if(DNINT(ln2) != 1.0)
		DBG_WARN(HERE, cbuf, "", TRUE);

#else

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, pi  = %20.3f,  DNINT(pi ) = %20.3f\n", RND_A, tpi, DNINT(tpi));
	ASSERT(HERE, DNINT(tpi) == 3.0, cbuf);

	sprintf(cbuf,"in check_nbits_in_types: RND_A = %20.3f, ln2 = %20.3f,  DNINT(ln2) = %20.3f\n", RND_A, ln2, DNINT(ln2));
	ASSERT(HERE, DNINT(ln2) == 1.0, cbuf);

#endif

	/* We typically need more information re. the FFT-mul params before being
	able to inteligently set the anti-thrashing array-padding params, so set = -1
	(which is taken to mean uninitialized) here:
	*/
	DAT_BITS = PAD_BITS = (int32)0xffffffff;

	FFT_MUL_BASE = (double)((uint64)1 << FFT_MUL_BITS);

/* Intend to relax this later to allow powers of 2 as large as 2^54: */
ASSERT(HERE, ((uint64)FFT_MUL_BASE >> 16) == 1, "util.c: FFT_MUL_BASE != 2^16");

	ASSERT(HERE, trailz64((uint64)FFT_MUL_BASE) == FFT_MUL_BITS, "mi64_cvt_double_uint64: trailz64((uint64)FFT_MUL_BASE) != FFT_MUL_BITS");
	ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_double_uint64: FFT_MUL_BASE not pure-integer!");
	ASSERT(HERE, FFT_MUL_BASE < 1.0*0x8000000*0x8000000, "mi64_cvt_double_uint64: FFT_MUL_BASE >= maximum allowed value of 2^54!");

	FFT_MUL_BASE_INV = 1.0/FFT_MUL_BASE;

	/* Test approximate 1/x and 1/sqrt(x) routines: */
	ftmp = finvest(1.5,  8);	ASSERT(HERE, fabs(ftmp - 0.666666666666667) < 4e-03, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(1.5, 53);	ASSERT(HERE, fabs(ftmp - 0.666666666666667) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(1.0, 53);	ASSERT(HERE, fabs(ftmp - 1.000000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(2.0, 53);	ASSERT(HERE, fabs(ftmp - 0.500000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(0.5, 53);	ASSERT(HERE, fabs(ftmp - 2.000000000000000) < 1e-14, "Unacceptable level of error in finvest() call!");
	ftmp = finvest(.75, 53);	ASSERT(HERE, fabs(ftmp - 1.333333333333333) < 1e-14, "Unacceptable level of error in finvest() call!");
	/* Try some large and small inputs: */
	ftmp = finvest(3.141592653589793e+15, 53);	ASSERT(HERE, fabs(ftmp - 3.183098861837907e-16) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = finvest(3.183098861837907e-16, 53);	ASSERT(HERE, fabs(ftmp - 3.141592653589793e+15) < 1e+00, "Unacceptable level of error in fisqrtest() call!");

	ftmp = fisqrtest(1.5,  8);	ASSERT(HERE, fabs(ftmp - 0.816496580927726) < 1e-3 , "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(1.5, 53);	ASSERT(HERE, fabs(ftmp - 0.816496580927726) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(1.0, 53);	ASSERT(HERE, fabs(ftmp - 1.000000000000000) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(2.0, 53);	ASSERT(HERE, fabs(ftmp - 0.707106781186548) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(0.5, 53);	ASSERT(HERE, fabs(ftmp - 1.414213562373095) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(0.3, 53);	ASSERT(HERE, fabs(ftmp - 1.825741858350554) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(.25, 53);	ASSERT(HERE, fabs(ftmp - 2.000000000000000) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(.75, 53);	ASSERT(HERE, fabs(ftmp - 1.154700538379251) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(3.0, 53);	ASSERT(HERE, fabs(ftmp - 0.577350269189626) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
	/* Try some large and small inputs: */
	ftmp = fisqrtest(3.141592653589793e+15, 53);	ASSERT(HERE, fabs(ftmp - 1.784124116152771e-08) < 1e-22, "Unacceptable level of error in fisqrtest() call!");
	ftmp = fisqrtest(3.183098861837907e-16, 53);	ASSERT(HERE, fabs(ftmp - 5.604991216397928e+07) < 1e-07, "Unacceptable level of error in fisqrtest() call!");

	/* Now do a whole mess of 'em: */
	rng_isaac_init(TRUE);
	for(i = 0; i < 100000; i++)
	{
		fran = rng_isaac_rand_double();
		fran = fabs(fran);
		if(fran > 0.0)
		{
			ftmp = finvest  (fran, 53);
			finv = 1.0/fran;
			ferr = (ftmp - finv)/(ftmp + finv);
			ASSERT(HERE, fabs(ferr) < 1e-14, "Unacceptable level of error in finvest  () call!");

			ftmp = fisqrtest(fran, 53);
			fsrt = 1.0/sqrt(fran);
			ferr = (ftmp - fsrt)/(ftmp + fsrt);
			ASSERT(HERE, fabs(ferr) < 1e-14, "Unacceptable level of error in fisqrtest() call!");
		}

		fran = rng_isaac_rand_double_norm_pos();
		if(fran < 0.0 || fran >= 1.0)
		{
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pos returns illegal value outside [0, 1): i = %d, %e\n", i,fran);
			ASSERT(HERE, 0, cbuf);
		}

		fran = rng_isaac_rand_double_norm_pm1();
		if(fabs(fran) >= 1.0)
		{
			sprintf(cbuf, "check_nbits_in_types: rng_isaac_rand_double_norm_pm1 returns illegal value outside (-1,+1): i = %d, %e\n", i, fran);
			ASSERT(HERE, 0, cbuf);
		}
	}

#ifdef USE_FMADD
	/* Test any FMADD-based routines, if def'd: */
	printf("Testing MUL50x50 routines...\n");
	test_mul50x50();
#endif

	return;
}

/********************/

/*
Originally wrote this as part of a workaround for the the DEC Unix V4.0 real*16 sincos bug,
which caused incorrect real*16 sincos results when the argument theta = m*pi/4 +- pi/512, m integer.
The more general purpose is to compare a pair of computed sin(theta), cos(theta) values
(either truncated-to-double extended-real ones, or fast-math-library doubles, or whatever)
to the results returned by calls to the standard math library sincos on the host system.
If a discrepancy between the input sincos data and the standard-library ones is detected
which exceeds some threshold (we use 1e-10 at present), print a warning and replace the
inputs with the latter. (If you instead suspect the latter of being the problematic ones,
you'll need to make the appropriate modifications).

Returns the largest absolute difference between real*8 and real*16 sin(theta) and cos(theta)
for the calling values, e.g. if the user wants to save the maximum discrepancy over a series of calls.
*/
#define	USE_DOUBLE_DEFAULT	0

double 	errprint_sincos(double *x, double *y, double theta)
{
	double tmp, adiff, maxdiff = 0.0;

	tmp = cos(theta);
	adiff = fabs(*x - tmp);
	if(adiff > maxdiff);
		maxdiff = adiff;
	if(adiff > 1e-10)
	{
		fprintf(stderr, "WARNING: real*16 sine error : theta = %20.15f, long double = %20.15f, double = %20.15f", theta,*x,tmp);
	#if USE_DOUBLE_DEFAULT
		fprintf(stderr, " ... using double-precision values instead.\n");
		*x = tmp;
	#else
		fprintf(stderr, "\n");
	#endif
	}

	tmp = sin(theta);
	adiff = fabs(*y - tmp);
	if(adiff > maxdiff);
		maxdiff = adiff;
	if(adiff > 1e-10)
	{
		fprintf(stderr, "WARNING: real*16 sine error : theta = %20.15f, long double = %20.15f, double = %20.15f", theta,*y,tmp);
	#if USE_DOUBLE_DEFAULT
		fprintf(stderr, " ... using double-precision values instead.\n");
		*y = tmp;
	#else
		fprintf(stderr, "\n");
	#endif
	}

	return adiff;
}

/**********************************/
/******* INFO, WARN ASSERT ********/
/**********************************/

void DEBUG_INFO(long line, char*file, char*info_string, char*info_file, int copy2stderr)
{
	FILE *fp = 0x0;

	if(STRNEQ(info_file, ""))
	{
		fp = fopen(info_file,"a");
		if(!fp)
			fprintf(stderr,"WARNING: unable to open file %s in call to DEBUG_INFO.\n", info_file);
	}

	if(fp)
	{
		fprintf(fp,"INFO: At line %lu of file %s:\n", line, file);
		fprintf(fp,"%s\n", info_string);
		fclose(fp); fp = 0x0;
	}

	if(copy2stderr || !fp)
	{
		fprintf(stderr,"INFO: At line %lu of file %s:\n", line, file);
		fprintf(stderr,"%s\n", info_string);
		fflush(stderr);
	}
}

void DEBUG_WARN(long line, char*file, char*warn_string, char*warn_file, int copy2stderr)
{
	FILE *fp = 0x0;

	if(STRNEQ(warn_file, ""))
	{
		fp = fopen(warn_file,"a");
		if(!fp)
			fprintf(stderr,"WARNING: unable to open file %s in call to DBG_WARN.\n", warn_file);
	}

	if(fp)
	{
		fprintf(fp,"WARN: At line %lu of file %s:\n", line, file);
		fprintf(fp,"%s\n", warn_string);
		fclose(fp); fp = 0x0;
	}

	if(copy2stderr || !fp)
	{
		fprintf(stderr,"WARN: At line %lu of file %s:\n", line, file);
		fprintf(stderr,"%s\n", warn_string);
		fflush(stderr);
	}
}

#ifdef	USE_C99
void ASSERT(char*func, long line, char*file, int expr, char*assert_string)
{
	/* Define a convenient spot to set a breakpoint: */
	if(!expr)
	{
		fprintf(stderr,"ERROR: Function %s, at line %lu of file %s\n", func, line, file);
		fprintf(stderr,"Assertion failed: %s\n", assert_string);
		/* Flush all output streams prior to asserting.
		We replace the original assert(0) call with an exit(EXIT_FAILURE),
		since some compilers seem to like to optimize away assertions. */
		fflush(NULL);
		exit(EXIT_FAILURE);
	}
}
#else
void ASSERT(long line, char*file, int expr, char*assert_string)
{
	/* Define a convenient spot to set a breakpoint: */
	if(!expr)
	{
		fprintf(stderr,"ERROR: at line %lu of file %s\n", line, file);
		fprintf(stderr,"Assertion failed: %s\n", assert_string);
		/* Flush all output streams prior to asserting.
		We replace the original assert(0) call with an exit(EXIT_FAILURE),
		since some compilers seem to like to optimize away assertions. */
		fflush(NULL);
		exit(EXIT_FAILURE);
	}
}
#endif

/***************/

void	WARN(char *typelist, ...)
{
	char *c;
	 int32 ival;
	uint32 uval;
	double dval;

	va_list varargs;
	va_start(varargs, typelist);
	/* Define a convenient spot to set a breakpoint: */
	for(c = typelist; *c; c++)
	{
		switch(*c)
		{
			case 'i':
				ival = va_arg(varargs, int32);
				break;
			case 'u':
				uval = va_arg(varargs,uint32);
				break;
			case 'd':
				dval = va_arg(varargs,double);
				break;
			default :
				ASSERT(HERE, 0,"0");
				break;
		}
	}
	va_end(varargs);
}

/***************/

/* Slow, default-int versions of leading and trailing-zero-counting algorithms: */
/*
int leadz(int x)
{
	uint32 i;
	for(i = 1; i < 32; i++)
	{
		if(x < 0)break;
		x <<= 1;
	}
	return(i-1);
}
*/

uint32 trailz32(uint32 x)
{
	uint32 i;
	for(i = 1; i < 32; i++)
	{
		if(x & 1)break;
		x >>= 1;
	}
	return(i-1);
}

uint32 trailz64(uint64 x)
{
	uint32 i;
	for(i = 1; i < 64; i++)
	{
		if(x & 0x0000000000000001)break;
		x >>= 1;
	}
	return(i-1);
}


/***********************************************************************************/
/* Fast 32 and 64-bit-specific versions of leading-zero-counting algorithms: */
/*
An illustrated view of the leading-zeros algorithm:
i1 =
00010111001010101001001110101010, has 3 lzeros.
shift := 16
-1>>16 =
00000000000000001111111111111111 < i1, so i1 has  < 16 lzeros, ly stays = 0, shift -= 8
-1>> 8 =
00000000111111111111111111111111 < i1, so i1 has  <  8 lzeros, ly stays = 0, shift -= 4
-1>> 4 =
00001111111111111111111111111111 < i1, so i1 has  <  4 lzeros, ly stays = 0, shift -= 2
-1>> 2 =
00111111111111111111111111111111 >=i1, so i1 has  >= 2 lzeros, ly += 2  = 2, shift += 1
-1>> 3 =
00011111111111111111111111111111 >=i1, so i1 has  >= 3 lzeros, ly += 1  = 3, shift += 0, done.


i2 =
00000000000000000000101110010101, has 20 lzeros.
shift := 16
-1>>16 =
00000000000000001111111111111111 >=i2, so i2 has  >=16 lzeros, ly +=16  =16, shift += 8
-1>>24 =
00000000000000000000000011111111 < i2, so i2 has  < 24 lzeros, ly stays =16, shift -= 4
-1>>20 =
00000000000000000000111111111111 >=i2, so i2 has  >=20 lzeros, ly +=    =20, shift += 2
-1>>22 =
00000000000000000000001111111111 < i2, so i2 has  < 22 lzeros, ly stays =20, shift -= 1
-1>>21 =
00000000000000000000011111111111 < i2, so i2 has  < 21 lzeros, ly stays =20, shift -= 0, done.
*/

/***************/

uint32 leadz32(uint32 i)
{
	uint32 lz, k, shift, ones_mask = 0xFFFFFFFF;
	if(i == 0) return 32;
	if((int32)i < 0) return 0;

	lz    =  0;
	shift = 16;
	k     = 16;
	while(k > 0)
	{
		if( (ones_mask >> shift) < i )
		{
			k >>= 1;
			shift -= k;
		}
		else
		{
			lz += k;
			k >>= 1;
			shift += k;
		}
	}

	DBG_ASSERT(HERE, ( (i << lz) >> lz == i ),"ERROR A in leadz32");
	DBG_ASSERT(HERE, ( i >> (32-lz)    == 0 ),"ERROR B in leadz32");

	return lz;
}

/***************/

uint32 leadz64(uint64 i)
{
	uint32 lz, k, shift;
	uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	if(i == 0) return 64;
	if(( int64)i < 0) return 0;

	lz    =  0;
	shift = 32;
	k     = 32;
	while(k > 0)
	{
		if( (ones_mask >> shift) < i )
		{
			k >>= 1;
			shift -= k;
		}
		else
		{
			lz += k;
			k >>= 1;
			shift += k;
		}
	}

	DBG_ASSERT(HERE, ( (i << lz) >> lz == i ),"ERROR A in leadz64");
	DBG_ASSERT(HERE, ( i >> (64-lz)    == 0 ),"ERROR B in leadz64");

	return lz;
}

uint32	leadz128(uint128 i)
{
	if(i.d1)
	{
		return leadz64(i.d1);
	}
	else
	{
		return leadz64(i.d0) + 64;
	}
}

uint32	leadz192(uint192 i)
{
	if(i.d2)
	{
		return leadz64(i.d2);
	}
	else if(i.d1)
	{
		return leadz64(i.d1) + 64;
	}
	else
	{
		return leadz64(i.d0) + 128;
	}
}

uint32	leadz256(uint256 i)
{
	if(i.d3)
	{
		return leadz64(i.d3);
	}
	else if(i.d2)
	{
		return leadz64(i.d2) + 64;
	}
	else if(i.d1)
	{
		return leadz64(i.d1) + 128;
	}
	else
	{
		return leadz64(i.d0) + 192;
	}
}

/***************/

uint32 isPow2(uint32 i32)
{
/*   If the input == 2^p, returns 1; otherwise returns 0. */
	uint32 i;
	ASSERT(HERE, i32, "isPow2: input zero!");	/* Check for zero in debug mode */
	i = trailz32(i32);
	if((i32 >> i) != 1)
		return 0;
	else
		return 1;
}

/***************/

uint32 isPow4(uint32 i32)
{
/*   If the input == 4^p, returns 1; otherwise returns 0. */
	uint32 i;
	ASSERT(HERE, i32, "isPow4: input zero!");	/* Check for zero in debug mode */
	i = trailz32(i32);
	if((i & 1) || (i32 >> i) != 1)
		return 0;
	else
		return 1;
}

/***************/

/* Population count: Returns the number of set bits in a 32-bit int.
   Optimized for the case where the low-order bits are those most-likely to be set.
*/
uint32 popCount(uint32 num)
{
	uint32 numones = 0;
	while(num)
	{
		if(num & 1)
			numones++;
		num >>= 1;
	}
	return numones;
}

/***************/

/* Extract (nbits) bits beginning at position (beg) */

uint32 ibits32(uint32 i, uint32 beg, uint32 nbits)
{
	uint32 ones_mask = 0xFFFFFFFF;
	return ( (i >> beg) & ~(ones_mask << nbits) );
}

uint64 ibits64(uint64 i, uint32 beg, uint32 nbits)
{
	uint64 ib;
	uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	ib = (i >> beg) & ~(ones_mask << nbits);
	return ( ib );
}

/***************/

/* Return result of copying (nbits) bits of a 64-bit integer x, starting at bit
(src_bit_start) into a target 64-bit integer y, starting at bit (tgt_bit_start).
The syntax mirrors that of the Fortran-90 MVBITS library function.
If bit-index parameters are not in [0,63], result will be unpredictable.
*                Example:         23                    41                      0 */
void	mvbits64(uint64 x, uint32 src_bit_start, uint32 nbits, uint64*y, uint32 tgt_bit_start)
{
	uint64 ones_mask = 0xFFFFFFFFFFFFFFFFull;
	uint64 mask;
	DBG_ASSERT(HERE, nbits <= 64, "util.c/mvbits64: nbits out or range\n");
	DBG_ASSERT(HERE, src_bit_start < 64, "util.c/mvbits64: src_bit_start out or range\n");
	DBG_ASSERT(HERE, tgt_bit_start < 64, "util.c/mvbits64: tgt_bit_start out or range\n");
	if(nbits == 0) return;
	mask = (ones_mask >> ((64-nbits)&63));	/* mask = 0x000001ffffffffff */
	/* Zero out the target bits: */
	*y &= ~(mask << tgt_bit_start);			/*   y &= 0xfffffe0000000000 */
	/* Copy the source bits into the gap: */
	*y += ((x >> src_bit_start) & mask) << tgt_bit_start;
	return;
}

/***************/

int pprimeF(uint32 p, uint32 base)
{
/* returns 1 if p is a base-z Fermat pseudoprime, 0 otherwise. */

	uint64 y = 1, n = p-1, flag;
	uint64 z = base;	/* Need a 64-bit to store intermediate products without overflow */

	while(n)
	{
		flag = n & 1;
		n >>= 1;
		if(flag) y = (y*z)%p;
		z = (z*z)%p;
		if(!z) return 0;
	}
	return((int)(y==1));
}

/***************/

int isPRP(uint32 p)
{
	/* TODO: replace/supplement this with a rigorous trial-divide test for p < 2^32 */
	return(pprimeF(p,2) && pprimeF(p,3) && pprimeF(p,5) && pprimeF(p,7) && pprimeF(p,11) && pprimeF(p,13));
}

/*******************/

/* Calculate 2^p mod q for p, q 32-bit unsigned ints. This can be used (among
other things) to effect a fast Fermat base-2 pseudoprime test, by calling with q = p-1.
*/
int twopmodq32(uint32 p, uint32 q)
{
	 int32 j;
	uint32 lead5, pshift, qhalf, qinv, zshift, start_index, x, lo, hi;

	ASSERT(HERE, (q&1) == 1, "twopmodq32: even modulus!");
	qhalf = q >> 1;	/* = (q-1)/2, since q odd. */

	pshift = p + 32;
	if(pshift < p)	/* Need special-casing for p just below 2^32  - the primes 2^32-(5,17) are good testcases here. */
	{
		j = -1;	/* leadz32(pshift) for 33-bit pshift goes negative */
		/* Extract leftmost 5 bits of pshift: */
		lead5 = 16 + (pshift >> 28);
	}
	else
	{
		/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
		j = leadz32(pshift);
		/* Extract leftmost 5 bits of pshift: */
		lead5 = ((pshift<<j) >> 27);
	}

	start_index = 32-j-5;	/* Leftward bit at which to start the l-r binary powering, assuming
							the leftmost 5 bits have already been processed via a shift (see next). */

	zshift = 31 - lead5;
	zshift <<= 1;		/* Doubling the shift count here takes cares of the first SQR_LOHI */
	pshift = ~pshift;	/* Overflow doesn't matter here, as long as we got the leading 5 bits of pshift right. */

	/*
	!    Find modular inverse (mod 2^32) of q in preparation for modular multiply.
	!    We use the simple and elegant iterative inversion method of Montgomery,
	!    which amounts to a modular analogue of Newton's method for iterative inversion:
	!
	!    0)   Zinv = Z                   ! Z*Zinv == 1 (mod 2^3)
	!    1)   Zinv = Zinv*(2 - Z*Zinv)   ! Z*Zinv == 1 (mod 2^6), etc.
	!
	!    where the number of correct bits (at the low end) doubles at each step,
	!    all arithmetic is modulo 2^32 and we repeat step (1) until we have the needed 32 bits.
	!
	!    We choose a different starting value of Zinv, XOR(3*Z, 2),
	!    so the first congruence holds modulo 2^4, thus requiring just 3 iterations.
	*/
	qinv = (q+q+q) ^ (uint32)2;	/* Overflow doesn't matter here, since we only care about the low 2 bits of 3*q. */

	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);
	qinv = qinv*((uint32)2 - q*qinv);

	/* Since zstart is a power of two < 2^32, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* For 64-bit hardware, Make sure we get a 32-bit shift result here by ANDing with 2^32-1: */
	lo = (qinv << zshift) & (uint32)0xffffffff;
	/* Emulate MULH64 here by getting full 64-bit product and right-shifting: */
	lo = (uint32)(((uint64)q * (uint64)lo) >> 32);
	x  = q - lo;

	if((pshift >> j) & (uint32)1)
	{
		DBG_ASSERT(HERE, x < q,"util.c: x < q");
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(x > qhalf)
		{
			x += x;
			x -= q;
		}
		else
		{
			x += x;
		}
	}

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI64(x,lo,hi): */
	#if 0
		uint64 t = (uint64)x * (uint64)x;	hi = (uint32)(t >> 32);	lo = (uint32)(t) & (uint32)0xffffffff;
		lo *= qinv;
		/* MULH32(q,lo,lo): */
		lo = (uint32)(((uint64)q * (uint64)lo) >> 32);
	#else
	  #ifdef MUL_LOHI32_SUBROUTINE
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		lo = __MULH32(q, lo);
	  #else
		MUL_LOHI32(x,x, lo,hi);
		lo *= qinv;
		MULH32(q,lo, lo);
	  #endif
	#endif

		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
	#if 0
		x = hi - lo;
		if(x > hi)
		{
			x += q;	/* had a borrow */
		}
	#else
		x = hi - lo + ((-(hi < lo)) & q);
	#endif

		if((pshift >> j) & (uint32)1)
		{
		/* Branchless version is much faster, but less readable, so give the branched one inside a #if 0: */
		#if 0
			if(x > qhalf)	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			{
				x = x + x;
				x -= q;
			}
			else
			{
				x = x + x;
			}
		#else
			x = x + x - ((-(x > qhalf)) & q);
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	return((int)((x + x - q) == 1));
}

/* Does an 8-fold base-2 PRP test on the prime candidates q0-7. */
int twopmodq32_x8(uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7)
{
	int retval = 0;
	 int32 j;
	uint64 t;
	uint32 start_index;
	uint32 lead0, pshift0, qinv0, zshift0, x0, lo0, hi0, qhalf0;
	uint32 lead1, pshift1, qinv1, zshift1, x1, lo1, hi1, qhalf1;
	uint32 lead2, pshift2, qinv2, zshift2, x2, lo2, hi2, qhalf2;
	uint32 lead3, pshift3, qinv3, zshift3, x3, lo3, hi3, qhalf3;
	uint32 lead4, pshift4, qinv4, zshift4, x4, lo4, hi4, qhalf4;
	uint32 lead5, pshift5, qinv5, zshift5, x5, lo5, hi5, qhalf5;
	uint32 lead6, pshift6, qinv6, zshift6, x6, lo6, hi6, qhalf6;
	uint32 lead7, pshift7, qinv7, zshift7, x7, lo7, hi7, qhalf7;

	DBG_ASSERT(HERE, (q0 < q1) && (q1 < q2) && (q2 < q3) && (q3 < q4) && (q4 < q5) && (q5 < q6) && (q6 < q7), "twopmodq32_x8: Inputs nonmonotone!");

	qhalf0 = q0 >> 1;	/* = (q-1)/2, since q odd. */
	qhalf1 = q1 >> 1;
	qhalf2 = q2 >> 1;
	qhalf3 = q3 >> 1;
	qhalf4 = q4 >> 1;
	qhalf5 = q5 >> 1;
	qhalf6 = q6 >> 1;
	qhalf7 = q7 >> 1;

	/* (p[i]-1)+32 = p + [31,33,37,39,41,49,54,60]: */
	pshift0 = q0 + 31;
	pshift1 = q1 + 31;
	pshift2 = q2 + 31;
	pshift3 = q3 + 31;
	pshift4 = q4 + 31;
	pshift5 = q5 + 31;
	pshift6 = q6 + 31;
	pshift7 = q7 + 31;

	/* Find number of leading zeros in p, use it to find the position of the leftmost ones bit: */
	j = leadz32(pshift0);
	if( leadz32(pshift7) != j )	/* Fused 8-fold algo needs all p's to have same bitlength */
	{
		retval  = (uint32)twopmodq32(q0-1, q0);
		retval += (uint32)twopmodq32(q1-1, q1) << 1;
		retval += (uint32)twopmodq32(q2-1, q2) << 2;
		retval += (uint32)twopmodq32(q3-1, q3) << 3;
		retval += (uint32)twopmodq32(q4-1, q4) << 4;
		retval += (uint32)twopmodq32(q5-1, q5) << 5;
		retval += (uint32)twopmodq32(q6-1, q6) << 6;
		retval += (uint32)twopmodq32(q7-1, q7) << 7;
		return retval;
	}

	if(pshift0 < q0)	/* Need special-casing for p just below 2^32  - the primes 2^32-(5,17) are good testcases here. */
	{
		j = -1;	/* leadz32(pshift) for 33-bit pshift goes negative */
		/* Extract leftmost 5 bits of pshift: */
		lead0 = 16 + (pshift0 >> 28);
		lead1 = 16 + (pshift1 >> 28);
		lead2 = 16 + (pshift2 >> 28);
		lead3 = 16 + (pshift3 >> 28);
		lead4 = 16 + (pshift4 >> 28);
		lead5 = 16 + (pshift5 >> 28);
		lead6 = 16 + (pshift6 >> 28);
		lead7 = 16 + (pshift7 >> 28);
	}
	else
	{
		/* Extract leftmost 5 bits of pshift and subtract from 32: */
		lead0 = ((pshift0<<j) >> 27);
		lead1 = ((pshift1<<j) >> 27);
		lead2 = ((pshift2<<j) >> 27);
		lead3 = ((pshift3<<j) >> 27);
		lead4 = ((pshift4<<j) >> 27);
		lead5 = ((pshift5<<j) >> 27);
		lead6 = ((pshift6<<j) >> 27);
		lead7 = ((pshift7<<j) >> 27);
	}

	start_index = 32-j-5;	/* Leftward bit at which to start the l-r binary powering, assuming
				 the leftmost 5 bits have already been processed via a shift (see next). */

	/* Doubling the shift count here takes cares of the first SQR_LOHI */
	zshift0 = 31 - lead0;	zshift0 <<= 1;	pshift0 = ~pshift0;
	zshift1 = 31 - lead1;	zshift1 <<= 1;	pshift1 = ~pshift1;
	zshift2 = 31 - lead2;	zshift2 <<= 1;	pshift2 = ~pshift2;
	zshift3 = 31 - lead3;	zshift3 <<= 1;	pshift3 = ~pshift3;
	zshift4 = 31 - lead4;	zshift4 <<= 1;	pshift4 = ~pshift4;
	zshift5 = 31 - lead5;	zshift5 <<= 1;	pshift5 = ~pshift5;
	zshift6 = 31 - lead6;	zshift6 <<= 1;	pshift6 = ~pshift6;
	zshift7 = 31 - lead7;	zshift7 <<= 1;	pshift7 = ~pshift7;

	/*
	Find modular inverse (mod 2^32) of q in preparation for modular multiply.
	*/
	qinv0 = (q0+q0+q0) ^ (uint32)2;
	qinv1 = (q1+q1+q1) ^ (uint32)2;
	qinv2 = (q2+q2+q2) ^ (uint32)2;
	qinv3 = (q3+q3+q3) ^ (uint32)2;
	qinv4 = (q4+q4+q4) ^ (uint32)2;
	qinv5 = (q5+q5+q5) ^ (uint32)2;
	qinv6 = (q6+q6+q6) ^ (uint32)2;
	qinv7 = (q7+q7+q7) ^ (uint32)2;
	for(j = 0; j < 3; ++j)
	{
		qinv0 = qinv0*((uint32)2 - q0*qinv0);
		qinv1 = qinv1*((uint32)2 - q1*qinv1);
		qinv2 = qinv2*((uint32)2 - q2*qinv2);
		qinv3 = qinv3*((uint32)2 - q3*qinv3);
		qinv4 = qinv4*((uint32)2 - q4*qinv4);
		qinv5 = qinv5*((uint32)2 - q5*qinv5);
		qinv6 = qinv6*((uint32)2 - q6*qinv6);
		qinv7 = qinv7*((uint32)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^32, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* For 64-bit hardware, Make sure we get a 32-bit shift result here by ANDing with 2^32-1: */
	lo0 = (qinv0 << zshift0) & (uint32)0xffffffff;
	lo1 = (qinv1 << zshift1) & (uint32)0xffffffff;
	lo2 = (qinv2 << zshift2) & (uint32)0xffffffff;
	lo3 = (qinv3 << zshift3) & (uint32)0xffffffff;
	lo4 = (qinv4 << zshift4) & (uint32)0xffffffff;
	lo5 = (qinv5 << zshift5) & (uint32)0xffffffff;
	lo6 = (qinv6 << zshift6) & (uint32)0xffffffff;
	lo7 = (qinv7 << zshift7) & (uint32)0xffffffff;

	/* lo = MULH64(q, lo): */
#ifdef MUL_LOHI32_SUBROUTINE
	lo0 = __MULH32(q0, lo0);
	lo1 = __MULH32(q1, lo1);
	lo2 = __MULH32(q2, lo2);
	lo3 = __MULH32(q3, lo3);
	lo4 = __MULH32(q4, lo4);
	lo5 = __MULH32(q5, lo5);
	lo6 = __MULH32(q6, lo6);
	lo7 = __MULH32(q7, lo7);
#else
	MULH32(q0,lo0, lo0);
	MULH32(q1,lo1, lo1);
	MULH32(q2,lo2, lo2);
	MULH32(q3,lo3, lo3);
	MULH32(q4,lo4, lo4);
	MULH32(q5,lo5, lo5);
	MULH32(q6,lo6, lo6);
	MULH32(q7,lo7, lo7);
#endif

	x0  = q0 - lo0;
	x1  = q1 - lo1;
	x2  = q2 - lo2;
	x3  = q3 - lo3;
	x4  = q4 - lo4;
	x5  = q5 - lo5;
	x6  = q6 - lo6;
	x7  = q7 - lo7;

	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
	if((pshift0 >> j) & (uint32)1){ DBG_ASSERT(HERE, x0 < q0,"util.c: x0 < q0"); x0 = x0 + x0 - ((-(x0 > qhalf0)) & q0); }
	if((pshift1 >> j) & (uint32)1){ DBG_ASSERT(HERE, x1 < q1,"util.c: x1 < q1"); x1 = x1 + x1 - ((-(x1 > qhalf1)) & q1); }
	if((pshift2 >> j) & (uint32)1){ DBG_ASSERT(HERE, x2 < q2,"util.c: x2 < q2"); x2 = x2 + x2 - ((-(x2 > qhalf2)) & q2); }
	if((pshift3 >> j) & (uint32)1){ DBG_ASSERT(HERE, x3 < q3,"util.c: x3 < q3"); x3 = x3 + x3 - ((-(x3 > qhalf3)) & q3); }
	if((pshift4 >> j) & (uint32)1){ DBG_ASSERT(HERE, x4 < q4,"util.c: x4 < q4"); x4 = x4 + x4 - ((-(x4 > qhalf4)) & q4); }
	if((pshift5 >> j) & (uint32)1){ DBG_ASSERT(HERE, x5 < q5,"util.c: x5 < q5"); x5 = x5 + x5 - ((-(x5 > qhalf5)) & q5); }
	if((pshift6 >> j) & (uint32)1){ DBG_ASSERT(HERE, x6 < q6,"util.c: x6 < q6"); x6 = x6 + x6 - ((-(x6 > qhalf6)) & q6); }
	if((pshift7 >> j) & (uint32)1){ DBG_ASSERT(HERE, x7 < q7,"util.c: x7 < q7"); x7 = x7 + x7 - ((-(x7 > qhalf7)) & q7); }

	for(j = start_index-2; j >= 0; j--)
	{
		/* SQR_LOHI64(x,lo,hi): */
	#ifdef MUL_LOHI32_SUBROUTINE
		MUL_LOHI32(x0,x0, lo0,hi0);
		MUL_LOHI32(x1,x1, lo1,hi1);
		MUL_LOHI32(x2,x2, lo2,hi2);
		MUL_LOHI32(x3,x3, lo3,hi3);
		MUL_LOHI32(x4,x4, lo4,hi4);
		MUL_LOHI32(x5,x5, lo5,hi5);
		MUL_LOHI32(x6,x6, lo6,hi6);
		MUL_LOHI32(x7,x7, lo7,hi7);
		lo0 *= qinv0;
		lo1 *= qinv1;
		lo2 *= qinv2;
		lo3 *= qinv3;
		lo4 *= qinv4;
		lo5 *= qinv5;
		lo6 *= qinv6;
		lo7 *= qinv7;
		lo0 = __MULH32(q0, lo0);
		lo1 = __MULH32(q1, lo1);
		lo2 = __MULH32(q2, lo2);
		lo3 = __MULH32(q3, lo3);
		lo4 = __MULH32(q4, lo4);
		lo5 = __MULH32(q5, lo5);
		lo6 = __MULH32(q6, lo6);
		lo7 = __MULH32(q7, lo7);
	#else
		MUL_LOHI32(x0,x0, lo0,hi0);
		MUL_LOHI32(x1,x1, lo1,hi1);
		MUL_LOHI32(x2,x2, lo2,hi2);
		MUL_LOHI32(x3,x3, lo3,hi3);
		MUL_LOHI32(x4,x4, lo4,hi4);
		MUL_LOHI32(x5,x5, lo5,hi5);
		MUL_LOHI32(x6,x6, lo6,hi6);
		MUL_LOHI32(x7,x7, lo7,hi7);
		lo0 *= qinv0;
		lo1 *= qinv1;
		lo2 *= qinv2;
		lo3 *= qinv3;
		lo4 *= qinv4;
		lo5 *= qinv5;
		lo6 *= qinv6;
		lo7 *= qinv7;
		MULH32(q0,lo0, lo0);
		MULH32(q1,lo1, lo1);
		MULH32(q2,lo2, lo2);
		MULH32(q3,lo3, lo3);
		MULH32(q4,lo4, lo4);
		MULH32(q5,lo5, lo5);
		MULH32(q6,lo6, lo6);
		MULH32(q7,lo7, lo7);
	#endif

		/* if(x < 0) x += q; */
		x0 = hi0 - lo0 + ((-(hi0 < lo0)) & q0);
		x1 = hi1 - lo1 + ((-(hi1 < lo1)) & q1);
		x2 = hi2 - lo2 + ((-(hi2 < lo2)) & q2);
		x3 = hi3 - lo3 + ((-(hi3 < lo3)) & q3);
		x4 = hi4 - lo4 + ((-(hi4 < lo4)) & q4);
		x5 = hi5 - lo5 + ((-(hi5 < lo5)) & q5);
		x6 = hi6 - lo6 + ((-(hi6 < lo6)) & q6);
		x7 = hi7 - lo7 + ((-(hi7 < lo7)) & q7);

		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if((pshift0 >> j) & (uint32)1){ DBG_ASSERT(HERE, x0 < q0,"util.c: x0 < q0"); x0 = x0 + x0 - ((-(x0 > qhalf0)) & q0); }
		if((pshift1 >> j) & (uint32)1){ DBG_ASSERT(HERE, x1 < q1,"util.c: x1 < q1"); x1 = x1 + x1 - ((-(x1 > qhalf1)) & q1); }
		if((pshift2 >> j) & (uint32)1){ DBG_ASSERT(HERE, x2 < q2,"util.c: x2 < q2"); x2 = x2 + x2 - ((-(x2 > qhalf2)) & q2); }
		if((pshift3 >> j) & (uint32)1){ DBG_ASSERT(HERE, x3 < q3,"util.c: x3 < q3"); x3 = x3 + x3 - ((-(x3 > qhalf3)) & q3); }
		if((pshift4 >> j) & (uint32)1){ DBG_ASSERT(HERE, x4 < q4,"util.c: x4 < q4"); x4 = x4 + x4 - ((-(x4 > qhalf4)) & q4); }
		if((pshift5 >> j) & (uint32)1){ DBG_ASSERT(HERE, x5 < q5,"util.c: x5 < q5"); x5 = x5 + x5 - ((-(x5 > qhalf5)) & q5); }
		if((pshift6 >> j) & (uint32)1){ DBG_ASSERT(HERE, x6 < q6,"util.c: x6 < q6"); x6 = x6 + x6 - ((-(x6 > qhalf6)) & q6); }
		if((pshift7 >> j) & (uint32)1){ DBG_ASSERT(HERE, x7 < q7,"util.c: x7 < q7"); x7 = x7 + x7 - ((-(x7 > qhalf7)) & q7); }
	}

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */
	retval += ((x0 + x0 - q0) == 1)     ;
	retval += ((x1 + x1 - q1) == 1) << 1;
	retval += ((x2 + x2 - q2) == 1) << 2;
	retval += ((x3 + x3 - q3) == 1) << 3;
	retval += ((x4 + x4 - q4) == 1) << 4;
	retval += ((x5 + x5 - q5) == 1) << 5;
	retval += ((x6 + x6 - q6) == 1) << 6;
	retval += ((x7 + x7 - q7) == 1) << 7;
	return retval;
}

/*******************/

/* Simple Euclidean GCD for 32-bit unsigned inputs. (Cf. Algorithm 2.1.4 in Crandall/Pomerance.)
For integers x, y with x, y > 0, returns GCD(x,y). If x or y = 0, returns max(x,y).
*/
uint32 gcd32(uint32 x, uint32 y)
{
	uint32 q, f;

	if(y == 0)
	{
		return x;
	}

	while(y)
	{
		/* Find quotient of current x/y and round toward zero: */
		q = x/y;
		/* Find y' and store in temporary: */
		f = x - q*y;
		/* Find x', i.e. move the old value of y into the slots for x: */
		x = y;
		/* New value of y: */
		y = f;
	}

	return(x);
}

uint64 gcd64(uint64 x, uint64 y)
{
	uint64 q, f;

	if(y == 0)
	{
		return x;
	}

	while(y)
	{
		q = x/y;
		f = x - q*y;
		x = y;
		y = f;
	}

	return(x);
}

/*******************/

/* Simple extended Euclidean GCD for 32-bit unsigned inputs. (Cf. Algorithm 2.1.4 in Crandall/Pomerance.)
For integers x, y with x, y > 0, returns integers {a,b,g} such that a*x + b*y = g = GCD(x,y).

When g = 1 and y > 0, the residues a and b are the inverses of x (mod y) and y (mod x), respectively.

The GCD g is the return value of the function; note that the multipliers a and b
overwrite the inputs x and y, so if the original inputs are needed subsequently,
they must be copied prior to calling the function.
*/
uint32 egcd32(uint32 *x, uint32 *y)
{
	uint32 g = *x, w = *y, q;
	int    a = 1, b = 0, u = 0, v = 1;
	/* Sign of these 3 doesn't matter since they're just temporaries: */
	uint32 d, e, f;

	if(*x == *y)
	{
		fprintf(stderr,"ERROR: eGCD of identical arguments x = y = %u is illegal!\n", *x);
		ASSERT(HERE, 0,"0");
	}
	else if(*x == 0)
	{
		fprintf(stderr,"ERROR: eGCD called with zero input: x = %u, y = %u\n", *x, *y);
		ASSERT(HERE, 0,"0");
	}
	else if(*y == 0)
	{
		fprintf(stderr,"ERROR: eGCD called with zero input: x = %u, y = %u\n", *x, *y);
		ASSERT(HERE, 0,"0");
	}

	while(w)
	{
		/* Find quotient of current x/y and round toward zero: */
		q = g/w;
		/* Find (u', v', w') and store in 3 temporaries: */
		d = a - q*u;
		e = b - q*v;
		f = g - q*w;
		/* Find (a', b', g'), i.e. move the old values of (u,v,w) into the slots for (a,b,g): */
		a = u;
		b = v;
		g = w;
		/* Recover new values of (u, v, w) from the temporaries: */
		u = d;
		v = e;
		w = f;
	}

	*x = a;
	*y = b;
	return(g);
}

/*********************/
/*
Finds multiplicative inverse of z (mod n).
*/
int modinv32(uint32 z, uint32 n)
{
	uint32 x = z, y = n;
	uint32 gcd;

	/* 01/26/04: Turns out this sign check isn't necessary, since the eGCD
				routine automatically handles the case x < y:

	if(x < y)
		gcd = egcd32(&y, &x);
	else if(x > y)
	*/

	gcd = egcd32(&x, &y);
	DBG_ASSERT(HERE, gcd == 1,"gcd in modinv21 is non-unity!");

/*printf("modinv(%u, %u) = %u\n", z, n, x);
	if(n > 100)ASSERT(HERE, 0,"0");
*/
	return x;
}

/********************/

/* Complex multiplication */
struct complex cmul(struct complex *a, struct complex *b)
{
	struct complex cout;
	cout.re = (*a).re*(*b).re - (*a).im*(*b).im;
	cout.im = (*a).re*(*b).im + (*a).im*(*b).re;
	return cout;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 128-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
*/
uint128 xmody128(uint128 x, uint128 y)
{
	uint32 lzx, lzy, nshiftl;
	uint128 t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	     if(x.d1)
		lzx = leadz64(x.d1);
	else
		lzx = leadz64(x.d0) + 64;

	     if(y.d1)
		lzy = leadz64(y.d1);
	else
		lzy = leadz64(y.d0) + 64;

	/* X < Y: return unmodified X. */
	if(lzx > lzy)
		return x;

	nshiftl = lzy - lzx;	/* nshiftlr = 64-nshiftl; */

	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT128(y, nshiftl, t);

		if(CMPULT128(t, x))
			SUB128(x, t, x);

		/* Right-shift t one place: */
		--nshiftl;
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT128(y, x))
		SUB128(x, y, x);

	return x;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 192-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
*/

uint192 xmody192(uint192 x, uint192 y)
{
	uint32 lzx, lzy, nshiftl;
	uint192 t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	     if(x.d2)
		lzx = leadz64(x.d2);
	else if(x.d1)
		lzx = leadz64(x.d1) + 64;
	else
		lzx = leadz64(x.d0) + 128;

	     if(y.d2)
		lzy = leadz64(y.d2);
	else if(y.d1)
		lzy = leadz64(y.d1) + 64;
	else
		lzy = leadz64(y.d0) + 128;

	/* X < Y: return unmodified X. */
	if(lzx > lzy)
		return x;

	nshiftl = lzy - lzx;

/*
printf("x =%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", x.d2, x.d1, x.d0);
printf("y =%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", y.d2, y.d1, y.d0);
printf("nshiftl = %u\n", nshiftl);
*/
	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT192(y, nshiftl, t);
/*printf("y<<%u=" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", nshiftl, t.d2, t.d1, t.d0); */

		if(CMPULT192(t, x))
		{
			SUB192(x, t, x);
/*printf("x*=%20" LLU "*2^128 + %20" LLU "*2^64 + %20" LLU "\n", x.d2, x.d1, x.d0); */
		}

		/* Right-shift t one place: */
		--nshiftl;
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT192(y, x))
		SUB192(x, y, x);

	return x;
}

/***********************************************************************************/
/*
Function to reduce x modulo y, where x and y are both 256-bit unsigned integers.
Algorithm is simple-but-slow bitwise shift-and-subtract scheme.
*/

uint256 xmody256(uint256 x, uint256 y)
{
	uint32 lzx, lzy, nshiftl;
	uint256 t;

	/* In preparation for x%y, Find the # of leading zeros in x and y. */
	lzx = LEADZ256(x);
	lzy = LEADZ256(y);

	nshiftl = lzy - lzx;
	while(nshiftl)
	{
		/* Use t to store the left-shifted versions of y: */
		LSHIFT256(y, nshiftl, t);
		if(CMPULT256(t, x))
		{
			SUB256(x, t, x);
		}
		/* Right-shift t one place: */
		--nshiftl;
	}
	/* Must ensure that this gets done once even if lzx == lzy: */
	if(CMPULT256(y, x))
		SUB256(x, y, x);

	return x;
}


/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 128-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x128_mod_y32 instead.
*/
uint32 x128_div_y32(uint128 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
/*printf("INIT: two64divy, two64mody = %20llu %20llu\n\n", two64divy, two64mody); */
	}

	/* Divide high digit by y, storing remainder in cy: */
	cy = (x->d1)%y;
	(x->d1) /= y;

	/* Remainder (must calculate this before modifying (x->d0), obviously): */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;

	/* Low digit of result: we must separately divide (x->d0) by y
	(making sure to add (x->d0)%y to  cy*two64mody first, so as not to drop a digit)
	because x->d0 may be as large as 2^64-1, and adding cy*two64mody
	prior to dividing risks unsigned integer overflow:
	*/
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;
/*printf("%20llu %20llu %2llu %2llu\n", x->d1, x->d0, cy, rem); */
	return (uint32)rem;
}

uint32 x128_mod_y32(uint128 x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Divide high digit by y, storing remainder in cy: */
	cy = (x.d1)%y;
	/* Remainder: */
	xlomody = (x.d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;

	return (uint32)rem;
}

/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 192-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x192_mod_y32 instead.
*/
uint32 x192_div_y32(uint192 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint128 t128;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 2 digits into a local uint128, then call 128-bit divide
	on those, return value is the carry into the low digit: */
	t128.d1 = (x->d2);
	t128.d0 = (x->d1);
	cy = x128_div_y32(&t128, y);
	(x->d2) = t128.d1;
	(x->d1) = t128.d0;

	/* Low digit: */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;

	return (uint32)rem;
}

uint32 x192_mod_y32(uint192 x, uint32 y)
{
	uint64 cy, rem;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint128 t128;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 2 digits into a local uint128, then call 128-bit divide
	on those, return value is the carry into the low digit: */
	t128.d1 = (x.d2);
	t128.d0 = (x.d1);
	cy = x128_div_y32(&t128, y);

	/* Low digit: */
	rem = (cy*two64mody + ((x.d0))%y)%y;

	return (uint32)rem;
}

/***********************************************************************************/
/*
Divide-with-Remainder of x by y, where x is a 256-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
Returns (x - x%y)/y in x, 32-bit remainder in the function result.

If you only want the remainder, not to perform the divide, call x256_mod_y32 instead.
*/
uint32 x256_div_y32(uint256 *x, uint32 y)
{
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint192 t192;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 3 digits into a local uint192, then call 192-bit divide
	on those, return value is the carry into the low digit: */
	t192.d2 = (x->d3);
	t192.d1 = (x->d2);
	t192.d0 = (x->d1);
	cy = x192_div_y32(&t192, y);
	(x->d3) = t192.d2;
	(x->d2) = t192.d1;
	(x->d1) = t192.d0;

	/* Low digit: */
	xlomody = (x->d0)%y;
	tsum = cy*two64mody + xlomody;
	rem = tsum%y;
	(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;

	return (uint32)rem;
}

uint32 x256_mod_y32(uint256 x, uint32 y)
{
	uint64 cy, rem;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;
	uint192 t192;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	/* Copy the upper 3 digits into a local uint192, then call 192-bit divide
	on those, return value is the carry into the low digit: */
	t192.d2 = (x.d3);
	t192.d1 = (x.d2);
	t192.d0 = (x.d1);
	cy = x192_div_y32(&t192, y);

	/* Low digit: */
	rem = (cy*two64mody + ((x.d0))%y)%y;

	return (uint32)rem;
}

/***********************************************************************************/

/* Need the uint64 ones of these because some compilers (e.g. MSVC, a.k.a .NET)
don't properly print 64-bit ints. */

/*
Returns decimal character representation of a 64-bit unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint64_base10_char(char char_buf[], uint64 q)
{
	uint32 i, n_dec_digits = 0, curr_digit;
	char c;
	/* 2^64 has 20 decimal digits - assume the user has allocated at least 20+1 for char_buf: */
	uint32 MAX_DIGITS = 20;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Needed to cast modulus 10 to uint32 here for result to come out correct: */
		curr_digit = q%(uint32)10;

		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		if(q != 0 || n_dec_digits == 0)
		{
			c = curr_digit + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;

		q /= 10;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/********************/

/*
Returns all-caps hexadecimal character representation of a uint64 in char_buf.
*/
int	convert_uint64_base16_char(char char_buf[], uint64 q)
{
	int i;
	const int hex_chars[16] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

	for(i=15; i >= 0; i--)
	{
		char_buf[i] = hex_chars[q & 15];
		q >>= 4;
	}
	char_buf[16] = '\0';
	return 0;
}

/*
For really large inputs we'll want to use base-10^19 for our mod, thus processing nearly one 64-bit
chunk at a time and cutting the number of expensive % operations by 19. But this will also require
us to count leading zeros in the leading (leftmost) base-10^19 word, which isn't worth it for small inputs.
*/

/*
Returns decimal character representation of a base-2^64 2-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint128_base10_char(char char_buf[], uint128 q128)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^128 has 39 decimal digits: */
	uint32 MAX_DIGITS = 39;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q128.d0 || q128.d1) || n_dec_digits == 0)
		{
			c = x128_div_y32(&q128, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/*
Returns decimal character representation of a base-2^64 3-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint192_base10_char(char char_buf[], uint192 q192)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^192 has 58 decimal digits: */
	uint32 MAX_DIGITS = 58;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q192.d0 || q192.d1 || q192.d2) || n_dec_digits == 0)
		{
			c = x192_div_y32(&q192, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/*
Returns decimal character representation of a base-2^64 4-word unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_uint256_base10_char(char char_buf[], uint256 q256)
{
	uint32 i, n_dec_digits = 0;
	char c;
	/* 2^256 has 78 decimal digits: */
	uint32 MAX_DIGITS = 78;

	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if((q256.d0 || q256.d1 || q256.d2 || q256.d3) || n_dec_digits == 0)
		{
			c = x256_div_y32(&q256, (uint32)10) + CHAROFFSET;
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}

	return (int)MAX_DIGITS-n_dec_digits;
}

/********************/
/* Basically a specialized version of the <stdlib.h> strtod function: */
double	convert_base10_char_double (char*char_buf)
{
	uint64 curr_sum = (uint64)0;
	double curr_mul = 0.0;
	uint32 i;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit, hi;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered. If we encounter a decimal point,
	the curr_mul multiplier is set = 1.0 and multiplied by 0.1 for every
	numeric digit found to the right of the DP.
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '.')	/* Found a decimal point */
			{
				ASSERT(HERE, curr_mul == 0.0,"curr_mul == 0.0");	/* Make sure this is the first . we've encountered */
				curr_mul = 1.0;
				continue;
			}
			else if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_double: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, curr_mul == 0.0,"curr_mul == 0.0");
			}
		}
		curr_mul *= 0.1;	/* Only has an effect if we're to the right of the DP */
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"convert_base10_char_double: curr_digit < 10");
		/* Store 10*currsum in a 128-bit product, so can check for overflow: */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64((uint64)10,curr_sum,&curr_sum,&hi);
	#else
		MUL_LOHI64((uint64)10,curr_sum, curr_sum, hi);
	#endif
		if(hi != 0)
		{
			fprintf(stderr, "ERROR: Mul-by-10 overflows in convert_base10_char_double: Offending input string = %s\n", char_buf);
			ASSERT(HERE, 0,"0");
		}
		curr_sum += curr_digit;	/* Since currsum now a multiple of 10, adding a single digit at the low end can't overflow */
	}

	/* If we encountered no DP we simply convert the pure-integer curr_sum to double
	and return that; otherwise we return (double)curr_sum*curr_mul .
	*/
#if 0
	printf("convert_base10_char_double: char_buf = %s, curr_sum = %llu, curr_mul = %lf\n",char_buf, curr_sum, curr_mul);
#endif
	if(curr_mul == 0.0)
	{
		curr_mul = (double)curr_sum;
	}
	else
	{
		curr_mul *= (double)curr_sum;
	}

	return curr_mul;
}

/********************/
/* Basically a 64-bit version of the <stdlib.h> strtoul function: */
uint64 convert_base10_char_uint64 (char*char_buf)
{
	uint64 curr_sum = (uint64)0;
	uint32 i;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit, hi;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint64: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"convert_base10_char_uint64: curr_digit < 10");
		/* Store 10*currsum in a 128-bit product, so can check for overflow: */
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64((uint64)10,curr_sum,&curr_sum,&hi);
	#else
		MUL_LOHI64((uint64)10,curr_sum, curr_sum, hi);
	#endif
		if(hi != 0)
		{
			fprintf(stderr, "ERROR: Mul-by-10 overflows in convert_base10_char_uint64: Offending input string = %s\n", char_buf);
			ASSERT(HERE, 0,"0");
		}
		curr_sum += curr_digit;	/* Since currsum now a multiple of 10, adding a single digit at the low end can't overflow */
	}

	return curr_sum;
}

uint128	convert_base10_char_uint128(char*char_buf)
{
	const uint32 LEN_MAX = 2;
	uint64 curr_sum[2] = {(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint128 x128;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint128: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT128: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x128.d0 = curr_sum[0];
	x128.d1 = curr_sum[1];
	return x128;
}

uint192	convert_base10_char_uint192(char*char_buf)
{
	const uint32 LEN_MAX = 3;
	uint64 curr_sum[3] = {(uint64)0,(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint192 x192;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint192: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT192: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x192.d0 = curr_sum[0];
	x192.d1 = curr_sum[1];
	x192.d2 = curr_sum[2];
	return x192;
}

uint256	convert_base10_char_uint256(char*char_buf)
{
	const uint32 LEN_MAX = 4;
	uint64 curr_sum[4] = {(uint64)0,(uint64)0,(uint64)0,(uint64)0};
	uint64 tmp = 0;
	uint256 x256;
	uint32 i, len = 1;
	int done_with_leading_whitespace = FALSE;
	char c;
	uint64 curr_digit;

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	for(i=0; i != 0xffffffff; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			if(isspace(c))
			{
				if(done_with_leading_whitespace)
					break;
				else
					continue;
			}

			done_with_leading_whitespace = TRUE;

			if(c == '\n' || c == '\0')
			{
				break;
			}
			else
			{
				fprintf(stderr,"convert_base10_char_uint256: isdigit(c) fails, s = %s, i = %u, c = %c\n", char_buf, i, c);
				ASSERT(HERE, 0,"0");
			}
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(curr_sum, (uint64)10, curr_sum, len);
		if(tmp != 0)
		{
			if(len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in CONVERT_BASE10_CHAR_UINT256: Offending input string = %s\n", char_buf);
				ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
			}
			curr_sum[len++] = tmp;
		}

		len += mi64_add_scalar(curr_sum, curr_digit, curr_sum, len);
		ASSERT(HERE, len <= LEN_MAX,"len <= LEN_MAX");
	}

	x256.d0 = curr_sum[0];
	x256.d1 = curr_sum[1];
	x256.d2 = curr_sum[2];
	x256.d3 = curr_sum[3];
	return x256;
}

/***********************/

/* Functions for 96/128/160/192-bit unsigned integer selected-bit testing: */
uint64	TEST_BIT96 (uint96 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT96(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT128(uint128 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT128(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT160(uint160 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT160(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT192(uint192 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT192(__x, __bit, __x);
	return (__x.d0 & 1);
}

uint64	TEST_BIT256(uint256 __x, uint32 __bit)
{
	/* Since call by value, can overwrite __x here: */
	RSHIFT256(__x, __bit, __x);
	return (__x.d0 & 1);
}

/***********************/

/* Given an IEEE-compliant normalized 64-bit float x, generates an approximate
floating-point inverse accurate to at least (numbits) bits of precision. */
double	finvest(double x, uint32 numbits)
{
	/* Used to store MS 8 non-hidden mantissa bits. We'd need to use a 16-bit int
	to allow for the possibility of a carryout (i.e. result = 256) from rounding
	the 9th-most-significant NHB into the upper 8 (which would involve
	additional logic to handle), we instead deal with the issue of rounding
	by assuming the midpoint - e.g. if truncating to the MS 8 NHBs yields
	a certain integer in [0,255], we assume the resulting roundoff error
	is always 0.5, i.e. our precomputed 1/x values are approximations to
	the resulting midpoints. This also avoids our having to treat an input
	of 1.00000000 as a special case, since we munge that to 1.000000001,
	whose inverse is < 1.0: */
	uint32 byteval;
	int ediff;
	uint32 nacc;
	uint64 itmp, mant, exp;
	double ftmp0, ftmp, err_num, err_den;

	/* Max. precision is 53 bits: */
	if(numbits > 53)
	{
		numbits = 53;
	}

	/* Unpack double into a uint64: */
	itmp = *(uint64 *)&x;
	/* Separate upper part of the significand from the sign/exponent fields: */
	exp  = (itmp >> 52) & MASK_EXP;
	mant =  itmp        & MASK_MANT;
	/* Make sure number is normalized: */
	ASSERT(HERE, exp != 0,"finvest: denormalized inputs illegal!");

	/* Store most-significant 8 non-hidden bits: */
	byteval = (mant >> 44) & 0x000000ff;

	/* Munge the exponent to get the inverse's exponent: double-precision
	1.0 has exponent 1023 and is its own inverse, so that is the corner case:
	Numbers in (1.0, 2.0) have exp = 1023 and inverses with exp = 1022, but
	1.0 *exactly* has inverse with exp = 1023. However, our approximate-midpoint
	scheme obviates the need for extra logic to handle this case - 1.0 gets
	approximated as 1.0 + {small}: */
	ediff = (int)exp - 1023;
	exp = (uint64)(1022 - ediff);

	/* Now get the approx-inverse byte and stick it into the most-significant
	8 non-hidden bits of the mantissa field: */
	mant = (uint64)byte_lookup_finvest[byteval] << 44;

	itmp = (itmp & MASK_SIGN) + (exp << 52) + mant;
	ftmp = *(double *)&itmp;

	/* Do as many Newton iterations as required - number of correct
	bits approximately doubles each iteration. The iteration we use
	for y = 1/x is
					y_{n+1} = y_n*[2 - x*y_n] ,
	which is nice, as it involves no divisions.
	*/
	/* Starting # of correct bits from table lookup = 8: */
	nacc = 8;
ftmp0 = ftmp;
	while(nacc < numbits)
	{
		ftmp = ftmp*(2.0 - x*ftmp);
		nacc += nacc;
	}
	err_num = ftmp - ftmp0;
	err_den = ftmp + ftmp0;
	if(fabs(err_num)/fabs(err_den) >= 2e-3)
	{
		sprintf(cbuf, "finvtest: ftmp0 too inaccurate! ftmp = %e, ftmp0 = %e, relerr = %e\n", ftmp, ftmp0,fabs(err_num)/fabs(err_den));
		ASSERT(HERE, 0, cbuf);
	}

	return ftmp;
}

/* Given an IEEE-compliant normalized 64-bit float x, generates an approximate
floating-point inverse square root accurate to at least (numbits) bits of precision.
This routine is very similar to finvest, so see the comments there for details. */
double	fisqrtest(double x, uint32 numbits)
{
	uint32 byteval;
	int ediff;
	uint32 nacc;
	uint64 itmp, mant, exp;
	double ftmp0, ftmp, err_num, err_den;

	/* Max. precision is 53 bits: */
	if(numbits > 53)
	{
		numbits = 53;
	}

	/* Unpack double into a uint64: */
	itmp = *(uint64 *)&x;
	/* Separate upper part of the significand from the sign/exponent fields: */
	exp  = (itmp >> 52) & MASK_EXP;
	mant =  itmp        & MASK_MANT;
	/* Make sure number is normalized: */
	ASSERT(HERE, exp != 0,"finvest: denormalized inputs illegal!");

	/* Store most-significant 9 non-hidden bits - we'll use either all
	or the high 8 of these, depending on the parity of the exponent: */
	byteval = (mant >> 43) & 0x000001ff;

	/* Munge the exponent to get the inverse square root's exponent: double-precision
	1.0 has exponent 1023 and is its own inverse, so that is the corner case:
	Numbers in (1.0, 4.0) have exp in [1023,1024] and ISQRTs with exp = 1022, but
	1.0 *exactly* has ISQRT with exp = 1023. However, our approximate-midpoint
	scheme obviates the need for extra logic to handle this case - 1.0 gets
	approximated as 1.0 + {small}. However, one additional twist in the 1/sqrt
	case is the asymmetry in the handling of ediff: e.g. 2.0 has ediff = +1 but
	maps to 0.707... with exp = 1022 (i.e. we want 1022 - ediff/2 for inputs > 1),
	but e.g. 0.5 has ediff = -1 but maps to 1/sqrt(0.5) = 1.414... with exp = 1023,
	and 0.3 has ediff = -2 and maps to 1/sqrt(0.3) = 1.825... also with exp = 1023,
	and .25 has ediff = -2 and maps to 1/sqrt(.25) = 2.000..., with exp = 1024,
	i.e. we want 1022 - (ediff-1)/2 for inputs < 1.
	*/
	ediff = (int)exp - 1023;	/* 1023 = 0x3ff */
	if(ediff >= 0)
	{
		exp = (uint64)(1022 - ediff/2);

		/* Since we need to handle mantissas in [1, 4), we differentiate via
		inputs in [1,2) and in [2,4) by examining ediff - if it's even it's the
		former interval and we need do nothing; if odd it's the latter and we
		need to add 2 to the floating version of the mantissa, i.e. 0x100 to byteval: */
		if(ediff & 0x1)
		{
			byteval += 0x100;	/* I realize "byteval" is a misnomer in this case... */
		}
		else
			byteval >>= 1;
	}
	else
	{
		exp = (uint64)(1022 - (ediff-1)/2);

		if(ediff & 0x1)
		{
			byteval += 0x100;	/* I realize "byteval" is a misnomer in this case... */
		}
		else
			byteval >>= 1;
	}

	/* Now get the approx-inverse byte and stick it into the most-significant
	8 non-hidden bits of the mantissa field: */

	mant = (uint64)byte_lookup_fisqrtest[byteval] << 44;

	itmp = (itmp & MASK_SIGN) + (exp << 52) + mant;
	ftmp = *(double *)&itmp;

	/* Do as many Newton iterations as required - number of correct
	bits approximately doubles each iteration. The iteration we use
	for y = 1/sqrt(x) is
					y_{n+1} = y_n*[3 - x*(y_n)^2]/2 ,
	which is nice, as it involves no divisions.
	*/
	/* Starting # of correct bits from table lookup = 8: */
	nacc = 8;
ftmp0 = ftmp;
	while(nacc < numbits)
	{
		ftmp = 0.5*ftmp*(3.0 - x*ftmp*ftmp);
		nacc += nacc;
	}
	err_num = ftmp - ftmp0;
	err_den = ftmp + ftmp0;
	if(fabs(err_num)/fabs(err_den) >= 2e-3)
	{
		sprintf(cbuf, "fisqrtest: ftmp0 too inaccurate! ftmp = %e, ftmp0 = %e, relerr = %e\n", ftmp, ftmp0,fabs(err_num)/fabs(err_den));
		ASSERT(HERE, 0, cbuf);
	}

	return ftmp;
}

/***********************/

/***********************/

