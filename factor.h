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
#ifndef factor_h_included
#define factor_h_included

#include "util.h"

/********
PLEASE REFER TO FACTOR.C FOR A DESCRIPTION OF THE APPLICABLE #DEFINES

	(except for FACTOR_PASS_MAX and MAX_BITS_P|Q, which are defined in Mdata.h,
	 and MAX_IMUL_WORDS and MUL_LOHI64_SUBROUTINE, defined in imul_macro.h).
********/

#undef P1WORD
#if(!defined(P2WORD) && !defined(P3WORD) && !defined(P4WORD) && !defined(NWORD))
	#define P1WORD
#endif

#if(defined(NWORD) || defined(P4WORD))
	#undef TRYQ
	#define TRYQ	1
#endif

#ifdef USE_FLOAT
	/* This will need to be made platform-dependent at some point: */
  #ifndef TRYQ
	#define TRYQ	4
  #elif(TRYQ != 1 && TRYQ != 2 && TRYQ != 4)
	#error USE_FLOAT option requires TRYQ = 1, 2 or 4
  #endif

	/* FP-based modmul currently only supported for one-word p's: */
	#ifdef P2WORD
		#error P2WORD may not be used together with USE_FLOAT!
	#endif
	#ifdef P3WORD
		#error P3WORD may not be used together with USE_FLOAT!
	#endif
	#ifdef P4WORD
		#error P4WORD may not be used together with USE_FLOAT!
	#endif

	#ifdef USE_65BIT
		#error USE_65BIT may not be used together with USE_FLOAT!
	#endif
	#ifdef USE_95BIT
		#error USE_65BIT may not be used together with USE_FLOAT!
	#endif
#elif defined(MUL_LOHI64_SUBROUTINE)
  #ifndef TRYQ
	#define TRYQ	4
  #elif(TRYQ != 1 && TRYQ != 4 && TRYQ != 4)
	#error MUL_LOHI64_SUBROUTINE option requires TRYQ = 1, 4 or 8
  #endif
#endif

#ifdef USE_FMADD
	/* This will need to be made platform-dependent at some point: */
  #ifndef TRYQ
	#define TRYQ	2
  #elif(TRYQ != 1 && TRYQ != 2 && TRYQ != 4)
	#error USE_FMADD option requires TRYQ = 1, 2 or 4
  #endif

	/* FMADD-based modmul currently only supported for one-word p's: */
	#ifdef P2WORD
		#error P2WORD may not be used together with USE_FMADD!
	#endif
	#ifdef P3WORD
		#error P3WORD may not be used together with USE_FMADD!
	#endif
	#ifdef P4WORD
		#error P4WORD may not be used together with USE_FMADD!
	#endif

	#ifdef USE_FLOAT
		#error USE_FLOAT may not be used together with USE_FMADD!
	#endif

	#ifdef USE_65BIT
		#error USE_65BIT may not be used together with USE_FMADD!
	#endif
	#ifdef USE_95BIT
		#error USE_65BIT may not be used together with USE_FMADD!
	#endif
#endif

/* Special-handling #define for 96-bit factor candidates: */
#ifdef USE_128x96
	/* Currently only values of 0,1,2 supported: */
	#if(USE_128x96 > 2)
		#error Unrecognized value of USE_128x96!
	#endif

	/* If 2-word-or-greater factoring not enabled, make sure factors < 2^96: */
	#if((MAX_BITS_Q > 96) && defined(P1WORD))
		#error Factoring > 96 bits requires PxWORD with x > 2 to be defined!
	#endif
#endif

/* Key debug #define off value of EWM_DEBUG (set in masterdefs.h) and FACTOR_STANDALONE, set at compile time: */
#undef FAC_DEBUG
#if EWM_DEBUG && defined(FACTOR_STANDALONE)
	#define FAC_DEBUG	1
#else
	#define FAC_DEBUG	0	/* Set == 1 to enable more-extensive self-checking in test_fac() */
#endif

#undef DBG_SIEVE
#if FAC_DEBUG
	#define	DBG_SIEVE	0	/* Set == 1 to enable sieve debugging (Note that no actual trial-divides will occur in this case) */
#endif

#ifndef TRYQ

	#if DBG_SIEVE
		#define TRYQ	0
	#elif defined(INTEGER_MUL_32)
		#define TRYQ	4
	#else
		#define TRYQ	8	/* # of candidates at a time to try. Must be of form 2^k, or zero to skip trial division step
					(e.g. to test sieve-related timings.) A value of 8 seems optimal on the Alpha 21264, which is
					reasonable - TRYQ < 8 may not allow all the integer MULs to complete by the time their results
					are needed, whereas TRYQ = 16 needs more registers than are available and at the same time
					doesn't pipeline significantly better than TRYQ = 8. */
	#endif

#endif

/* Don't have 16-input versions of the twopmodq routines for q > 2^63, so only allow TRYQ up to 8. */
#ifndef TRYQ
	#error TRYQ undefined!
#endif
#if(TRYQ != 0 && TRYQ != 1 && TRYQ != 2 && TRYQ != 4 && TRYQ != 8)
	#error	Illegal value of TRYQ
#endif
#if(TRYQ == 2 && !defined(USE_FLOAT) && !defined(USE_FMADD))
	#error	TRYQ = 2 only allowed if USE_FLOAT is defined
#endif

/* Make sure the TRYQ = 4 fused macros are only used on the PPC32: */
#ifdef MOD_INI_Q4
	#ifndef CPU_SUBTYPE_PPC32
		#error TRYQ = 4 fused macros should only used on the PPC!
	#endif
#endif

/* Factoring-only globals: */
extern int restart;
extern uint64 checksum1;	/* Sum (mod 2^64) of all q's tried during a predetermined interval */
extern uint64 checksum2;	/* Sum (mod 2^64) of all 2^p mod q's for  a predetermined interval */

#ifdef FACTOR_STANDALONE
	/* Declare a blank STATFILE string to ease program logic: */
	extern char STATFILE[];
#endif

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* factor.c: */

#ifndef FACTOR_STANDALONE
	/* exponents > 64 bits require standalone-mode build: */
	int factor(char *pstring, double log2_min_factor, double log2_max_factor);
#endif

int		test_fac(void);

uint64	test_modsqr64    (uint64  x, uint64  q);
uint96	test_modsqr96    (uint96  x, uint96  q);
uint128	test_modsqr128_96(uint128 x, uint128 q);
uint128	test_modsqr128   (uint128 x, uint128 q);

uint64	twopmodq63    (uint64 p, uint64 q);
uint64	twopmodq63_q4 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3);
uint64	twopmodq63_q8 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);
uint64	twopmodq63_x8 (          uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);
uint64	twopmodq64    (uint64 p, uint64 q);
uint64	twopmodq64_q4 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3);
uint64	twopmodq64_q8 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);
uint64	twopmodq65    (uint64 p, uint64 q);
uint64	twopmodq65_q4 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3);
uint64	twopmodq65_q8 (uint64 p, uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);

void cvt_uint78_3word_double	(uint96 __x, double*__fword0, double*__fword1, double*__fword2);
void cvt78_3word_double_uint96	(double __fword0, double __fword1, double __fword2, uint96*__x);
void mulq78_3word_double		(double __x0,double __x1,double __x2, double __y0,double __y1,double __y2, double*__lo0,double*__lo1,double*__lo2);

uint64	twopmodq78_3WORD_DOUBLE   (uint64 p, uint96 q);
uint64	twopmodq78_3WORD_DOUBLE_q2(uint64 p, uint96 q0, uint96 q1);
uint64	twopmodq78_3WORD_DOUBLE_q4(uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3);

uint64	twopmodq100_2WORD_DOUBLE   (uint64 p, uint96 q);
uint64	twopmodq100_2WORD_DOUBLE_q2(uint64 p, uint96 q0, uint96 q1);
uint64	twopmodq100_2WORD_DOUBLE_q4(uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3);

uint96	twopmodq96    (uint64 p, uint96 q);
uint64	twopmodq96_q4 (uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3);
uint64	twopmodq96_q8 (uint64 p, uint96 q0, uint96 q1, uint96 q2, uint96 q3, uint96 q4, uint96 q5, uint96 q6, uint96 q7);
uint64	twopmodq96_q8_const_qhi(uint64 p, uint64 qhi, uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);

uint128	twopmodq128_96   (uint64 p, uint128 q);
uint64	twopmodq128_96_q4(uint64 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3);
uint64	twopmodq128_96_q8(uint64 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3, uint128 q4, uint128 q5, uint128 q6, uint128 q7);
uint64	twopmodq128_96_q8_const_qhi(uint64 p, uint64 qhi, uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint64 q4, uint64 q5, uint64 q6, uint64 q7);

uint128	twopmodq128   (uint64 p, uint128 q);

uint128	twopmodq128x2	(uint64 *p, uint128 q);
uint64	twopmodq128_q4	(uint128 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3);
uint64	twopmodq128_q8	(uint128 p, uint128 q0, uint128 q1, uint128 q2, uint128 q3, uint128 q4, uint128 q5, uint128 q6, uint128 q7);

uint160	twopmodq160   (uint160 p, uint160 q);
uint64	twopmodq160_q4(uint160 p, uint160 q0, uint160 q1, uint160 q2, uint160 q3);
uint64	twopmodq160_q8(uint160 p, uint160 q0, uint160 q1, uint160 q2, uint160 q3, uint160 q4, uint160 q5, uint160 q6, uint160 q7);

uint192	twopmodq192   (uint192 p, uint192 q);
uint64	twopmodq192_q4(uint192 p, uint192 q0, uint192 q1, uint192 q2, uint192 q3);
uint64	twopmodq192_q8(uint192 p, uint192 q0, uint192 q1, uint192 q2, uint192 q3, uint192 q4, uint192 q5, uint192 q6, uint192 q7);

uint256	twopmodq256   (uint256 p, uint256 q);
uint64	twopmodq256_q4(uint256 p, uint256 q0, uint256 q1, uint256 q2, uint256 q3);
uint64	twopmodq256_q8(uint256 p, uint256 q0, uint256 q1, uint256 q2, uint256 q3, uint256 q4, uint256 q5, uint256 q6, uint256 q7);

uint32	CHECK_PKMOD60(uint32 pmod60, uint32 kmod60);

#endif	/* factor_h_included */

