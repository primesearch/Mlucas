/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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
/****************************************************************************
 * Nonstandard floating-point functionality (e.g. fused MUL/ADD)
   is accessed or emulated via this header. Specifically, we currently
   support the following math functions via macros - per our usual macro
   naming convention, any macro whose name starts with __ returns a value
   (of type double unless specifically noted), and all arguments are also
   of type double unless otherwise indicated. Single-precision versions
   of the same macros (if we ever decide they are needed) would append 'S'
   to the same set of macro names:
 *   __FNABS	- Floating Negative Absolute Value
 *   __FNINT	- Floating Convert to Integer Word
 *   __FINT		- Floating Convert to Integer Word with Round toward Zero
 *   __FMADD	- Floating Multiply-Add (Double-Precision)
 *   __FMSUB	- Floating Multiply-Subract (Double-Precision)
 *   __FNMADD	- Floating Negative Multiply-Add (Double-Precision)
 *   __FNMSUB	- Floating Negative Multiply-Subtract (Double-Precision)
 *   __FINVEST	- Floating Reciprocal Estimate
 *   __FISQRTEST- Floating Reciprocal Square Root Estimate
 *   __FSEL		- Floating Select
 ****************************************************************************/
#ifndef float_intrin_h_included
#define float_intrin_h_included

/* This file encapsulates platform-specific-and-generic inline functions
(or macros) for the following double-precision floating-point operations:
*/
/* Macro Name:		Input(s)	Output		Description
	__FABS			a			 |a|		Floating-Point Absolute Value
	__FNABS			a			-|a|		Floating Negative Absolute Value
	__FNINT			a			nint(a)		Floating-Point Round-to-Nearest
	__FINT			a			 int(a)		Floating-Point Round-toward-Zero (a.k.a. integer truncation)
	__FMADD			a,b,c		 a*b + c	Floating Multiply-Add
	__FMSUB			a,b,c		 a*b - c	Floating Multiply-Subtract
	__FNMADD		a,b,c		-a*b - c	Floating Negative Multiply-Add
	__FNMSUB		a,b,c		-a*b + c	Floating Negative Multiply-Subtract
	__FSEL			a,b,c		a>=0?b:c	Floating Select

For now we prefer to define the following 2 in non-macro form in util.c:

	__FINVEST		cf. finvest in util.c	Floating Reciprocal Estimate
	__FISQRTEST		cf. fisqrtest in util.c	Floating Reciprocal Square Root Estimate
*/

/* For now, always prefer non-ASM-based versions of NINT and INT: */
#define	__FNINT(a)		DNINT(a)	/* cf. util.h */
#define	__FINT(a)		(double)((int64)(a))

/* Specific include files and/or ASM macros for various platforms: */
/* PowerPC: */
#if(defined(CPU_IS_PPC))

  #if 0	// ppc_intrinsics.h not currently used in the code - this is test-include code for possible future use only.

	/* Gnu C compiler: */
	#if(defined(COMPILER_TYPE_GCC))
		/* Get the hardware intrinsics file, usually at
		/usr/include/gcc/darwin/default/ppc_intrinsics.h .
		The XLC compiler includes the same intrinsics automatically.
		*/
		#include <ppc_intrinsics.h>
	/* Metrowerks CodeWarrior C compiler: */
	#elif(defined(COMPILER_TYPE_MWERKS))
		#include "ppc_intrinsics.h"
	#endif

	/* Make sure the ppc_intrinsics header file was successfully included: */
	#ifndef _PPC_INTRINSICS_H_
		#error Failed to find ppc_intrinsics.h file!
	#endif

	/* Now alias the needed subset of intrinsics to their global (platform-independent) names */
	#define	__FABS		__fabs
	#define	__FNABS		__fnabs
	#define	__FMADD		__fmadd
	#define	__FMSUB		__fmsub
	#define	__FNMADD	__fnmadd
	#define	__FNMSUB	__fnmsub
	#define	__FSEL		__fsel
	/* Use functions finvest and fisqrtest in util.c instead:
	#define	__FINVEST	__fres
	#define	__FISQRTEST	__frsqrte
	*/

  #else

	/* Use ASM macros to access the needed hardware functions: */
	/* double-precision |a| : */
	#define __FABS(a)\
	({\
	  double result;\
	  __asm__ ("fabs %0, %1" : "=f" (result) : "f" (a)); result;\
	})
	/* double-precision -|a| : */
	#define __FNABS(a)\
	({\
	  double result;\
	  __asm__ ("fnabs %0, %1" : "=f" (result) : "f" (a)); result;\
	})
	/* double-precision (a*b + c) : */
	#define __FMADD(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fmadd %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision (a*b - c) : */
	#define __FMSUB(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fmsub %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision -(a*b + c) : */
	#define __FNMADD(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fnmadd %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision -(a*b - c) : */
	#define __FNMSUB(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fnmsub %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision select: a>=0?b:c : */
	#define __FSEL(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fsel %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})

  #endif

/* Itanium: */
#elif(defined(CPU_IS_IA64))

	/* Use ASM macros to access the needed hardware functions: */
	/* double-precision |a| : */
	#define __FABS(a)\
	({\
	  double result;\
	  __asm__ ("fabs %0, %1" : "=f" (result) : "f" (a)); result;\
	})
	/* double-precision -|a| : */
	#define __FNABS(a)\
	({\
	  double result;\
	  __asm__ ("fnegabs %0, %1" : "=f" (result) : "f" (a)); result;\
	})
	/* double-precision (a*b + c) : */
	#define __FMADD(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fma.pc.s0 %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision (a*b - c) : */
	#define __FMSUB(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fms.pc.s0 %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision -(a*b + c) : */
	#define __FNMADD(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fnma.pc.s0 %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); result;\
	})
	/* double-precision -(a*b - c) : Note ia64 has no hardware fnms instruction, so we use -fms */
	#define __FNMSUB(a,b,c)\
	({\
	  double result;\
	  __asm__ ("fnmsub %0, %1, %2, %3" : "=f" (result) : "f" (a), "f" (b), "f" (c)); -result;\
	})
	/* double-precision select: a>=0?b:c : */
	#define	__FSEL(a,b,c)	((a) >= 0 ? (b) : (c))

/* For all other platforms, substitute generic equivalents: */
#else

	#define	__FABS(a)		((a) >= 0 ? (a) : (-a))
	#define	__FNABS(a)		((a) >= 0 ? (-a) : (a))
	// Fused-multiply-add|sub emulation:
	#define	__FMADD(a,b,c)	 (a)*(b) + (c)
	#define	__FMSUB(a,b,c)	 (a)*(b) - (c)
	// Negated-product versions of FMA/FMS:
	#define	__FNMADD(a,b,c)	 (c) - (a)*(b)
	#define	__FNMSUB(a,b,c)	-(a)*(b) - (c)
	#define	__FSEL(a,b,c)	((a) >= 0 ? (b) : (c))

	// Name aliases for FMA instructions:

	// 0th quartet suitable for rvalues in 4-operand FMA4 syntax, as supported for SIMD math on SSE5-supporting AMD CPUs:
	#define  FMA4(a,b,c)	__FMADD(a,b,c)
	#define  FMS4(a,b,c)	__FMSUB(a,b,c)
	#define FNMA4(a,b,c)	__FNMADD(a,b,c)
	#define FNMS4(a,b,c)	__FNMSUB(a,b,c)

	// Next 3 quartets emulate 3-operand FMA syntax available for SIMD math on AVX2-supporting Intel CPUs.
	// Here we emulate GCC inline-asm syntax, i.e. FMA3(src2/mem2,src1,src0) : src0 = [result of FMA operation]
	// thus the rightmost of the 3 inputs is overwritten with the result. Each quartet corr. to one of the 3
	// distinct "flavors" of FMA, which differ in which 2 operands are multiplied together and which added.

// 1. VFMADD132 src0, src1, src2/mem2 : src0 = src0*[src2/mem2] + src1	[note how 132-mnemonic = src-operand index order in RHS]
//	i.e. uses the src0 multiplicand as both a source and as the destination, src1 is addend.
/*** Note this differs from the Intel FMA3 ISA, which defines FMA132(c,a,b) as b = +- b.c +- a.
	Thus when translating code prototyped using the 132-macros below, must swap rightmost 2 args [AT&T/GCC syntax]. ***/
	#define  FMA132(c,a,b)	a =  __FMADD(a,c,b)	// a = + a.c + b
	#define  FMS132(c,a,b)	a =  __FMSUB(a,c,b)	// a = + a.c - b
	#define FNMA132(c,a,b)	a = __FNMADD(a,c,b)	// a = - a.c + b
	#define FNMS132(c,a,b)	a = __FNMSUB(a,c,b)	// a = - a.c - b

// 2. VFMADD213 src0, src1, src2/mem2 : src0 = src0*src1 + [src2/mem2]	[why not call this one ...123?]
//	i.e. uses the src0 multiplicand as both a source and as the destination, [src2/mem2] is addend:
	#define  FMA213(a,b,c)	a =  __FMADD(a,b,c)	// a = + a.b + c
	#define  FMS213(a,b,c)	a =  __FMSUB(a,b,c)	// a = + a.b - c
	#define FNMA213(a,b,c)	a = __FNMADD(a,b,c)	// a = - a.b + c
	#define FNMS213(a,b,c)	a = __FNMSUB(a,b,c)	// a = - a.b - c

// 3. VFMADD231 src0, src1, src2/mem2 : src0 = src1*[src2/mem2] + src0
//	i.e. uses the src0 addend as both a source and as the destination:
	#define  FMA231(c,b,a)	a =  __FMADD(b,c,a)	// a = + b.c + a
	#define  FMS231(c,b,a)	a =  __FMSUB(b,c,a)	// a = + b.c - a
	#define FNMA231(c,b,a)	a = __FNMADD(b,c,a)	// a = - b.c + a
	#define FNMS231(c,b,a)	a = __FNMSUB(b,c,a)	// a = - b.c - a

// In each case, the first 2 of the 3 input args are the multiplicands, and the rightmost digit in the 3-digit
// instruction mnemonic indicates which of the 3 input args is overwritten, in right-to-left "Nth arg from the right"
// fashion. In our notation above, it is always the 'a'-input which is overwritten.

#endif	/* if(defined("your architecture here")) */

#endif	/* float_intrin_h_included */

