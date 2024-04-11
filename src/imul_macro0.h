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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef imul_macro0_h_included
#define imul_macro0_h_included

#include "util.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
Wide integer multiply macros, with ASM fragments to access efficient non-standard-C
machine instructions wherever feasible. Also includes macros for wide-integer add,
subtract, < > <= == comparison, and bitwise logical shift.

There are multiple levels of related header files here:

	imul_macro0.h - "Tier 0" macros, i.e. macros which contain no local variables and call no other macros or only other Tier 0 macros

	imul_macro1.h - "Tier 1" macros, i.e. macros which contain 1 nesting of local variables, either by having local variables themselves or by calling another Tier 1 macro (but the resulting nested macro must genuinely have one level of local variables)

	imul_macro2.h - "Tier 2" macros, i.e. macros which contain 2 nestings of local variables, etc.

The reason we do this is that some compilers (e.g. MSVC) screw up on variable-scoping
when 2 nested macros both use a same-named local variable, and the outer macro passes
that local variable as an argument to the inner macro, for instance the following code
will cause MSVC to give a "warning C4700: uninitialized variable '_t' used" message
and a likely runtime error because the compiler does not properly treat the 2 instances
of '_t' in the respective macros as distinct local variable.

	// This is the tier 0 (inner) macro:
	#define FOO(x,y)\
	{\
		int _t;\
		_t = x*x;\
		y = _t;\
	}

	// This is the tier 1 (outer) macro:
	#define BAR(a,b,c)\
	{\
		int _t;\
		_t = a+b;\
		FOO(_t,c);\
	}

By using Tiered macros, we can adopt a simple local-variable naming convention
that avoids these problems, e.g. prepend a unique tier-dependent number of _ to
macro-local variable names. In our case we prepend _ for Tier 0, __ for Tier 1, etc,
so the above becomes

	// This is the tier 0 (inner) macro:
	#define FOO(x,y)\
	{\
		int _t;\
		_t = x*x;\
		y = _t;\
	}

	// This is the tier 1 (outer) macro:
	#define BAR(a,b,c)\
	{\
		int __t;\
		__t = a+b;\
		FOO(__t,c);\
	}

and now there is no longer any potential variable-name collision, either within the macros
or the with functions using them (if we declare no _-prepended variables local to the functions).
*/

#undef	MAX_IMUL_WORDS
#define MAX_IMUL_WORDS	3	/* Max. # of 64-bit words allowed for multiplicands;
							   64 x this should match the bit count of the largest MUL macro. */

/* Integer multiply macros. We define these to take advantage of specialized
   hardware instructions wherever possible. All operations are designed to be
   overwrite-safe, i.e. to yield correct results in the case where one or more
   of the input and output addresses coincide. For example, in the macro
   MUL_LOHI64(_x,_y,_lo,_hi), this requires us to store the _lo output in a temp
   before calculating the _hi output, since the former operation potentially
   overwrites _x and/or _y.

	Arguments of form _arg are presumed to be unsigned 64-bit ints.
	Arguments of form ***32 are presumed to be unsigned 32-bit ints.

	__MULL32	(x32,y32     )	    Lower 32 bits of  64-bit product of x32 and y32.
	__MULH32	(x32,y32     )	    Upper 32 bits of  64-bit product of x32 and y32.
	MULL32		(x32,y32,lo32)		Lower 32 bits of  64-bit product of x32 and y32 returned in lo32. (Only need in cases where the compiler doesn't properly handle (uint32)(x32 * y32).)
	MULH32		(x32,y32,hi32)		Upper 32 bits of  64-bit product of x32 and y32 returned in hi32.
	MUL_LOHI32	(_x,_y,_lo,_hi)		Lower/Upper 32 bits of  64-bit product of _x and _y returned in _lo and _hi, respectively.
	SQR_LOHI32	(_x,   _lo,_hi)		Lower/Upper 32 bits of  64-bit square  of _x        returned in _lo and _hi, respectively.

	__MULL64	(_x,_y       )	Lower       64 bits of 128-bit product of _x and _y.
	__MULH64	(_x,_y       )	      Upper 64 bits of 128-bit product of _x and _y.
	MULL64		(_x,_y,_lo   )	Lower       64 bits of 128-bit product of _x and _y returned in _lo.
	MULH64		(_x,_y,   _hi)	      Upper 64 bits of 128-bit product of _x and _y returned in          _hi.
	MUL64x32	(_x,_y,_lo,_hi)		Lower 64 and Upper 32 bits of 96-bit product of _x < 2^64 and _y < 2^32.
	MUL_LOHI64	(_x,_y,_lo,_hi)		Lower/Upper 64 bits of 128-bit product of _x and _y returned in _lo and _hi, respectively.
	SQR_LOHI64	(_x,   _lo,_hi)		Lower/Upper 64 bits of 128-bit square  of _x        returned in _lo and _hi, respectively.
*/
#undef __MULL32
#undef __MULH32
#undef MULL32
#undef MULH32
#undef MUL_LOHI32
#undef SQR_LOHI32

#undef __MULL64
#undef __MULH64
#undef MULL64
#undef MULH64
#undef MUL64x32
#undef MUL_LOHI64
#undef SQR_LOHI64

/* Workaround for a SunStudio-for-AMD64 compiler bug: */
#if(defined(CPU_IS_X86_64) && (defined(COMPILER_TYPE_SUNC) || defined(COMPILER_TYPE_ICC)))
	#define MUL_LOHI64_SUBROUTINE
#endif

#if(defined(CPU_IS_X86)  && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_ICC)))
//	#define MUL_LOHI64_SUBROUTINE
#endif

/* Subroutine form, in case compiler has a problem with the macro form:. */
#ifdef MUL_LOHI64_SUBROUTINE
  #if defined(VERBOSE_HEADERS) && defined(COMPILER_TYPE_GCC)
	#warning Defined MUL_LOHI64_SUBROUTINE
	#warning Using C/32-bit version of MUL_LOHI32 macro
  #endif

	#define MUL_LOHI32(_x32,_y32,_lo,_hi)\
	{\
		uint64 _x = (uint64)(_x32), _y = (uint64)(_y32);\
		uint64 _tt = _x * _y;\
		_lo = (uint32) _tt;\
		_hi = (uint32)(_tt >> 32);\
	}
	#define SQR_LOHI32(_x,   _lo,_hi)	MUL_LOHI32(_x,_x,_lo,_hi)
	#define __MULL32(x32,y32     )	 ((uint32)(x32)*(uint32)(y32))
	#define __MULH32(x32,y32     )	(((uint32)(x32)*(uint64)(y32)) >> 32)
	#define MULL32(  x32,y32,lo32)	lo32 = (uint32)((uint32)(x32)*(uint32)(y32))
	#define MULH32(  x32,y32,hi32)	hi32 = __MULH32((uint32)(x32),(uint64)(y32))
	/*
	Even though the high output for 64x32-bit is always < 2^32,
	declare y and hi here as 64-bit ints to allow flexibility for caller:
	*/
	void	MUL64x32(  uint64 x, uint64 y, uint64 *lo, uint64 *hi);
	void	MUL_LOHI64(uint64 x, uint64 y, uint64 *lo, uint64 *hi);
	void	SQR_LOHI64(uint64 x,           uint64 *lo, uint64 *hi);
	uint64	__MULL64(  uint64 x, uint64 y);
	uint64	__MULH64(  uint64 x, uint64 y);

	#define	MULL64(_x, _y, _lo)	_lo = __MULL64((_x), (_y))
	#define	MULH64(_x, _y, _hi)	_hi = __MULH64((_x), (_y))

	/* The Apple-style fused macros fail on x86/MSVC, so since their
	use keys off whether MOD_INI_Q4 is #defined, just comment this out for now:
	#define MOD_INI_Q4(\
	 q0,qinv0\
	,q1,qinv1\
	,q2,qinv2\
	,q3,qinv3\
	)\
	{\
		ql0  = (uint32) q0   ;\
		ql1  = (uint32) q1   ;\
		ql2  = (uint32) q2   ;\
		ql3  = (uint32) q3   ;\
		qil0 = (uint32) qinv0;\
		qil1 = (uint32) qinv1;\
		qil2 = (uint32) qinv2;\
		qil3 = (uint32) qinv3;\
		\
		qh0  = (uint32)(q0    >> 32);\
		qh1  = (uint32)(q1    >> 32);\
		qh2  = (uint32)(q2    >> 32);\
		qh3  = (uint32)(q3    >> 32);\
		qih0 = (uint32)(qinv0 >> 32);\
		qih1 = (uint32)(qinv1 >> 32);\
		qih2 = (uint32)(qinv2 >> 32);\
		qih3 = (uint32)(qinv3 >> 32);\
	}
	*/

	/* For each input xj, calculates the following sequence:

		SQR_LOHI64(xj,loj,hij);
		loj = MULL64(loj,qinvj);
		yj  = MULH64(loj,qj);
	*/
	#define MOD_SQR_Q4(\
	 x0,hi0,y0\
	,x1,hi1,y1\
	,x2,hi2,y2\
	,x3,hi3,y3\
	)\
	{\
		uint32 ah0,al0,bh0,bl0,ah1,al1,bh1,bl1,ah2,al2,bh2,bl2,ah3,al3,bh3,bl3;\
		uint32 midh,midl,mjdh,mjdl,ss,tt;	/* Use s0-3, t0-3 here and below to avoid artificial execution serialization? */\
		/*uint64 a,b,c,d;*/\
	\
		/* SQR_LOHI64(xj,loj,hij), loj in (alj,ahj), hij in (blj,bhj) : */\
		al0 = (uint32) x0;\
		al1 = (uint32) x1;\
		al2 = (uint32) x2;\
		al3 = (uint32) x3;\
	\
		ah0 = (uint32)(x0 >> 32);\
		ah1 = (uint32)(x1 >> 32);\
		ah2 = (uint32)(x2 >> 32);\
		ah3 = (uint32)(x3 >> 32);\
	\
		MULH32(ah0,ah0,bh0); bl0 = ah0*ah0; MULH32(ah0,al0,midh); midl = ah0*al0; MULH32(al0,al0,ah0); al0 = al0*al0; ss = (midl << 1); ah0 += ss; tt = (midh << 1) | (midl >> 31); bl0 += (ah0 < ss) + tt; bh0 += (bl0 < tt) + (midh >> 31);\
		MULH32(ah1,ah1,bh1); bl1 = ah1*ah1; MULH32(ah1,al1,midh); midl = ah1*al1; MULH32(al1,al1,ah1); al1 = al1*al1; ss = (midl << 1); ah1 += ss; tt = (midh << 1) | (midl >> 31); bl1 += (ah1 < ss) + tt; bh1 += (bl1 < tt) + (midh >> 31);\
		MULH32(ah2,ah2,bh2); bl2 = ah2*ah2; MULH32(ah2,al2,midh); midl = ah2*al2; MULH32(al2,al2,ah2); al2 = al2*al2; ss = (midl << 1); ah2 += ss; tt = (midh << 1) | (midl >> 31); bl2 += (ah2 < ss) + tt; bh2 += (bl2 < tt) + (midh >> 31);\
		MULH32(ah3,ah3,bh3); bl3 = ah3*ah3; MULH32(ah3,al3,midh); midl = ah3*al3; MULH32(al3,al3,ah3); al3 = al3*al3; ss = (midl << 1); ah3 += ss; tt = (midh << 1) | (midl >> 31); bl3 += (ah3 < ss) + tt; bh3 += (bl3 < tt) + (midh >> 31);\
	\
		/* Use | rather than + here (and below) to prevent compiler from using a 32-bit add-with-carry: */\
		hi0 = (uint64)bl0 | ((uint64)bh0 << 32);\
		hi1 = (uint64)bl1 | ((uint64)bh1 << 32);\
		hi2 = (uint64)bl2 | ((uint64)bh2 << 32);\
		hi3 = (uint64)bl3 | ((uint64)bh3 << 32);\
	\
		/*lo0 = (uint64)al0 + ((uint64)ah0 << 32);	SQR_LOHI64(x0,&a,&b);	if(a != lo0) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,a,lo0);	if(b != hi0) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,b,hi0);*/\
		/*lo1 = (uint64)al1 + ((uint64)ah1 << 32);	SQR_LOHI64(x1,&a,&b);	if(a != lo1) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,a,lo1);	if(b != hi1) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,b,hi1);*/\
		/*lo2 = (uint64)al2 + ((uint64)ah2 << 32);	SQR_LOHI64(x2,&a,&b);	if(a != lo2) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,a,lo2);	if(b != hi2) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,b,hi2);*/\
		/*lo3 = (uint64)al3 + ((uint64)ah3 << 32);	SQR_LOHI64(x3,&a,&b);	if(a != lo3) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,a,lo3);	if(b != hi3) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,b,hi3);*/\
	\
		/* loj = MULL64(loj,qinvj) : */\
	\
		/*a = ((uint64)al0 + ((uint64)ah0 << 32))*qinv0;*/\
		/*b = ((uint64)al1 + ((uint64)ah1 << 32))*qinv1;*/\
		/*c = ((uint64)al2 + ((uint64)ah2 << 32))*qinv2;*/\
		/*d = ((uint64)al3 + ((uint64)ah3 << 32))*qinv3;*/\
	\
		ss = ah0*qil0; tt = al0*qih0; MULH32(al0,qil0,ah0); al0 = al0*qil0; ah0 += ss + tt;\
		ss = ah1*qil1; tt = al1*qih1; MULH32(al1,qil1,ah1); al1 = al1*qil1; ah1 += ss + tt;\
		ss = ah2*qil2; tt = al2*qih2; MULH32(al2,qil2,ah2); al2 = al2*qil2; ah2 += ss + tt;\
		ss = ah3*qil3; tt = al3*qih3; MULH32(al3,qil3,ah3); al3 = al3*qil3; ah3 += ss + tt;\
	\
		/* We really only need to save loj for the 65-bit case: */\
		lo0 = (uint64)al0 | ((uint64)ah0 << 32);	/*if(a != lo0) printf("o0\n");*/\
		lo1 = (uint64)al1 | ((uint64)ah1 << 32);	/*if(b != lo1) printf("o1\n");*/\
		lo2 = (uint64)al2 | ((uint64)ah2 << 32);	/*if(c != lo2) printf("o2\n");*/\
		lo3 = (uint64)al3 | ((uint64)ah3 << 32);	/*if(d != lo3) printf("o3\n");*/\
	\
		/* yj = MULH64(loj,qj) is here, store outputs in 64-bit calling argument : */\
		MULH32(ah0,qh0,bh0); bl0 = ah0*qh0; MULH32(ah0,ql0,midh); midl = ah0*ql0; MULH32(qh0,al0,mjdh); mjdl = qh0*al0; MULH32(ah0,al0,ql0); ss = midl + mjdl; ah0 += ss; tt = (ss < midl) + midh + mjdh; bl0 += (ah0 < ss) + tt; bh0 += (tt < midh) + (bl0 < tt);\
		MULH32(ah1,qh1,bh1); bl1 = ah1*qh1; MULH32(ah1,ql1,midh); midl = ah1*ql1; MULH32(qh1,al1,mjdh); mjdl = qh1*al1; MULH32(ah1,al1,ql1); ss = midl + mjdl; ah1 += ss; tt = (ss < midl) + midh + mjdh; bl1 += (ah1 < ss) + tt; bh1 += (tt < midh) + (bl1 < tt);\
		MULH32(ah2,qh2,bh2); bl2 = ah2*qh2; MULH32(ah2,ql2,midh); midl = ah2*ql2; MULH32(qh2,al2,mjdh); mjdl = qh2*al2; MULH32(ah2,al2,ql2); ss = midl + mjdl; ah2 += ss; tt = (ss < midl) + midh + mjdh; bl2 += (ah2 < ss) + tt; bh2 += (tt < midh) + (bl2 < tt);\
		MULH32(ah3,qh3,bh3); bl3 = ah3*qh3; MULH32(ah3,ql3,midh); midl = ah3*ql3; MULH32(qh3,al3,mjdh); mjdl = qh3*al3; MULH32(ah3,al3,ql3); ss = midl + mjdl; ah3 += ss; tt = (ss < midl) + midh + mjdh; bl3 += (ah3 < ss) + tt; bh3 += (tt < midh) + (bl3 < tt);\
	\
		y0 = (uint64)bl0 | ((uint64)bh0 << 32);\
		y1 = (uint64)bl1 | ((uint64)bh1 << 32);\
		y2 = (uint64)bl2 | ((uint64)bh2 << 32);\
		y3 = (uint64)bl3 | ((uint64)bh3 << 32);\
	\
	/*a = MULH64(lo0,q0);	if(a != y0) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,a,lo0);*/\
	/*a = MULH64(lo1,q1);	if(a != y1) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,a,lo1);*/\
	/*a = MULH64(lo2,q2);	if(a != y2) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,a,lo2);*/\
	/*a = MULH64(lo3,q3);	if(a != y3) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,a,lo3);*/\
	}

/********************************************************************************/
/*** MACHINE-SPECIFIC ASM MACROS FOR 64X64=>128-BIT INTEGER UNSIGNED MULTIPLY ***/
/********************************************************************************/

/* Alpha: */
#elif(defined(CPU_IS_ALFA))

	/* Assume __UMULH has already been #defined in platform.h for this arch: */
	#define MUL_LOHI64(_x,_y,_lo,_hi){uint64 _t = (_x)*(_y); _hi = __UMULH((_x), (_y));	_lo = _t;}
	#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
	#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
	#define	__MULL64(	 _x,_y      )	        (_x)*(_y)
	#define	__MULH64(	 _x,_y      )	__UMULH((_x),(_y))
	#define MULL64(	 _x,_y, _lo     )	_lo =         (_x)*(_y)
	#define MULH64(  _x,_y,      _hi)	_hi = __UMULH((_x),(_y))

/* 64-bit x86 (Sun Studio has a bug related to inlines, so we skip it for now) */
#elif(defined(CPU_IS_X86_64))
  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning X86_64-type CPU detected
  #endif
	/* Intel C: */
	#if(defined(COMPILER_TYPE_ICC))
		#error No 64-bit MUL intrinsics supported for amd64/em64t under Intel C
	/* Gnu C: */
	#elif(defined(COMPILER_TYPE_GCC))
		#define __MULH64(_x, _y)	\
		  ({ uint64 _lo, _hi;		\
		  __asm__("mulq %3" : "=a" (_lo), "=d" (_hi) : "%0" (_x), "rm" (_y)); _hi; })
		#define MUL_LOHI64(_x,_y,_lo,_hi)	__asm__("mulq %3" : "=a" (_lo), "=d" (_hi) : "%0" (_x), "rm" (_y) );
		#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
		#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
		#define MULH64(  _x,_y,     _hi) _hi = __MULH64((_x), (_y))

	/* Sun C: */
	#elif(defined(COMPILER_TYPE_SUNC))

		extern uint64_t __MULH64(uint64_t, uint64_t);
		extern void   MUL_LOHI64(uint64_t, uint64_t, uint64_t, uint64_t);
		#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
		#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
		#define MULH64(  _x,_y,     _hi)	_hi = __MULH64((_x), (_y))

	/* MSVC: */
	#elif(defined(COMPILER_TYPE_MSVC))

		#define MUL_LOHI64(_xin,_yin,_lo,_hi)	\
		{\
			uint64 _x, _y;\
			_x = _xin; _y = _yin;	/* In case either of the inputs is a literal */\
			__asm	movq	rax, _x		\
			__asm	movq	r15, _y		\
			__asm	mulq	r15,rax 	/* Result of rax*_y stored in rdx:rax as _hi:_lo	*/	\
			__asm	movq	_lo, rax	\
			__asm	movq	_hi, rdx	\
		}

		#define MULH64(_x,_y,_hi)	\
		{\
			uint64 _lo;	\
			__asm	movq	rax, _x		\
			__asm	mulq	_y	/* Result of rax*_y stored in rdx:rax as _hi:_lo	*/	\
			__asm	movq	_lo, rax	\
			__asm	movq	_hi, rdx	\
		}

		#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
		#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)

	#else
		#error unknown compiler for AMD64.
	#endif

	#define   MULL64(_x,_y,_lo     )	_lo = (_x)*(_y)
	#define	__MULL64(_x, _y        )	      (_x)*(_y)

/* 32-bit X86, Gnu C or MSVC compiler: */
#elif(defined(CPU_IS_X86))
  #if defined(VERBOSE_HEADERS) && (defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_NVCC))
	#warning X86_32-type CPU detected
  #endif
	#if(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

		/* The 32x32->64 MUL and 32-bit pairwise ADD macros are adapted from GMP longlong.h: */

		/* (GMP's add_ssaaaa, with reordered arguments)
			ADD32PAIR(low_addend_1, high_addend_1, low_addend_2, high_addend_2, low_sum, high_sum)
			adds two 64-bit unsigned integers, each composed of two 32-bit pieces
			HIGH_ADDEND_1 and LOW_ADDEND_1, and HIGH_ADDEND_2 and LOW_ADDEND_2
			respectively.  The result is placed in HIGH_SUM and LOW_SUM.  Overflow
			(i.e. carry out) is not stored anywhere, and is lost.
		*/
		#define ADD32PAIR(al, ah, bl, bh, sl, sh)	\
		  __asm__ (	"addl %5,%1\n\t"					\
		  			"adcl %3,%0\n\t"					\
			   : "=r" ((uint32)(sh)), "=&r" ((uint32)(sl))		\
			   : "%0" ((uint32)(ah)), "g" ((uint32)(bh)),		\
				 "%1" ((uint32)(al)), "g" ((uint32)(bl)))

		/* (GMP's umul_ppmm, with reordered arguments):
			MUL_LOHI32(input1, input2, low_prod, high_prod)
			multiplies two 32-bit unsigned integers MULTIPLER and MULTIPLICAND,
			and generates a two-word product with lower and upper halves in LOW_PROD and HIGH_PROD, respectively.
		*/
		#define MUL_LOHI32(_x,_y,_lo,_hi)\
		  __asm__ ("mull %3"							\
			   : "=a" (_lo), "=d" (_hi)					\
			   : "%0" (_x), "rm" (_y))

		#define MULH32(_x,_y,_hi)\
		  ({ uint32 _lo;	\
		  __asm__ ("mull %3"							\
			   : "=a" (_lo), "=d" (_hi)					\
			   : "%0" (_x), "rm" (_y)); })

		#define __MULH32(_x, _y)	\
		  ({ uint32 _lo, _hi;		\
		  __asm__("mull %3" : "=a" (_lo), "=d" (_hi) : "%0" (_x), "rm" (_y)); _hi; })

		#define __MULL32(x32,y32     )	                ((uint32)(x32)*(uint32)(y32))
		#define MULL32(  x32,y32,lo32)	lo32 = (uint32) ((uint32)(x32)*(uint32)(y32))
		#define SQR_LOHI32(_x,   _lo,_hi)	MUL_LOHI32(_x,_x,_lo,_hi)

		/* Low 64 bits of product of uint64 inputs _x and _y returned in uint64 _lo */
		#define	MULL64(_x,_y,_lo)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _a,_b,_c,_d;\
			uint32 _lo32,_hi32;\
		\
			_ux.u64 = _x;\
			_a = _ux.u32[0];	/* x_lo */\
			_b = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_c = _uy.u32[0];	/* y_lo */\
			_d = _uy.u32[1];	/* y_hi */\
			MUL_LOHI32(_a,_c,_lo32,_hi32);	\
			_hi32 += _a*_d + _b*_c;	\
		\
			_ux.u32[0] = _lo32;\
			_ux.u32[1] = _hi32;\
			_lo= _ux.u64;\
		}

		/* 64x32=>96-bit product algorithm:
		represent the inputs as x = a + b*2^32, y = c ( < 2^32), then do 4
		32-bit MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .

		On the IA32, we do 32x64-bit via 2 32x32-bit MULs (2 mull).

		Even though the high output for 64x32-bit is always < 2^32,
		assume _y and _hi here are 64-bit ints to allow flexibility for caller.
		*/
		#define MUL64x32(_x, _y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _a,_b,_c;\
			uint32 _lo32,_md32,_hi32;\
			uint32 _bclo;\
		\
			_ux.u64 = _x;\
			_a = _ux.u32[0];	/* x_lo */\
			_b = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_c = _uy.u32[0];	/* y_lo */\
		\
			MUL_LOHI32(_a,_c,_lo32,_md32);	\
			MUL_LOHI32(_b,_c,_bclo,_hi32);	\
		\
		{\
		  __asm__ volatile (\
			"movl	%[_bclo], %%eax		\n\t"\
			"addl	%%eax	,%[_md32]	\n\t"/* (0) md32 + (b*c)_lo                     , result (middle 32 bits) in md32, carryout in CF bit. */\
			"movl	%[_hi32], %%eax		\n\t"\
			"adcl	$0	 	, %%eax		\n\t"/* (1)    0 + (b*c)_hi + (carryin from (0)), result ( upper 32 bits) in hi32, carryout in CF bit - should be zero! */\
			"movl	%%eax	,%[_hi32]	\n\t"/* Move high result to hi32 */	\
			: /* outputs: none */\
			: [_bclo] "m" (_bclo)	/* All inputs from memory/register here */\
			 ,[_md32] "m" (_md32)	\
			 ,[_hi32] "m" (_hi32)	\
			: "cc","memory","eax"	/* Clobbered registers */\
			);\
		}\
		\
			_ux.u32[0] = _lo32;\
			_ux.u32[1] = _md32;\
			_lo= _ux.u64;\
			_uy.u32[0] = _hi32;\
			_uy.u32[1] = 0;\
			_hi= _uy.u64;\
		}

		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_,bh_,bl_;\
			uint32 hahbh_,hahbl_,halbh_,halbl_,lahbh_,lahbl_,lalbh_,lalbl_;\
			uint32 sumhh_,sumhl_,sumlh_;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			bl_ = _uy.u32[0];	/* y_lo */\
			bh_ = _uy.u32[1];	/* y_hi */\
		\
		/* EWM: for 32-bit unsigned inputs these have a nominal latency of ?? cycles. */\
		/* For each MUL output we include the order in which it is used in the add/carry */\
		/* code below (6 steps, labeled (0)-(5)), in order to properly schedule the MUL: */\
			MUL_LOHI32(ah_,bh_,lahbh_,hahbh_);	/* (x_hi*y_hi)l,h 4,2 */\
			MUL_LOHI32(al_,bh_,lalbh_,halbh_);	/* (x_lo*y_hi)l,h 0,1 */\
			MUL_LOHI32(ah_,bl_,lahbl_,hahbl_);	/* (x_hi*y_lo)l,h 3,1 */\
			MUL_LOHI32(al_,bl_,lalbl_,halbl_);	/* (x_lo*y_lo)l,h -,0 */\
		/* Bits   0- 31 : lalbl                           */\
		/* Bits  32- 63 :*halbl + lahbl +*lalbh   (sumlh) */\
		/* Bits  64- 95 :*hahbl +*halbh + lahbh   (sumhl) */\
		/* Bits  96-127 : hahbh                   (sumhh) */\
		{\
		  __asm__ volatile (\
			"movl	%[halbl_], %%eax	\n\t"\
			"addl	%[lalbh_], %%eax	\n\t"\
			"movl	%%eax	,%[sumlh_]	\n\t"/* Move high result to hi32 */	\
			"movl	%[hahbl_], %%eax	\n\t"\
			"adcl	%[halbh_], %%eax	\n\t"\
			"movl	%%eax	,%[sumhl_]	\n\t"/* Move high result to hi32 */	\
			"movl	%[hahbh_], %%eax	\n\t"\
			"adcl	$0	 	, %%eax		\n\t"/* (1)    0 + (b*c)_hi + (carryin from (0)), result ( upper 32 bits) in hi32, carryout in CF bit - should be zero! */\
			"movl	%%eax	,%[sumhh_]	\n\t"/* Move high result to hi32 */	\
			"movl	%[lahbl_], %%eax	\n\t"\
			"addl	%%eax	,%[sumlh_]	\n\t"\
		/*	"addl	%[sumlh_], %%eax	\n\t"*/\
		/*	"movl	%%eax	,%[sumlh_]	\n\t"*/\
			"movl	%[lahbh_], %%eax	\n\t"\
			"adcl	%%eax	,%[sumhl_]	\n\t"\
		/*	"adcl	%[sumhl_], %%eax	\n\t"*/\
		/*	"movl	%%eax	,%[sumhl_]	\n\t"*/\
			"adcl	$0	,%[sumhh_]	\n\t"\
		/*	"movl	$0, %%eax	\n\t"*/\
		/*	"adcl	%[sumhh_], %%eax	\n\t"*/\
		/*	"movl	%%eax	,%[sumhh_]	\n\t"*/\
			: /* outputs: none */\
			: [hahbh_] "m" (hahbh_)	/* All inputs from memory/register here */\
			 ,[hahbl_] "m" (hahbl_)	\
			 ,[halbh_] "m" (halbh_)	\
			 ,[halbl_] "m" (halbl_)	\
			 ,[lahbh_] "m" (lahbh_)	\
			 ,[lahbl_] "m" (lahbl_)	\
			 ,[lalbh_] "m" (lalbh_)	\
			 ,[sumhh_] "m" (sumhh_)	\
			 ,[sumhl_] "m" (sumhl_)	\
			 ,[sumlh_] "m" (sumlh_)	\
			: "cc","memory","eax"	/* Clobbered registers */\
			);\
		}\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi= _ux.u64;\
			_uy.u32[0] = lalbl_;\
			_uy.u32[1] = sumlh_;\
			_lo= _uy.u64;\
		}

		#define SQR_LOHI64(_x,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_;\
			uint32 hahah_,lahah_,hahal_,lahal_,halal_,lalal_;\
			uint32 sumhh_,sumhl_,sumlh_;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
		\
			MUL_LOHI32(ah_,ah_,lahah_,hahah_);	/* <64:127> */\
			MUL_LOHI32(ah_,al_,lahal_,hahal_);	/* <32: 95> */\
			MUL_LOHI32(al_,al_,lalal_,halal_);	/* < 0: 63> */\
			/* 2*<lahal,hahal>: */\
			sumlh_ = lahal_ << 1;\
			sumhl_ =(hahal_ << 1) | (lahal_>> 31);\
			sumhh_ = hahah_ + (hahal_>> 31);\
		{\
		  __asm__ volatile (\
			"movl	%[halal_], %%eax	\n\t"\
			"addl	%[sumlh_], %%eax	\n\t"\
			"movl	%%eax	,%[sumlh_]	\n\t"\
			"movl	%[lahah_], %%eax	\n\t"\
			"adcl	%[sumhl_], %%eax	\n\t"\
			"movl	%%eax	,%[sumhl_]	\n\t"\
			"adcl	$0	,%[sumhh_]	\n\t"\
		/*	"movl	$0, %%eax	\n\t"*/\
		/*	"adcl	%[sumhh_], %%eax	\n\t"*/\
		/*	"movl	%%eax	,%[sumhh_]	\n\t"*/\
			: /* outputs: none */\
			: [halal_] "m" (halal_)	/* All inputs from memory/register here */\
			 ,[lahah_] "m" (lahah_)	\
			 ,[sumhh_] "m" (sumhh_)	\
			 ,[sumhl_] "m" (sumhl_)	\
			 ,[sumlh_] "m" (sumlh_)	\
			: "cc","memory","eax"	/* Clobbered registers */\
			);\
		}\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi = _ux.u64;\
			_uy.u32[0] = lalal_;\
			_uy.u32[1] = sumlh_;\
			_lo = _uy.u64;\
		}

	#elif(defined(COMPILER_TYPE_MSVC))

		/* 32-bit implementation of 64-bit unsigned add-sans-carry */
		/******** MSVC asm doesn't like varnames ending in h and l ... e.g. treats al as (long)a ********/
		#define ADD32PAIR(al_, ah_, bl_, bh_, sl_, sh_)	\
		{\
			do {\
			__asm	mov	eax, al_		\
			__asm	add	eax, bl_	/* Add low halves, result in eax, CF bit set to indicate carryout */	\
			__asm	mov	sl_, eax	/* Move low result to sl */	\
			__asm	mov	eax, ah_		\
			__asm	adc	eax, bh_	/* Add high halves + carryin, result in eax, CF bit will be discarded */	\
			__asm	mov	sh_, eax	/* Move high result to sh */	\
			} while(0);\
		}

		#define __MULL32(x32,y32     )	                ((uint32)(x32)*(uint32)(y32))
		#define MULL32(  x32,y32,lo32)	lo32 = (uint32) ((uint32)(x32)*(uint32)(y32))

		/* Multiplies two 32-bit unsigned integers _x and _y,
			and generates a two-word (64-bit) product in _lo and _hi.
		*/
		#define MUL_LOHI32(_x,_y,_lo,_hi)\
		{\
			__asm	mov	eax, _x		\
			__asm	mul	_y	/* Result of eax*_y stored in edx:eax as _hi:_lo	*/	\
			__asm	mov	_lo, eax		\
			__asm	mov	_hi, edx		\
		}
		#define SQR_LOHI32(_x,   _lo,_hi)	MUL_LOHI32(_x,_x,_lo,_hi)

		#define MULH32(_x,_y,_hi)\
		{\
			__asm	mov	eax, _x		\
			__asm	mul	_y			\
			__asm	mov	_hi, edx	\
		}
		#define __MULH32(_x,_y)	(( uint32 _hi; MULH32(_x,_y,_hi); _hi ))

		/* Low 64 bits of product of uint64 inputs _x and _y returned in uint64 _lo */
		#define	MULL64(_x,_y,_lo)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _a,_b,_c,_d;\
			uint32 _lo32,_hi32;\
		\
			_ux.u64 = _x;\
			_a = _ux.u32[0];	/* x_lo */\
			_b = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_c = _uy.u32[0];	/* y_lo */\
			_d = _uy.u32[1];	/* y_hi */\
			MUL_LOHI32(_a,_c,_lo32,_hi32);	\
			_hi32 += _a*_d + _b*_c;	\
		\
			_ux.u32[0] = _lo32;\
			_ux.u32[1] = _hi32;\
			_lo= _ux.u64;\
		}

		/* 64x32=>96-bit product algorithm:
		represent the inputs as x = a + b*2^32, y = c ( < 2^32), then do 4
		32-bit MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .

		On the IA32, we do 32x64-bit via 2 32x32-bit MULs (2 mull).

		Even though the high output for 64x32-bit is always < 2^32,
		assume _y and _hi here are 64-bit ints to allow flexibility for caller.
		*/
		#define MUL64x32(_x, _y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _a,_b,_c;\
			uint32 _lo32,_md32,_hi32;\
			uint32 _bclo;\
		\
			_ux.u64 = _x;\
			_a = _ux.u32[0];	/* x_lo */\
			_b = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_c = _uy.u32[0];	/* y_lo */\
		\
			MUL_LOHI32(_a,_c,_lo32,_md32);	\
			MUL_LOHI32(_b,_c,_bclo,_hi32);	\
		\
			{\
			__asm	mov	eax, _bclo		\
			__asm	add	_md32, eax	/* (0) md32 + (b*c)_lo                     , result (middle 32 bits) in md32, carryout in CF bit. */\
			__asm	mov	eax, _hi32		\
			__asm	adc	eax, 0		/* (1)    0 + (b*c)_hi + (carryin from (0)), result ( upper 32 bits) in hi32, carryout in CF bit - should be zero! */\
			__asm	mov	_hi32, eax	/* Move high result to hi32 */	\
		}\
		\
			_ux.u32[0] = _lo32;\
			_ux.u32[1] = _md32;\
			_lo= _ux.u64;\
			_uy.u32[0] = _hi32;\
			_uy.u32[1] = 0;\
			_hi= _uy.u64;\
		}

		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* IA32 is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_,bh_,bl_;\
			uint32 hahbh_,hahbl_,halbh_,halbl_,lahbh_,lahbl_,lalbh_,lalbl_;\
			uint32 sumhh_,sumhl_,sumlh_;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			bl_ = _uy.u32[0];	/* y_lo */\
			bh_ = _uy.u32[1];	/* y_hi */\
		\
		/* EWM: for 32-bit unsigned inputs these have a nominal latency of ?? cycles. */\
		/* For each MUL output we include the order in which it is used in the add/carry */\
		/* code below (6 steps, labeled (0)-(5)), in order to properly schedule the MUL: */\
			MUL_LOHI32(ah_,bh_,lahbh_,hahbh_);	/* (x_hi*y_hi)l,h 4,2 */\
			MUL_LOHI32(al_,bh_,lalbh_,halbh_);	/* (x_lo*y_hi)l,h 0,1 */\
			MUL_LOHI32(ah_,bl_,lahbl_,hahbl_);	/* (x_hi*y_lo)l,h 3,1 */\
			MUL_LOHI32(al_,bl_,lalbl_,halbl_);	/* (x_lo*y_lo)l,h -,0 */\
		/* Bits   0- 31 : lalbl                           */\
		/* Bits  32- 63 :*halbl + lahbl +*lalbh   (sumlh) */\
		/* Bits  64- 95 :*hahbl +*halbh + lahbh   (sumhl) */\
		/* Bits  96-127 : hahbh                   (sumhh) */\
		\
			{\
			__asm	mov	eax, halbl_	\
			__asm	add	eax, lalbh_	/* (0) halbl + lalbh                       , result (partial sum) in sumlh,carryout in CF bit. */\
			__asm	mov	sumlh_, eax	\
			__asm	mov	eax, hahbl_	\
			__asm	adc	eax, halbh_	/* (1) halbh + hahbl + (carryin from sumlh), result (partial sum) in sumhl,carryout in CF bit. */\
			__asm	mov	sumhl_, eax	\
			__asm	mov	eax, hahbh_	\
			__asm	adc	eax, 0		/* (2)     0 + hahbh + (carryin from sumhl), result (partial sum) in sumhh,carryout in CF) bit - should be zero! */\
			__asm	mov	sumhh_, eax	\
			__asm	mov	eax, lahbl_	\
			__asm	add	sumlh_, eax	/* (3) sumlh + lahbl                       , result (sumlh), carryout in CF bit. */\
			__asm	mov	eax, lahbh_	\
			__asm	adc	sumhl_, eax	/* (4) sumhl + lahbh + (carryin from sumlh), result (sumhl), carryout in CF bit. */\
			__asm	mov	eax, 0		\
			__asm	adc	sumhh_, eax	/* (5)     0 + sumhh + (carryin from sumhl), result (sumhh), carryout in CF bit - should be zero! */\
			}\
		\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi= _ux.u64;\
			_uy.u32[0] = lalbl_;\
			_uy.u32[1] = sumlh_;\
			_lo= _uy.u64;\
		}

		#define SQR_LOHI64(_x,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_;\
			uint32 hahah_,lahah_,hahal_,lahal_,halal_,lalal_;\
			uint32 sumhh_,sumhl_,sumlh_;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
		\
			MUL_LOHI32(ah_,ah_,lahah_,hahah_);	/* <64:127> */\
			MUL_LOHI32(ah_,al_,lahal_,hahal_);	/* <32: 95> */\
			MUL_LOHI32(al_,al_,lalal_,halal_);	/* < 0: 63> */\
			/* 2*<lahal,hahal>: */\
			sumlh_ = lahal_ << 1;\
			sumhl_ =(hahal_ << 1) | (lahal_>> 31);\
			sumhh_ = hahah_ + (hahal_>> 31);\
		\
			{\
 			__asm	mov	eax, halal_	\
			__asm	add	sumlh_, eax	/* (0) sumlh + halal                     , result in sumlh, carryout in CF */\
 			__asm	mov	eax, lahah_	\
			__asm	adc	sumhl_, eax	/* (1) sumhl + lahah + (carryin from (0)), result in sumhl, carryout in CF */\
			__asm	mov	eax, 0		\
			__asm	adc	sumhh_, eax	/* (2)     0 + sumhh + (carryin from (1)), result in sumhh, carryout in CF - should be zero! */\
			}\
		\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi = _ux.u64;\
			_uy.u32[0] = lalal_;\
			_uy.u32[1] = sumlh_;\
			_lo = _uy.u64;\
		}

	#else
		#error unknown compiler for IA32!
	#endif

		/* Generic 128-bit multiply algorithm, returning only the upper 64 bits (high part) of the result.
		Actual calling arguments are assumed to be 64-bit ints - user must make sure this is true.
		Result written into hi.
		*/
		#define MULH64(_x, _y,_hi)\
		{\
			uint64 _tt;					\
			MUL_LOHI64((_x), (_y), _tt, _hi);	\
		}
		#define	__MULL64(_x, _y)	({ uint64 _lo; MULL64((_x), (_y), _lo); _lo; })
		#define	__MULH64(_x, _y)	({ uint64 _lo; MULH64((_x), (_y), _lo); _lo; })

/* Itanium: */
#elif(defined( CPU_IS_IA64 ))
	#if(defined(COMPILER_TYPE_ICC))

		/* Itanium system under Linux. Compile using icc, not cc. */
		#if(defined(OS_TYPE_LINUX))
			/* Functions to get lower and upper 64 bits of 64x64 -> 128-bit unsigned product.
			in the multiply/add versions of these, the add is always into the LOWER 64 bits
			of the product - the result of __MULH64_ADD is only modified by this if the add
			produces a carry out of the low 64 bits of the result, in which case the result
			of __MULH64_ADD(_x, _y, _add) is one larger than that of __MULH64(_x, _y).
			*/
			#define  __MULL64(_x, _y)				_m64_xmalu((_x), (_y),      0 )
			#define  __MULH64(_x, _y)				_m64_xmahu((_x), (_y),      0 )
			#define  __MULL64_ADD(_x, _y, _add)	_m64_xmalu((_x), (_y), (_add))
			#define  __MULH64_ADD(_x, _y, _add)	_m64_xmahu((_x), (_y), (_add))
		#else
			#error unknown OS type for Itanium.
		#endif

	/* IA-64 (Itanium), HP C or C++ compiler for HPUX: Eddie Gornish of HP writes:

	"In order to get the correct behavior, you need to use _Asm_setf and _Asm_getf, such as:

	#define  __MULL64(_x, _y)        (uint64)_Asm_getf(_FR_SIG, _Asm_xma(_XM_LU, _Asm_setf(_FR_SIG, (_x)), _Asm_setf(_FR_SIG, (_y)), (__fpreg)0))

	Casting to an __fpreg gives you a floating-point number on the FP-side (i.e. set/fcvt combination).
	You need a fixed-point number on the FP-side (i.e. just setf).  There is really no type that represents a
	fixed-point number on the FP-side, so simple C-style casting will not work."
	*/
	#elif(defined(COMPILER_TYPE_HPC))
		#define  __MULL64(_x, _y)				(uint64)_Asm_getf(_FR_SIG, _Asm_xma(_XM_LU, _Asm_setf(_FR_SIG, (_x)), _Asm_setf(_FR_SIG, (_y)), (__fpreg)0))
		#define  __MULH64(_x, _y)				(uint64)_Asm_getf(_FR_SIG, _Asm_xma(_XM_HU, _Asm_setf(_FR_SIG, (_x)), _Asm_setf(_FR_SIG, (_y)), (__fpreg)0))
		#define  __MULL64_ADD(_x, _y, _add)	(uint64)_Asm_getf(_FR_SIG, _Asm_xma(_XM_LU, _Asm_setf(_FR_SIG, (_x)), _Asm_setf(_FR_SIG, (_y)), _Asm_setf(_FR_SIG, (_add)))
		#define  __MULH64_ADD(_x, _y, _add)	(uint64)_Asm_getf(_FR_SIG, _Asm_xma(_XM_HU, _Asm_setf(_FR_SIG, (_x)), _Asm_setf(_FR_SIG, (_y)), _Asm_setf(_FR_SIG, (_add)))

	#elif(defined(COMPILER_TYPE_GCC))
		/* Functions to get lower 64 bits of 64x64 -> 128-bit unsigned product: */
		#define __MULL64(_x, _y)	\
		  ({ unsigned long rslt;\
		  __asm__("xma.lu %0 = %1, %2, f0" : "=f"(rslt) : "%f"(_x), "f"(_y)); rslt; })
		#define __MULL64_ADD(_x, _y, add)	\
		  ({ unsigned long rslt;\
		  __asm__("xma.lu %0 = %1, %2, %3" : "=f"(rslt) : "%f"(_x), "f"(_y), "f"(addr)); rslt; })
		/* Functions to get upper 64 bits of 64x64 -> 128-bit unsigned product: */
		#define __MULH64(_x, _y)	\
		  ({ unsigned long rslt;\
		  __asm__("xma.hu %0 = %1, %2, f0" : "=f"(rslt) : "%f"(_x), "f"(_y)); rslt; })
		#define __MULH64_ADD(_x, _y, add)	\
		  ({ unsigned long rslt;\
		  __asm__("xma.hu %0 = %1, %2, %3" : "=f"(rslt) : "%f"(_x), "f"(_y), "f"(addr)); rslt; })
	#else
		#error unknown compiler for Itanium.
	#endif

	/* Prepend an extra _ onto the t-temp here, since temp "_t" is also frequently used in bigger macros: */
	#define MUL_LOHI64(    _x, _y,       _lo,_hi) {uint64 _t = (_x)*(_y); _hi = __MULH64((_x), (_y));	_lo = _t;}
	#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
	#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
	#define MUL_LOHI64_ADD(_x, _y, _add,_lo,_hi) {uint64 _t = __MULL64_ADD((_x), (_y), (_add)); _hi = __MULH64_ADD((_x), (_y), (_add));	_lo = _t;}
	#define SQR_LOHI64_ADD(_x,      _add,_lo,_hi) {uint64 _t = __MULL64_ADD((_x), (_x), (_add)); _hi = __MULH64_ADD((_x), (_x), (_add));	_lo = _t;}

	#define MULL64(_x, _y,_lo     )	_lo = __MULL64((_x), (_y))
	#define MULH64(_x, _y,     _hi)	_hi = __MULH64((_x), (_y))

	#define MULL64_LOW_ADD(_x, _y, _add,_lo     )	_lo = __MULL64_ADD((_x), (_y), (_add))
	#define MULH64_ADD(    _x, _y, _add,     _hi)	_hi = __MULH64_ADD((_x), (_y), (_add))

/* PowerPC: */
#elif(defined(CPU_IS_PPC))
	/* 64-bit (e.g. G5 and later): */
	#if(defined(CPU_SUBTYPE_PPC64))

		/* High 64 bits of unsigned 64-bit (double-word) product (analog of Alpha UMULH): */

		#if(defined(COMPILER_TYPE_GCC))
			#define __MULL64(_x,_y)	\
				({ uint64 _lo;			\
				__asm__("mulld  %0,%1,%2" : "=r"(_lo) : "%r"(_x), "r"(_y));	_lo; })
			#define __MULH64(_x,_y)	\
				({ uint64 _hi;			\
				__asm__("mulhdu %0,%1,%2" : "=r"(_hi) : "%r"(_x), "r"(_y));	_hi; })
		#elif(defined(COMPILER_TYPE_XLC))
			#define	__MULL64(_x,_y)	(_x)*(_y)
			#define __MULH64(_x,_y)	__mulhdu((_x),(_y))
		#else
			#error unknown compiler for PPC64.
		#endif

		/* Prepend an extra _ onto the t-temp here, since temp "_t" is also frequently used in bigger macros: */
		#define MUL_LOHI64(_x,_y,_lo,_hi)	{uint64 _t = (_x)*(_y); _hi = __MULH64((_x), (_y));	_lo = _t;}
		#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
		#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
		#define MULL64(    _x,_y,_lo     )	_lo = (_x)*(_y)
		#define MULH64(    _x,_y,     _hi)	_hi = __MULH64((_x), (_y))

	/* 32-bit (G3 AND G4): */
	#elif(defined(CPU_SUBTYPE_PPC32))

		/* EWM: v18: Needed to prepend _ to names in e.g. MUL_LOHI64 to eliminate compiler errors related to
		var-name collisions, e.g. qfloat.c declaring uint64 a,b,c,d,lo,hi and calling MUL_LOHI64(a,b,lo,hi)
		collided with unnamed union '{ uint32 u32[2]; uint64 u64 } a;' decl in original form of said macro.
		I normally make a habit of using _ and __ local--varname-prefixing in macros, but the ppc32 wide-mul
		code was gifted to me by a ppc-expert coder at the time, and since it built and ran, aside from a few
		comments and other tweaks, I never changed that.

		The source of the problem took some ferreting out, as GCC was flagging the error rather unhelpfully as
			error: assigning to 'uint64' (aka 'unsigned long long') from incompatible type
				  'union <anonymous at [file,line]>'
				...
				  note: instantiated from:
					a.u64 = _x;\
		Well, _x is dummy arg-name standing for a uint64 calling argument, so where's the problem? The problem
		is that if the calling function uses a var named 'a' for the 1st arg of the macro call, that translates as
					a.u64 = a;\
		and the compiler assumes th RHS 'a' refers to the aforementioned same-named locally declared union.
		Now ppc32 builds are no longer a thing as of this writing, but still nice to maintain buildability if possible.
		Also switched from the aforementioned union to 32-bit pointer-based access by way of initial workaround, which
		also still used 'a' a local varname, but now of a local uint64. Once the reason for the original GCC error
		became clear, it also became clear why the pointer-access workaround worked, it effectively resulted in this
		assignment:
			[local uint64 var]a = [calling function uint64]a;\
		So kept that but still prepended _ to all macro-local varnames just by way of it being good standard practice.
		*/
	    /* The following optimized Macros require TRYQ = 4 in factor.h
		   (Thanks to Klaus Kastens for initial versions of these.)
		*/

		/* KK: Returns the high 32 bits of a 32x32-bit unsigned integer product: */
		/* NOTE: This only works reliable if the output isn't the same expression */
		/* as one of the inputs! i.e. MULH32(a,a,b) might not work */

		#if(defined(COMPILER_TYPE_GCC))

			#define MULL32(	_x32,_y32,_lo32) __asm__("mullw  %0,%1,%2" : "=r"(_lo32) : "%r"((uint32)(_x32)), "r"((uint32)(_y32)));
			#define MULH32(	_x32,_y32,_hi32) __asm__("mulhwu %0,%1,%2" : "=r"(_hi32) : "%r"((uint32)(_x32)), "r"((uint32)(_y32)));
			#define __MULL32(	_x32,_y32     )	\
				({ uint32 _lo32;					\
				__asm__("mullw  %0,%1,%2" : "=r"(_lo32) : "%r"((uint32)(_x32)), "r"((uint32)(_y32)));	_lo32; })
			#define __MULH32(	_x32,_y32     )	\
				({ uint32 _hi32;					\
				__asm__("mulhwu %0,%1,%2" : "=r"(_hi32) : "%r"((uint32)(_x32)), "r"((uint32)(_y32)));	_hi32; })

		#elif(defined(COMPILER_TYPE_XLC))
			/* If XLC intrinsic available, prefer that over Gnu-style ASM syntax
			since intrinisc form may allow compiler to do better optimization: */
			/* XLC has no sppecial intrinsics for low-half MUL: */
			#define MULL32(	 _x32,_y32,_lo32) _lo32 = (uint32)((uint32)(_x32)*(uint32)(_y32));
			#define MULH32(	 _x32,_y32,_hi32) _hi32 = __mulhwu((uint32)(_x32),(uint32)(_y32));
			#define __MULL32(_x32,_y32      )         (uint32)((uint32)(_x32)*(uint32)(_y32));
			#define __MULH32(_x32,_y32      )         __mulhwu((uint32)(_x32),(uint32)(_y32));

		#else
			#error unknown compiler for PPC32.
		#endif

		#if 0
			#define	MULL64(_x,_y,_lo)	({\
				uint32 _xlo = (uint32)(_x);			\
				uint32 _ylo = (uint32)(_y);			\
				uint32 _xhi = (uint32)(_x >> 32);	\
				uint32 _yhi = (uint32)(_y >> 32);	\
				\
				_lo = (uint64)__MULL32(_xlo,_ylo)	\
					+(((uint64)__MULH32(_xlo,_ylo)) << 32)	\
					+(((uint64)__MULL32(_xlo,_yhi) + (uint64)__MULL32(_xhi,_ylo)) << 32);	\
				})
		#else
			#define MULL64(_x,_y,_lo)\
			{\
				uint32 _ah,_al,_bh,_bl;\
				uint32 _hahbh,_hahbl,_halbh,_halbl,_lahbh,_lahbl,_lalbh,_lalbl;\
				uint32 _sumhh,_sumhl,_sumlh;\
				uint64 _a = _x,_b = _y;\
				uint32 *_aptr = (uint32 *)&_a, *_bptr = (uint32 *)&_b;\
				/* PPC is big-endian, hence hi32 at aptr, lo32 at aptr+1: */\
				_ah = * _aptr;		/* x_hi */\
				_al = *(_aptr+1);	/* x_lo */\
				_bh = * _bptr;		/* y_hi */\
				_bl = *(_bptr+1);	/* y_lo */\
			\
				MULL32(_al,_bl,_lalbl);	/* (x_lo*y_lo)_lo */\
				MULH32(_al,_bl,_halbl);	/* (x_lo*y_lo)_hi */\
				MULL32(_al,_bh,_lalbh);	/* (x_lo*y_hi)_lo */\
				MULL32(_ah,_bl,_lahbl);	/* (x_hi*y_lo)_lo */\
			/* Bits   0- 31 : lalbl                           */\
			/* Bits  32- 63 : halbl + lahbl + lalbh   (no carryouts needed) */\
				*(_bptr+1) = _lalbl;\
				* _bptr    = _halbl + _lahbl + _lalbh;\
				_lo = _b;\
			}
		#endif

		/* 64x32=>96-bit product algorithm:
		represent the inputs as x = a + b*2^32, y = c ( < 2^32), then do 4
		32-bit MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .

		On the G4, we do 32x64-bit via 4 32x32-bit MULs (2 mullw, 2 mulhw). Taking
		advantage of any potential asymmetries in input sizes is especially important
		because on the G4 only IU2 can do MUL, and each MUL must finish before the
		next can start, so shaving a cycle off the nominal 4-cycle MUL latency for
		via the early-exit multiply feature when the B-operand has high half
		zero (more precisely, high 15 bits all 0 or all 1) is useful.
			In our case _y may be smaller than 32 bits, so we order the 32-bit MUL
		inputs to take advantage of that, by making the smaller input the 2nd input
		argument to the mullw and mulhw calls.

		Even though the high output for 64x32-bit is always < 2^32,
		assume _y and _hi here are 64-bit ints to allow flexibility for caller.
		*/
		#define MUL64x32(_x, _y,_lo,_hi)\
		{\
			uint32 _a,_b,_c;\
			uint32 _lo32,_md32,_hi32;\
			uint32 _bclo;\
			uint64 _a = _x,_b = _y;\
			uint32 *_aptr = (uint32 *)&_a, *_bptr = (uint32 *)&_b;\
			/* PPC is big-endian, hence hi32 at aptr, lo32 at aptr+1: */\
			_b = * _aptr;		/* x_hi */\
			_a = *(_aptr+1);	/* x_lo */\
			_c = *(_bptr+1);	/* y_lo */\
		\
			MULL32(_a,_c,_lo32);	/* (a*c)_lo */\
			MULH32(_a,_c,_md32);	/* (a*c)_hi */\
			MULL32(_b,_c,_bclo);	/* (b*c)_lo */\
			MULH32(_b,_c,_hi32);	/* (b*c)_hi */\
		\
		  __asm__(\
			"addc  %0,%2,%3\n\t"	/* (0) md32 + (b*c)_lo                     , result (middle 32 bits) in %0, carryout in XER(CA) bit. */\
			"addze %1,%4"			/* (1)    0 + (b*c)_hi + (carryin from (0)), result ( upper 32 bits) in %1, carryout in XER(CA) bit - should be zero! */\
			: "=r"(_md32), "=r"(_hi32)\
			: "%0"(_md32), "r"(_bclo), "1"(_hi32));\
		\
			* _aptr    = _md32;\
			*(_aptr+1) = _lo32;\
			* _bptr    = 0;\
			*(_bptr+1) = _hi32;\
			_lo = _a;\
			_hi = _b;\
		}

		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			uint32 _ah,_al,_bh,_bl;\
			uint32 _hahbh,_hahbl,_halbh,_halbl,_lahbh,_lahbl,_lalbh,_lalbl;\
			uint32 _sumhh,_sumhl,_sumlh;\
			uint64 _a = _x,_b = _y;\
			uint32 *_aptr = (uint32 *)&_a, *_bptr = (uint32 *)&_b;\
			/* PPC is big-endian, hence hi32 at aptr, lo32 at aptr+1: */\
			_ah = * _aptr;		/* x_hi */\
			_al = *(_aptr+1);	/* x_lo */\
			_bh = * _bptr;		/* y_hi */\
			_bl = *(_bptr+1);	/* y_lo */\
		\
		/* EWM: for 32-bit unsigned inputs these have a nominal latency of 10 cycles. */\
		/* For each MUL output we include the order in which it is used in the add/carry */\
		/* code below (6 steps, labeled (0)-(5)), in order to properly schedule the MUL: */\
			MULH32(_ah,_bh,_hahbh);	/* (x_hi*y_hi)_hi (2) */\
			MULL32(_ah,_bh,_lahbh);	/* (x_hi*y_hi)_lo (4) */\
			MULH32(_al,_bh,_halbh);	/* (x_lo*y_hi)_hi (1) */\
			MULL32(_al,_bh,_lalbh);	/* (x_lo*y_hi)_lo (0) */\
			MULH32(_ah,_bl,_hahbl);	/* (x_hi*y_lo)_hi (1) */\
			MULL32(_ah,_bl,_lahbl);	/* (x_hi*y_lo)_lo (3) */\
			MULH32(_al,_bl,_halbl);	/* (x_lo*y_lo)_hi (0) */\
			MULL32(_al,_bl,_lalbl);	/* (x_lo*y_lo)_lo (-) */\
		/* Bits   0- 31 : lalbl                           */\
		/* Bits  32- 63 :*halbl + lahbl +*lalbh   (sumlh) */\
		/* Bits  64- 95 :*hahbl +*halbh + lahbh   (sumhl) */\
		/* Bits  96-127 : hahbh                   (sumhh) */\
		\
			/* KK:  addc/adde/addze forces execution serialization! Other solution for carry propagation? */\
			/* EWM: these all have latency 1, so serialization not a problem... */\
		  __asm__(\
			"addc  %0,%3,%4\n\t"	/* (0) halbl + lalbh                       , result (partial sumlh) in %0, carryout in XER(CA) bit. */\
			"adde  %1,%5,%6\n\t"	/* (1) halbh + hahbl + (carryin from sumlh), result (partial sumhl) in %1, carryout in XER(CA) bit. */\
			"addze %2,%7"			/* (2)     0 + hahbh + (carryin from sumhl), result (partial sumhh) in %2, carryout in XER(CA) bit - should be zero! */\
			: "=r"(_sumlh), "=r"(_sumhl), "=r"(_sumhh)\
			: "%r"(_halbl), "r"(_lalbh), "%r"(_halbh), "r"(_hahbl), "r"(_hahbh));\
		\
		  __asm__(\
			"addc  %0,%3,%4\n\t"	/* (3) sumlh + lahbl                       , result (sumlh) in %0, carryout in XER(CA) bit. */\
			"adde  %1,%5,%6\n\t"	/* (4) sumhl + lahbh + (carryin from sumlh), result (sumhl) in %1, carryout in XER(CA) bit. */\
			"addze %2,%7"			/* (5)     0 + sumhh + (carryin from sumhl), result (sumhh) in %2, carryout in XER(CA) bit - should be zero! */\
			: "=r"(_sumlh), "=r"(_sumhl), "=r"(_sumhh)\
			: "%0"(_sumlh), "r"(_lahbl), "%1"(_sumhl), "r"(_lahbh), "2"(_sumhh));\
		\
			* _aptr    = _sumhh;\
			*(_aptr+1) = _sumhl;\
			* _bptr    = _sumlh;\
			*(_bptr+1) = _lalbl;\
			_hi = _a;\
			_lo = _b;\
		}

		#define SQR_LOHI64(_x,_lo,_hi)\
		{\
			uint32 _ah,_al;\
			uint32 _hahah,_lahah,_hahal,_lahal,_halal,_lalal;\
			uint32 _sumhh,_sumhl,_sumlh;\
			uint64 _a = _x,_b;\
			uint32 *_aptr = (uint32 *)&_a, *_bptr = (uint32 *)&_b;\
			/* PPC is big-endian, hence hi32 at aptr, lo32 at aptr+1: */\
			_ah = * _aptr;		/* x_hi */\
			_al = *(_aptr+1);	/* x_lo */\
		\
			MULH32(_ah,_ah,_hahah);\
			MULL32(_ah,_ah,_lahah);\
			MULH32(_ah,_al,_hahal);\
			MULL32(_ah,_al,_lahal);\
			MULH32(_al,_al,_halal);\
			MULL32(_al,_al,_lalal);\
		\
			_sumlh =  _lahal << 1;\
			_sumhl = (_hahal << 1) | (_lahal >> 31);\
			_sumhh =  _hahah + (_hahal >> 31);\
		\
			/* adde/addze forces execution serialization! Other solution for carry propagation? */\
		  __asm__(\
			"addc  %0,%3,%4\n\t"\
			"adde  %1,%5,%6\n\t"\
			"addze %2,%7"\
			: "=r"(_sumlh), "=r"(_sumhl), "=r"(_sumhh)\
			: "%0"(_sumlh), "r"(_halal), "%1"(_sumhl), "r"(_lahah), "2"(_sumhh));\
		\
			* _aptr    = _sumhh;\
			*(_aptr+1) = _sumhl;\
			* _bptr    = _sumlh;\
			*(_bptr+1) = _lalal;\
			_hi = _a;\
			_lo = _b;\
		}

		#define MULH64(_x,_y,_hi)\
		{\
			uint32 _ah,_al,_bh,_bl;\
			uint32 _hahbh,_hahbl,_halbh,_halbl,_lahbh,_lahbl,_lalbh;\
			uint32 _sumhh,_sumhl,_sumlh;\
			uint64 _a = _x,_b = _y;\
			uint32 *_aptr = (uint32 *)&_a, *_bptr = (uint32 *)&_b;\
			/* PPC is big-endian, hence hi32 at aptr, lo32 at aptr+1: */\
			ah = * _aptr;		/* x_hi */\
			al = *(_aptr+1);	/* x_lo */\
			bh = * _bptr;		/* y_hi */\
			bl = *(_bptr+1);	/* y_lo */\
		\
			MULH32(_ah,_bh,_hahbh);\
			MULL32(_ah,_bh,_lahbh);\
			MULH32(_al,_bh,_halbh);\
			MULL32(_al,_bh,_lalbh);\
			MULH32(_ah,_bl,_hahbl);\
			MULL32(_ah,_bl,_lahbl);\
			MULH32(_al,_bl,_halbl);\
		\
			/* adde/addze forces execution serialization! Other solution for carry propagation? */\
		  __asm__(\
			"addc  %0,%3,%4\n\t"\
			"adde  %1,%5,%6\n\t"\
			"addze %2,%7"\
			: "=r"(_sumlh), "=r"(_sumhl), "=r"(_sumhh)\
			: "%r"(_halbl), "r"(_lalbh), "%r"(_halbh), "r"(_hahbl), "r"(_hahbh));\
		\
		  __asm__(\
			"addc  %0,%3,%4\n\t"\
			"adde  %1,%5,%6\n\t"\
			"addze %2,%7"\
			: "=r"(_sumlh), "=r"(_sumhl), "=r"(_sumhh)\
			: "%0"(_sumlh), "r"(_lahbl), "%1"(_sumhl), "r"(_lahbh), "2"(_sumhh));\
		\
			* _aptr    = _sumhh;\
			*(_aptr+1) = _sumhl;\
			_hi = a;	/* return value */\
		}

		#define	__MULL64(_x, _y)	({ uint64 _lo; MULL64((_x), (_y), _lo); _lo; })
		#define	__MULH64(_x, _y)	({ uint64 _lo; MULH64((_x), (_y), _lo); _lo; })

		#define MOD_INI_Q4(\
		 q0,qinv0\
		,q1,qinv1\
		,q2,qinv2\
		,q3,qinv3\
		)\
		{\
			ql0  = (uint32) q0   ;\
			ql1  = (uint32) q1   ;\
			ql2  = (uint32) q2   ;\
			ql3  = (uint32) q3   ;\
			qil0 = (uint32) qinv0;\
			qil1 = (uint32) qinv1;\
			qil2 = (uint32) qinv2;\
			qil3 = (uint32) qinv3;\
			\
			DBG_ASSERT((ql0  >> 32) == 0,"MOD_INI_Q4: (ql0  >> 32) == 0");\
			DBG_ASSERT((ql1  >> 32) == 0,"MOD_INI_Q4: (ql1  >> 32) == 0");\
			DBG_ASSERT((ql2  >> 32) == 0,"MOD_INI_Q4: (ql2  >> 32) == 0");\
			DBG_ASSERT((ql3  >> 32) == 0,"MOD_INI_Q4: (ql3  >> 32) == 0");\
			DBG_ASSERT((qil0 >> 32) == 0,"MOD_INI_Q4: (qil0 >> 32) == 0");\
			DBG_ASSERT((qil1 >> 32) == 0,"MOD_INI_Q4: (qil1 >> 32) == 0");\
			DBG_ASSERT((qil2 >> 32) == 0,"MOD_INI_Q4: (qil2 >> 32) == 0");\
			DBG_ASSERT((qil3 >> 32) == 0,"MOD_INI_Q4: (qil3 >> 32) == 0");\
			\
			qh0  = (uint32)(q0    >> 32);\
			qh1  = (uint32)(q1    >> 32);\
			qh2  = (uint32)(q2    >> 32);\
			qh3  = (uint32)(q3    >> 32);\
			qih0 = (uint32)(qinv0 >> 32);\
			qih1 = (uint32)(qinv1 >> 32);\
			qih2 = (uint32)(qinv2 >> 32);\
			qih3 = (uint32)(qinv3 >> 32);\
			\
			DBG_ASSERT((qh0  >> 32) == 0,"MOD_INI_Q4: (qh0  >> 32) == 0");\
			DBG_ASSERT((qh1  >> 32) == 0,"MOD_INI_Q4: (qh1  >> 32) == 0");\
			DBG_ASSERT((qh2  >> 32) == 0,"MOD_INI_Q4: (qh2  >> 32) == 0");\
			DBG_ASSERT((qh3  >> 32) == 0,"MOD_INI_Q4: (qh3  >> 32) == 0");\
			DBG_ASSERT((qih0 >> 32) == 0,"MOD_INI_Q4: (qih0 >> 32) == 0");\
			DBG_ASSERT((qih1 >> 32) == 0,"MOD_INI_Q4: (qih1 >> 32) == 0");\
			DBG_ASSERT((qih2 >> 32) == 0,"MOD_INI_Q4: (qih2 >> 32) == 0");\
			DBG_ASSERT((qih3 >> 32) == 0,"MOD_INI_Q4: (qih3 >> 32) == 0");\
		}

		/* For each input xj, calculates the following sequence:

			SQR_LOHI64(xj,loj,hij);
			loj = MULL64(loj,qinvj);
			yj  = MULH64(loj,qj);
		*/
		#define MOD_SQR_Q4(\
		 x0,hi0,y0\
		,x1,hi1,y1\
		,x2,hi2,y2\
		,x3,hi3,y3\
		)\
		{\
			uint32 ah0,al0,bh0,bl0,ah1,al1,bh1,bl1,ah2,al2,bh2,bl2,ah3,al3,bh3,bl3;\
			uint32 midh,midl,mjdh,mjdl,ss,tt;	/* Use s0-3, t0-3 here and below to avoid artificial execution serialization? */\
			uint64 a,b,c,d;	/* DEBUG */\
		\
			/* SQR_LOHI64(xj,loj,hij), loj in (alj,ahj), hij in (blj,bhj) : */\
			al0 = (uint32) x0;\
			al1 = (uint32) x1;\
			al2 = (uint32) x2;\
			al3 = (uint32) x3;\
		\
			ah0 = (uint32)(x0 >> 32);\
			ah1 = (uint32)(x1 >> 32);\
			ah2 = (uint32)(x2 >> 32);\
			ah3 = (uint32)(x3 >> 32);\
		\
			MULH32(ah0,ah0,bh0); bl0 = ah0*ah0; MULH32(ah0,al0,midh); midl = ah0*al0; MULH32(al0,al0,ah0); al0 = al0*al0; ss = (midl << 1); ah0 += ss; tt = (midh << 1) | (midl >> 31); bl0 += (ah0 < ss) + tt; bh0 += (bl0 < tt) + (midh >> 31);\
			MULH32(ah1,ah1,bh1); bl1 = ah1*ah1; MULH32(ah1,al1,midh); midl = ah1*al1; MULH32(al1,al1,ah1); al1 = al1*al1; ss = (midl << 1); ah1 += ss; tt = (midh << 1) | (midl >> 31); bl1 += (ah1 < ss) + tt; bh1 += (bl1 < tt) + (midh >> 31);\
			MULH32(ah2,ah2,bh2); bl2 = ah2*ah2; MULH32(ah2,al2,midh); midl = ah2*al2; MULH32(al2,al2,ah2); al2 = al2*al2; ss = (midl << 1); ah2 += ss; tt = (midh << 1) | (midl >> 31); bl2 += (ah2 < ss) + tt; bh2 += (bl2 < tt) + (midh >> 31);\
			MULH32(ah3,ah3,bh3); bl3 = ah3*ah3; MULH32(ah3,al3,midh); midl = ah3*al3; MULH32(al3,al3,ah3); al3 = al3*al3; ss = (midl << 1); ah3 += ss; tt = (midh << 1) | (midl >> 31); bl3 += (ah3 < ss) + tt; bh3 += (bl3 < tt) + (midh >> 31);\
		\
			/* Use | rather than + here (and below) to prevent compiler from using a 32-bit add-with-carry: */\
			hi0 = (uint64)bl0 | ((uint64)bh0 << 32);\
			hi1 = (uint64)bl1 | ((uint64)bh1 << 32);\
			hi2 = (uint64)bl2 | ((uint64)bh2 << 32);\
			hi3 = (uint64)bl3 | ((uint64)bh3 << 32);\
		/* DEBUG:\
		lo0 = (uint64)al0 + ((uint64)ah0 << 32);	SQR_LOHI64(x0,&a,&b);	if(a != lo0) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,a,lo0);	if(b != hi0) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,b,hi0);	\
		lo1 = (uint64)al1 + ((uint64)ah1 << 32);	SQR_LOHI64(x1,&a,&b);	if(a != lo1) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,a,lo1);	if(b != hi1) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,b,hi1);	\
		lo2 = (uint64)al2 + ((uint64)ah2 << 32);	SQR_LOHI64(x2,&a,&b);	if(a != lo2) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,a,lo2);	if(b != hi2) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,b,hi2);	\
		lo3 = (uint64)al3 + ((uint64)ah3 << 32);	SQR_LOHI64(x3,&a,&b);	if(a != lo3) printf("x,a,lo = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,a,lo3);	if(b != hi3) printf("x,b,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,b,hi3);	\
		*/\
			/* loj = MULL64(loj,qinvj) : */\
		/* DEBUG:\
		a = ((uint64)al0 + ((uint64)ah0 << 32))*qinv0;	\
		b = ((uint64)al1 + ((uint64)ah1 << 32))*qinv1;	\
		c = ((uint64)al2 + ((uint64)ah2 << 32))*qinv2;	\
		d = ((uint64)al3 + ((uint64)ah3 << 32))*qinv3;	\
		*/\
			ss = ah0*qil0; tt = al0*qih0; MULH32(al0,qil0,ah0); al0 = al0*qil0; ah0 += ss + tt;\
			ss = ah1*qil1; tt = al1*qih1; MULH32(al1,qil1,ah1); al1 = al1*qil1; ah1 += ss + tt;\
			ss = ah2*qil2; tt = al2*qih2; MULH32(al2,qil2,ah2); al2 = al2*qil2; ah2 += ss + tt;\
			ss = ah3*qil3; tt = al3*qih3; MULH32(al3,qil3,ah3); al3 = al3*qil3; ah3 += ss + tt;\
		\
			/* We really only need to save loj for the 65-bit case: */\
			lo0 = (uint64)al0 | ((uint64)ah0 << 32);	/*if(a != lo0) printf("o0\n");*/\
			lo1 = (uint64)al1 | ((uint64)ah1 << 32);	/*if(b != lo1) printf("o1\n");*/\
			lo2 = (uint64)al2 | ((uint64)ah2 << 32);	/*if(c != lo2) printf("o2\n");*/\
			lo3 = (uint64)al3 | ((uint64)ah3 << 32);	/*if(d != lo3) printf("o3\n");*/\
		\
			/* yj = MULH64(loj,qj) is here, store outputs in 64-bit calling argument : */\
			MULH32(ah0,qh0,bh0); bl0 = ah0*qh0; MULH32(ah0,ql0,midh); midl = ah0*ql0; MULH32(qh0,al0,mjdh); mjdl = qh0*al0; MULH32(ah0,al0,ql0); ss = midl + mjdl; ah0 += ss; tt = (ss < midl) + midh + mjdh; bl0 += (ah0 < ss) + tt; bh0 += (tt < midh) + (bl0 < tt);\
			MULH32(ah1,qh1,bh1); bl1 = ah1*qh1; MULH32(ah1,ql1,midh); midl = ah1*ql1; MULH32(qh1,al1,mjdh); mjdl = qh1*al1; MULH32(ah1,al1,ql1); ss = midl + mjdl; ah1 += ss; tt = (ss < midl) + midh + mjdh; bl1 += (ah1 < ss) + tt; bh1 += (tt < midh) + (bl1 < tt);\
			MULH32(ah2,qh2,bh2); bl2 = ah2*qh2; MULH32(ah2,ql2,midh); midl = ah2*ql2; MULH32(qh2,al2,mjdh); mjdl = qh2*al2; MULH32(ah2,al2,ql2); ss = midl + mjdl; ah2 += ss; tt = (ss < midl) + midh + mjdh; bl2 += (ah2 < ss) + tt; bh2 += (tt < midh) + (bl2 < tt);\
			MULH32(ah3,qh3,bh3); bl3 = ah3*qh3; MULH32(ah3,ql3,midh); midl = ah3*ql3; MULH32(qh3,al3,mjdh); mjdl = qh3*al3; MULH32(ah3,al3,ql3); ss = midl + mjdl; ah3 += ss; tt = (ss < midl) + midh + mjdh; bl3 += (ah3 < ss) + tt; bh3 += (tt < midh) + (bl3 < tt);\
		\
			y0 = (uint64)bl0 | ((uint64)bh0 << 32);\
			y1 = (uint64)bl1 | ((uint64)bh1 << 32);\
			y2 = (uint64)bl2 | ((uint64)bh2 << 32);\
			y3 = (uint64)bl3 | ((uint64)bh3 << 32);\
		/* DEBUG:\
			a = MULH64(lo0,q0);	if(a != y0) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x0,a,lo0);	\
			a = MULH64(lo1,q1);	if(a != y1) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x1,a,lo1);	\
			a = MULH64(lo2,q2);	if(a != y2) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x2,a,lo2);	\
			a = MULH64(lo3,q3);	if(a != y3) printf("lo,q,hi = %20" PRIu64 " %20" PRIu64 " %20" PRIu64 "\n",x3,a,lo3);	\
		*/\
		}

	#else
		#error Unrecognized or undefined CPU_SUBTYPE value for PPC!
    #endif

/* Generic macro form: */

#elif !defined(MUL_LOHI64_SUBROUTINE)	/* User compile-time override flag to force subroutine form rather than macro: */

  #ifdef COMPILER_TYPE_GCC
   #ifdef VERBOSE_HEADERS
	#warning Using generic C/32-bit version of MUL64 macros
   #endif
	#define __MULL32(x32,y32     )                   ((uint32)(x32)*(uint32)(y32))
	#define __MULH32(x32,y32     )                  (((uint32)(x32)*(uint64)(y32)) >> 32)
	#define MULL32(  x32,y32,lo32)	lo32 = (uint32) ((uint32)(x32)*(uint32)(y32))
	#define MULH32(  x32,y32,hi32)	hi32 = __MULH32((uint32)(x32),(uint64)(y32))

	#define	MULL64(_x,_y,_lo)\
	{\
			uint32 xlo = (uint32)(_x);			\
			uint32 ylo = (uint32)(_y);			\
			uint32 xhi = (uint32)(_x >> 32);	\
			uint32 yhi = (uint32)(_y >> 32);	\
			\
			_lo = (uint64)__MULL32( xlo,ylo)	\
				+(((uint64)__MULH32(xlo,ylo)) << 32)	\
				+(((uint64)__MULL32( xlo,yhi) + (uint64)__MULL32( xhi,ylo)) << 32);	\
	}

	/* 64x32=>96-bit product algorithm:
	represent the inputs as x = a + b*2^32, y = c ( < 2^32), then do 4
	32-bit MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .

	Even though the high output for 64x32-bit is always < 2^32,
	assume _y and _hi here are 64-bit ints to allow flexibility for caller.
	*/
   #if EWM_DEBUG
	#define  MUL64x32(_x, _y,_lo,_hi)\
	{\
		char s0[21],s1[21];\
		uint64 _t,_a,_b;\
		\
		ASSERT(((uint64)(_y) >> 32) == 0,"MUL64x32: ((_y) >> 32) == 0");\
		MUL_LOHI64((_x), (uint64)(_y), _a, _b);\
		\
		_lo = ((uint32)((_x) & 0x00000000ffffffff)) * (_y);	/* a*c */\
		_t  = ((uint32)((_x) >> 32)) * (_y);				/* b*c */\
		_hi = (_t >> 32);\
		_t <<= 32;\
		_lo +=  _t;\
		_hi += (_lo < _t);\
		\
		if(_a != _lo || _b != _hi)\
		{\
			printf("x = %s, y = %s\n", s0[convert_uint64_base10_char(s0,_x )], s1[convert_uint64_base10_char(s1,_y)]);\
			printf("LO= %s, A = %s\n", s0[convert_uint64_base10_char(s0,_lo)], s1[convert_uint64_base10_char(s1,_a)]);\
			printf("HI= %s, B = %s\n", s0[convert_uint64_base10_char(s0,_hi)], s1[convert_uint64_base10_char(s0,_b)]);\
			ASSERT(0,"0");\
		}\
	}
   #else
	#define  MUL64x32(_x, _y,_lo,_hi)\
	{\
		uint64 _t;\
		\
		_lo = ((uint32)((_x) & 0x00000000ffffffff)) * (_y);	/* a*c */\
		_t  = ((uint32)((_x) >> 32)) * (_y);					/* b*c */\
		_hi = (_t >> 32);\
		_t <<= 32;\
		_lo +=  _t;\
		_hi += (_lo < _t);\
	}
   #endif

	/* Generic 128-bit product macros: represent the inputs as
	x = a + b*2^32, y = c + d*2^32, and then do 4 MULs and a bunch of
	adds-with-carry to get x*y = b*d*2^64 + (a*d + b*c)*2^32 + a*c .
	Actual calling arguments are assumed to be 64-bit ints - user must
	make sure this is true. Result written into lo, hi.
	*/
	#define MUL_LOHI64(_x,_y,_lo,_hi)\
	{\
		uint64 _a,_b,_c,_d,_ac,_ad,_bc;				\
		/*_a = (_x) & (uint64)0x00000000ffffffff;*/		\
		_a = ((_x)<<32)>>32;								\
		_b =  (_x)>>32;									\
		/*_c = (_y) & (uint64)0x00000000ffffffff;*/		\
		_c = ((_y)<<32)>>32;								\
		_d =  (_y)>>32;									\
		/* Calculate 4 subproducts in order in which they are first used */\
		_ac = _a*_c;										\
		_bc = _b*_c;										\
		_hi = _b*_d;										\
		_ad = _a*_d;										\
		_lo  =  _ac;		/* use _lo to store copy of _ac */\
		_ac +=         (_bc<<32);	_hi += (_ac < _lo);	\
		_lo  =  _ac + (_ad<<32);	_hi += (_lo < _ac);	\
						_hi += (_bc>>32) + (_ad>>32);	\
	}

	/* Generic 128-bit squaring algorithm: represent the input as
	x = a + b*2^32, and then do 3 MULs and a bunch of add-with-carries
	to get x^2 = b^2*2^64 + a*b*2^33 + a^2 . Actual calling arguments
	are assumed to be 64-bit ints - user must make sure this is true.
	Result written into lo, hi.
	*/
	#define SQR_LOHI64(_x,_lo,_hi)\
	{\
		uint64 _a,_b,_aa,_ab;			\
		/*_a = (_x) & (uint64)0x00000000ffffffff;*/\
		_a = ((_x)<<32)>>32;				\
		_b =  (_x)>>32;					\
		_aa = _a*_a;						\
		_ab = _a*_b;						\
		_hi  = _b*_b;					\
		_lo  = _aa + (_ab<<33);			\
		_hi += (_ab>>31) + (_lo < _aa);	\
	}

	/* Generic 128-bit multiply algorithm, returning only the upper 64 bits (high part) of the result.
	Actual calling arguments are assumed to be 64-bit ints - user must make sure this is true.
	Result written into hi.
	*/
	#define MULH64(_x, _y,_hi)\
	{\
		uint64 _tt;					\
		MUL_LOHI64((_x), (_y), _tt, _hi);	\
	}

  #elif defined(COMPILER_TYPE_NVCC)

	#define MULL32(	 _x,_y, _lo     )	_lo =          (_x)*(_y)
	#define MULH32(  _x,_y,      _hi)	_hi = __umulhi((_x),(_y))
	#define __MULL32(_x, _y)  ({ uint32 _lo =          (_x)*(_y) ; _lo; })
	#define __MULH32(_x, _y)  ({ uint32 _hi = __umulhi((_x),(_y)); _hi; })
	#define MUL_LOHI32(_x,_y,_lo,_hi)\
	{\
		_lo =          (_x)*(_y) ;\
		_hi = __umulhi((_x),(_y));\
	}

   #if 0	// My hand-rolled 64-bit-int emulations below are faster than NVCC's emulations-via-intrinsics:

	#warning Using CUDA intrinsic MUL64 macros
	#define MUL_LOHI64(_x,_y,_lo,_hi){uint64 _t = (_x)*(_y); _hi = __umul64hi((_x), (_y));	_lo = _t;}
	#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)
	#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
	#define	__MULL64(	 _x,_y      )	           (_x)*(_y)
	#define	__MULH64(	 _x,_y      )	__umul64hi((_x),(_y))
	#define MULL64(	 _x,_y, _lo     )	_lo =            (_x)*(_y)
	#define MULH64(  _x,_y,      _hi)	_hi = __umul64hi((_x),(_y))

   #else

	// These single-PTX-instructions inline-ASMs are taken staright from Oliver Weihe's MFAKTC my_intrinsics.h file:
	// R = A + B: carryin =  no; carryout = yes:
	__device__ static unsigned int __add_cc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("add.cc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}
	// R = A + B: carryin = yes; carryout = yes:
	__device__ static unsigned int __addc_cc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("addc.cc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}
	// R = A + B: carryin = yes; carryout =  no:
	__device__ static unsigned int __addc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("addc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}

	// R = A - B: carryin =  no; carryout = yes:
	__device__ static unsigned int __sub_cc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("sub.cc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}
	// R = A - B: carryin = yes; carryout = yes:
	__device__ static unsigned int __subc_cc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("subc.cc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}
	// R = A - B: carryin = yes; carryout =  no:
	__device__ static unsigned int __subc(unsigned int a, unsigned int b)
	{
	  unsigned int r;
	  asm("subc.u32 %0, %1, %2;" : "=r" (r) : "r" (a) , "r" (b));
	  return r;
	}

	// R = (A * B).lo + C: carryin =  no; carryout =  no:
	__device__ static unsigned int __umad32(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("mad.lo.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).lo + C: carryin =  no; carryout = yes:
	__device__ static unsigned int __umad32_cc(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).lo + C: carryin = yes; carryout =  no:
	__device__ static unsigned int __umad32c(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("madc.lo.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).lo + C: carryin = yes; carryout = yes:
	__device__ static unsigned int __umad32c_cc(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("madc.lo.cc.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}

	// R = (A * B).hi + C: carryin =  no; carryout =  no:
	__device__ static unsigned int __umad32hi(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("mad.hi.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).hi + C: carryin =  no; carryout = yes:
	__device__ static unsigned int __umad32hi_cc(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("mad.hi.cc.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).hi + C: carryin = yes; carryout =  no:
	__device__ static unsigned int __umad32hic(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("madc.hi.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}
	// R = (A * B).hi + C: carryin = yes; carryout = yes:
	__device__ static unsigned int __umad32hic_cc(unsigned int a, unsigned int b, unsigned int c)
	{
	  unsigned int r;
	  asm("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (r) : "r" (a) , "r" (b), "r" (c));
	  return r;
	}

	/* Low 64 bits of product of uint64 inputs _x and _y returned in uint64 _lo */
	#define	MULL64(_x,_y,_lo)\
	{\
		union {\
			uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
			uint64 u64;\
		} _ux,_uy;\
	\
		uint32 _a,_b,_c,_d;\
		uint32 _lo32,_hi32;\
	\
		_ux.u64 = _x;\
		_a = _ux.u32[0];	/* x_lo */\
		_b = _ux.u32[1];	/* x_hi */\
		_uy.u64 = _y;\
		_c = _uy.u32[0];	/* y_lo */\
		_d = _uy.u32[1];	/* y_hi */\
		MUL_LOHI32(_a,_c,_lo32,_hi32);\
		_hi32 += _a*_d + _b*_c;	\
	\
		_ux.u32[0] = _lo32;\
		_ux.u32[1] = _hi32;\
		_lo= _ux.u64;\
	}

	/* 64x32=>96-bit product algorithm:
	represent the inputs as x = a + b*2^32, y = c ( < 2^32), then do 4
	32-bit MULs and a bunch of add-with-carries to get x*y = b*c*2^32 + a*c .

	On the IA32, we do 32x64-bit via 2 32x32-bit MULs (2 mull).

	Even though the high output for 64x32-bit is always < 2^32,
	assume _y and _hi here are 64-bit ints to allow flexibility for caller.
	*/
	#define MUL64x32(_x, _y,_lo,_hi)\
	{\
		union {\
			uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
			uint64 u64;\
		} _ux,_uy;\
	\
		uint32 _a,_b,_c;\
		uint32 _lo32,_md32,_hi32;\
		uint32 _bclo,cy;\
	\
		_ux.u64 = _x;\
		_a = _ux.u32[0];	/* x_lo */\
		_b = _ux.u32[1];	/* x_hi */\
		_uy.u64 = _y;\
		_c = _uy.u32[0];	/* y_lo */\
	\
		MUL_LOHI32(_a,_c,_lo32,_md32);	\
		MUL_LOHI32(_b,_c,_bclo,_hi32);	\
	\
		_md32 += _bclo;	\
		cy = _md32 < _bclo;	/* (0) md32 + (b*c)_lo                     , result (middle 32 bits) in md32, carryout in CY. */\
		_hi32 += cy;		/* (1)    0 + (b*c)_hi + (carryin from (0)), result ( upper 32 bits) in hi32, carryout in CY - should be zero! [In fact assume so.] */\
	\
		_ux.u32[0] = _lo32;\
		_ux.u32[1] = _md32;\
		_lo= _ux.u64;\
		_uy.u32[0] = _hi32;\
		_uy.u32[1] = 0;\
		_hi= _uy.u64;\
	}

	#if 1	// PTX-code version of MUL_LOHI64:

		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _al,_ah,_bl,_bh;\
			uint32 _r0,_r1,_r2,_r3;\
		\
			_ux.u64 = _x;\
			_al = _ux.u32[0];	/* x_lo */\
			_ah = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_bl = _uy.u32[0];	/* y_lo */\
			_bh = _uy.u32[1];	/* y_hi */\
		/* Old 1-line-at-a-time inline-ASM version: */\
		/*	MULL32(_al,_bl,_r0);				// r0 =(_al*_bl).lo   , no carry-out */\
		/*	MULH32(_al,_bl,_r1);				// r1 =(_al*_bl).hi   , no carry-out */\
		/*	_r1 = __umad32_cc(_ah,_bl,_r1);		// r1+=(_ah*_bl).lo   , may carry-out */\
		/*	_r2 = __umad32hic(_ah,_bl,  0);		// r2 =(_ah*_bl).hi+cy, no carry-out */\
		/*	_r1 = __umad32_cc(_al,_bh,_r1);		// r1+=(_al*_bh).lo   , may carry-out */\
		/*	_r2 = __umad32hic_cc(_al,_bh,_r2);	// r2+=(_al*_bh).hi+cy, may carry-out */\
		/*	_r3 = __addc(  0,  0);				// r3 = cy, no carry-out */\
		/*	_r2 = __umad32_cc(_ah,_bh,_r2);		// r2+=(_ah*_bh).lo   , may carry-out */\
		/*	_r3 = __umad32hic(_ah,_bh,_r3);		// r3+=(_ah*_bh).hi+cy */\
		/* Using straight PTX inline-ASM: */\
		asm volatile (\
			"mul.lo.u32		%0,%4,%6;		\n\t"/* r0 =(r4*r6).lo   , no carry-out */\
			"mul.hi.u32		%1,%4,%6;		\n\t"/* r1 =(r4*r6).hi   , no carry-out */\
			"mad.lo.cc.u32	%1,%5,%6,%1;	\n\t"/* r1+=(r5*r6).lo   , may carry-out */\
			"madc.hi.u32	%2,%5,%6, 0;	\n\t"/* r2 =(r5*r6).hi+cy, no carry-out */\
			"mad.lo.cc.u32	%1,%4,%7,%1;	\n\t"/* r1+=(r4*r7).lo   , may carry-out */\
			"madc.hi.cc.u32	%2,%4,%7,%2;	\n\t"/* r2+=(r4*r7).hi+cy, may carry-out */\
			"addc.u32		%3, 0, 0;		\n\t"/* r3 = cy, no carry-out */\
			"mad.lo.cc.u32	%2,%5,%7,%2;	\n\t"/* r2+=(r5*r7).lo   , may carry-out */\
			"madc.hi.u32	%3,%5,%7,%3;	\n\t"/* r3+=(r5*r7).hi+cy */\
		  : "=r" (_r0) , "=r" (_r1) "=r" (_r2) , "=r" (_r3) : "r" (_al), "r" (_ah) , "r" (_bl), "r" (_bh));\
		\
			_ux.u32[0] = _r0;\
			_ux.u32[1] = _r1;\
			_lo= _ux.u64;\
			_uy.u32[0] = _r2;\
			_uy.u32[1] = _r3;\
			_hi= _uy.u64;\
		}

		#define SQR_LOHI64(_x,   _lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _al,_ah,_bl,_bh;\
			uint32 _r0,_r1,_r2,_r3;\
		\
			_ux.u64 = _x;\
			_al = _ux.u32[0];	/* x_lo */\
			_ah = _ux.u32[1];	/* x_hi */\
		\
		asm volatile (\
			"mul.lo.u32		%0,%4,%4;		\n\t"/* r0 =(r4*r4).lo   ,  no carry-out */\
			"mul.hi.u32		%1,%4,%4;		\n\t"/* r1 =(r4*r4).hi   ,  no carry-out */\
			"mad.lo.cc.u32	%1,%4,%5,%1;	\n\t"/* r1+=(r4*r5).lo   , may carry-out */\
			"madc.hi.u32	%2,%4,%5, 0;	\n\t"/* r2 =(r4*r5).hi+cy,  no carry-out */\
			"mad.lo.cc.u32	%1,%4,%5,%1;	\n\t"/* r1+=(r4*r5).lo   , may carry-out */\
			"madc.hi.cc.u32	%2,%4,%5,%2;	\n\t"/* r2+=(r4*r5).hi+cy, may carry-out */\
			"addc.u32		%3, 0, 0;		\n\t"/* r3 = cy, no carry-out */\
			"mad.lo.cc.u32	%2,%5,%5,%2;	\n\t"/* r2+=(r5*r5).lo   , may carry-out */\
			"madc.hi.u32	%3,%5,%5,%3;	\n\t"/* r3+=(r5*r5).hi+cy */\
		  : "=r" (_r0) , "=r" (_r1) "=r" (_r2) , "=r" (_r3) : "r" (_al), "r" (_ah));\
		\
			_ux.u32[0] = _r0;\
			_ux.u32[1] = _r1;\
			_lo= _ux.u64;\
			_uy.u32[0] = _r2;\
			_uy.u32[1] = _r3;\
			_hi= _uy.u64;\
		}

		#define MULH64_TEST(_x, _y,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _al,_ah,_bl,_bh;\
			uint32 _r0,_r1,_r2,_r3;\
		\
			_ux.u64 = _x;\
			_al = _ux.u32[0];	/* x_lo */\
			_ah = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_bl = _uy.u32[0];	/* y_lo */\
			_bh = _uy.u32[1];	/* y_hi */\
		\
			MULL32(_al,_bl,_r0);				/*printf("[1]: %10u\n",_r0);*/	/* r0 =(_al*_bl).lo   ,  no carry-out */\
			MULH32(_al,_bl,_r1);				/*printf("[2]: %10u\n",_r1);*/	/* r1 =(_al*_bl).hi   ,  no carry-out */\
			_r1 = __umad32_cc(_ah,_bl,_r1);		/*printf("[3]: %10u\n",_r1);*/	/* r1+=(_ah*_bl).lo   , may carry-out */\
			_r2 = __umad32hic(_ah,_bl,  0);		/*printf("[4]: %10u\n",_r2);*/	/* r2 =(_ah*_bl).hi+cy,  no carry-out */\
			_r1 = __umad32_cc(_al,_bh,_r1);		/*printf("[5]: %10u\n",_r1);*/	/* r1+=(_al*_bh).lo   , may carry-out */\
			_r2 = __umad32hic_cc(_al,_bh,_r2);	/*printf("[6]: %10u\n",_r2);*/	/* r2+=(_al*_bh).hi+cy, may carry-out */\
			_r3 = __addc(  0,  0);				/*printf("[7]: %10u\n",_r3);*/	/* r3 = cy, no carry-out */\
			_r2 = __umad32_cc(_ah,_bh,_r2);		/*printf("[8]: %10u\n",_r2);*/	/* r2+=(_ah*_bh).lo   , may carry-out */\
			_r3 = __umad32hic(_ah,_bh,_r3);		/*printf("[9]: %10u\n",_r3);*/	/* r3+=(_ah*_bh).hi+cy */\
		\
			_uy.u32[0] = _r2;\
			_uy.u32[1] = _r3;\
			_hi= _uy.u64;\
		}
		// Above test version of the MULH64 macro reveals that CY flag preservation across successive 1-line-inline-ASM calls
		// is not guaranteed - in my tests, the crucial carry (in the sense that its loss was hosing my Mfactor self-tests)
		// out of [5] was getting lost as written, but suddenly preserved when I added the (now-commented-out) printf statements.
		// Having the printing on not being a practical option, we now do things the right way, namely via a single 9-instruction
		// block of inline-ASM for the crucial MAC sequence. (Subsequently did similarly with the other MUL macros):
		#define MULH64(_x, _y,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 _al,_ah,_bl,_bh;\
			uint32 _r2,_r3;\
		\
			_ux.u64 = _x;\
			_al = _ux.u32[0];	/* x_lo */\
			_ah = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			_bl = _uy.u32[0];	/* y_lo */\
			_bh = _uy.u32[1];	/* y_hi */\
		\
		/* Don't use %0 value but need for output-enum, so stick value into _al to get rid of 'unused	\
		variable' compiler warning resulting from use of dedicated var 'r1' to hold %o output value: */\
		asm volatile (\
			"mul.hi.u32		%0,%3,%5;		\n\t"/* r1 =(r4*r6).hi   ,  no carry-out */\
			"mad.lo.cc.u32	%0,%4,%5,%0;	\n\t"/* r1+=(r5*r6).lo   , may carry-out */\
			"madc.hi.u32	%1,%4,%5, 0;	\n\t"/* r2 =(r5*r6).hi+cy,  no carry-out */\
			"mad.lo.cc.u32	%0,%3,%6,%0;	\n\t"/* r1+=(r4*r7).lo   , may carry-out */\
			"madc.hi.cc.u32	%1,%3,%6,%1;	\n\t"/* r2+=(r4*r7).hi+cy, may carry-out */\
			"addc.u32		%2, 0, 0;		\n\t"/* r3 = cy, no carry-out */\
			"mad.lo.cc.u32	%1,%4,%6,%1;	\n\t"/* r2+=(r5*r7).lo   , may carry-out */\
			"madc.hi.u32	%2,%4,%6,%2;	\n\t"/* r3+=(r5*r7).hi+cy */\
		  : "=r" (_al) "=r" (_r2) , "=r" (_r3) : "r" (_al), "r" (_ah) , "r" (_bl), "r" (_bh));\
		\
			_uy.u32[0] = _r2;\
			_uy.u32[1] = _r3;\
			_hi = _uy.u64;\
		}

	#else	// Generic (non-PTX-code) version of above:

		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_,bh_,bl_;\
			uint32 hahbh_,hahbl_,halbh_,halbl_,lahbh_,lahbl_,lalbh_,lalbl_;\
			uint32 sumhh_,sumhl_,sumlh_,cy;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			bl_ = _uy.u32[0];	/* y_lo */\
			bh_ = _uy.u32[1];	/* y_hi */\
		\
		/* EWM: for 32-bit unsigned inputs these have a nominal latency of ?? cycles. */\
		/* For each MUL output we include the order in which it is used in the add/carry */\
		/* code below (6 steps, labeled (0)-(5)), in order to properly schedule the MUL: */\
			MUL_LOHI32(ah_,bh_,lahbh_,hahbh_);	/* (x_hi*y_hi)l,h 4,2 */\
			MUL_LOHI32(al_,bh_,lalbh_,halbh_);	/* (x_lo*y_hi)l,h 0,1 */\
			MUL_LOHI32(ah_,bl_,lahbl_,hahbl_);	/* (x_hi*y_lo)l,h 3,1 */\
			MUL_LOHI32(al_,bl_,lalbl_,halbl_);	/* (x_lo*y_lo)l,h -,0 */\
		/* Bits   0- 31 : lalbl                           */\
		/* Bits  32- 63 :*halbl + lahbl +*lalbh   (sumlh) */\
		/* Bits  64- 95 :*hahbl +*halbh + lahbh   (sumhl) */\
		/* Bits  96-127 : hahbh                   (sumhh) */\
		\
			sumlh_  = halbl_ + lalbh_;	/* (0) halbl + lalbh                       , result (partial sum) in sumlh,carryout in CY. */\
			cy  = sumlh_ < halbl_;\
			halbh_ += cy;\
			cy  = halbh_ < cy;\
			sumhl_  = hahbl_ + halbh_;	/* (1) halbh + hahbl + (carryin from sumlh), result (partial sum) in sumhl,carryout in CY. */\
			cy += sumhl_ < hahbl_;\
			sumhh_	= hahbh_ + cy;		/* (2)     0 + hahbh + (carryin from sumhl), result (partial sum) in sumhh,carryout in CF) bit - should be zero! */\
		\
			sumlh_ += lahbl_;\
			cy  = sumlh_ < lahbl_;		/* (3) sumlh + lahbl                       , result (sumlh), carryout in CY. */\
			lahbh_ += cy;\
			cy  = lahbh_ < cy;\
			sumhl_ += lahbh_;\
			cy += sumhl_ < lahbh_;		/* (4) sumhl + lahbh + (carryin from sumlh), result (sumhl), carryout in CY. */\
			sumhh_ += cy;				/* (5)     0 + sumhh + (carryin from sumhl), result (sumhh), carryout in CY - should be zero! [In fact assume so.] */\
		\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi= _ux.u64;\
			_uy.u32[0] = lalbl_;\
			_uy.u32[1] = sumlh_;\
			_lo= _uy.u64;\
		}

		// Compute square of uint64 x via 32-bit x = x.lo + x.hi<<32 decomposition and 32-bit arithmetic:
		#define SQR_LOHI64(_x,_lo,_hi)\
		{\
			union {\
			uint32 u32[2];\
			uint64 u64;\
			} a,b;\
		\
			uint32 ah,al,cy;\
			uint32 hahah,lahah,hahal,lahal,halal,lalal;\
			uint32 sumhh,sumhl,sumlh;\
		\
			a.u64 = _x;\
			ah = a.u32[1];	/* NVCC is little-endian, LS 32 bits in [0]. */\
			al = a.u32[0];\
		/* Compute 64-bit subproducts in (2^32.ah + al)^2 = 2^64.[ah^2] + 2^32.2.[ah.al] + [al^2]: */\
			MULH32(ah,ah,hahah);	/* ah^2, hi half */\
			 MULL32(ah,ah,lahah);	/* ah^2, lo half */\
			MULH32(ah,al,hahal);	/* ah.al,hi half */\
			 MULL32(ah,al,lahal);	/* ah.al,lo half */\
			MULH32(al,al,halal);	/* al^2, hi half */\
			 MULL32(al,al,lalal);	/* al^2, lo half */\
		/* 2*<lahal,hahal>: */\
			sumlh = lahal + lahal;\
			sumhl = (hahal + hahal) | ((int)lahal < 0);	/* Shifts slow on most CUDA devices, so replace << 1 by add and >> 31 by signed-less-than-0 */\
			sumhh = hahah + ((int)hahal < 0);\
		\
			sumlh += halal;/* (0) sumlh + halal                     , result in sumlh, carryout in CY */\
			cy  = sumlh < halal;\
			lahah +=    cy;/* (1) sumhl + lahah + (carryin from (0)), result in sumhl, carryout in CY */\
			cy  = lahah < cy;\
			sumhl += lahah;\
			cy += sumhl < lahah;\
			sumhh += cy;		/* (2)     0 + sumhh + (carryin from (1)), result in sumhh, carryout in CY - should be zero! (In fact assume so.) */\
		\
			a.u32[1] = sumhh;\
			a.u32[0] = sumhl;\
			_hi = a.u64;\
			b.u32[1] = sumlh;\
			b.u32[0] = lalal;\
			_lo = b.u64;\
		}

		/* Generic 128-bit multiply algorithm, returning only the upper 64 bits (high part) of the result.
		Actual calling arguments are assumed to be 64-bit ints - user must make sure this is true.
		Result written into hi.
		EWM: *to-do* Currently this is just MUL_LOHI64, with the _lo-arg-related stuff cut out.
					Need to optimize, e.g. use optimizations that generate the carry into the high half
					without requiring exact computation of the entire lower-half partial products.
		*/
		#define MULH64(_x, _y,_hi)\
		{\
			union {\
				uint32 u32[2];	/* NVCC is little-endian, LS 32 bits in [0]. */\
				uint64 u64;\
			} _ux,_uy;\
		\
			uint32 ah_,al_,bh_,bl_;\
			uint32 hahbh_,hahbl_,halbh_,halbl_,lahbh_,lahbl_,lalbh_;\
			uint32 sumhh_,sumhl_,sumlh_,cy;\
		\
			_ux.u64 = _x;\
			al_ = _ux.u32[0];	/* x_lo */\
			ah_ = _ux.u32[1];	/* x_hi */\
			_uy.u64 = _y;\
			bl_ = _uy.u32[0];	/* y_lo */\
			bh_ = _uy.u32[1];	/* y_hi */\
		\
		/* EWM: for 32-bit unsigned inputs these have a nominal latency of ?? cycles. */\
		/* For each MUL output we include the order in which it is used in the add/carry */\
		/* code below (6 steps, labeled (0)-(5)), in order to properly schedule the MUL: */\
			MUL_LOHI32(ah_,bh_,lahbh_,hahbh_);	/* (x_hi*y_hi)l,h 4,2 */\
			MUL_LOHI32(al_,bh_,lalbh_,halbh_);	/* (x_lo*y_hi)l,h 0,1 */\
			MUL_LOHI32(ah_,bl_,lahbl_,hahbl_);	/* (x_hi*y_lo)l,h 3,1 */\
			MULH32    (al_,bl_,       halbl_);	/* (x_lo*y_lo)l,h -,0 */\
		/* Bits   0- 31 : lalbl                           */\
		/* Bits  32- 63 :*halbl + lahbl +*lalbh   (sumlh) */\
		/* Bits  64- 95 :*hahbl +*halbh + lahbh   (sumhl) */\
		/* Bits  96-127 : hahbh                   (sumhh) */\
		\
			sumlh_  = halbl_ + lalbh_;	/* (0) halbl + lalbh                       , result (partial sum) in sumlh,carryout in CY. */\
			cy  = sumlh_ < halbl_;\
			halbh_ += cy;\
			cy  = halbh_ < cy;\
			sumhl_  = hahbl_ + halbh_;	/* (1) halbh + hahbl + (carryin from sumlh), result (partial sum) in sumhl,carryout in CY. */\
			cy += sumhl_ < hahbl_;\
			sumhh_	= hahbh_ + cy;		/* (2)     0 + hahbh + (carryin from sumhl), result (partial sum) in sumhh,carryout in CF) bit - should be zero! */\
		\
			sumlh_ += lahbl_;\
			cy  = sumlh_ < lahbl_;		/* (3) sumlh + lahbl                       , result (sumlh), carryout in CY. */\
			lahbh_ += cy;\
			cy  = lahbh_ < cy;\
			sumhl_ += lahbh_;\
			cy += sumhl_ < lahbh_;		/* (4) sumhl + lahbh + (carryin from sumlh), result (sumhl), carryout in CY. */\
			sumhh_ += cy;				/* (5)     0 + sumhh + (carryin from sumhl), result (sumhh), carryout in CY - should be zero! [In fact assume so.] */\
		\
			_ux.u32[0] = sumhl_;\
			_ux.u32[1] = sumhh_;\
			_hi= _ux.u64;\
		}

	#endif	// endif(PTX-code version of MUL_LOHI64)?

   #endif	// 0 [NVCC's emulations-via-intrinsics]] | 1 [hand-rolled 64-bit-int emulations] ?

  #endif	// GCC | NVCC ?

	#define __MULL64(_x,_y)	({ uint64 _lo; MULL64(_x,_y,_lo); _lo; })
	#define __MULH64(_x,_y)	({ uint64 _hi; MULH64(_x,_y,_hi); _hi; })

#endif


/* On systems where no fast 128-bit hardware integer multiply is available
   and 32x32 => 64-bit MUL is slow, it may pay to use a floating emulation of MULH64.
   This requires the smaller of the two multiplicands (we assume this is passed
   in x) to be no larger than 52 bits in size  - y may be as large as 64 bits.
*/
#if USE_FLOATING_MULH64

	static const double scale = 1.0/(131072.0 * 131072.0 * 131072.0);

	#define MULT_LOHI_52x64(x, y, lo, hi)								\
	{																	\
		/* We expect the float <==> int conversions here will be slow, so	*/\
		/* if we're doing a lot of these a multi-operand version should		*/\
		/* alleviate the type-conversion latencies:							*/\
		const double dhi = (double)(x) * scale * (double)((y) >> 11);	\
		const uint64 hi4est = (uint64)(dhi) - 1;						\
		lo = (x) * (y);													\
		hi = (hi4est + (3 & ((lo >> 62) - hi4est))) >> 2;				\
	}

#endif

#ifndef __MULL32
	#define __MULL32(x32,y32     )	 ((uint32)(x32)*(uint32)(y32))
#endif
#ifndef __MULH32
	#define __MULH32(x32,y32     )	(((uint32)(x32)*(uint64)(y32)) >> 32)
#endif
#ifndef MULL32
	#define MULL32(  x32,y32,lo32)	lo32 = (uint32)((uint32)(x32)*(uint32)(y32))
#endif
#ifndef MULH32
	#define MULH32(  x32,y32,hi32)	hi32 = __MULH32((uint32)(x32),(uint64)(y32))
#endif
#ifndef MUL_LOHI32
	#define MUL_LOHI32(_x32,_y32,_lo,_hi)\
	{\
		uint64 _x = (uint64)(_x32) & 0x00000000FFFFFFFFull, _y = (uint64)(_y32) & 0x00000000FFFFFFFFull;\
		uint64 _tt = _x * _y;\
		_lo = (uint32) _tt;\
		_hi = (uint32)(_tt >> 32);\
	}
#endif
#ifndef SQR_LOHI32
	#define SQR_LOHI32(_x,   _lo,_hi)	MUL_LOHI32(_x,_x,_lo,_hi)
#endif

#ifdef MUL_LOHI64_SUBROUTINE
	/* 64-bit stuff below defined as functions (not macros) in imul_macro.c */
#else
	#ifndef __MULL64
		#error __MULL64 Undefined!
	#endif
	#ifndef __MULH64
		#error __MULH64 Undefined!
	#endif
	#ifndef MULL64
		#error MULL64 Undefined!
	#endif
	#ifndef MULH64
		#error MULH64 Undefined!
	#endif
	#ifndef MUL64x32
		#error MUL64x32 Undefined!
	#endif
	#ifndef MUL_LOHI64
		#error MUL_LOHI64 Undefined!
	#endif
	#ifndef SQR_LOHI64
		#error SQR_LOHI64 Undefined!
  #endif
#endif

#ifdef __cplusplus
}
#endif

#endif	/* imul_macro0_h_included */

