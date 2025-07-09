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
#ifndef align_h_included
#define align_h_included

#include "types.h"

/* These are for basic memory allocation, and to force alignment of array data on desired-byte boundaries.
We use the normally-not-recommended immediate-overwrite-of-pointer form of realloc() because if the returned
pointer is null we exit immediately, thus the resulting memory leak is never an issue.

In the Align macros we cast pointers to longto accommodate architectures which use 64-bit address arithmetic.
Note that rather than simply assuming sizeof(void *) <= sizeof(long), we check this at program invocation, in
util.c::check_nbits_in_types()>
*/

#define ALLOC_INT(_p,_n)	(int           *)realloc(_p,(_n)*sizeof(int           )+256)
#define ALIGN_INT(_p)		(int           *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_UINT(_p,_n)	(uint32        *)realloc(_p,(_n)*sizeof(uint32        )+256)
#define ALIGN_UINT(_p)		(uint32        *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_INT64(_p,_n)	(int64         *)realloc(_p,(_n)*sizeof(int64         )+256)
#define ALIGN_INT64(_p)		(int64         *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_UINT64(_p,_n)	(uint64        *)realloc(_p,(_n)*sizeof(uint64        )+256)
#define ALIGN_UINT64(_p)	(uint64        *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_UINT128(_p,_n)(uint128       *)realloc(_p,(_n+_n)*sizeof(uint64     )+256)
#define ALIGN_UINT128(_p)	(uint128       *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_FLOAT(_p,_n)	(float         *)realloc(_p,(_n)*sizeof(float         )+256)
#define ALIGN_FLOAT(_p)		(float         *)(((intptr_t)(_p) | 63)+1)

#define ALLOC_DOUBLE(_p,_n)	(double        *)realloc(_p,(_n)*sizeof(double        )+512)
#define ALIGN_DOUBLE(_p)	(double        *)(((intptr_t)(_p) | 127)+1)

#define ALLOC_f128(_p,_n)	(__float128    *)realloc(_p,(_n)*sizeof(__float128    )+512)
#define ALIGN_f128(_p)		(__float128    *)(((intptr_t)(_p) | 127)+1)

#define ALLOC_COMPLEX(_p,_n)(struct complex*)realloc(_p,(_n)*sizeof(struct complex)+512)
#define ALIGN_COMPLEX(_p)	(struct complex*)(((intptr_t)(_p) | 127)+1)

// Vector-double|uint64-alloc used by SIMD builds; register size difference between YMM and XMM taken care of by def of vec_dbl in types.h:
#ifdef USE_SSE2

	#define ALLOC_VEC_DBL(_p,_n)(vec_dbl*)realloc(_p,(_n)*sizeof(vec_dbl)+512)
	#define ALIGN_VEC_DBL(_p)	(vec_dbl*)(((intptr_t)(_p) | 127)+1)

	#define ALLOC_VEC_U64(_p,_n)(vec_u64*)realloc(_p,(_n)*sizeof(vec_u64)+512)
	#define ALIGN_VEC_U64(_p)	(vec_u64*)(((intptr_t)(_p) | 127)+1)

#else	// In scalar-mode simply use the above double|uint64 macros:

	#define ALLOC_VEC_DBL(_p,_n)	ALLOC_DOUBLE(_p,_n)
	#define ALIGN_VEC_DBL(_p)		ALIGN_DOUBLE(_p)

	#define ALLOC_VEC_U64(_p,_n)	ALLOC_UINT64(_p,_n)
	#define ALIGN_VEC_U64(_p)		ALIGN_UINT64(_p)

#endif

#define ALLOC_POINTER(_p,_ptr_type,_n)(_ptr_type*)realloc(_p,(_n)*sizeof(_ptr_type)+64)
#define ALIGN_POINTER(_p,_ptr_type)	  (_ptr_type*)(((intptr_t)(_p) | 63)+1)

#define ALLOC_QFLOAT(_p,_n)	ALLOC_UINT128(_p,_n)
#define ALIGN_QFLOAT(_p)	ALIGN_UINT128(_p)

/*
 On the x86 family, alignment of the stack is very important
 This uses the GNU gcc  __builtin_alloca function to align doubles properly
 This is taken from GNU/FFTW package
*/
#ifdef COMPILER_TYPE_GCC
#	if (defined(__i386))
#		define HACK_ALIGN_STACK_EVEN(){					\
		if( (((uint64) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);\
		}

#		define HACK_ALIGN_STACK_ODD() {					\
		if(!(((uint64) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);\
		}
#	else
#		define HACK_ALIGN_STACK_EVEN() /* */
#		define HACK_ALIGN_STACK_ODD() /* */
#	endif
#else
#	define HACK_ALIGN_STACK_EVEN() /* */
#	define HACK_ALIGN_STACK_ODD() /* */
#endif


#endif	/* align_h_included */
