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
Note that rather than simply assuming sizeof(void *) == sizeof(long), we check this at program invocation, in
util.c::check_nbits_in_types()>
*/

#define ALLOC_INT(_p,_n)	(int           *)realloc(_p,(_n)*sizeof(int           )+1024)
#define ALIGN_INT(_p)		(int           *)(((long)(_p) | 63)+1)

#define ALLOC_UINT(_p,_n)	(uint          *)realloc(_p,(_n)*sizeof(uint          )+1024)
#define ALIGN_UINT(_p)		(uint          *)(((long)(_p) | 63)+1)

#define ALLOC_INT64(_p,_n)	(int64         *)realloc(_p,(_n)*sizeof(int64         )+1024)
#define ALIGN_INT64(_p)		(int64         *)(((long)(_p) | 63)+1)

#define ALLOC_UINT64(_p,_n)	(uint64        *)realloc(_p,(_n)*sizeof(uint64        )+1024)
#define ALIGN_UINT64(_p)	(uint64        *)(((long)(_p) | 63)+1)

#define ALLOC_DOUBLE(_p,_n)	(double        *)realloc(_p,(_n)*sizeof(double        )+1024)
#define ALIGN_DOUBLE(_p)	(double        *)(((long)(_p) | 127)+1)

#define ALLOC_COMPLEX(_p,_n)(struct complex*)realloc(_p,(_n)*sizeof(struct complex)+1024)
#define ALIGN_COMPLEX(_p)	(struct complex*)(((long)(_p) | 127)+1)


/*
 On the x86 family, alignment of the stack is very important
 This uses the GNU gcc  __builtin_alloca function to align doubles properly
 This is taken from GNU/FFTW package
*/
#ifdef COMPILER_TYPE_GCC
#	if (defined(__i386))
#		define HACK_ALIGN_STACK_EVEN(){					\
		if( (((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);\
		}

#		define HACK_ALIGN_STACK_ODD() {					\
		if(!(((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);\
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
