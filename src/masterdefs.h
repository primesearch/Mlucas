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
#ifndef masterdefs_h_included
#define masterdefs_h_included

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>	// Nov 2021: Add to provide POSIX case-insensitive string compare string compare strcasecmp() and strncasecmp();
						// cf. https://stackoverflow.com/questions/5820810/case-insensitive-string-comparison-in-c
#include <time.h>

#ifdef macintosh
	#include <console.h>	/* Macintosh CW */
#endif

#undef  EWM_DEBUG
#define EWM_DEBUG		0	/* Set = 1 to turn on various debugging diagnostics, especially DBG_ASSERT, defined in util.c . */

/* cf. util.h|c : If debug enabled, alias DBG_ASSERT to ASSERT (a function defined
in util.c), otherwise alias the entire 4-argument DBG_ASSERT invocation to "Bolivian"
(to paraphrase ex-heavyweight boxing champ Mike Tyson.) */
#if EWM_DEBUG
	#define DBG_ASSERT ASSERT
	#define DBG_WARN   WARN
	#define DBG_INFO   INFO
#else	/* Bolivian - lump both the FILE and LINE args together as a single __here, that's why it looks like these take 1 less arg than the underlying functions: */
	#define DBG_ASSERT(__arg1, __arg2)	/* */
	#define DBG_WARN(__here, __arg2, __arg3, __arg4)	/* */
	#define DBG_INFO(__here, __arg2, __arg3, __arg4)	/* */
#endif

/*******************************************************************************
   Mlucas-specific master #defines:
*******************************************************************************/

/* Set = 1 to do a simple FFT/IFFT-returns-original-inputs test
(sans weighting and dyadic squaring) using pseudorandom inputs:
*/
#undef  FFT_DEBUG
#define FFT_DEBUG	0

#undef  NOBRANCH
#define NOBRANCH	1	/* Switch between branched and branchless versions of various key sequences. */

#ifndef	LO_ADD
	#define	LO_ADD		1	/* TRUE = use algorithm with more mul and fewer add */
#endif

#undef	N_LEADING_RADICES
#define	N_LEADING_RADICES	8	/* # of intervals we split adjacent power-of-2 transform lengths into */

#endif	/* masterdefs_h_included */
