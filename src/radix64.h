/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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
#ifndef radix64_included
#define radix64_included

#include "radix32.h"

	#define c64_1  ((double)0.99518472667219688624)
	#define s64_1  ((double)0.09801714032956060199)	/* exp(1*I*twopi/64) */		
	#define c64_3  ((double)0.95694033573220886494)
	#define s64_3  ((double)0.29028467725446236764)	/* exp(3*I*twopi/64) */		
	#define c64_5  ((double)0.88192126434835502971)
	#define s64_5  ((double)0.47139673682599764856)	/* exp(5*I*twopi/64) */		
	#define c64_7  ((double)0.77301045336273696081)
	#define s64_7  ((double)0.63439328416364549822)	/* exp(7*I*twopi/64) */		

#endif	/* #ifndef radix64_included */
