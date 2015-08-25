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
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef radix128_included
#define radix128_included

#include "radix64.h"

	#define c128_1 ((double)0.99879545620517239271)
	#define s128_1 ((double)0.04906767432741801425)	/* exp(1*I*twopi/128) */		
	#define c128_3 ((double)0.98917650996478097345)
	#define s128_3 ((double)0.14673047445536175165)	/* exp(3*I*twopi/128) */		
	#define c128_5 ((double)0.97003125319454399260)
	#define s128_5 ((double)0.24298017990326388994)	/* exp(5*I*twopi/128) */		
	#define c128_7 ((double)0.94154406518302077841)
	#define s128_7 ((double)0.33688985339222005068)	/* exp(7*I*twopi/128) */		
	#define c128_9 ((double)0.90398929312344333158)
	#define s128_9 ((double)0.42755509343028209431)	/* exp(9*I*twopi/128) */		
	#define c128_b ((double)0.85772861000027206990)
	#define s128_b ((double)0.51410274419322172658)	/* exp(b*I*twopi/128) */		
	#define c128_d ((double)0.80320753148064490981)
	#define s128_d ((double)0.59569930449243334345)	/* exp(d*I*twopi/128) */		
	#define c128_f ((double)0.74095112535495909118)
	#define s128_f ((double)0.67155895484701840061)	/* exp(f*I*twopi/128) */		

#endif	/* #ifndef radix128_included */
