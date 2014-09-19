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
#ifndef radix256_included
#define radix256_included

#include "radix128.h"

	#define c256_01 ((double)0.99969881869620422011)
	#define s256_01 ((double)0.02454122852291228802)	/* exp(01*I*twopi/256) */
	#define c256_03 ((double)0.99729045667869021613)
	#define s256_03 ((double)0.07356456359966742351)	/* exp(03*I*twopi/256) */
	#define c256_05 ((double)0.99247953459870999816)
	#define s256_05 ((double)0.12241067519921619847)	/* exp(05*I*twopi/256) */
	#define c256_07 ((double)0.98527764238894124478)
	#define s256_07 ((double)0.17096188876030122632)	/* exp(07*I*twopi/256) */
	#define c256_09 ((double)0.97570213003852854447)
	#define s256_09 ((double)0.21910124015686979717)	/* exp(09*I*twopi/256) */
	#define c256_0b ((double)0.96377606579543986670)
	#define s256_0b ((double)0.26671275747489838626)	/* exp(0b*I*twopi/256) */
	#define c256_0d ((double)0.94952818059303666721)
	#define s256_0d ((double)0.31368174039889147658)	/* exp(0d*I*twopi/256) */
	#define c256_0f ((double)0.93299279883473888774)
	#define s256_0f ((double)0.35989503653498814869)	/* exp(0f*I*twopi/256) */
	#define c256_11 ((double)0.91420975570353065467)
	#define s256_11 ((double)0.40524131400498987082)	/* exp(11*I*twopi/256) */
	#define c256_13 ((double)0.89322430119551532038)
	#define s256_13 ((double)0.44961132965460659995)	/* exp(13*I*twopi/256) */
	#define c256_15 ((double)0.87008699110871141870)
	#define s256_15 ((double)0.49289819222978403677)	/* exp(15*I*twopi/256) */
	#define c256_17 ((double)0.84485356524970707332)
	#define s256_17 ((double)0.53499761988709721055)	/* exp(17*I*twopi/256) */
	#define c256_19 ((double)0.81758481315158369658)
	#define s256_19 ((double)0.57580819141784530063)	/* exp(19*I*twopi/256) */
	#define c256_1b ((double)0.78834642762660626210)
	#define s256_1b ((double)0.61523159058062684536)	/* exp(1b*I*twopi/256) */
	#define c256_1d ((double)0.75720884650648454767)
	#define s256_1d ((double)0.65317284295377676396)	/* exp(1d*I*twopi/256) */
	#define c256_1f ((double)0.72424708295146692105)
	#define s256_1f ((double)0.68954054473706692449)	/* exp(1f*I*twopi/256) */

#endif	/* #ifndef radix256_included */
