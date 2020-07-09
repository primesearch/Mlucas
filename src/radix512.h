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
#ifndef radix512_included
#define radix512_included

#include "radix256.h"

	#define c512_01 ((double)0.99992470183914454092)
	#define s512_01 ((double)0.01227153828571992607)	/* exp(01*I*twopi/512) */
	#define c512_03 ((double)0.99932238458834950089)
	#define s512_03 ((double)0.03680722294135883230)	/* exp(03*I*twopi/512) */		
	#define c512_05 ((double)0.99811811290014920712)
	#define s512_05 ((double)0.06132073630220857774)	/* exp(05*I*twopi/512) */
	#define c512_07 ((double)0.99631261218277801263)
	#define s512_07 ((double)0.08579731234443989040)	/* exp(07*I*twopi/512) */		
	#define c512_09 ((double)0.99390697000235604155)
	#define s512_09 ((double)0.11022220729388305873)	/* exp(09*I*twopi/512) */
	#define c512_0b ((double)0.99090263542778002511)
	#define s512_0b ((double)0.13458070850712618623)	/* exp(0b*I*twopi/512) */		
	#define c512_0d ((double)0.98730141815785838241)
	#define s512_0d ((double)0.15885814333386144158)	/* exp(0d*I*twopi/512) */
	#define c512_0f ((double)0.98310548743121632720)
	#define s512_0f ((double)0.18303988795514095840)	/* exp(0f*I*twopi/512) */		
	#define c512_11 ((double)0.97831737071962763313)
	#define s512_11 ((double)0.20711137619221854957)	/* exp(11*I*twopi/512) */
	#define c512_13 ((double)0.97293995220556014550)
	#define s512_13 ((double)0.23105810828067111950)	/* exp(13*I*twopi/512) */		
	#define c512_15 ((double)0.96697647104485210912)
	#define s512_15 ((double)0.25486565960451457139)	/* exp(15*I*twopi/512) */
	#define c512_17 ((double)0.96043051941556581124)
	#define s512_17 ((double)0.27851968938505310503)	/* exp(17*I*twopi/512) */		
	#define c512_19 ((double)0.95330604035419383697)
	#define s512_19 ((double)0.30200594931922806681)	/* exp(19*I*twopi/512) */
	#define c512_1b ((double)0.94560732538052132579)
	#define s512_1b ((double)0.32531029216226293393)	/* exp(1b*I*twopi/512) */		
	#define c512_1d ((double)0.93733901191257492328)
	#define s512_1d ((double)0.34841868024943456820)	/* exp(1d*I*twopi/512) */
	#define c512_1f ((double)0.92850608047321556602)
	#define s512_1f ((double)0.37131719395183754318)	/* exp(1f*I*twopi/512) */		
	#define c512_21 ((double)0.91911385169005774400)
	#define s512_21 ((double)0.39399204006104810836)	/* exp(21*I*twopi/512) */
	#define c512_23 ((double)0.90916798309052237667)
	#define s512_23 ((double)0.41642956009763718231)	/* exp(23*I*twopi/512) */		
	#define c512_25 ((double)0.89867446569395384316)
	#define s512_25 ((double)0.43861623853852763738)	/* exp(25*I*twopi/512) */
	#define c512_27 ((double)0.88763962040285394789)
	#define s512_27 ((double)0.46053871095824002336)	/* exp(27*I*twopi/512) */		
	#define c512_29 ((double)0.87607009419540660724)
	#define s512_29 ((double)0.48218377207912274823)	/* exp(29*I*twopi/512) */
	#define c512_2b ((double)0.86397285612158673808)
	#define s512_2b ((double)0.50353838372571755840)	/* exp(2b*I*twopi/512) */		
	#define c512_2d ((double)0.85135519310526514244)
	#define s512_2d ((double)0.52458968267846890591)	/* exp(2d*I*twopi/512) */
	#define c512_2f ((double)0.83822470555483804338)
	#define s512_2f ((double)0.54532498842204642200)	/* exp(2f*I*twopi/512) */		
	#define c512_31 ((double)0.82458930278502526468)
	#define s512_31 ((double)0.56573181078361319707)	/* exp(31*I*twopi/512) */
	#define c512_33 ((double)0.81045719825259479195)
	#define s512_33 ((double)0.58579785745643886000)	/* exp(33*I*twopi/512) */		
	#define c512_35 ((double)0.79583690460888353651)
	#define s512_35 ((double)0.60551104140432551359)	/* exp(35*I*twopi/512) */
	#define c512_37 ((double)0.78073722857209447856)
	#define s512_37 ((double)0.62485948814238637675)	/* exp(37*I*twopi/512) */		
	#define c512_39 ((double)0.76516726562245892617)
	#define s512_39 ((double)0.64383154288979146473)	/* exp(39*I*twopi/512) */
	#define c512_3b ((double)0.74913639452345932577)
	#define s512_3b ((double)0.66241577759017176077)	/* exp(3b*I*twopi/512) */		
	#define c512_3d ((double)0.73265427167241283493)
	#define s512_3d ((double)0.68060099779545305024)	/* exp(3d*I*twopi/512) */
	#define c512_3f ((double)0.71573082528381865446)
	#define s512_3f ((double)0.69837624940897285320)	/* exp(3f*I*twopi/512) */

#endif	/* #ifndef radix512_included */
