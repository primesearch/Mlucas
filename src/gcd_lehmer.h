/*******************************************************************************
*                                                                              *
*   (C) 1997-2015 by Ernst W. Mayer.                                           *
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
#ifndef gcd_lehmer_h_included
#define gcd_lehmer_h_included

#include "Mlucas.h"
#include "genFFT_mul.h"

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* gcd_lehmer.c: */
uint32	mi64_gcd(
	uint64 u[], uint64 v[], uint32 const ndim,
	const uint32 EGCD, uint64 Ap[], uint64 Bp[], uint32 *len_AB, uint32 *sign_AB,
	const uint32 HALF, uint64 Cp[], uint64 Dp[], uint32 *len_CD, uint32 *sign_CD, const uint32 len_targ);

uint32	matrix_vector_product_sub(uint64c abmul[], uint64c cdmul[], uint64 *uv_ptr[], uint32 len);
uint32	matrix_vector_product_add(uint64c abmul[], uint64c cdmul[], uint64 *uv_ptr[], uint32 len);

int		CMP_LT_PROD192	(uint64 a, uint64 xlo, uint64 xhi, uint64 b, uint64 ylo, uint64 yhi);
int		pprime192		(uint192 p, uint64 z);
uint192	bitwise_mod192	(uint192 x, uint192 y);
/*
void	mv_dwtvarbase_to_int64	(x,p,m,u,ndim);
*/
void	gcd_init();
int		test_gcd();

#ifdef __cplusplus
}
#endif

#endif	/* gcd_lehmer_h_included */

