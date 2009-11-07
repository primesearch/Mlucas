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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef mi64_h_included
#define mi64_h_included

#include "types.h"
#include "imul_macro.h"
#include "util.h"

/*******************************************************************************
   Function prototypes. The corresponding function definitions will either
   be in a {function name}.c file or (for cases where a .c file contains
   multiple function definitions) in the given .c file:
*******************************************************************************/

/* Set-X-equal-to-Y, Set-X-equal-to-scalar-A, Set-X-equal-to-0, set selected bit: */
void	mi64_set_eq				(uint64 x[], const uint64 y[], uint32 len);
void	mi64_set_eq_scalar		(uint64 x[], const uint64   a, uint32 len);
void	mi64_clear				(uint64 x[], uint32 len);
void	mi64_set_bit			(uint64 x[], uint32 bit);
int		mi64_test_bit			(const uint64 x[], uint32 bit);

/* bitwise shift y = (x <<|>> nbits): */
void	mi64_shl 				(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);
void	mi64_shrl				(const uint64 x[], uint64 y[], uint32 nbits, uint32 len);

/* unsigned compare: these are all built from just two elemental compares: < and == */
uint32	mi64_cmpult				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmpeq				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmplt_scalar		(const uint64 x[], uint64 a, uint32 len);
uint32	mi64_cmpgt_scalar		(const uint64 x[], uint64 a, uint32 len);
uint32	mi64_cmpeq_scalar		(const uint64 x[], uint64 a, uint32 len);
/* 11/14/2008: replace with macros:
uint32	mi64_cmpugt				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmpule				(const uint64 x[], const uint64 y[], uint32 len);
uint32	mi64_cmpuge				(const uint64 x[], const uint64 y[], uint32 len);
*/
#define	mi64_cmpugt(x, y, len)	 mi64_cmpult(y, x, len)
#define	mi64_cmpule(x, y, len)	!mi64_cmpugt(x, y, len)
#define mi64_cmpuge(x, y, len)	!mi64_cmpult(x, y, len)

/* leading & trailing binary zeros: */
uint32	mi64_trailz				(const uint64 x[], uint32 len);
uint32	mi64_leadz				(const uint64 x[], uint32 len);

/* Extract leading 64 significant bits, or 128-most-significant-bits-starting-at-a-specified-bit-position. */
uint32	mi64_extract_lead64		(const uint64 x[], uint32 len, uint64*result);
double	mi64_cvt_double			(const uint64 x[], uint32 len);
void 	mi64_extract_lead128	(const uint64 x[], uint32 len, uint32 nshift, uint64 lead_x[]);

/* Compare-to-zero and find-leading-nonzero-element: */
uint32	mi64_iszero				(const uint64 x[], uint32 len);
uint32	mi64_getlen				(const uint64 x[], uint32 len);
void	mi64_setlen				(uint64 x[], uint32 oldlen, uint32 newlen);

/* add/subtract (vector +- vector) : */
uint64	mi64_add				(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
uint64	mi64_sub				(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);

/* arithmetic and logical negation: */
void	mi64_nega				(const uint64 x[], uint64 y[], uint32 len);
void	mi64_negl				(const uint64 x[], uint64 y[], uint32 len);

/* add/subtract (vector +- scalar) : */
uint64	mi64_add_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);
uint64	mi64_sub_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);

/* unsigned multiply (vector * scalar) : */
uint64	mi64_mul_scalar			(const uint64 x[], uint64 a, uint64 y[], uint32 len);

/* unsigned multiply (vector * vector) : */
void	mi64_mul_vector			(const uint64 x[], uint32 lenX, const uint64 y[], uint32 lenY, uint64 z[], uint32 *lenZ);
void	mi64_mul_vector_hi_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);
void	mi64_mul_vector_lo_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len);

/* Routines needed for multiword uint64<==>double conversion and FFT-based mul: */
uint32	mi64_cvt_uint64_double	(const uint64 x[], const uint64 y[], uint32 len, double a[]);
uint32	mi64_cvt_double_uint64	(const double a[], uint32   n, uint64 x[], uint64 y[]);

/* returns 1 if the multiword unsigned integer p is a base-z Fermat pseudoprime, 0 otherwise: */
uint32	mi64_pprimeF			(const uint64 p[], uint64 z, uint32 len);

/* Slow bit-at-a-time method to get [optionally] quotient q = x/y and [optionally] remainder r = x%y: */
void	mi64_div				(const uint64 x[], uint64 y[], uint64 q[], uint64 r[], uint32 len);
/* Fast is-divisible-by-scalar, div-by-scalar y]: */
int		mi64_is_div_by_scalar32 (const uint32 x[], uint32 a, uint32 len);
uint32	mi64_is_div_by_scalar32_x4(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 len);
uint32	mi64_is_div_by_scalar32_x8(const uint32 x[], uint32 tryq0, uint32 tryq1, uint32 tryq2, uint32 tryq3, uint32 tryq4, uint32 tryq5, uint32 tryq6, uint32 tryq7, uint32 len);
int		mi64_is_div_by_scalar64 (const uint64 x[], uint64 a, uint32 len);
uint32	mi64_div_y32			(uint64 x[], uint32 y, uint64 q[], uint32 len);

/* Basic I/O routines: */
int		convert_mi64_base10_char(char char_buf[], const uint64 x[], uint32 len);
uint64 *convert_base10_char_mi64(const char*char_buf, uint32 *len);

/* Modular powering: returns 1 if 2^p == 1 (mod q), and (optionally) returns the full residue in res[]: */
uint32	mi64_twopmodq			(const uint64 p[], uint32 len_p, uint64 q[], uint32 len_q, uint64*res);


#endif	/* mi64_h_included */
