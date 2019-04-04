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

#include "util.h"

/* Set the value of the round constant used for fast NINT emulation: */
void	get_fp_rnd_const(double*RND_A, double*RND_B)
{
#if(FP_MANTISSA_BITS_DOUBLE == 64)	/* X86 64-mantissa-bit register doubles: */
	*RND_A = 3.0*0x4000000*0x2000000*0x800;
	*RND_B =12.0*0x2000000*0x1000000*0x800;
	fprintf(stderr,"INFO: using 64-bit-significand form of floating-double rounding constant for scalar-mode DNINT emulation.\n");
#else	/* These assume IEEE64-compliant double-precision hardware arithmetic: */
	*RND_A = 3.0*0x4000000*0x2000000;
	*RND_B =12.0*0x2000000*0x1000000;
	fprintf(stderr,"INFO: using 53-bit-significand form of floating-double rounding constant for scalar-mode DNINT emulation. \n");
#endif
}

