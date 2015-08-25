/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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
#ifndef gen_fft_h_included
#define gen_fft_h_included

#ifdef __cplusplus
extern "C" {
#endif

/* Enumeration constant of the various supported values for the MODE argument to genFFT_mul().
As the use of an enum implies, these modes are mutually exclusive:

	Mode				Description
	----------------	-------------------------
	INIT_ARRAYS			Init FFT-related bit-reversal-index and roots-of-unity data, using input x-array for scratch storage

	The rest assume the function has been previously called in INIT_ARRAYS mode for the FFT length in question:

	FORWARD_FFT_ONLY	The fFFT of the input X-array is computed and stored in-place
	AUTO_SQUARE			The fFFT of the input X-array is computed, followed by a wrapper/dyadic-square step and an iFFT, all in-place.
	MUL_PRECOMPUTED		The X-array is assumed to contain an untransformed input vector, and the Y-array to contain a data vector which was previously-transformed by calling this routine in FORWARD_FFT_ONLY mode. The fFFT of the input X-array is computed, followed by a wrapper/dyadic-mul-with-Y-transform step and an iFFT. The result is returned in X; Y is unaffected. (I.e. this is designed for the common case where we have a constant vector which will be used to multiply many sets of inouts).

*/
enum mode {INIT_ARRAYS, FORWARD_FFT_ONLY, AUTO_SQUARE, MUL_PRECOMPUTED};

/* genFFT_mul.c: */
void  genFFT_mul(double x[], double y[], int n, int INIT_ARRAYS, int MODE);
void  genFFT_mul_process_chunk(double a[], double ab_mul[], double cd_mul[], int n, struct complex rt0[], struct complex rt1[], int index[], int ii, int nradices_prim, int radix_prim[], int MODE);

void pairFFT_mul(double x[], double y[], int n, int INIT_ARRAYS, int FORWARD_FFT_ONLY);
void pairFFT_mul_process_chunk(double a[], double ab_mul[], double cd_mul[], int n, struct complex rt0[], struct complex rt1[], int index[], int ii, int nradices_prim, int radix_prim[], int FORWARD_FFT_ONLY, int skip_square);
void radix16_pairFFT_mul(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int INIT_ARRAYS, int FORWARD_FFT_ONLY, int skip_square);


/* The complex/rel wrapper and dyadic-mul step, combined with the final-fFFt/initial-iFFT radix pass: */
void	radix16_genFFT_wrapper_mul(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int MODE);
void	radix32_genFFT_wrapper_mul(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int MODE);

#ifdef __cplusplus
}
#endif

#endif	/* gen_fft_h_included */

