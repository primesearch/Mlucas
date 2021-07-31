/*******************************************************************************
*                                                                              *
*   (C) 1997-2020 by Ernst W. Mayer.                                           *
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

#include "Mlucas.h"

/***************/

/*
Macro versions of these are in pair_square.h, since radix32_wrapper_square.c also needs to inline those;
SSE2 macros for this are in sse2_macro_gcc64.h.
*/
void pair_square(double *x1, double *y1, double *x2, double *y2, double c, double s)
{
/*
!   Given complex scalars H[j] = (x1,y1) and H[N-j] = (x2,y2) along with complex exponential E = (c,s),
!   calculates I[j] = H[j]^2 + {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4 and its complex conjugate I~,
!   returns the former in H[j] and the latter in H[N-j].
*/
	// Use that (H[j] - H~[N-j])^2 = H(j)^2 - 2*H(j)*H~(N-j) + H~(N-j)^2 to efficiently compute both (H[j]-H~[N-j])^2 and H[j]^2:
#if 0
	double rt0,rt1,rt2,rt3,it1,it2,it3;
	// H[j] = (r1,i1); H[N-j] = (r2,i2):
	rt1 = *x1;	it1 = *y1;	rt2 = *x2;	it2 = *y2;	// H[j]-H~[N-j] = (r1-r2,i1+i2); ()^2 = [(r1-r2)^2-(i1+i2)^2] + 2.I.[(r1-r2).(i1+i2)]
												// = [(r1^2-i1^2) + (r2^2-i2^2) - 2.(r1.r2+i1.i2)] + 2.I.[(r1.i1-r2.i2) - (i1.r2-r1.i2)]
	// Calculate cross product terms:
	rt3 = rt1*rt2 + it1*it2; rt3 = rt3 + rt3;	// 2.(r1.r2 + i1.i2)
	it3 = it1*rt2 - rt1*it2; it3 = it3 + it3;	// 2.(i1.r2 - r1.i2)
	// Now calculate square terms and store back in the same temporaries:
	rt0 = (rt1 + it1)*(rt1 - it1); it1 = rt1*it1; it1 = it1 + it1; rt1 = rt0;	// rt1,it1 = (r1^2-i1^2); 2.r1.i1
	rt0 = (rt2 + it2)*(rt2 - it2); it2 = rt2*it2; it2 = it2 + it2; rt2 = rt0;	// rt2,it2 = (r2^2-i2^2); 2.r2.i2
	// {1 + exp(2*pi*I*j/N)}*{H[j]-H~[N-j]}^2/4 :
	rt3 = rt1 + rt2 - rt3;	// Re(H[j]-H~[N-j])
	it3 = it1 - it2 - it3;	// Im(H[j]-H~[N-j])
	rt0 = ((c + 1.0)*rt3 - s*it3)*0.25;
	it3 = (s*rt3 + (c + 1.0)*it3)*0.25;
	// And now complete and store the results:
	*x1 = (rt1 - rt0);	// Re(I[j])
	*y1 = (it1 - it3);	// Im(I[j])
	// N-j terms are as above, but with the replacements: rt1<-->rt2, it1<-->it2, it3|-->-it3:
	*x2 = (rt2 - rt0);
	*y2 = (it2 + it3);
// Cost: [22 add, 12 mul], compared to [18 add, 18 mul] for generic-mul version ... seems too add-heavy.
#elif 0	// Quick test of mul version of this function, using square inputs:
	double re,im,tt;
/*...gather the 4 complex elements which are to be combined...*/
		//	Re{H[j]}	Im{H[j]}	Re{I[j]}	Im{I[j]}	Re{H[N-j]}	Im{H[N-j]}	Re{I[N-j]}	Im{I[N-j]}
	double r1 = *x1,	i1 = *y1,	r2 = *x1,	i2 = *y1,	r3 = *x2,	i3 = *y2,	r4 = *x2,	i4 = *y2;
// calculate 2nd square-like term and store in temp...
	re = r3*r4 - i3*i4;	// re := Re{H(n2-j)*I(n2-j)}
	im = r3*i4 + i3*r4;	// im := Im{H(n2-j)*I(n2-j)}
// calculate difference terms...
	r3 = r1 - r3;		// r3 := Re{H(j)-H~(n2-j)}
	i3 = i1 + i3;		// i3 := Im{H(j)-H~(n2-j)}
	r4 = r2 - r4;		// r4 := Re{I(j)-I~(n2-j)}
	i4 = i2 + i4;		// i4 := Im{I(j)-I~(n2-j)}
// now calculate 1st square-like term and store back in H(j) slot...
	tt = r1*r2 - i1*i2;			// r1 := Re{H(j)*I(j)}
	i1 = r1*i2 + i1*r2; r1 = tt;// i1 := Im{H(j)*I(j)}
// calculate the complex products to build the second term...
	tt = r3*r4 - i3*i4;			// Re{(H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])}
	i3 = r3*i4 + i3*r4; r3 = tt;// Im{(H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])}
	tt = ((c + 1.0)*r3 - s*i3)*0.25;	// Re{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])/4}
	i3 = (s*r3 + (c + 1.0)*i3)*0.25;	// Im{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])/4}
// and now complete and store the results.
	*x1 = (r1-tt);	// Re{M(j)}
	*y1 = (i1-i3);	// Im{M(j)}
// N-j terms are as above, but with the replacements: r1<-->r2, i1<-->i2, i3|-->-i3.
	*x2 = (re-tt);	// Re{M(N-j)}
	*y2 = (im+i3);	// Im{M(N-j)}
#else
	double re,im,tt, r1 = *x1, i1 = *y1, r2 = *x2, i2 = *y2, cc = (c + 1.0)*0.25, ss = s*0.25;
	// H[j]-H~[N-j] = (r1-r2,i1+i2); ()^2 = [(r1-r2)^2-(i1+i2)^2] + 2.I.[(r1-r2).(i1+i2)]
// calculate 2nd square-like term and store in temp...
	re = (r2+i2)*(r2-i2);	// re := Re{H(n2-j)^2}
	im = r2*i2 + i2*r2;		// im := Im{H(n2-j)^2}
// calculate difference terms...
	r2 = r1 - r2;			// r2 := Re{H(j)-H~(n2-j)}
	i2 = i1 + i2;			// i2 := Im{H(j)-H~(n2-j)}
// now calculate 1st square-like term and store back in H(j) slot...
	tt = (r1+i1)*(r1-i1);		// r1 := Re{H(j)^2}
	i1 = r1*i1 + i1*r1; r1 = tt;// i1 := Im{H(j)^2}
// calculate the complex products to build the second term...
	tt = (r2+i2)*(r2-i2);		// Re{(H[j] - H~[N/2-j])^2}
	i2 = r2*i2 + i2*r2; r2 = tt;// Im{(H[j] - H~[N/2-j])^2}
	tt = (cc*r2 - ss*i2);	// Re{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])^2/4}
	i2 = (ss*r2 + cc*i2);	// Im{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])^2/4}
// and now complete and store the results.
	*x1 = (r1-tt);	// Re{M(j)}
	*y1 = (i1-i2);	// Im{M(j)}
// N-j terms are as above, but with the replacements: r1<-->r2, i1<-->i2, i3|-->-i3.
	*x2 = (re-tt);	// Re{M(N-j)}
	*y2 = (im+i2);	// Im{M(N-j)}
// Cost: [19 add, 15 mul] ... or [16 add, 18 mul] if replace re-part-of-cmuls (r+i)*(r-i) with r^2-i^2.
// Can save another [2 add, 2 mul] by precomputing cc = (c + 1.0)/4 and ss = s/4.
#endif
}

// Jul 2019: This routine adapted from my vintage 1999 mersenne_pm1.f90 code, with input-indec swap 2 <--> 3:
void pair_mul(
	double *x1, double *y1, double *x2, double *y2, const double sx3, const double sy3, const double sx4, const double sy4,
	const double c, const double s)
{
/*
!   Given complex scalars H[j] = (x1,y1), H[N-j] = (x2,y2) and (const)I[j] = (x3,y3), I[N-j] = (x4,y4)
!   along with complex exponential E = (c,s),
!   calculates M[j] = H[j]*I[j] + {1 + exp(4*pi*I*j/N)}*{H[j]-H~[N-j]}*{I[j]-I~[N-j]}/4 and its complex conjugate M~,
!   returns the former in H[j] and the latter in H[N-j], thus overwriting those non-const inputs.
*/
	double re,im,tt, cc = (c + 1.0)*0.25, ss = s*0.25;
/*...gather the 4 complex elements which are to be combined...*/
		//	Re{H[j]}	Im{H[j]}	Re{H[N-j]}	Im{H[N-j]}	Re{I[j]}	Im{I[j]}	Re{I[N-j]}	Im{I[N-j]}
	double r1 = *x1,	i1 = *y1,	r2 = *x2,	i2 = *y2,	r3 = sx3,	i3 = sy3,	r4 = sx4,	i4 = sy4;

/*...Have: H, H~, I, I~	need: H*I, H~*I~, H - H~, I - I~. Use the sequence:
	Find H~I~, store in tmp
	Find H-H~, store in H~
	Find I-I~, store in I~
	Find HI, store in H
	Store H~I~ in I
*/
// calculate 2nd square-like term and store in temp...
	re = r2*r4 - i2*i4;	// re := Re{H(n2-j)*I(n2-j)}
	im = r2*i4 + i2*r4;	// im := Im{H(n2-j)*I(n2-j)}
// calculate difference terms...
	r2 = r1 - r2;		// r2 := Re{H(j)-H~(n2-j)}
	i2 = i1 + i2;		// i2 := Im{H(j)-H~(n2-j)}
	r4 = r3 - r4;		// r4 := Re{I(j)-I~(n2-j)}
	i4 = i3 + i4;		// i4 := Im{I(j)-I~(n2-j)}
// now calculate 1st square-like term and store back in H(j) slot...
	tt = r1*r3 - i1*i3;			// r1 := Re{H(j)*I(j)}
	i1 = r1*i3 + i1*r3; r1 = tt;// i1 := Im{H(j)*I(j)}
// calculate the complex products to build the second term...
	tt = r2*r4 - i2*i4;			// Re{(H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])}
	i2 = r2*i4 + i2*r4; r2 = tt;// Im{(H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])}
	tt = (cc*r2 - ss*i2);	// Re{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])/4}
	i2 = (ss*r2 + cc*i2);	// Im{(1 + exp(4*pi*I*j/N)) * (H[j] - H~[N/2-j])*(I[j] - I~[N/2-j])/4}
// and now complete and store the results.
	*x1 = (r1-tt);	// Re{M(j)}
	*y1 = (i1-i2);	// Im{M(j)}
// N-j terms are as above, but with the replacements: r1<-->r3, i1<-->i3, i2|-->-i2.
	*x2 = (re-tt);	// Re{M(N-j)}
	*y2 = (im+i2);	// Im{M(N-j)}
// Cost: 16 add, 16 mul [Ignoring the (1 add, 2 mul) cost of the cc,ss precomputation]
}

