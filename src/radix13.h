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
#ifndef radix13_included
#define radix13_included

/* The C macro in this file implements a scalar-doubles 13-DFT bas on Van Buskirk's algorithm, but here
first a copy of my original early-2000s notes (cf. e.g. radix13_difN_cy_dit1.c in the 2009 version of Mlucas)
re. the construction of a reasonably efficient 13-DFT via length-6 subconvolutions: a length-6 cyclic for the
cosine terms, and a length-6 acyclic for the sine terms.

Radix-13 DFT using CRT for 6x6 subconvos:

The cosine terms of the output need the following combos:

c1*(x1+x12)+c2*(x2+x11)+c3*(x3+x10)+c4*(x4+x9)+c5*(x5+x8)+c6*(x6+x7)
c2*(x1+x12)+c4*(x2+x11)+c6*(x3+x10)+c5*(x4+x9)+c3*(x5+x8)+c1*(x6+x7)
c3*(x1+x12)+c6*(x2+x11)+c4*(x3+x10)+c1*(x4+x9)+c2*(x5+x8)+c5*(x6+x7)
c4*(x1+x12)+c5*(x2+x11)+c1*(x3+x10)+c3*(x4+x9)+c6*(x5+x8)+c2*(x6+x7)
c5*(x1+x12)+c3*(x2+x11)+c2*(x3+x10)+c6*(x4+x9)+c1*(x5+x8)+c4*(x6+x7)
c6*(x1+x12)+c1*(x2+x11)+c5*(x3+x10)+c2*(x4+x9)+c4*(x5+x8)+c3*(x6+x7)

where cX := cos(X*2*pi/13). These give no hint as to what permutation of input and/or trig. term indices we'll
need in order to put this in the form of a convolution. For the sine terms things seem even more dicy, since there
we not only have to put the index pattern in the form of a convolution (and the same perutation that accomplishes
this for the cosine terms also does the same for the sine terms), but we also need to get the + and - signs
in the form of an acayclic convo. This is because (unlike radix-11) there doesn't appear to be a pattern of
sign flips on the inputs ((xi-xj) and sine multipliers) which makes the sine part look like a cyclic convo.
But first let's deal just with the index permutation: my program permute.c gives the permutation (and perhaps
it's not unique, since the program stops as soon as it finds the first successful permutation) (0,5,2,4,3,1)
as one which yields a convolution index pattern. Let's apply it to the sine terms, defined by

+s1*(x1-x12)+s2*(x2-x11)+s3*(x3-x10)+s4*(x4-x9)+s5*(x5-x8)+s6*(x6-x7) = S1
+s2*(x1-x12)+s4*(x2-x11)+s6*(x3-x10)-s5*(x4-x9)-s3*(x5-x8)-s1*(x6-x7) = S2
+s3*(x1-x12)+s6*(x2-x11)-s4*(x3-x10)-s1*(x4-x9)+s2*(x5-x8)+s5*(x6-x7) = S3
+s4*(x1-x12)-s5*(x2-x11)-s1*(x3-x10)+s3*(x4-x9)-s6*(x5-x8)-s2*(x6-x7) = S4
+s5*(x1-x12)-s3*(x2-x11)+s2*(x3-x10)-s6*(x4-x9)-s1*(x5-x8)+s4*(x6-x7) = S5
+s6*(x1-x12)-s1*(x2-x11)+s5*(x3-x10)-s2*(x4-x9)+s4*(x5-x8)-s3*(x6-x7) = S6 .

Letting s(1,2,3,4,5,6)
      = b(0,5,2,4,3,1) gives

+b0*(x1-x12)+b5*(x2-x11)+b2*(x3-x10)+b4*(x4-x9)+b3*(x5-x8)+b1*(x6-x7) = S1
+b5*(x1-x12)+b4*(x2-x11)+b1*(x3-x10)-b3*(x4-x9)-b2*(x5-x8)-b0*(x6-x7) = S2
+b2*(x1-x12)+b1*(x2-x11)-b4*(x3-x10)-b0*(x4-x9)+b5*(x5-x8)+b3*(x6-x7) = S3
+b4*(x1-x12)-b3*(x2-x11)-b0*(x3-x10)+b2*(x4-x9)-b1*(x5-x8)-b5*(x6-x7) = S4
+b3*(x1-x12)-b2*(x2-x11)+b5*(x3-x10)-b1*(x4-x9)-b0*(x5-x8)+b4*(x6-x7) = S5
+b1*(x1-x12)-b0*(x2-x11)+b3*(x3-x10)-b5*(x4-x9)+b4*(x5-x8)-b2*(x6-x7) = S6 ,

and rearranging the rows to reorder the b-terms of the leftmost column yields

+b0*(x1-x12)+b5*(x2-x11)+b2*(x3-x10)+b4*(x4-x9)+b3*(x5-x8)+b1*(x6-x7) = S1
+b1*(x1-x12)-b0*(x2-x11)+b3*(x3-x10)-b5*(x4-x9)+b4*(x5-x8)-b2*(x6-x7) = S6
+b2*(x1-x12)+b1*(x2-x11)-b4*(x3-x10)-b0*(x4-x9)+b5*(x5-x8)+b3*(x6-x7) = S3
+b3*(x1-x12)-b2*(x2-x11)+b5*(x3-x10)-b1*(x4-x9)-b0*(x5-x8)+b4*(x6-x7) = S5
+b4*(x1-x12)-b3*(x2-x11)-b0*(x3-x10)+b2*(x4-x9)-b1*(x5-x8)-b5*(x6-x7) = S4
+b5*(x1-x12)+b4*(x2-x11)+b1*(x3-x10)-b3*(x4-x9)-b2*(x5-x8)-b0*(x6-x7) = S2

We can see that the columns have convolution index patterns, i.e. the six distinct
circular shifts of the vector (0,1,2,3,4,5). Next, letting

a0 = (x1-x12), a1 = (x2-x11), a4 = (x3-x10), a2 = (x4-x9), a3 = (x5-x8), a5 = (x6-x7),
(and doing analogously for the cosine terms, except there we don't need to worry about signs) we get

+a0.b0+a1.b5+a2.b4+a3.b3+a4.b2+a5.b1 = S1
+a0.b1-a1.b0-a2.b5+a3.b4+a4.b3-a5.b2 = S6
+a0.b2+a1.b1-a2.b0+a3.b5-a4.b4+a5.b3 = S3
+a0.b3-a1.b2-a2.b1-a3.b0+a4.b5+a5.b4 = S5
+a0.b4-a1.b3+a2.b2-a3.b1-a4.b0-a5.b5 = S4
+a0.b5+a1.b4-a2.b3-a3.b2+a4.b1-a5.b0 = S2 .

Now we need to generate the proper sign pattern for an acyclic.

Flipping the sign of b0 (i.e. s1) and then flipping the sign of the entire first row gives:

+a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 =-S1;	a0, b0 frozen
+a0.b1+a1.b0-a2.b5+a3.b4+a4.b3-a5.b2 = S6;	need to flip a3.b4,a4.b3
+a0.b2+a1.b1+a2.b0+a3.b5-a4.b4+a5.b3 = S3;	need to flip a5.b3,a3.b5
+a0.b3-a1.b2-a2.b1+a3.b0+a4.b5+a5.b4 = S5;	need to flip a0.b3,a3.b0, then entire row
+a0.b4-a1.b3+a2.b2-a3.b1+a4.b0-a5.b5 = S4;	need to flip a1.b3,a3.b1
+a0.b5+a1.b4-a2.b3-a3.b2+a4.b1+a5.b0 = S2;	need to flip a2.b3,a3.b2

Note that every one of the needed 2-term sign flips involves a3 in one term and b3 in the other, and we can flip both
without affecting the already-finished first row, so flipping signs of a3 and b3 (i.e. s5) and then of row 4:

+a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 =-S1
+a0.b1+a1.b0-a2.b5-a3.b4-a4.b3-a5.b2 = S6
+a0.b2+a1.b1+a2.b0-a3.b5-a4.b4-a5.b3 = S3
+a0.b3+a1.b2+a2.b1+a3.b0-a4.b5-a5.b4 =-S5
+a0.b4+a1.b3+a2.b2+a3.b1+a4.b0-a5.b5 = S4
+a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 = S2 .

The LHS is now in the form of a length-6 acyclic convolution. We next turn to the issue of how to efficiently perform the
two needed convolutions (length-6 cyclic for cosine terms, length-6 acyclic for sine terms), which involves factoring the
two underlying order-6 polynomials into their irreducible components and then applying the Chinese Remainder Theorem (CRT).


Given variable A = [a0,a1,a2,a3,a4,a5], constant B = [b0,b1,b2,b3,b4,b5],
compute convo(A,B) modulo P := x^6 - 1 (cyclic) and modulo Q := x^6 + 1 (acyclic).
As our trial vectors for debugging, we use A = [3,-1,4,-1,-5,-9] and B = [2,-7,-1,8,2,8],
and (for checksum purposes) assume these digits are with respect to a base-10 expansion,
i.e.
	A = 3+10*(-1+10*( 4+10*(-1+10*(-5+10*(-9))))) = -950607
	B = 2+10*(-7+10*(-1+10*( 8+10*( 2+10*( 8))))) = +827832 .

********* CYCLIC CONVO: **********

Outputs are, in terms of coefficients of various powers of x:

x^0: a0.b0+a1.b5+a2.b4+a3.b3+a4.b2+a5.b1 =  66
x^1: a0.b1+a1.b0+a2.b5+a3.b4+a4.b3+a5.b2 = -24
x^2: a0.b2+a1.b1+a2.b0+a3.b5+a4.b4+a5.b3 = -78
x^3: a0.b3+a1.b2+a2.b1+a3.b0+a4.b5+a5.b4 = -63
x^4: a0.b4+a1.b3+a2.b2+a3.b1+a4.b0+a5.b5 = -81
x^5: a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 =  72


P factors into x^6-1 = (x^3-1)*(x^3+1) = (x-1)*(x+1)*(x^2-x+1)*(x^2+x+1) := P0*P1*P2*P3,
so first reduce both input polys modulo the factors, using that:

modulo P0, x^n == +1 ;
modulo P1, x^n == +1 if n even, -1 if n odd ;
modulo P2, x^2 == +x-1, x^3 == -1, x^4 == -x, x^5 == -x^2 == -x+1 ;
modulo P3, x^2 == -x-1, x^3 == +1, x^4 == +x, x^5 == +x^2 == -x-1 .

											COST:
A mod P0 :=     x-1 = a0+a1+a2+a3+a4+a5			-9	5 ADD (form 0+2+4, 1+3+5, add)
      P1 :=     x+1 = a0-a1+a2-a3+a4-a5			13	1 ADD (subtract the two 3-term sums found above)
      P2 := x^2-x+1 = [a0-a2-a3+a5] + [a1+a2-a4-a5]*x	-9+17x	8 ADD (form d=0-2, e=3-5, f=1-5, g=2-4,
      P3 := x^2+x+1 = [a0-a2+a3-a5] + [a1-a2+a4-a5]*x	 7-  x		 then A mod P2 = (d-e)+(f+g)*x
												A mod P3 = (d+e)+(f-g)*x .

similar for B, but we can precompute the latter, so they're basically free.	12, -6, 3-18x, 3-12x

Check: for x := 10,
A%(x-1) = A%9 = 0 == 9, A%(x+1) = A%11 = 2 == 13, A%(x^2-x+1) = A%91 = -21 == -9+17*10, A%(x^2+x+1) = A%111 = -3 = 7-10, ok.
B%(x-1) = B%9 = 3 ==12, B%(x+1) = B%11 = 5 == -6, B%(x^2-x+1) = B%91 =   5 ==  3-18*10, B%(x^2+x+1) = B%111 =105 = 3-12*10, ok.

				COST:
polypro(A,B) mod P0:	0 ADD, 1 MUL, output = p00		-108 ==   0 mod   9, ok
                 P1:	0 ADD, 1 MUL, output = p10		-78  ==  10 mod  11, ok
                 P2:	3 ADD, 3 MUL, output = p20 + x*p21	(-9.3) + (9.18+17.3)*x + (-17.18).x^2 = -27 + 213.x - 306.x^2
                                                             == -27 + 213.x - 306.(x-1) modulo P2
                                                              =  279 - 93.x == -14 mod  91, ok
                 P3:	3 ADD, 3 MUL, output = p30 + x*p31 .	( 7.3) + (-7.12-1.3)*x + (  1.12).x^2 =  21 -  87.x +  12.x^2
                                                             ==  21 -  87.x -  12.(x+1) modulo P3
                                                              =   9 -  99.x ==  18 mod 111, ok
LENGTH-2 SUBCONVOS: full-length polypro output is

(a+b.x)*(c+d.x) = a.c + (b.c+a.d).x + b.d.x^2, cyclic convo is modulo x^2-1,
 i.e. folds the x^2 coeff. in with the x^0, giving (a.c+b.d) + (b.c+a.d).x .

Now,

(a+b).(c+d) = a.c+a.d+b.c+b.d
(a-b).(c-d) = a.c-a.d-b.c+b.d

so to get convo using just 2 muls, precompute e = (c+d)/2, f = (c-d)/2 (free), then calculate

y0 = (a+b).e
y1 = (a-b).f
z0 = y0+y1	(x^0 output)
z1 = y0-y1	(x^1 output).

Now, modulo P2 = x^2-x+1, our convo output is (a.c-b.d) + (b.c+a.d + b.d).x,
i.e. x^0 has 2.b.d subtracted w.r.to convo output, x^1 has extra b.d added in.
In fact this looks like a complex multiply (x = I) with the extra b.d. added to the imaginary output.
So we use a modified Karatsuba: precompute (c+d), then calculate
					a=-9, b=17, c=3, d=-18
m0 = a.c				-27		(d-e).(m-n)
m1 = b.d				-306		(f+g).(o+p)	a+b = d-e+f+g = (d+f) + (g-e)
y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d	8.-15 = -120	NB: a+b = [a0-a2-a3+a5] + [a1+a2-a4-a5] = a0+a1-a3-a4 = (a0-a4) + (a1-a3) = u + v
							similarly, c+d = x + y, so the product = u.x+v.x+u.y+v.y
z0 = m0-m1	(x^0 output)		279
z1 = y0-m	(x^1 output).		-93	COST: 3 mul, 3 add	(same # total ops as length-2 cyclic convo, but equal # of mul and add.)

Modulo P3 = x^2+x+1, our convo output is (a.c-b.d) + (b.c+a.d - b.d).x,
i.e. x^0 has 2.b.d subtracted w.r.to convo output, x^1 has b.d subtracted.
To get this, we precompute (d-c), then calculate
					a=7, b=-1, c=3, d=-12
m0 = a.c				21		(d+e).(m+n)
m1 = b.d				12		(g-f).(p-o)	a-b = d+e+f-g = (d+f) - (g-e)
y0 = (a-b).(d-c) =-a.c+a.d+b.c-b.d	8.-15 = -120	NB: a-b = [a0-a2+a3-a5] - [a1-a2+a4-a5] = a0-a1+a3-a4 = (a0-a4) - (a1-a3) = u - v
							similarly, d-c = y - x, so the product =-u.x+v.x+u.y-v.y
z0 = m0-m1	(x^0 output)		9
z1 = y0+m0	(x^1 output).		-99	COST: 3 mul, 3 add


Lastly, need to combine the 4 modular polynomial product outputs according to CRT.
First, find the four inverse polys, t0-3:

t0*P1*P2*P3 == 1 modulo P0.
Now:
  P1*P2*P3 = (x+1)*(x^4+x^2+1) = x^5+x^4+x^3+x^2+x+1 == +6 modulo x-1, so t0 = +1/6.

t1*P0*P2*P3 == 1 modulo P1.
Now:
  P0*P2*P3 = (x-1)*(x^4+x^2+1) = x^5-x^4+x^3-x^2+x-1 == -6 modulo x+1, so t1 = -1/6.

t2*P0*P1*P3 == 1 modulo P2.
Now:
  P0*P1*P3 = (x^2-1)*(x^2+x+1) =     x^4+x^3    -x-1 == -2x-2  modulo x^2-x+1. Try multiplying by x:
x*P0*P1*P3                     = x^5+x^4    -x^2-x   == (-x+1)-x-(x-1)-x = -4x+2,
so (x-2)*P0*P1*P3 == -4x+2+4x+4 = 6, so t2 = (x-2)/6.

t3*P0*P1*P2 == 1 modulo P3.
Now:
  P0*P1*P2 = (x^2-1)*(x^2-x+1) =     x^4-x^3    +x-1 == +2x-2  modulo x^2+x+1. Try multiplying by x:
x*P0*P1*P2                     = x^5-x^4    +x^2-x   == -4x-2,
so (x+2)*P0*P1*P3 == -6, so t3 = -(x+2)/6.


Next, find the s-terms:

s0 := t0*P1*P2*P3 =  (x^5+x^4+x^3+x^2+x+1)/6 ,
s1 := t1*P0*P2*P3 = -(x^5-x^4+x^3-x^2+x-1)/6 ,
s2 := t2*P0*P1*P3 = (x-2)*(x^4+x^3    -x-1)/6 = (x^5-x^4-2.x^3-x^2+x+2)/6 ,
s3 := t3*P0*P1*P2 =-(x+2)*(x^4-x^3    +x-1)/6 =-(x^5+x^4-2.x^3+x^2+x-2)/6 ,

Then, the polynomial product mod (x^6 - 1) = p00*s0 + p10*s1 + (p20 + x*p21)*s2 + (p30 + x*p31)*s3,
where the signature of each p-term times the respective s-polynomial is

p00: s0*6    =  1  +x+x^2  +x^3  +x^4+x^5
p10: s1*6    =  1  -x+x^2  -x^3  +x^4-x^5
p20: s2*6    =  2  +x-x^2-2.x^3  -x^4+x^5
p21: s2*6*x ==  1+2.x+x^2  -x^3-2.x^4-x^5 mod (x^6 - 1)
p30: s3*6    =  2  -x-x^2+2.x^3  -x^4-x^5
p31: s3*6*x == -1+2.x-x^2  -x^3+2.x^4-x^5 mod (x^6 - 1) .

So, crunching all the numbers and dividing by 6 (the divide can be absorbed into the b-constants,
so is free) gives the following output coefficients, where we group p00&p10, p20&p30, p21&p31 by letting
a,b = (p00+-p10)/6, c,d = (p20+-p30)/6, e,f = (p21+-p31)/6                   p00 = -108, p10 = -78		COST: 6 ADD
                                                                             p20 =  279, p30 =   9
                                                                             p21 = - 93, p31 = -99, so
											 a = -31, b= -5, c = +48, d = +45, e = -32, f = +1

x^0: a + 2.c +   f					66
x^1: b +   d + 2.e				 -24
x^2: a -   c +   f = (x^0 coeff) - 3.(c)	 -78
x^3: b - 2.d -   e = (x^5 coeff) - 3.(d)	 -63
x^4: a -   c - 2.f = (x^2 coeff) - 3.(f)	 -81
x^5: b +   d -   e = (x^1 coeff) - 3.(e)		72, looks good!

so the total cost of the reconstruction phase = 16 add, and we can save one add by calculating (a-c)
and using it in both the x^2 and x^4 terms, and can save another add by calculating (b+d)
and using it in both the x^1 and x^5 terms, or using (b-e) in both the x^3 and x^5 terms.			COST: 14 ADD
We can also use the following sequence, which trades 4 multiplies for 6 adds, and thus
might be preferable for a floating-point implementation, where the muls could be done
alongside the adds:

x^2 = a - c + f		2 add
x^0 = x^2 + 3.c		1 add, 1 mul
x^4 = x^2 - 3.f		1 add, 1 mul
x^5 = b + d - e		2 add
x^1 = x^5 + 3.e		1 add, 1 mul
x^3 = x^5 - 3.d		1 add, 1 mul

GRAND TOTAL: 40 add, 8 mul (or 34 add, 12 mul), compared to 30 add, 36 mul for the naive scheme.
We can save some further adds by using that the X0 output of the radix-13 DFT = x0 + (A mod P0), i.e. needs just 1 more add,
and by rolling the x0 term into the above calculation (i.e. it needs just
1 add to include the real part of x0 in the p00 output, which then propagates it to all 6 x^k coeffs),
and lastly by using the fact that P1 also appears in the X0 output computation (X0 = x0 + A mod P1).


********* ACYCLIC CONVO: **********

For A = [3,-1,4,-1,-5,-9] and B = [2,-7,-1,8,2,8], outputs are, in terms of coefficients of various powers of x:

x^0: a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 = -54
x^1: a0.b1+a1.b0-a2.b5-a3.b4-a4.b3-a5.b2 = -22
x^2: a0.b2+a1.b1+a2.b0-a3.b5-a4.b4-a5.b3 = 102
x^3: a0.b3+a1.b2+a2.b1+a3.b0-a4.b5-a5.b4 =  53
x^4: a0.b4+a1.b3+a2.b2+a3.b1+a4.b0-a5.b5 =  63
x^5: a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 =  72


Q factors into x^6+1 = (x^2+1)*(x^4-x^2+1) := Q0*Q1,
so first reduce both input polys modulo the factors, using that:

modulo Q0, x^2 == -1, x^3 == -x, x^4 == +1, x^5 == +x ;
modulo Q1, x^4 == x^2-1, x^5 == x^3-x .

Check: for x := 10,
A%(x^2+1) = A%101 =  5, A%(x^4-x^2+1) = A%9901 =  -111
B%(x^2+1) = B%101 = 36, B%(x^4-x^2+1) = B%9901 = +6049
A*B%(x^2+1) = 5.36%101 = 79, A*B%(x^4-x^2+1) = 1829.
													COST:
A mod Q0 :=     x^2+1 = [a0-a2+a4] + [a1-a3+a5]*x				-6 -9.x			4 ADD
      Q1 := x^4-x^2+1 = [a0-a4] + [a1-a5]*x + [a2+a4]*x^2 + [a3+a5]*x^3		8+8.x-x^2-10.x^3	4 ADD

similar for B, but we can precompute the latter, so they're basically free.		B mod Q0 = 5-7.x, B mod Q1 = 0-15.x+x^2+16.x^3

				COST:
convo(A,B) mod Q0:	 3 MUL, 4 ADD, output = q00 + x*q01
               Q1:	12 MUL,15 ADD, output = q10 + x*q11 + x^2*q21 + x^3*q31.

LENGTH-2 SUBCONVO: full-length polypro output is

(a+b.x)*(c+d.x) = a.c + (b.c+a.d).x + b.d.x^2, acyclic convo is modulo x^2+1,
 i.e. folds the -x^2 coeff. in with the x^0, giving (a.c-b.d) + (b.c+a.d).x , which is equivalent to complex multiply with constant b,c.
So we use Karatsuba: precompute (c+d), then calculate
						a=-6, b=-9, c=5, d=-7
m0 = a.c					-30
m1 = b.d					63
y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d		30
z0 = m0-m1	(x^0 output)			q00 = -93
z1 = y0-m0-m1	(x^1 output).			COST: 3 mul, 4 add	q01 =  -3, q00+q01.10 = -93-30 == 79 mod 101, ok.

LENGTH-4 SUBCONVO: full-length polypro output is

(a+b.x+c.x^2+d.x^3)*(e+f.x+g.x^2+h.x^3) = a.e + (a.f+b.e).x + (a.g+b.f+c.e).x^2 + (a.h+b.g+c.f+d.e).x^3 + (b.h+c.g+d.f).x^4 + (c.h+d.g).x^5 + d.h.x^6,

Length-4 cyclic convo is modulo x^4-1, i.e. folds this in such a fashion as to give the following output coefficients:

x^0: a.e+b.h+c.g+d.f
x^1: a.f+b.e+c.h+d.g
x^2: a.g+b.f+c.e+d.h
x^3: a.h+b.g+c.f+d.e

Length-4 acyclic convo is modulo x^4+1, i.e. folds this in such a fashion as to give the following output coefficients:

x^0: a.e-b.h-c.g-d.f
x^1: a.f+b.e-c.h-d.g
x^2: a.g+b.f+c.e-d.h
x^3: a.h+b.g+c.f+d.e

modulo Q1 = x^4-x^2+1, the polypro outputs are (using that x^4 == x^2-1, x^5 == x^3-x, x^6 == x^4-x^2 == -1 mod Q1)
						a,b,c,d = 8,8,-1,-10	e,f,g,h = 0,-15,+1,16
x^0: a.e-b.h-c.g-d.f - d.h			8.  0 - 8. 16 - -1. +1 - -10.-15 - (-10.16) =    0-128+  1-150 + 160 = -117
x^1: a.f+b.e-c.h-d.g				8.-15 + 8.  0 - -1. 16 - -10. +1            = -120+  0+ 16+ 10       =  -94
x^2: a.g+b.f+c.e + (b.h+c.g+d.f)		8. +1 + 8.-15 + -1.  0 + (8.16-1.1+10.15)   =   +8-120+  0 + 279     = +165
x^3: a.h+b.g+c.f+d.e + (c.h+d.g)		8. 16 + 8. +1 + -1.-15 + -10.  0 + (-16-10) =  128+  8+ 15+  0 -  26 = +125

(a+b*x+c*x^2+d*x^3)*(e+f*x+g*x^2+h*x^3)%9901 = 1829, ok.

****************
08/26/02: Try the following sequence:

1) Do length-4 cyclic convo (5 mul, 15 add)
2) Via side calculation, obtain x = (d.h), y = (c.h+d.g), z = (b.h+c.g+d.f); (naive opcount is 6 mul, 3 add)
3) Obtain polypro modulo Q1 by modifying the cyclic convo outputs as follows:

	(x^0 term) - 2.z - x;
	(x^1 term) - 2.y
	(x^2 term) - 2.x
	(x^3 term)       + y, which needs an additional (3 mul, 5 add) or (8 add).

Thus, via this route the mod-Q1 polypro costs (14 mul, 23 add) or (11 mul, 26 add). Both are much add-heavier than we'd like.

****************

...doesn't look very promising. Unless can find an algorithm in Nuss. for polypro modulo Q1,
perhaps should try separating even and odd-order terms and then doing 4 subconvos modulo x^2-x+1:

									intermediates needed:			total ops needed:
(a+c.x)*(e+g.x) == (a.e-c.g) + (c.e+a.g + c.g).x mod x^2-x+1 = m0+n0.x	a.e, c.g, a+c plus 1 mul, 2 add		3 mul, 3 add		+1   +7.x
(a+c.x)*(f+h.x) == (a.f-c.h) + (c.f+a.h + c.h).x mod x^2-x+1 = m1+n1.x	a.f, c.h      plus 1 mul, 2 add		3 mul, 2 add	-104 +127.x
(b+d.x)*(e+g.x) == (b.e-d.g) + (d.e+b.g + d.g).x mod x^2-x+1 = m2+n2.x	b.e, d.g, b+d plus 1 mul, 2 add		3 mul, 3 add	 +10   -2.x
(b+d.x)*(f+h.x) == (b.f-d.h) + (d.f+b.h + d.h).x mod x^2-x+1 = m3+n3.x	b.f, d.h,     plus 1 mul, 2 add		3 mul, 2 add	 +40 +118.x

NB: a+c = a0+a2, b+d = a1+a3, more accurate to get directly from the a-inputs this way.

%%%%%%%%%%%%%% STUFF to TRY %%%%%%%%%%%%%%%%

Rewriting the mod-Q1 convo in terms of indexed vectors (a,b,c,d) = a0-3 and (e,f,g,h) = b0-3
may make it easier to look for structure in the computation:

x^0: a0.b0-a1.b3-a2.b2-a3.b1 - a3.b3
x^1: a0.b1+a1.b0-a2.b3-a3.b2
x^2: a0.b2+a1.b1+a2.b0       + (a1.b3+a2.b2+a3.b1)
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + (a2.b3+a3.b2)

Note: all elts except b1 appear either 4 or 6 times; b1 appears 5 times.

Now, regular length-4 acyclic (i.e. convo modulo x^4+1) gives

x^0: a0.b0-a1.b3-a2.b2-a3.b1		need to sub a3.b3
x^1: a0.b1+a1.b0-a2.b3-a3.b2		no change needed
x^2: a0.b2+a1.b1+a2.b0-a3.b3		need to add a3.b3 + (a1.b3+a2.b2+a3.b1), which = -[(x^0 term) - a0.b0 - a3.b3]
x^3: a0.b3+a1.b2+a2.b1+a3.b0		need to add (a2.b3+a3.b2)
					TOTAL: 4 mul, 6 add to effect corrections

If replace inputs b0,b1,b2,b3 with b0+b2,b1+b3,b2,b3, then get

x^0: a0.b0-a1.b3-a2.b2-a3.b1 + a0.b2-a3.b3	need to sub a0.b2
x^1: a0.b1+a1.b0-a2.b3-a3.b2 + a0.b3+a1.b2	need to sub a0.b3 and a1.b2
x^2: a0.b2+a1.b1+a2.b0-a3.b3 + a1.b3+a2.b2	need to add a3.b1
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + a2.b3+a3.b2	no change needed
					TOTAL: 4 mul, 4 add to effect corrections. Can we do even better?

replace b1 input with b1+b3, a2 input with a0+a2:

x^0: a0.b0-a1.b3-a2.b2-a3.b1 - a3.b3-a0.b2	need to add a0.b2
x^1: a0.b1+a1.b0-a2.b3-a3.b2			no change needed
x^2: a0.b2+a1.b1+a2.b0-a3.b3 + a1.b3+a0.b0	need to add a3.b3-a0.b0+a2.b2+a3.b1
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + a0.b1+a2.b3	need to add (a2.b3+a3.b2)

Try using sign flips to make mod-Q1 convo look more like a length-4 cyclic:

x^0: a0.b0+a1.b3+a2.b2+a3.b1 + a3.b3		chg sign of a0, then of entire row; a0, b0 frozen
x^1:-a0.b1+a1.b0-a2.b3-a3.b2
x^2:-a0.b2+a1.b1+a2.b0       + (a1.b3+a2.b2+a3.b1)
x^3:-a0.b3+a1.b2+a2.b1+a3.b0 + (a2.b3+a3.b2)


Aside:
Nussbaumer's length-3 cyclic convo algorithm needs 4 mul, 11 add; his length-3 acyclic, OTOH, needs 5 mul and a whopping 20 add.
This is insane, since we can very easily convert a length-3 acyclic into a length-3 cyclic, as follows:

length-3 acyclic:

x^0: a0.b0-a1.b2-a2.b1
x^1: a0.b1+a1.b0-a2.b2
x^2: a0.b2+a1.b1+a2.b0

x^0: a0.b0+a1.b2+a2.b1		flip sign of a0, and then of entire row; a0, b0 frozen
x^1: a0.b1+a1.b0+a2.b2		flip sign of a2 and b1 (which leaves first row unchanged)
x^2: a0.b2+a1.b1+a2.b0		flip sign of entire row, done!

For length-2 acyclic, this fails, since we can only flip an even number of signs this way, and the acyclic has an odd # of minuses.

x^0: a0.b0-a1.b1
x^1: a0.b1+a1.b0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We then get the final x^0-3 coeffs via

x^0: m0-n3		1 add		-117
x^1: m1+m2		1 add		 -94
x^2: m3+n0+n3		2 add		+165
x^3: n1+n2		1 add		+125

So the total cost for the polypro mod Q1 is 12 MUL, 15 ADD, compared to 5 mul, 15 add for simple length-4 cyclic.
(Multiply count seems high - can reduce to 9 MUL via length-4 cyclic convo with certain input sign changes + clever
combinations of correction terms, but only at the cost of ~10 more adds.)

Lastly, need to combine the 2 modular polynomial product outputs according to CRT.

First, find the two inverse polys, t0 and t1:

t0*Q1 == 1 modulo Q0.
Now:
  1*Q1 = x^4-x^2+1 == 3, so t0 = 1/3.

t1*Q0 == 1 modulo Q1.
Now:
  1*Q0 = x^2+  1,
  x*Q0 = x^3+  x,
x^2*Q0 = x^4+x^2 == 2.x^2-1 modulo Q1, so (x^2-2  )*Q0 == -3 modulo Q1, so t1 = (2-x^2)/3 .

Next, find the s-terms:

s0 := t0*Q1 mod Q = Q1/3              = (1-x^2+x^4)/3 ;
s1 := t1*Q0 mod Q = (2-x^2)*(1+x^2)/3 = (2+x^2-x^4)/3 .

product mod (x^6 + 1) = (q00 + x*q01)*s0 + (q10 + x*q11 + x^2*q12 + x^3*q13)*s1,
where the signature of each q-term times the respective t-polynomial is

q00: s0*3      =  1      -x^2      +x^4
q01: s0*3*x   ==      x        -x^3    +x^5
q10: s1*3      =  2      +x^2      -x^4
q11: s1*3*x   ==    2.x        +x^3    -x^5
q12: s1*3*x^2 ==  1    +2.x^2      +x^4     mod (x^6 + 1)
q13: s1*3*x^3 ==      x      +2.x^3    +x^5 mod (x^6 + 1) .

So, crunching all the numbers and dividing by 3 (the divide can be absorbed into the b-constants,
so is free) gives the following output coefficients, where we define a,b = q00+-q10, c = q10+q12, d,e = q01+-q11, f = q11+q13:
						q00,01 = -93,-3	q10,11,12,13 = -117,-94,+165,+125	a,b,c = -210,+24,+48	d,e,f = -97,+91,+31
x^0:   q00      +2.q10        +q12       = a+c					-162/3 = -54, ok
x^1:         q01      +2.q11        +q13 = d+f					 +66/3 = +22, ok
x^2:  -q00        +q10      +2.q12       = c+c-a				+306/3 = 102, ok
x^3:        -q01        +q11      +2.q13 = f+f-d				+159/3 = +53, ok
x^4:   q00        -q10        +q12       = b+q12 (don't really need b)	+189/3 = +63, ok
x^5:         q01        -q11        +q13 = e+q13 (don't really need e)	+216/3 = +72, ok. It works!

So the total cost of the reconstruction = 14 add.
We can also use the following sequence, which trades 4 multiplies for 6 adds, and thus
might be preferable for a floating-point implementation, where the muls could be done
alongside the adds:

x^4 = q00 - q10 + q12	2 add
x^0 = x^4 + 3.q10	1 add, 1 mul
x^2 = 3.q12 - x^4	1 add, 1 mul
x^5 = q01 - q11 + q13	2 add
x^1 = x^5 + 3.q11	1 add, 1 mul
x^3 = 3.q13 - x^5	1 add, 1 mul


GRAND TOTAL: 41 add, 15 mul (or 35 add, 19 mul), compared to 30 add, 36 mul for the naive scheme.


To do the complete radix-13 DFT, we then need to do as follows:

!   We refer to the terms C1,2,3,4,5,6 (which do not explicitly involving the imaginary constant I)
!   as the "cosine part" of the output, and S1,2,3,4,5,6 (those multiplied by I) as the "sine part."
!                                                                                   opcount for general odd-prime radix R:
!   Form      (x1+-x12),  (x2+-x11),  (x3+-x10),  (x4+-x9),  (x5+-x8),  (x6+-x7) :  0 FMUL, 24 FADD
!   Form X0                                                                      :           2 FADD
!   Form x0+c1*(x1+x12)+c2*(x2+x11)+c3*(x3+x10)+c4*(x4+x9)+c5*(x5+x8)+c6*(x6+x7) :
!   Form x0+c2*(x1+x12)+c4*(x2+x11)+c6*(x3+x10)+c5*(x4+x9)+c3*(x5+x8)+c1*(x6+x7) :
!   Form x0+c3*(x1+x12)+c6*(x2+x11)+c4*(x3+x10)+c1*(x4+x9)+c2*(x5+x8)+c5*(x6+x7) :
!   Form x0+c4*(x1+x12)+c5*(x2+x11)+c1*(x3+x10)+c3*(x4+x9)+c6*(x5+x8)+c2*(x6+x7) :
!   Form x0+c5*(x1+x12)+c3*(x2+x11)+c2*(x3+x10)+c6*(x4+x9)+c1*(x5+x8)+c4*(x6+x7) :
!   Form x0+c6*(x1+x12)+c1*(x2+x11)+c5*(x3+x10)+c2*(x4+x9)+c4*(x5+x8)+c3*(x6+x7) : 16 FMUL, 82 FADD (2 real length-6 cyclics, plus 2 adds for x0+...)
!
!   Form    s1*(x1-x12)+s2*(x2-x11)+s3*(x3-x10)+s4*(x4-x9)+s5*(x5-x8)+s6*(x6-x7) :
!   Form    s2*(x1-x12)+s4*(x2-x11)+s6*(x3-x10)-s5*(x4-x9)-s3*(x5-x8)-s1*(x6-x7) :
!   Form    s3*(x1-x12)+s6*(x2-x11)-s4*(x3-x10)-s1*(x4-x9)+s2*(x5-x8)+s5*(x6-x7) :
!   Form    s4*(x1-x12)-s5*(x2-x11)-s1*(x3-x10)+s3*(x4-x9)-s6*(x5-x8)-s2*(x6-x7) :
!   Form    s5*(x1-x12)-s3*(x2-x11)+s2*(x3-x10)-s6*(x4-x9)-s1*(x5-x8)+s4*(x6-x7) :
!   Form    s6*(x1-x12)-s1*(x2-x11)+s5*(x3-x10)-s2*(x4-x9)+s4*(x5-x8)-s3*(x6-x7) : 30 FMUL, 82 FADD (2 real length-6 acyclics)
!   Form X1,2,3,4,5,6,7,8,9,10,11,12                                             :  0 FMUL, 24 FADD
!
!   Totals :                                                                       46 FMUL, 214 FADD (62 FMUL, 190 FADD for LO_ADD version)
!                                                                     compared to 144 FMUL, 192 FADD for naive scheme,
!                                                                        and just  16 FMUL,  98 FADD for radix-12,
!                                                                        and just  32 FMUL, 160 FADD for radix-14.

UPSHOT: Definitely better than naive radix-13; cut multiply count by more than a factor of 3, with ~10% more adds.
        Still relatively less efficient than radix-12, but not terribly worse than the next-higher composite radix,
        R = 14 (1.55x as many muls per point, 1.43x as many adds per point.)
*/

#define DC1 ((double)-1.66277335223429382656331008444690382)
#define DC2 ((double) 0.73124599097534822519618254560377760)
#define DC3 ((double) 1.0070740657275332544937477077369340)
#define DC4 ((double)-0.30816846519175820059367219184688543)
#define DC5 ((double) 0.81698338691215549726750306085822509)
#define DC6 ((double) 0.22376403323791637458136993583748652)
#define DS1 ((double) 0.57514072947400312136838554745545335)
#define DS2 ((double)-0.17413860115213590500566079492926474)
#define DS3 ((double)-0.33582506518644535421963182119524142)
#define DS4 ((double) 4.5240494294812713569277280991401412E-0002)
#define DS5 ((double) 1.1543953381323634420147226757584967)
#define DS6 ((double)-8.7981928766792081008399945211312915E-0002)
#define DS7 ((double) 0.90655220171271016880349079977456841)
#define DS8 ((double) 1.1971367726043428094538453399784083)
#define DS9 ((double)-0.24784313641965327321123187598392850)
#define DSa ((double)-0.86131170741789745523421351878316690)
#define DSb ((double)-4.2741434471979367439122664219911502E-0002)

/*
The basic tangent-based radix-13 DFT macro is defined in radix13.h::RADIX_13_DFT(...).
Here is an 8-register-optimized version which moreover mimics x86-style destructive (2-input, one overwritten with output)
register arithmetic. This version uses the copies-into-the-__t-temps to allow for an in-place (inputs overlap outputs) calling
convention; if only an out-of-place version is needed (as in the carry routine proper, delete the __t*r = __A*r assignments at
the beginning, and replace the __t*r right-hand-side operands in the opening the yi-terms of the imaginary parts with __A*r.

This version increases the number of multiplies somewhat relative to the original temporary-heavy implemetation (from 68 to 88)
since it uses muls-by-2 to effect the in-place doubling needed by the very frequent radix-2 butterfly operation sequence which
takes inputs x,y and outputs x+y and x-y:

	x -= y
	y *= 2
	y += x

These extra constants are needed for this version:
	DC23 =  DC2/DC3	(Note that DC6/DC4 = -DC2/DC3)
	DC54 =  DC5/DC4
	DC65 =  DC6/DC5
	DS63 = -DS6/DS3
	DS74 = -DS7/DS4
	DS85 = -DS8/DS5
	DS93 = -DS9/DS3
	DSa4 = -DSa/DS4
	DSb5 = -DSb/DS5
Note we do not need the above consts DC2,DC5,DC6,DS6,DS7,DS8,DS9,DSa,DSb, i.e. this versions needs the same number of precomputed constants, 17:
*/
#define DC23 ((double) 0.7261094450357824054685101554)
#define DC54 ((double) -2.651093408937175306253240338)
#define DC65 ((double)  0.273890554964217594531489845)
#define DS63 ((double) -0.2619873794052440703141563891)
#define DS74 ((double)-20.03851230725071170999037075)
#define DS85 ((double) -1.037024954155764805303318302)
#define DS93 ((double) -0.7380126205947559296858436109)
#define DSa4 ((double) 19.03851230725071170999037075)
#define DSb5 ((double)  0.03702495415576480530331830225)

#define RADIX_13_XYZ(\
	__A0r,__A0i,__A1r,__A1i,__A2r,__A2i,__A3r,__A3i,__A4r,__A4i,__A5r,__A5i,__A6r,__A6i,__A7r,__A7i,__A8r,__A8i,__A9r,__A9i,__Aar,__Aai,__Abr,__Abi,__Acr,__Aci,\
	__B0r,__B0i,__B1r,__B1i,__B2r,__B2i,__B3r,__B3i,__B4r,__B4i,__B5r,__B5i,__B6r,__B6i,__B7r,__B7i,__B8r,__B8i,__B9r,__B9i,__Bar,__Bai,__Bbr,__Bbi,__Bcr,__Bci\
)\
{\
	double xr0, xr1, xr2, xr3, xr4, xr5, xr6, xr7;\
	double __t1r,__t2r,__t3r,__t4r,__t5r,__t6r,__t7r,__t8r,__t9r,__tar,__tbr,__tcr;\
		__t1r = __A1r;\
		__t2r = __A2r;\
		__t3r = __A3r;\
		__t4r = __A4r;\
		__t5r = __A5r;\
		__t6r = __A6r;\
		__t7r = __A7r;\
		__t8r = __A8r;\
		__t9r = __A9r;\
		__tar = __Aar;\
		__tbr = __Abr;\
		__tcr = __Acr;\
	/***************/\
	/* REAL PARTS: */\
	/***************/\
	/* xr-terms need 8 registers for each side: */\
		xr7 = __A6r;\
		xr5 = __A5r;\
		xr4 = __A4r;\
		xr6 = __A3r;\
		xr3 = __A2r;\
		xr1 = __A1r;\
		xr7 += __A7r;\
		xr5 += __A8r;\
		xr4 += __A9r;\
		xr6 += __Aar;\
		xr3 += __Abr;\
		xr1 += __Acr;\
		xr1 -= xr5;	xr5 *= 2.0;\
		xr3 -= xr6;	xr6 *= 2.0;\
		xr4 -= xr7;	xr7 *= 2.0;\
		xr5 += xr1;\
		xr6 += xr3;\
		xr7 += xr4;\
		xr2 = xr5;\
		xr2 += xr6;\
		xr2 += xr7;\
		xr0 = __A0r;\
		xr0 += xr2;\
		xr2 *= DC1;\
		__B0r = xr0;\
		xr0 += xr2;\
		xr2 = DC3;\
		xr6 *= xr2;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		__B1r = xr6;\
		__B2r = xr5;\
		__B3r = xr7;\
		xr2 = DC23;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		xr6 *= xr2;\
		xr5 += xr0;\
		xr7 += xr0;\
		xr6 += xr0;\
		xr5 += __B1r;\
		xr7 += __B2r;\
		xr6 += __B3r;\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
		xr2 = DC65;\
		xr0 = DC54;\
		__B1r = xr4;\
		__B2r = xr1;\
		__B3r = xr3;\
		xr4 *= xr2;  /* *= DC65 */\
		xr1 *= xr2;  /* *= DC65 */\
		xr3 *= xr2;  /* *= DC65 */\
		xr4 += __B3r;\
		xr1 -= __B1r;\
		xr3 += __B2r;\
		xr4 *= xr0;  /* *= DC54 */\
		xr1 *= xr0;  /* *= DC54 */\
		xr3 *= xr0;  /* *= DC54 */\
		xr2 = DC4;\
		xr4 += __B2r;\
		xr1 -= __B3r;\
		xr3 -= __B1r;\
		xr4 *= xr2;  /* *= DC4 */\
		xr1 *= xr2;  /* *= DC4 */\
		xr3 *= xr2;  /* *= DC4 */\
	/* Spill into destination outputs: */\
		xr6 -= xr4;	xr4 *= 2.0;\
		xr7 -= xr1;	xr1 *= 2.0;\
		xr5 -= xr3;	xr3 *= 2.0;\
		__B3r = xr7;\
		__B8r = xr5;\
		xr4 += xr6;\
		xr7 += xr1;\
		xr5 += xr3;\
		__B1r = xr5;\
		__B2r = xr7;\
		__B4r = xr6;\
		__B6r = xr4;\
	/* yi-terms: */\
		xr1 = __A1i;\
		xr2 = __A2i;\
		xr5 = __A3i;\
		xr3 = __A4i;\
		xr4 = __A5i;\
		xr6 = __A6i;\
		xr1 -= __Aci;\
		xr2 -= __Abi;\
		xr5 -= __Aai;\
		xr3 -= __A9i;\
		xr4 -= __A8i;\
		xr6 -= __A7i;\
		xr7 = xr1;\
		xr0 = xr2;\
		xr7 -= xr3;\
		xr0 += xr4;\
		xr7 += xr5;\
		xr0 += xr6;\
		xr1 += xr3;\
		xr5 += xr3;\
		xr2 -= xr4;\
		xr6 -= xr4;\
		xr4 = xr0;\
		xr3 = xr7;\
		xr0 *= DS2;\
		xr7 *= DS2;\
		xr4 *= DS1;\
		xr3 *= DS1;\
		xr4 -= xr7;\
		xr3 += xr0;\
		__Bcr = xr4;	/* tmp-store in Bcr */\
		xr0 = xr1;\
		xr4 = xr5;\
		xr0 += xr2;\
		xr4 += xr6;\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
		__Bbr = xr0;	/* tmp-store in Bbr */\
		  xr7 = xr4;\
		xr4 *= DS63;\
		xr0 *= DS93;\
		xr4 += __Bbr;\
		xr0 +=   xr7;\
		xr4 *= DS3;\
		xr0 *= DS3;\
		/*\
		xr7 = xr4+DS4*xr2-DS7*xr6;\
		xr2 = xr0+DS4*xr6-DSa*xr2;\
		*/\
		__Bbr = xr2;\
		  xr7 = xr6;\
		xr6 *= DS74;\
		xr2 *= DSa4;\
		xr6 += __Bbr;\
		xr2 +=   xr7;\
		xr6 *= DS4;\
		xr2 *= DS4;\
		xr6 += xr4;\
		xr2 += xr0;\
		/*\
		xr7 = xr4+DS5*xr1-DS8*xr5;\
		xr0 = xr0+DS5*xr5-DSb*xr1;\
		*/\
		__Bbr = xr1;\
		  xr7 = xr5;\
		xr5 *= DS85;\
		xr1 *= DSb5;\
		xr5 += __Bbr;\
		xr1 +=   xr7;\
		xr5 *= DS5;\
		xr1 *= DS5;\
		xr5 += xr4;\
		xr1 += xr0;\
		xr4 = __Bcr;\
		/*\
		xr7 = xr3+xr6;\
		xr6 = xr6-xr3+xr2;\
		xr2 = xr3+xr2;\
		xr0 = xr4+xr5;\
		xr5 = xr5-xr4+xr1;\
		xr4 = xr4+xr1;\
		*/\
		xr7 = xr3;\
		xr0 = xr4;\
		xr7+= xr6;\
		xr0+= xr5;\
		xr6-= xr3;\
		xr5-= xr4;\
		xr6+= xr2;\
		xr5+= xr1;\
		xr2+= xr3;\
		xr4+= xr1;\
	/* Combine xr and yi-terms to get real parts of outputs: */\
		xr1 = 2.0;\
		xr3 = __B6r;	xr3 -= xr7;	__B6r = xr3;	xr7 *= xr1;	xr7 += xr3;	__B7r = xr7;\
		xr3 = __B8r;	xr3 -= xr6;	__B8r = xr3;	xr6 *= xr1;	xr6 += xr3;	__B5r = xr6;\
		xr3 = __B2r;	xr3 -= xr2;	__B2r = xr3;	xr2 *= xr1;	xr2 += xr3;	__Bbr = xr2;\
		xr3 = __B3r;	xr3 -= xr0;	__B3r = xr3;	xr0 *= xr1;	xr0 += xr3;	__Bar = xr0;\
		xr3 = __B4r;	xr3 -= xr5;	__B4r = xr3;	xr5 *= xr1;	xr5 += xr3;	__B9r = xr5;\
		xr3 = __B1r;	xr3 -= xr4;	__B1r = xr3;	xr4 *= xr1;	xr4 += xr3;	__Bcr = xr4;\
		/*\
		__B6r -= xr7;	xr7 *= xr1;	xr7 += __B6r;	__B7r = xr7;\
		__B8r -= xr6;	xr6 *= xr1;	xr6 += __B8r;	__B5r = xr6;\
		__B2r -= xr2;	xr2 *= xr1;	xr2 += __B2r;	__Bbr = xr2;\
		__B3r -= xr0;	xr0 *= xr1;	xr0 += __B3r;	__Bar = xr0;\
		__B4r -= xr5;	xr5 *= xr1;	xr5 += __B4r;	__B9r = xr5;\
		__B1r -= xr4;	xr4 *= xr1;	xr4 += __B1r;	__Bcr = xr4;\
		*/\
	/***************/\
	/* IMAG PARTS: Swap __**r <--> __**i, Use __t*r (instead of __A*r) as inputs for yr-terms, Replace __B[j] with __B[13-j] for j > 0: */\
	/***************/\
	/* xi-terms need 8 registers for each side: */\
		xr7 = __A6i;\
		xr5 = __A5i;\
		xr4 = __A4i;\
		xr6 = __A3i;\
		xr3 = __A2i;\
		xr1 = __A1i;\
		xr7 += __A7i;\
		xr5 += __A8i;\
		xr4 += __A9i;\
		xr6 += __Aai;\
		xr3 += __Abi;\
		xr1 += __Aci;\
		xr1 -= xr5;	xr5 *= 2.0;\
		xr3 -= xr6;	xr6 *= 2.0;\
		xr4 -= xr7;	xr7 *= 2.0;\
		xr5 += xr1;\
		xr6 += xr3;\
		xr7 += xr4;\
		xr2 = xr5;\
		xr2 += xr6;\
		xr2 += xr7;\
		xr0 = __A0i;\
		xr0 += xr2;\
		xr2 *= DC1;\
		__B0i = xr0;\
		xr0 += xr2;\
		xr2 = DC3;\
		xr6 *= xr2;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		__Bci = xr6;\
		__Bbi = xr5;\
		__Bai = xr7;\
		xr2 = DC23;\
		xr5 *= xr2;\
		xr7 *= xr2;\
		xr6 *= xr2;\
		xr5 += xr0;\
		xr7 += xr0;\
		xr6 += xr0;\
		xr5 += __Bci;\
		xr7 += __Bbi;\
		xr6 += __Bai;\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
		xr2 = DC65;\
		xr0 = DC54;\
		__Bci = xr4;\
		__Bbi = xr1;\
		__Bai = xr3;\
		xr4 *= xr2;  /* *= DC65 */\
		xr1 *= xr2;  /* *= DC65 */\
		xr3 *= xr2;  /* *= DC65 */\
		xr4 += __Bai;\
		xr1 -= __Bci;\
		xr3 += __Bbi;\
		xr4 *= xr0;  /* *= DC54 */\
		xr1 *= xr0;  /* *= DC54 */\
		xr3 *= xr0;  /* *= DC54 */\
		xr2 = DC4;\
		xr4 += __Bbi;\
		xr1 -= __Bai;\
		xr3 -= __Bci;\
		xr4 *= xr2;  /* *= DC4 */\
		xr1 *= xr2;  /* *= DC4 */\
		xr3 *= xr2;  /* *= DC4 */\
	/* Spill into destination outputs: */\
		xr6 -= xr4;	xr4 *= 2.0;\
		xr7 -= xr1;	xr1 *= 2.0;\
		xr5 -= xr3;	xr3 *= 2.0;\
		__Bai = xr7;\
		__B5i = xr5;\
		xr4 += xr6;\
		xr7 += xr1;\
		xr5 += xr3;\
		__Bci = xr5;\
		__Bbi = xr7;\
		__B9i = xr6;\
		__B7i = xr4;\
	/* yr-terms: */\
		xr1 = __t1r;\
		xr2 = __t2r;\
		xr5 = __t3r;\
		xr3 = __t4r;\
		xr4 = __t5r;\
		xr6 = __t6r;\
		xr1 -= __tcr;\
		xr2 -= __tbr;\
		xr5 -= __tar;\
		xr3 -= __t9r;\
		xr4 -= __t8r;\
		xr6 -= __t7r;\
		xr7 = xr1;\
		xr0 = xr2;\
		xr7 -= xr3;\
		xr0 += xr4;\
		xr7 += xr5;\
		xr0 += xr6;\
		xr1 += xr3;\
		xr5 += xr3;\
		xr2 -= xr4;\
		xr6 -= xr4;\
		xr4 = xr0;\
		xr3   = xr7;\
		xr0 *= DS2;\
		xr7 *= DS2;\
		xr4 *= DS1;\
		xr3   *= DS1;\
		xr4 -= xr7;\
		xr3   += xr0;\
		__B1i = xr4;	/* tmp-store in B1i */\
		xr0 = xr1;\
		xr4 = xr5;\
		xr0 += xr2;\
		xr4 += xr6;\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
		__B2i = xr0;	/* tmp-store in B2i */\
		  xr7 = xr4;\
		xr4 *= DS63;\
		xr0 *= DS93;\
		xr4 += __B2i;\
		xr0 +=   xr7;\
		xr4 *= DS3;\
		xr0 *= DS3;\
		/*\
		xr7 = xr4+DS4*xr2-DS7*xr6;\
		xr2 = xr0+DS4*xr6-DSa*xr2;\
		*/\
		__B2i = xr2;\
		  xr7 = xr6;\
		xr6 *= DS74;\
		xr2 *= DSa4;\
		xr6 += __B2i;\
		xr2 +=   xr7;\
		xr6 *= DS4;\
		xr2 *= DS4;\
		xr6 += xr4;\
		xr2 += xr0;\
		/*\
		xr7 = xr4+DS5*xr1-DS8*xr5;\
		xr0 = xr0+DS5*xr5-DSb*xr1;\
		*/\
		__B2i = xr1;\
		  xr7 = xr5;\
		xr5 *= DS85;\
		xr1 *= DSb5;\
		xr5 += __B2i;\
		xr1 +=   xr7;\
		xr5 *= DS5;\
		xr1 *= DS5;\
		xr5 += xr4;\
		xr1 += xr0;\
		xr4 = __B1i;\
		/*\
		xr7 = xr3+xr6;\
		xr6 = xr6-xr3+xr2;\
		xr2 = xr3+xr2;\
		xr0 = xr4+xr5;\
		xr5 = xr5-xr4+xr1;\
		xr4 = xr4+xr1;\
		*/\
		xr7 = xr3;\
		xr0 = xr4;\
		xr7+= xr6;\
		xr0+= xr5;\
		xr6-= xr3;\
		xr5-= xr4;\
		xr6+= xr2;\
		xr5+= xr1;\
		xr2+= xr3;\
		xr4+= xr1;\
	/* Combine xr and yi-terms to get real parts of outputs: */\
		xr1 = 2.0;\
		xr3 = __B7i;	xr3 -= xr7;	__B7i = xr3;	xr7 *= xr1;	xr7 += xr3;	__B6i = xr7;\
		xr3 = __B5i;	xr3 -= xr6;	__B5i = xr3;	xr6 *= xr1;	xr6 += xr3;	__B8i = xr6;\
		xr3 = __Bbi;	xr3 -= xr2;	__Bbi = xr3;	xr2 *= xr1;	xr2 += xr3;	__B2i = xr2;\
		xr3 = __Bai;	xr3 -= xr0;	__Bai = xr3;	xr0 *= xr1;	xr0 += xr3;	__B3i = xr0;\
		xr3 = __B9i;	xr3 -= xr5;	__B9i = xr3;	xr5 *= xr1;	xr5 += xr3;	__B4i = xr5;\
		xr3 = __Bci;	xr3 -= xr4;	__Bci = xr3;	xr4 *= xr1;	xr4 += xr3;	__B1i = xr4;\
/* Totals: 164 FADD, 88 FMUL. */\
}

#endif	/* #ifndef radix13_included */
