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
#ifndef radix17_included
#define radix17_included

// Radix-17 complex DFT using cyclic|acyclic length-8 subconvo for cosine terms, respectively.
#define NAIVE_CONVO	0	// Useful to impl. this first to debug overall DFT scheme, then drop in the optimized convos

#if NAIVE_CONVO
	// Length-8 cyclic convolution macro b = a_dc + (a * x), with a_dc being the DC signal component,
	// and any constant-inputs assumed to be in x, if such are applicable.
	// In our fast Nussbaumer-style version of this routine we propagate a_dc to all outputs with O(1)
	// operations rather than O(n) as here:
	#define cyclic_8( _a_dc,\
		_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,\
		_x0,_x1,_x2,_x3,_x4,_x5,_x6,_x7,\
		_b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7 \
	)\
	{\
		_b0 = _a_dc + _a0*_x0 + _a7*_x1 + _a6*_x2 + _a5*_x3 + _a4*_x4 + _a3*_x5 + _a2*_x6 + _a1*_x7;\
		_b1 = _a_dc + _a1*_x0 + _a0*_x1 + _a7*_x2 + _a6*_x3 + _a5*_x4 + _a4*_x5 + _a3*_x6 + _a2*_x7;\
		_b2 = _a_dc + _a2*_x0 + _a1*_x1 + _a0*_x2 + _a7*_x3 + _a6*_x4 + _a5*_x5 + _a4*_x6 + _a3*_x7;\
		_b3 = _a_dc + _a3*_x0 + _a2*_x1 + _a1*_x2 + _a0*_x3 + _a7*_x4 + _a6*_x5 + _a5*_x6 + _a4*_x7;\
		_b4 = _a_dc + _a4*_x0 + _a3*_x1 + _a2*_x2 + _a1*_x3 + _a0*_x4 + _a7*_x5 + _a6*_x6 + _a5*_x7;\
		_b5 = _a_dc + _a5*_x0 + _a4*_x1 + _a3*_x2 + _a2*_x3 + _a1*_x4 + _a0*_x5 + _a7*_x6 + _a6*_x7;\
		_b6 = _a_dc + _a6*_x0 + _a5*_x1 + _a4*_x2 + _a3*_x3 + _a2*_x4 + _a1*_x5 + _a0*_x6 + _a7*_x7;\
		_b7 = _a_dc + _a7*_x0 + _a6*_x1 + _a5*_x2 + _a4*_x3 + _a3*_x4 + _a2*_x5 + _a1*_x6 + _a0*_x7;\
	}

	// Length-8 acyclic convolution macro b = a * x,
	// with any constant-inputs assumed to be in x, if such are applicable.
  #if 0
	#define acyclic_8(\
		_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,\
		_x0,_x1,_x2,_x3,_x4,_x5,_x6,_x7,\
		_b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7 \
	)\
	{\
		_b0 = _a0*_x0 - _a7*_x1 - _a6*_x2 - _a5*_x3 - _a4*_x4 - _a3*_x5 - _a2*_x6 - _a1*_x7;\
		_b1 = _a1*_x0 + _a0*_x1 - _a7*_x2 - _a6*_x3 - _a5*_x4 - _a4*_x5 - _a3*_x6 - _a2*_x7;\
		_b2 = _a2*_x0 + _a1*_x1 + _a0*_x2 - _a7*_x3 - _a6*_x4 - _a5*_x5 - _a4*_x6 - _a3*_x7;\
		_b3 = _a3*_x0 + _a2*_x1 + _a1*_x2 + _a0*_x3 - _a7*_x4 - _a6*_x5 - _a5*_x6 - _a4*_x7;\
		_b4 = _a4*_x0 + _a3*_x1 + _a2*_x2 + _a1*_x3 + _a0*_x4 - _a7*_x5 - _a6*_x6 - _a5*_x7;\
		_b5 = _a5*_x0 + _a4*_x1 + _a3*_x2 + _a2*_x3 + _a1*_x4 + _a0*_x5 - _a7*_x6 - _a6*_x7;\
		_b6 = _a6*_x0 + _a5*_x1 + _a4*_x2 + _a3*_x3 + _a2*_x4 + _a1*_x5 + _a0*_x6 - _a7*_x7;\
		_b7 = _a7*_x0 + _a6*_x1 + _a5*_x2 + _a4*_x3 + _a3*_x4 + _a2*_x5 + _a1*_x6 + _a0*_x7;\
		/* Totals: [56 ADD, 64 MUL]. */\
	}
  #else	// Try a length-4 complex right-angle-transform variant ... I used this to debug the "Nussbaumerized" version below:
	#define acyclic_8(\
		_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,\
		_x0,_x1,_x2,_x3,_x4,_x5,_x6,_x7,\
		_b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7 \
	)\
	{\
		const double _isrt2 = .70710678118654752440, _c1 = .92387953251128675613, _s1 = .38268343236508977172;\
		double _t;\
		double _a0r,_a0i,_a1r,_a1i,_a2r,_a2i,_a3r,_a3i;\
		double _x0r,_x0i,_x1r,_x1i,_x2r,_x2i,_x3r,_x3i;\
		/* Apply right-angle weights to a0-7 ... the + I* is automatic due to our complex-pairing of
		stride-4-separated real inputs, e.g. (a0 + I.a4) is a no-op: */\
		_a0r = _a0;					_a0i = _a4;\
		_a1r = _c1*_a1 - _s1*_a5;	_a1i = _s1*_a1 + _c1*_a5;\
		_a2r = _isrt2*(_a2 - _a6);	_a2i = _isrt2*(_a2 + _a6);\
		_a3r = _s1*_a3 - _c1*_a7;	_a3i = _c1*_a3 + _s1*_a7;\
		/* Apply right-angle weights to x0-7: */\
		_x0r = _x0;					_x0i = _x4;\
		_x1r = _c1*_x1 - _s1*_x5;	_x1i = _s1*_x1 + _c1*_x5;\
		_x2r = _isrt2*(_x2 - _x6);	_x2i = _isrt2*(_x2 + _x6);\
		_x3r = _s1*_x3 - _c1*_x7;	_x3i = _c1*_x3 + _s1*_x7;\
		/* b0 = a0.x0+a3.x1+a2.x2+a1.x3: */\
		_b0  = _a0r*_x0r-_a0i*_x0i; _b4  = _a0r*_x0i+_a0i*_x0r;\
		_b0 += _a3r*_x1r-_a3i*_x1i; _b4 += _a3r*_x1i+_a3i*_x1r;\
		_b0 += _a2r*_x2r-_a2i*_x2i; _b4 += _a2r*_x2i+_a2i*_x2r;\
		_b0 += _a1r*_x3r-_a1i*_x3i; _b4 += _a1r*_x3i+_a1i*_x3r;\
		/* b1 = a1.x0+a0.x1+a3.x2+a2.x3: */\
		_b1  = _a1r*_x0r-_a1i*_x0i; _b5  = _a1r*_x0i+_a1i*_x0r;\
		_b1 += _a0r*_x1r-_a0i*_x1i; _b5 += _a0r*_x1i+_a0i*_x1r;\
		_b1 += _a3r*_x2r-_a3i*_x2i; _b5 += _a3r*_x2i+_a3i*_x2r;\
		_b1 += _a2r*_x3r-_a2i*_x3i; _b5 += _a2r*_x3i+_a2i*_x3r;\
		/* b2 = a2.x0+a1.x1+a0.x2+a3.x3: */\
		_b2  = _a2r*_x0r-_a2i*_x0i; _b6  = _a2r*_x0i+_a2i*_x0r;\
		_b2 += _a1r*_x1r-_a1i*_x1i; _b6 += _a1r*_x1i+_a1i*_x1r;\
		_b2 += _a0r*_x2r-_a0i*_x2i; _b6 += _a0r*_x2i+_a0i*_x2r;\
		_b2 += _a3r*_x3r-_a3i*_x3i; _b6 += _a3r*_x3i+_a3i*_x3r;\
		/* b3 = a3.x0+a2.x1+a1.x2+a0.x3: */\
		_b3  = _a3r*_x0r-_a3i*_x0i; _b7  = _a3r*_x0i+_a3i*_x0r;\
		_b3 += _a2r*_x1r-_a2i*_x1i; _b7 += _a2r*_x1i+_a2i*_x1r;\
		_b3 += _a1r*_x2r-_a1i*_x2i; _b7 += _a1r*_x2i+_a1i*_x2r;\
		_b3 += _a0r*_x3r-_a0i*_x3i; _b7 += _a0r*_x3i+_a0i*_x3r;\
		/* Undo right-angle weights for last 3 complex outputs: */\
		_t = _b1; _b1 = _c1*_t + _s1*_b5;	_b5 = _c1*_b5 - _s1*_t;\
		_t = _b2; _b2 = _isrt2*(_t + _b6);	_b6 = _isrt2*(_b6 - _t);\
		_t = _b3; _b3 = _s1*_t + _c1*_b7;	_b7 = _s1*_b7 - _c1*_t;\
		/* Totals: [68 ADD, 84 MUL], if we handle the x0-7 weighting-of-trig-consts in preprocessing;
		The right-angle weight/unweight steps account for [12 ADD, 20 MUL] of that total opcount.
		Compare to the [56 ADD, 64 MUL] of the above naive 8-acyclic, and [52 ADD, 40 MUL] of the
		"Nussbaumerized" version of the right-angle scheme. */\
	}
  #endif

#else	// Optimized 8-convos:

	// Nussbaumer's 8-convo with an added DC term, instruction reordering to better reveal
	// butterfly structure and FMA opportunities where they exist, and reduced temp-variable count.
	// Permits in-place (inputs same mem-addresses as outputs):
	#define cyclic_8( _a_dc,\
		_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,\
		_x0,_x1,_x2,_x3,_x4,_x5,_x6,_x7,_x8,_x9,_x10,_x11,_x12,_x13,\
		_b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7 \
	)\
	{\
	double _y0,_y1,_y2,_y3,_y4,_y5,_y6,_y7,_y8,_y9,_y10,_y11,_y12,_y13,_y14,_y15,_y16,_y17,_y18,_y19;\
		_y0  = _a0 + _a4;	_y6  = _a0 - _a4;\
		_y1  = _a1 + _a5;	_y7  = _a1 - _a5;\
		_y2  = _a2 + _a6;	_y8  = _a2 - _a6;\
		_y3  = _a3 + _a7;	_y9  = _a3 - _a7;\
		_y4  = _y0 + _y2;	_y10 = _y0 - _y2;\
		_y5  = _y1 + _y3;	_y11 = _y1 - _y3;\
		_y12 = _y4 + _y5;	_y13 = _y4 - _y5;\
		/* Done with y0-5 at this point; re-use those in the multiply section */\
		_y14 = _y7  + _y9 ;\
		_y15 = _y6  + _y8 ;\
		_y17 = _y8  - _y9 ;\
		_y18 = _y6  - _y7 ;\
		_y19 = _y10 + _y11;\
		_y16 = _y15 - _y14;										/* [20 add, 0 mul] */\
		/* Mpy by the trig consts _x and inject the DC term: */\
		_y12 = _y12 * _x6  + _a_dc;\
		_y13 = _y13 * _x7 ;\
		_y16 = _y16 * _x10;\
		_y17 = _y17 * _x11;\
		_y18 = _y18 * _x12;\
		_y19 = _y19 * _x13;										/* [1 add, 6 mul] or [0 add, 6 fma] */\
		/* Done with y14-19 at this point; re-use those in final-assembly section */\
		/* Begin final assembly, including 3-stage propagation of the DC term to 2/4/8 outputs: */\
		_y14 = _y14 * _x8  + _y16;\
		_y15 = _y15 * _x9  - _y16;\
		_y9  = _y17 + _y9  * _x3 ;\
		_y8  = _y17 - _y8  * _x2 ;\
		_y7  = _y7  * _x1  + _y18;\
		_y6  = _y6  * _x0  - _y18;\
		_y11 = _y19 - _y11 * _x5 ;\
		_y10 = _y19 + _y10 * _x4 ;\
		_y16 = _y13 + _y12;\
		_y12 = _y12 - _y13;										/* [10 add, 8 mul] or [2 add, 8 fma] */\
		/* y7,9 have DC term on entry: */\
		_y0 = _y14 - _y9 ;\
		_y1 = _y11 + _y16;\
		_y2 = _y15 + _y8 ;\
		_y3 = _y10 + _y12;\
		_y4 = _y14 + _y7 ;\
		_y5 = _y16 - _y11;\
		_y6 = _y15 + _y6 ;\
		_y7 = _y12 - _y10;\
		/* Outputs: u11,13,15,17 have DC term coming into this final butterfly: */\
		_b0 = _y1 + _y0;	_b4 = _y1 - _y0;\
		_b1 = _y3 + _y2;	_b5 = _y3 - _y2;\
		_b2 = _y5 + _y4;	_b6 = _y5 - _y4;\
		_b3 = _y7 + _y6;	_b7 = _y7 - _y6;					/* [16 add, 0 mul] */\
		/* Totals: [47 ADD, 14 MUL] or [38 ADD, 14 FMA]. Making use of the obvious FMA oppos cuts 9 ADD. */\
	}

	/* Nussbaumer's length-8 acyclic convolution has a whopping [21 mul, 77 add], so instead use the standard
	trick whereby we multiply the [0,...,n-1]th elts of each of our input vectors by the j = [0,...,n-1]th
	(2n) roots of unity exp(I.j.Pi/n) and then ffed them to a suitably complexified analog of the above 8-convo.
	Said 16th roots of unity, letting exp(I.Pi/8) = c1 + I.s1 and exp(I.Pi/4) = (1 + I)/sqrt(2) := (1 + I)*isrt2,
	are:
		w0 = 1					w1 =  c1 + I.s1
		w2 = (1 + I).isrt2		w3 =  s1 + I.c1
		w4 = I					w5 = -s1 + I.c1
		w6 = (-1 + I).isrt2		w7 = -c1 + I.s1 .
	Now look at the opening 8-convo butterfly sequence modifed with these acyclic-DWT multipliers:
		y0  = a0.w0 + a4.w4 = a0+I.a4 = [a0r-a4i,a0i+a4r], y0,y4 cost [0 mul, 4 add]
		y1  = a1.w1 + a5.w5 = a1.[c1+I.s1]+a5.[-s1+I.c1] = c1.[a1+I.a5]+s1.[I.a1-a5] = c1.[a1r-a5i,a1i+a5r]-s1.[a1i+a5r,-(a1r-a5i)]
			so let a1r5i = a1r-a5i and a1i5r = a1i+a5r and we have y1 = c1.[a1r5i,a1i5r]-s1.[a1i5r,-a1r5i]
		 = [c1.a1r5i-s1.a1i5r,c1.a1i5r+s1.a1r5i] and y5 = [same, but with a1i5r,a1r5i switched: same 4 products assembled differently] .
			Thus y1,y5 cost [4 mul, 4 add], which is half the cost of the naive [a1.w1, a5.w5]+butterfly computation.
		y2  = a2.w2 + a6.w6 = [w2.(1+I) + w6.(-1+I)].isrt2 = [w2r-w2i,w2r+w2i].isrt2 + [-w6r-w6i,w6r-w6i].isrt2; y2,y6 cost [2 mul, 6 add]
		y3  = a3.w3 + a7.w7; y3,y7 cost [4 mul, 4 add]
	So the acyclic-DWT-weighted y0-y7 computation costs [10 mul, 18 add], versus [0 mul, 16 add] for the 2 x real 8-convo version of same.
	The other part of the 8-convo where the complex acyclic-DWT-weighting increases the cost is the middle 14-muls section. The problem is
	that the trig-term combos used by Nussbaumer 8-convo obscure the simple complex-roots symmetries we took advantage of above, so applying
	the same acyclic-DWT multipliers to the trig terms means the ensuing 14-mul sequence will now cost 14 full cmul = [56 mul, 28 add]. Ouch!
	Thus this approach incurs the [28 mul, 92 add] cost of the 2 basic real-data 8-convos, plus a complex-mul penalty of [52 mul, 30 add],
	for a total opcount of [80 mul, 122 add]. Bad as that is, it's comparable to Nussbaumer's 8-acyclic, 2x of which needs [42 mul, 154 add].
	Compare to our naive 8-acyclic approach, 2x of which costs [128 mul, 112 add] - save 48 mul at a cost of 10 more add, but wreck the
	FMA-niceness of the naive approach with its [128 fma (16 mul), 0 add]. But we can do better...

	Q: What about instead using a pair of complex length-4 right-angle convolutions to effect the 8-acyclic convos of the Re,Im parts?
	Recall our 8-convo-effecting DWT weigts are
		w0 = 1					w1 =  c1 + I.s1
		w2 = (1 + I).isrt2		w3 =  s1 + I.c1
		w4 = I					w5 = -s1 + I.c1
		w6 = (-1 + I).isrt2		w7 = -c1 + I.s1 .
	That would mean a *complex* 4-convo with inputs
		A = [w0.(a0 + I.a4), w1.(a1 + I.a5), w2.(a2 + I.a6), w3.(a3 + I.a7)]
		X = [w0.(x0 + I.x4), w1.(x1 + I.x5), w2.(x2 + I.x6), w3.(x3 + I.x7)]
	Plugging in w0-7 and expressing in terms of the reals a0-7, x0-7, c1,s1,isrt2:
		A = [(a0 + I.a4), (c1 + I.s1).(a1 + I.a5), (1 + I).isrt2.(a2 + I.a6), (s1 + I.c1).(a3 + I.a7)]
		  = [(a0 + I.a4), (c1.a1-s1.a5) + I.(s1.a1+c1.a5), isrt2.(a2-a6) + I.isrt2.(a2+a6), (s1.a3-c1.a7) + I.(c1.a3+s1.a7)]
		X = [(x0 + I.x4), (c1.x1-s1.x5) + I.(s1.x1+c1.x5), isrt2.(x2-x6) + I.isrt2.(x2+x6), (s1.x3-c1.x7) + I.(c1.x3+s1.x7)] ,
		thus need [10 mul, 6 add] to apply DWT weights to each input-vec, the A-weighting is the only one which adds to our
		opcount since the X-terms assumed constant and thus precomputable.
	The ensuing Nussbaumer-style 4-convo costs [5 cmul, 15 cadd] = [20 mul, 40 add], but now must add the
	[20 mul, 12 add] cost of DWT weighting the A-inputs and unweighting the B-outputs.
	Total = [52 ADD, 40 MUL] or [30 ADD, 42 FMA], 5 more ADD than Nussbaumer 8-cyclic, and a whopping 26 more MUL.
	NOTE: As an alternative to the right-angle scheme can do an 8-cyclic followed by 16 FMA to effect sign-flips ==> [38 ADD, 30 FMA]
			That may be our best cycle-count-optimized option on hardware with similar ADD and FMA throughput.
	Still a better add-and-total-opcount and much more FMA-izable than the [77 ADD, 21 MUL] of Nussbaumer's pure-real 8-acyclic!
	So let's compute the associated trig-term constants:

		pi = 4*a(1); t = 2*pi/17; [then compute sj = s(j*t) for j = 1,...,8]
	For 8-acyclic, sine terms ordered as:
		h0 = s4;
		h1 = s3;
		h2 =-s2;
		h3 = s7;
		h4 = s1;
		h5 = s5;
		h6 = s8;
		h7 = s6;
	Applying acyclic-DWT weights:
		isrt2 = 1/sqrt(2); cc = c(pi/8); ss = s(pi/8)
		r = h1; h1 = cc*r - ss*h5;	h5 = ss*r + cc*h5;
		r = h2; h2 = isrt2*(r - h6);	h6 = isrt2*(r + h6);
		r = h3; h3 = ss*r - cc*h7;	h7 = cc*r + ss*h7;
	Leaves the following [re,im] pairs:
		h0 =  0.99573417629503452186	h4 =  0.36124166618715294873
		h1 =  0.45894830467224528388	h5 =  1.23117518643485793262
		h2 = -0.60630528816617163954	h6 = -0.34644422799046121289
		h3 = -0.53581491587833643377	h7 =  0.79174787216011104930
	Complexify:
		h0r = h0; h0i = h4;
		h1r = h1; h1i = h5;
		h2r = h2; h2i = h6;
		h3r = h3; h3i = h7;
	Then compute complex versions of Nussbaumer's 4-convo trig consts:
		x0r = h0r + h2r
		x1r = h1r + h3r
		x2r = (x0r + x1r)/4
		x3r = (x0r - x1r)/4
		x4r = ((h0r - h2r) - (h1r - h3r))/2
		x5r = ((h0r - h2r) + (h1r - h3r))/2
		x6r = (h0r - h2r)/2, and analogously for the imaginary parts.
	Gives
		x0r = 0.38942888812886288232;	x0i = 0.01479743819669173584;
		x1r = -.07686661120609114989;	x1i = 2.02292305859496898192;
		x2r = 0.07814056923069293310;	x2i = 0.50943012419791517944;
		x3r = 0.11657387483373850805;	x3i = -.50203140509956931152;
		x4r = 0.30363812195531222187;	x4i = 0.13412928995143363915;
		x5r = 1.29840134250589393952;	x5i = 0.57355660422618052247;
		x6r = 0.80101973223060308070;	x6i = 0.35384294708880708081;
	only the 5 pairs [x2-x5] are actually sent to the complex Nussbaumer 4-convo macro.
	*/
	// Assumes right-angle weights have been applied to trig terms x0-7 on entry:
	#define acyclic_8(\
		_a0,_a1,_a2,_a3,_a4,_a5,_a6,_a7,\
		_x0r,_x0i,_x1r,_x1i,_x2r,_x2i,_x3r,_x3i,_x4r,_x4i, _isrt2,_c1,_s1,\
		_b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7 \
	)\
	{\
		double _t;\
		double _a0r,_a0i,_a1r,_a1i,_a2r,_a2i,_a3r,_a3i;\
		double _y0r,_y0i,_y1r,_y1i,_y2r,_y2i,_y3r,_y3i,_y4r,_y4i,_y5r,_y5i,_y6r,_y6i;\
		double _m0r,_m0i,_m1r,_m1i,_m2r,_m2i,_m3r,_m3i,_m4r,_m4i;\
		double _u0r,_u0i,_u1r,_u1i,_u2r,_u2i,_u3r,_u3i;\
		/* Apply right-angle weights to a0-7 ... the + I* is automatic due to our complex-pairing of
		stride-4-separated real inputs, e.g. (a0 + I.a4) is a no-op: */\
		_a0r = _a0;					_a0i = _a4;\
		_a1r = _c1*_a1 - _s1*_a5;	_a1i = _s1*_a1 + _c1*_a5;\
		_a2r = _isrt2*(_a2 - _a6);	_a2i = _isrt2*(_a2 + _a6);\
		_a3r = _s1*_a3 - _c1*_a7;	_a3i = _c1*_a3 + _s1*_a7;					/* [6 add, 10 mul] or [0 add, 11 fma] */\
		/* Opening butterflies: */\
		_y0r = _a0r + _a2r;								_y0i = _a0i + _a2i;\
		_y1r = _a1r + _a3r;								_y1i = _a1i + _a3i;\
		_y4r = _a0r - _a2r;								_y4i = _a0i - _a2i;\
		_y5r = _a1r - _a3r;								_y5i = _a1i - _a3i;\
		_y2r = _y0r + _y1r;								_y2i = _y0i + _y1i;\
		_y3r = _y0r - _y1r;								_y3i = _y0i - _y1i;\
		_y6r = _y4r + _y5r;								_y6i = _y4i + _y5i;		/* [14 add, 0 mul] */\
		/* Multiply sequence is now complex due to right-angle transform: */\
		_m0r = _y2r * _x0r - _y2i * _x0i;	_m0i = _y2i * _x0r + _y2r * _x0i;\
		_m1r = _y3r * _x1r - _y3i * _x1i;	_m1i = _y3i * _x1r + _y3r * _x1i;\
		_m2r = _y4r * _x2r - _y4i * _x2i;	_m2i = _y4i * _x2r + _y4r * _x2i;\
		_m3r = _y5r * _x3r - _y5i * _x3i;	_m3i = _y5i * _x3r + _y5r * _x3i;\
		_m4r = _y6r * _x4r - _y6i * _x4i;	_m4i = _y6i * _x4r + _y6r * _x4i;	/* [10 add, 20 mul] or [0 add, 20 fma] */\
		/* Closing 2 levels of butterflies: */\
		_u0r = _m0r + _m1r;								_u0i = _m0i + _m1i;\
		_u1r = _m0r - _m1r;								_u1i = _m0i - _m1i;\
		_u2r = _m4r - _m3r;								_u2i = _m4i - _m3i;\
		_u3r = _m4r - _m2r;								_u3i = _m4i - _m2i;\
		/* Outputs: */\
		_b0 = _u0r + _u2r;								_b4 = _u0i + _u2i;\
		_b1 = _u1r + _u3r;								_b5 = _u1i + _u3i;\
		_b2 = _u0r - _u2r;								_b6 = _u0i - _u2i;\
		_b3 = _u1r - _u3r;								_b7 = _u1i - _u3i;		/* [16 add, 0 mul] */\
		/* Undo right-angle weights for last 3 complex outputs: */\
		_t = _b1; _b1 = _c1*_t + _s1*_b5;	_b5 = _c1*_b5 - _s1*_t;\
		_t = _b2; _b2 = _isrt2*(_t + _b6);	_b6 = _isrt2*(_b6 - _t);\
		_t = _b3; _b3 = _s1*_t + _c1*_b7;	_b7 = _s1*_b7 - _c1*_t;				/* [6 add, 10 mul] or [0 add, 11 fma] */\
		/* Totals: [52 ADD, 40 MUL] or [30 ADD, 42 FMA]; the right-angle weight/unweight ==> [12 ADD, 20 MUL] or [22 FMA] of that. (Ugh!)
		Compare to [77 ADD, 21 MUL] of the Nussbaumer pure-real 8-acyclic and [46 ADD, 14 MUL] of the 8-cyclic. */\
	}

#endif

void RADIX_17_DFT(
	double A0r,double A0i,double A1r,double A1i,double A2r,double A2i,double A3r,double A3i,double A4r,double A4i,double A5r,double A5i,double A6r,double A6i,double A7r,double A7i,double A8r,double A8i,double A9r,double A9i,double Aar,double Aai,double Abr,double Abi,double Acr,double Aci,double Adr,double Adi,double Aer,double Aei,double Afr,double Afi,double Agr,double Agi,
	double*B0r,double*B0i,double*B1r,double*B1i,double*B2r,double*B2i,double*B3r,double*B3i,double*B4r,double*B4i,double*B5r,double*B5i,double*B6r,double*B6i,double*B7r,double*B7i,double*B8r,double*B8i,double*B9r,double*B9i,double*Bar,double*Bai,double*Bbr,double*Bbi,double*Bcr,double*Bci,double*Bdr,double*Bdi,double*Ber,double*Bei,double*Bfr,double*Bfi,double*Bgr,double*Bgi
)
{
#if NAIVE_CONVO
	/* Trigonometric consts cJ = cos(2*J*pi/17), sJ = sin(2*J*pi/17)): */
	const double
	c1 =  .93247222940435580457, s1 =  .36124166618715294873,
	c2 =  .73900891722065911592, s2 =  .67369564364655721170,
	c3 =  .44573835577653826741, s3 =  .89516329135506232206,
	c4 =  .09226835946330199525, s4 =  .99573417629503452187,
	c5 = -.27366299007208286351, s5 =  .96182564317281907041,
	c6 = -.60263463637925638916, s6 =  .79801722728023950334,
	c7 = -.85021713572961415212, s7 =  .52643216287735580026,
	c8 = -.98297309968390177827, s8 =  .18374951781657033160;
#else
	#include "radix16.h"
	const double isrt2 = ISRT2, cc = c16, ss = s16;
	const double
	/* 20 consts needed by Nussbaumer 8-convo in terms of our permuted cosine terms, h0-7 = c4,c3,c2,c7,c1,c5,c8,c6,
	then form temps a0 = h0+h4, a1 = h1+h5, a2 = h2+h6, a3 = h3+h7, a4 = a0+a2, a5 = a1+a3, 14 terms needed for the convo are: */
	a6  =  0.79760102082331790481,	/* (((h2-h6)-(h0-h4)) - ((h1-h5)-(h3-h7)))/2 */
	a7  =  1.51700236667193903571,	/* (((h2-h6)-(h0-h4)) + ((h1-h5)+(h3-h7)))/2 */
	a8  =  0.67679849673088522641,	/* (((h2-h6)+(h0-h4)) + ((h1-h5)+(h3-h7)))/2 */
	a9  =  0.92438099608124298938,	/* (((h2-h6)+(h0-h4)) + ((h1-h5)-(h3-h7)))/2 */
	a10 =  0.08905559162060637074,	/* ((a1-a3)-(a0-a2))/4 */
	a11 =  0.72340797728605660183,	/* ((a1-a3)+(a0-a2))/4 */
	a12 = -0.06249999999999999997,	/* (a4+a5)/8 */
	a13 =  0.25769410160110378435,	/* (a4-a5)/8 */
	a14 = -0.29631068529534802316,	/* ((h0-h4) - (h3-h7))/2 */
	a15 = -0.06040126204621633920,	/* ((h0-h4) + (h1-h5))/2 */
	a16 = -0.42010193497052690465,	/* (h0-h4)/2 */
	a17 =  0.44088907348175354244,	/* ((h2-h6)+(h0-h4))/2 */
	a18 =  1.28109294342280735174,	/* ((h2-h6)-(h0-h4))/2 */
	a19 =  0.31717619283272511554,	/* (a0-a2)/4 */
	//b0r = 0.38942888812886288232,	b0i = 0.01479743819669173584,
	//b1r = -.07686661120609114989,	b1i = 2.02292305859496898192,
	b2r = 0.07814056923069293310,	b2i = 0.50943012419791517944,
	b3r = 0.11657387483373850805,	b3i = -.50203140509956931152,
	b4r = 0.30363812195531222187,	b4i = 0.13412928995143363915,
	b5r = 1.29840134250589393952,	b5i = 0.57355660422618052247,
	b6r = 0.80101973223060308070,	b6i = 0.35384294708880708081;
#endif
	double y0r,y1r,y2r,y3r,y4r,y5r,y6r,y7r,y0i,y1i,y2i,y3i,y4i,y5i,y6i,y7i;
	double C1r,C2r,C3r,C4r,C5r,C6r,C7r,C8r,C1i,C2i,C3i,C4i,C5i,C6i,C7i,C8i;
	double S1r,S2r,S3r,S4r,S5r,S6r,S7r,S8r,S1i,S2i,S3i,S4i,S5i,S6i,S7i,S8i;
	/* Cosine terms are:
		C1 = c1*(x1+xG)+c2*(x2+xF)+c3*(x3+xE)+c4*(x4+xD)+c5*(x5+xC)+c6*(x6+xB)+c7*(x7+xA)+c8*(x8+x9)
		C2 = c2*(x1+xG)+c4*(x2+xF)+c6*(x3+xE)+c8*(x4+xD)+c7*(x5+xC)+c5*(x6+xB)+c3*(x7+xA)+c1*(x8+x9)
		C3 = c3*(x1+xG)+c6*(x2+xF)+c8*(x3+xE)+c5*(x4+xD)+c2*(x5+xC)+c1*(x6+xB)+c4*(x7+xA)+c7*(x8+x9)
		C4 = c4*(x1+xG)+c8*(x2+xF)+c5*(x3+xE)+c1*(x4+xD)+c3*(x5+xC)+c7*(x6+xB)+c6*(x7+xA)+c2*(x8+x9)
		C5 = c5*(x1+xG)+c7*(x2+xF)+c2*(x3+xE)+c3*(x4+xD)+c8*(x5+xC)+c4*(x6+xB)+c1*(x7+xA)+c6*(x8+x9)
		C6 = c6*(x1+xG)+c5*(x2+xF)+c1*(x3+xE)+c7*(x4+xD)+c4*(x5+xC)+c2*(x6+xB)+c8*(x7+xA)+c3*(x8+x9)
		C7 = c7*(x1+xG)+c3*(x2+xF)+c4*(x3+xE)+c6*(x4+xD)+c1*(x5+xC)+c8*(x6+xB)+c2*(x7+xA)+c5*(x8+x9)
		C8 = c8*(x1+xG)+c1*(x2+xF)+c7*(x3+xE)+c2*(x4+xD)+c6*(x5+xC)+c3*(x6+xB)+c5*(x7+xA)+c4*(x8+x9)
	Letting
		b0-7 = c4,c3,c2,c7,c1,c5,c8,c6 and
		y0 = (x6+xB), y1 = (x8+x9), y2 = (x5+xC), y3 = (x1+xG), y4 = (x7+xA), y5 = +(x2+xF), y6 = (x3+xE), y7 = (x4+xD),
	Reordering the S-outputs as shown turns the above into an 8-cyclic in standard form:
		C5 = a0*y0+a7*y1+a6*y2+a5*y3+a4*y4+a3*y5+a2*y6+a1*y7
		C8 = a1*y0+a0*y1+a7*y2+a6*y3+a5*y4+a4*y5+a3*y6+a2*y7
		C6 = a2*y0+a1*y1+a0*y2+a7*y3+a6*y4+a5*y5+a4*y6+a3*y7
		C4 = a3*y0+a2*y1+a1*y2+a0*y3+a7*y4+a6*y5+a5*y6+a4*y7
		C3 = a4*y0+a3*y1+a2*y2+a1*y3+a0*y4+a7*y5+a6*y6+a5*y7
		C2 = a5*y0+a4*y1+a3*y2+a2*y3+a1*y4+a0*y5+a7*y6+a6*y7
		C7 = a6*y0+a5*y1+a4*y2+a3*y3+a2*y4+a1*y5+a0*y6+a7*y7
		C1 = a7*y0+a6*y1+a5*y2+a4*y3+a3*y4+a2*y5+a1*y6+a0*y7
	*/
	y0r = A6r+Abr;
	y1r = A8r+A9r;
	y2r = A5r+Acr;
	y3r = A1r+Agr;
	y4r = A7r+Aar;
	y5r = A2r+Afr;
	y6r = A3r+Aer;
	y7r = A4r+Adr;
	*B0r = A0r + y0r+y1r+y2r+y3r+y4r+y5r+y6r+y7r;
#if NAIVE_CONVO
	cyclic_8( A0r,
		y0r,y1r,y2r,y3r,y4r,y5r,y6r,y7r,
		c4,c3,c2,c7,c1,c5,c8,c6,
		C5r,C8r,C6r,C4r,C3r,C2r,C7r,C1r
	)
#else
	cyclic_8( A0r,
		y0r,y1r,y2r,y3r,y4r,y5r,y6r,y7r,
		a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,
		C5r,C8r,C6r,C4r,C3r,C2r,C7r,C1r
	)
#endif
	y0i = A6i+Abi;
	y1i = A8i+A9i;
	y2i = A5i+Aci;
	y3i = A1i+Agi;
	y4i = A7i+Aai;
	y5i = A2i+Afi;
	y6i = A3i+Aei;
	y7i = A4i+Adi;
	*B0i = A0i + y0i+y1i+y2i+y3i+y4i+y5i+y6i+y7i;
#if NAIVE_CONVO
	cyclic_8( A0i,
		y0i,y1i,y2i,y3i,y4i,y5i,y6i,y7i,
		c4,c3,c2,c7,c1,c5,c8,c6,
		C5i,C8i,C6i,C4i,C3i,C2i,C7i,C1i
	)
#else
	cyclic_8( A0i,
		y0i,y1i,y2i,y3i,y4i,y5i,y6i,y7i,
		a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,
		C5i,C8i,C6i,C4i,C3i,C2i,C7i,C1i
	)
#endif
	/*
	Sine terms are:
		S1 = s1*(x1-xG)+s2*(x2-xF)+s3*(x3-xE)+s4*(x4-xD)+s5*(x5-xC)+s6*(x6-xB)+s7*(x7-xA)+s8*(x8-x9)
		S2 = s2*(x1-xG)+s4*(x2-xF)+s6*(x3-xE)+s8*(x4-xD)-s7*(x5-xC)-s5*(x6-xB)-s3*(x7-xA)-s1*(x8-x9)
		S3 = s3*(x1-xG)+s6*(x2-xF)-s8*(x3-xE)-s5*(x4-xD)-s2*(x5-xC)+s1*(x6-xB)+s4*(x7-xA)+s7*(x8-x9)
		S4 = s4*(x1-xG)+s8*(x2-xF)-s5*(x3-xE)-s1*(x4-xD)+s3*(x5-xC)+s7*(x6-xB)-s6*(x7-xA)-s2*(x8-x9)
		S5 = s5*(x1-xG)-s7*(x2-xF)-s2*(x3-xE)+s3*(x4-xD)+s8*(x5-xC)-s4*(x6-xB)+s1*(x7-xA)+s6*(x8-x9)
		S6 = s6*(x1-xG)-s5*(x2-xF)+s1*(x3-xE)+s7*(x4-xD)-s4*(x5-xC)+s2*(x6-xB)+s8*(x7-xA)-s3*(x8-x9)
		S7 = s7*(x1-xG)-s3*(x2-xF)+s4*(x3-xE)-s6*(x4-xD)+s1*(x5-xC)+s8*(x6-xB)-s2*(x7-xA)+s5*(x8-x9)
		S8 = s8*(x1-xG)-s1*(x2-xF)+s7*(x3-xE)-s2*(x4-xD)+s6*(x5-xC)-s3*(x6-xB)+s5*(x7-xA)-s4*(x8-x9)
	See the N == 17 comments in dft_sine_term_opt.c.txt (we append 'txt' in Mlucas releases to prevent this main()-containing
	C source file from being compiled in release-build mode) for details of the reindexing/relabeling magic here, but defining:

	b0-7 = s4,s3,-s2,s7,s1,s5,s8,s6 and
	z0 = (x6-xB), z1 = (x8-x9), z2 = (x5-xC), z3 = (x1-xG), z4 = (x7-xA), z5 = -(x2-xF), z6 = (x3-xE), z7 = (x4-xD),
	flipping signs on and reordering the S-outputs as shown turns the above into an 8-acyclic in standard form:

	-	S5 = +b0*z0-b7*z1-b6*z2-b5*z3-b4*z4-b3*z5-b2*z6-b1*z7
	-	S8 = +b1*z0+b0*z1-b7*z2-b6*z3-b5*z4-b4*z5-b3*z6-b2*z7
	-	S6 = +b2*z0+b1*z1+b0*z2-b7*z3-b6*z4-b5*z5-b4*z6-b3*z7
		S4 = +b3*z0+b2*z1+b1*z2+b0*z3-b7*z4-b6*z5-b5*z6-b4*z7
		S3 = +b4*z0+b3*z1+b2*z2+b1*z3+b0*z4-b7*z5-b6*z6-b5*z7
	-	S2 = +b5*z0+b4*z1+b3*z2+b2*z3+b1*z4+b0*z5-b7*z6-b6*z7
		S7 = +b6*z0+b5*z1+b4*z2+b3*z3+b2*z4+b1*z5+b0*z6-b7*z7
		S1 = +b7*z0+b6*z1+b5*z2+b4*z3+b3*z4+b2*z5+b1*z6+b0*z7

	Alternatively, we can negating the first 3 rows above - here with resulting #terms-with-minus-sign tabulated in rcol:
		S5 = -b0*z0+b7*z1+b6*z2+b5*z3+b4*z4+b3*z5+b2*z6+b1*z7	1
		S8 = -b1*z0-b0*z1+b7*z2+b6*z3+b5*z4+b4*z5+b3*z6+b2*z7	1
		S6 = -b2*z0-b1*z1-b0*z2+b7*z3+b6*z4+b5*z5+b4*z6+b3*z7	3
		S4 = +b3*z0+b2*z1+b1*z2+b0*z3-b7*z4-b6*z5-b5*z6-b4*z7	4
		S3 = +b4*z0+b3*z1+b2*z2+b1*z3+b0*z4-b7*z5-b6*z6-b5*z7	3
		S2 = +b5*z0+b4*z1+b3*z2+b2*z3+b1*z4+b0*z5-b7*z6-b6*z7	2
		S7 = +b6*z0+b5*z1+b4*z2+b3*z3+b2*z4+b1*z5+b0*z6-b7*z7	1
		S1 = +b7*z0+b6*z1+b5*z2+b4*z3+b3*z4+b2*z5+b1*z6+b0*z7	0
	We then do a *cyclic* convo, which ignores all the minus signs, and follow that with 16 += 2.b*.z* adds
	to sign-flip the 16 terms with - signs in the above [0,1,1,2,2,3,3,4]-zeros-per-row formulation:
		S5 = -b0*z0
		S8 = -b1*z0-b0*z1
		S6 = -b2*z0-b1*z1-b0*z2
		S4 = -b7*z4-b6*z5-b5*z6-b4*z7
		S3 = -b7*z5-b6*z6-b5*z7
		S2 = -b7*z6-b6*z7
		S7 = -b7*z7
	*/
	y0r = A6r-Abr;
	y1r = A8r-A9r;
	y2r = A5r-Acr;
	y3r = A1r-Agr;
	y4r = A7r-Aar;
	y5r =-A2r+Afr;
	y6r = A3r-Aer;
	y7r = A4r-Adr;
#if NAIVE_CONVO
	acyclic_8(
		y0r,y1r,y2r,y3r,y4r,y5r,y6r,y7r,
		s4,s3,-s2,s7,s1,s5,s8,s6,
		S5r,S8r,S6r,S4r,S3r,S2r,S7r,S1r
	)
#else
	acyclic_8(
		y0r,y1r,y2r,y3r,y4r,y5r,y6r,y7r,
		b2r,b2i,b3r,b3i,b4r,b4i,b5r,b5i,b6r,b6i, isrt2,cc,ss,
		S5r,S8r,S6r,S4r,S3r,S2r,S7r,S1r
	)
#endif
printf("S1-8r = %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f\n",S1r,S2r,S3r,S4r,S5r,S6r,S7r,S8r);
	y0i = A6i-Abi;
	y1i = A8i-A9i;
	y2i = A5i-Aci;
	y3i = A1i-Agi;
	y4i = A7i-Aai;
	y5i =-A2i+Afi;
	y6i = A3i-Aei;
	y7i = A4i-Adi;
#if NAIVE_CONVO
	acyclic_8(
		y0i,y1i,y2i,y3i,y4i,y5i,y6i,y7i,
		s4,s3,-s2,s7,s1,s5,s8,s6,
		S5i,S8i,S6i,S4i,S3i,S2i,S7i,S1i
	)
#else
	acyclic_8(
		y0i,y1i,y2i,y3i,y4i,y5i,y6i,y7i,
		b2r,b2i,b3r,b3i,b4r,b4i,b5r,b5i,b6r,b6i, isrt2,cc,ss,
		S5i,S8i,S6i,S4i,S3i,S2i,S7i,S1i
	)
#endif
printf("S1-8i = %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f\n",S1i,S2i,S3i,S4i,S5i,S6i,S7i,S8i);
	// Output-terms butterfly. Recall S2,5,6,8 are negated:
	// Bk = Ck + I*Sk, k = 1,...,8:
	*B1r = C1r - S1i;
	*B2r = C2r + S2i;
	*B3r = C3r - S3i;
	*B4r = C4r - S4i;
	*B5r = C5r + S5i;
	*B6r = C6r + S6i;
	*B7r = C7r - S7i;
	*B8r = C8r + S8i;
	*B1i = C1i + S1r;
	*B2i = C2i - S2r;
	*B3i = C3i + S3r;
	*B4i = C4i + S4r;
	*B5i = C5i - S5r;
	*B6i = C6i - S6r;
	*B7i = C7i + S7r;
	*B8i = C8i - S8r;
	// B[17-k] = Ck - I*Sk, k = 8,...,1:
	*B9r = C8r - S8i;
	*Bar = C7r + S7i;
	*Bbr = C6r - S6i;
	*Bcr = C5r - S5i;
	*Bdr = C4r + S4i;
	*Ber = C3r + S3i;
	*Bfr = C2r - S2i;
	*Bgr = C1r + S1i;
	*B9i = C8i + S8r;
	*Bai = C7i - S7r;
	*Bbi = C6i + S6r;
	*Bci = C5i + S5r;
	*Bdi = C4i - S4r;
	*Bei = C3i - S3r;
	*Bfi = C2i + S2r;
	*Bgi = C1i - S1r;
/* Total Cost = [5*(N-1) ADD] + 2*[convos cost] .
	Naive (quadratic-convos) scheme: cyclic_8 = [64 ADD, 56 MUL], acyclic_8 = [56 ADD, 56 MUL] ==> Total = [320 ADD, 224 MUL].
	For the optimized-convos scheme: cyclic_8 = [47 ADD, 14 MUL], acyclic_8 = [52 ADD, 40 MUL] ==> Total = [278 ADD, 108 MUL],
	compared to the [274 ADD, 82 MUL] opcount for the Selesnick/Burrus 17-DFT.
	For FMA-optimized-convos scheme: cyclic_8 = [38 ADD, 14 FMA], acyclic_8 = [38 ADD, 30 FMA] ==> Total = [232 ADD, 88 FMA].
	Compare to FMAized 1-DFT which needs [120 ADD, 78 FMA], n*log(n) extrapolation to N=17 gives [173 ADD, 113 FMA], so our
	17-DFT has worse ADD count and better MUL count than predicted by simple extrapolation.
*/
}

#endif	/* #ifndef radix17_included */

#if 0
// bc code for Nussbaumer 8-convo consts:
pi = 4*a(1)
t = 2*pi/17
h0 = c(4*t)
h1 = c(3*t)
h2 = c(2*t)
h3 = c(7*t)
h4 = c(1*t)
h5 = c(5*t)
h6 = c(8*t)
h7 = c(6*t)
a0  = h0+h4
a1  = h1+h5
a2  = h2+h6
a3  = h3+h7
a4  = a0+a2
a5  = a1+a3
a6  = (((h2-h6)-(h0-h4)) - ((h1-h5)-(h3-h7)))/2
a7  = (((h2-h6)-(h0-h4)) + ((h1-h5)+(h3-h7)))/2
a8  = (((h2-h6)+(h0-h4)) + ((h1-h5)+(h3-h7)))/2
a9  = (((h2-h6)+(h0-h4)) + ((h1-h5)-(h3-h7)))/2
a10 = ((a1-a3)-(a0-a2))/4
a11 = ((a1-a3)+(a0-a2))/4
a12 = (a4+a5)/8
a13 = (a4-a5)/8
a14 = ((h0-h4) - (h3-h7))/2
a15 = ((h0-h4) + (h1-h5))/2
a16 = (h0-h4)/2
a17 = ((h2-h6)+(h0-h4))/2
a18 = ((h2-h6)-(h0-h4))/2
a19 = (a0-a2)/4
#endif
