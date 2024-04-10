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

#if 0
	[Last updated Mar 2022]
	This is a C-ified version of my 2003-dated F90-based FGT-exploration code.
	Most of the comments about op- and cycle counts reflect that dating, at which
	time the leding candidate architectures for such arithmetic were Alpha and MIPS,
	due to their decent (and improving, as embodied by e.g. tha third-gen Alpha 21264)
	64 x 64 ==> 128-bit integer multiply, a.k.a. MUL_LOHI64, support. (Alpha via the
	separate MULQ, UMULH instruction pair; MIPS with its single DMULTU instruction.)

	As of this date, x86_64 has taken over, with recent iterations having excellent
	MUL_LOHI64 support via its single MUL instruction, which appears as 'mulq' in
	GCC-style inline assembler, not to be confused with the aforementoned Alpha MULQ,
	returning the low half of the full 128 bit product.

	Thus x86_64 MUL is like a fully pipelined MIPS DMULTU with 2 key differences:

	1. MUL uses the legacy x86 2-operand MUL syntax, but note the resulting register
	constraints (i.e. need to move result out of the A:D register pair prior to the
	next MUL) is mitigated in terms of performance impact by hardware-level register
	aliasing, i.e. the user-required 'move result out of A:D' typically ends up being
	a no-op at the microcode-execution level.

	2. Haswell-and-beyond iterations of x86_64 add a 3-operand MULX instruction which
	further eases the register-shuffling burden on the assembly programmer.

	And the x86_64 ADD/SUB both set a carry flag, which can be used to speed mutiowrd
	arithmetic, e.g. via an ADD/ADC instruction pair in the 2-word case relevant to mod-q
	arithmetic here.

	The biggest problem with M61-mod FGT on x86_64 is lack of SIMD 128-bit integer product
	support, relative to the excellent floating-point SIMD MUL+FMA support. The latest-gen
	avx-512 SIMD offers the following 64-bit SIMD integer multiply instructions:

	o VPMUL[U]DQ: Multiply packed [un]signed doubleword integers in even-indexed 32-bit
	dwords of zmm2 by packed [un]signed doubleword integers in even-indexed 32-bit dwords
	of zmm3/m512/m64bcst, and store the 8 quadword results in zmm1 using writemask k1.

	o VPMULL[D|Q]: Multiply the packed [dword|qword] signed integers in zmm2 and zmm3/m512/
	m32bcst and store the low [32|64] bits of each product in zmm1 under writemask k1.
	(Equivalent to [16|8]-way SIMD version of IMULL[D|Q].)

	o [Needs AVX512IFMA instruction extensions]
	VPMADD52[L|H]UQ: Multiply unsigned 52-bit integers in zmm2 and zmm3/m128 and add the
	low 52 bits of the 104-bit product to the qword unsigned integer accumulators in zmm1
	using writemask k1.

	Q: Can we use some combination of the above, plus perhaps some bytewise or similar low-
	bitness MUL and/or LUT to synthesize a fast 512-bit SIMD analog of Alpha MULQ+UMULH ?
#endif

#ifdef USE_FGT61

// Prototypes for functions defined below are in util.h, which includes fgt_m61.h
#include "util.h"

#define q 	0x1FFFFFFFFFFFFFFFull	// 2^61 - 1
#define q2	0x3FFFFFFFFFFFFFFEull	// 2*q
#define q4	0x7FFFFFFFFFFFFFFCull	// 4*q
#define bb	0x2000000000000006ull	// q + 7; Used in bounds checking this is the max. possible output of QREDUCE.

/***************/

/* We desire a product of 2 64-bit integers x*y modulo q = 2^61-1.
If we scale x and y (e.g. premultiply one by 8 or premultiply one by 2
and the other by 4) so that the scaled product x~*y~ = 2^64*hi + lo = 8*x*y, then
hi + (lo/8) returns x*y, partially reduced. We assume such an input scaling has been done.

OUTPUT RANGE: Cf. the lengthy comment preceding CMUL_MODQ for details; under the above assumptions
*********** to-do: repeat AGM computation for MUL of form prodq8(2*x, 4*y) ***********
for 0 <= a <= 2^61 + sqrt(2^61),  0 <= b < 2^64, (b divisible by 8),
	        0 <= prodq8(a,b) <= 2q.

COST: one MUL_LOHI64, two IOPs.
*/
uint64 prodq8(const uint64 x, const uint64 y)
{
	uint64 lo,hi;

	// Find hi, lo so 2^64*hi + lo = x*y
#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(x,y, &lo,&hi);
#else
	MUL_LOHI64(x,y,  lo, hi);
#endif
	ASSERT((lo & 7) == 0, "ERROR: product not divisible by 8 in PRODQ8!");
//if(hi + (lo >> 3) > q2) fprintf(stderr, "PRODQ8 inputs: %llu,%llu, outputs: %llu,%llu, result = %llu\n",x,y,lo,hi,hi + (lo >> 3));
//	ASSERT(hi + (lo >> 3) <= q2, "ERROR: result out of range in PRODQ8!");
	return hi + (lo >> 3);		// hi + (lo/8)
}

/***************/

// Returns A*X, where X is in [0,2^61) and A is unsigned <= 8.
// [2015: Replace elaborate case-based impl of original with simple MULQ, which is fast on more or less all 64-bit arches.]
uint64 mul_by_3bit(const uint64 a, const uint64 x)
{
	ASSERT((x >> 61) == 0, "ERROR: x out of range in MUL_BY_3BIT!");
	return a * x;
}

/***************/

#if 1

// Real-inputs multiply. Input bounds:
// A may be as large as 2^63-1 = 4q + 3;
// B may be as large as 2^62-1 = 2q + 1.
// Output bounds: *********************** To-Do! *********************
uint64 rmul_modq(const uint64 a, const uint64 b)
{
	ASSERT(a < 0x8000000000000000ull && b < 0x4000000000000000ull, "Input(s) out of range!");
	return prodq8(a<<1, b<<2);
}

#else	// Hybrid Int/float version of above; strictly experimental,
		// payoff only on arches lacking UMULH and/or where MUL is slow.
/*
Combines a MULQ (64-bit integer multiply) with a floating multiply
to replace UMULH for calculation of the full 128-bit product of two
unsigned integer operands, x*y, and returns the product modulo q = 2^61-1.
The second operand (y) is assumed to be fully reduced modulo q, i.e. to be
in [0,q] (typically this will be the real or imaginary part of a precomputed
modular root of unity.) The first (x) is assumed to be partially reduced,
i.e. to be in [0,B], where B = q+7 is the upper bound on the result of a
modular reduction of a 64-bit unsigned operand u64 via x = (u64 >> 61) + (u64 & q).

Using IEEE 64-bit floating multiply for the high part of the integer product,
we can exactly get a product of at most (53 + 64 - 1) = 116 bits (53 for the
floating high product, 64 for the integer low product, subtract 1 overlap
bit used for error correction of the low bit of the floating product.) This
limits us to 58-bit operands. We can still efficiently obtain a 122-bit product
by splitting the operands:

x = a*2^58 + b,  y = c*2^58 + d,  where a in [0,8], c in [0,7], b,d in [0,2^58).

The full product is thus decomposed as

x*y = a*c*2^116 + (a*d + c*b)*2^58 + b*d,

where b*d has 116 bits and can be got via mixed-floating/integer multiply.
This seems to increase our total multiply count, but since a and c are small,
the products involving them can be obtained via simply shifts and adds, e.g.
using a small select_case(small_multiplicand) structure. It is easy to show
(using the AGM inequality) that a*c*2^58 + (a*d + c*b) < 2^64, so we can
simplify this part into a*y + c*b, and thus have just two multiplies-by-small
to perform. Separately reducing a*y + c*b and b*d modulo q, we obtain a result
in [0, 2q], the same bound as for the all-integer mult_lohi algorithm.

WAIT - must assume x < 2^58 also! Otherwise, can`t convert (x >> 5) exactly
to floating...

	1101001011110010111111010110000110111001010011010001

							    1101000000000101111100101000011111101100001010011111001001010111

	11010010111100101111110101100001101110010100110100101101000000000101111100101000011111101100001010011111001001010111
*/
uint64 rmul_modq(const uint64 x, const uint64 y)
{
	uint64 x58,y58,ay,cb,bd_lo,bd_hi,bd_modq,error_mod4, hi4est;
	const uint64 two58m1 = (1ull<<58)-1;
	double dhi;
	const double scaleinv = (double)0x20000000*0x10000000;	/* (double)2^57 = 2^29 * 2^28 */
	int dbg = FALSE;

	// b*d is here...
	x58 = x & two58m1;
	y58 = y & two58m1;

	bd_lo = x58*y58;	// lo = x58*y58 mod 2^64 = MULQ(x58, y58)

	dhi = (double)x58 * scaleinv * (double)(y58 >> 5);	// Since y precomputed, can eventually eliminate one FMUL and one conversion.

	// Error correction using two overlap bits is here.
	hi4est = (uint64)dhi - 1;
	error_mod4 = ( (bd_lo >> 62) - hi4est ) & 3ull; // Error mod 4
	bd_hi = (hi4est + error_mod4) >> 2;
	ASSERT(bd_lo <= (bd_hi << 3) + bd_lo, "ERROR: overflow of b*d(lo + hi>>3) summand!");
	bd_modq = qreduce((bd_hi << 3) + bd_lo);

	ay = mul_by_3bit((x >> 58),y);
	cb = mul_by_3bit((y >> 58), (x & two58m1));
	ASSERT(cb <= ay+cb, "ERROR: overflow of ay+cb summand!");
	bd_modq = qreduce((bd_hi << 3) + bd_lo);

	// Now form [(a*y + c*b)*2^58 + b*d] mod q.
	return qreduce(mul_pow2_modq( ay+cb, 58 ) + bd_modq);
}
#endif

/***************/

/* Complex multiplication of the Gaussian integer a0 + I*a1 by b0 + I*b1
modulo q = 2^61 - 1, with result xout + I*yout.

All arguments are treated as unsigned.

OPERAND SIZES: a0 and a1 are assumed to result from a modular reduction
of the form (x >> 61) + (x & q). b0 and b1 are assumed to be precomputed
constants (e.g. roots of unity) premultiplied by 8 (cf. below). That is,

	a0 and a1 are in [0, B] := [0, q + 7]  = [0, 2^61 + 6]
	b0 and b1 are in [0,8*q] = [0, 2^64 - 8].

NOTE: The constants b0 and b1 are assumed to have been premultiplied by 8
prior to calling this routine, which simplifies the modular reduction
of a 128-bit product significantly. Without the premultiplication, the
modular reduction of a 128-bit product (hi*2^64 + lo) requires 8 elemental
integer operations (IOPS; defined, using C syntax, as + - >> << & ~ | ^,
plus compares, CMOV and architecture-specific elemental operations such
as Alpha SxADD (scaled add, x=4 or 8), but excluding loads and stores):

	tmp = (hi << 3) + (lo >> 61) + (lo & q);
	return( (tmp >> 61) + (tmp & q) ).

Premultiplying b0 and b1 by 8 reduces the cost to just 5 IOPS:

	tmp = hi + (lo >> 3);				(normalize hi)
	return( (tmp >> 61) + (tmp & q) )	(normalize lo).

Since in general b0 + I*b1 is a precomputed root of unity, the multiply-
by-8 can be done in the precomputation phase, meaning we get it for free.

In some cases we can do just the 128-bit multiply and hi + lo/8 steps,
eschewing the modular reduction. To this end, we define the following
separate procedures for multiply-and-fold-hi-into-lo and modular reduction:

	  procedure prodq8  (x, y)
		  Find hi, lo so 2^64*hi + lo = x*y
		  return hi + (lo/8)
	  end prodq8

If x*y == 0 (mod 8), this returns (x*y/8) mod q, partially reduced.
The function has one 128-bit multiply (or two half-multiplies) and two IOPs.

	  procedure qreduce (x)
		  return( (x >> 61) + (x & q) )
	  end prodq8

Sometimes we`ll want a negated modular product
(e.g. when multiplying (a0 + I*a1) * (0 + I*b1) = (-a1*b1) + I*(a0*b1) )
while keeping the output in [0, q]. One can simplify

	q - [(tmp >> 61) + (tmp & q)]		(4 IOPs)

to (q & ~tmp) - (tmp >> 61), saving an IOP on Alpha (which has the BIC
machine instruction to effect and-not) and on SPARC (which has ANDN).
MIPS has NOR, and we can use (q & ~tmp) = NOR(~q, tmp) if we can afford a
register for ~q = 7 * 2^61.  The operation count on MIPS will be no worse
even if it needs an explicit NOT.

Since the FGT (especially power-of-2 radices) will involve sequential
pairwise adds and subtracts of CMUL_MOD outputs, we would like the modded
multiply outputs to be in [0, 2q], thus ensuring that multiply outputs can
be added/subtracted pairwise (with adds of k*q as needed to make subtract
outputs nonnegative) through two levels of a radix-8 or radix-16 transform
and still have results in [0,8*q], or really [0,2^64-1].
Peter Montgomery notes that in fact the a-operands can be somewhat
larger than q+7 and still satisfy this constraint, as follows:

	The value of lo will be exactly divisible by 8.
	We want to show tmp <= 2*q where tmp = hi + lo/8 = prodq8(a, b).
	There are two cases to consider:

	(1) When a + b < 2^62, the arithmetic-geometric mean inequality (AGM) gives

		a*b <= (a+b)^2/4 <= (2^62 - 1)^2/4 = 2^61 * q + 1/4.

	Here hi <= q and lo/8 <= q, so tmp <= 2q.

	(2) Otherwise a + b >= 2^62, i.e. a >= 2^62 - b > 2^61, so 2^61 <= a <= 2^61 + sqrt(2^61).
	Also b >= 2^62 - a > 2^61 - sqrt(2^61).
	So 0 < (a - 2^61)*(2^61 - b) < 2^61.  From a*b = 2^61*hi + lo/8
	and 0 <= lo/8 < 2^61 we conclude

		lo/8 = 2^61 - (a - 2^61)*(2^61 - b) = a*b - 2^61*(a + b - 2^61 - 1)

	and hence

	 hi = a + b - 2^61 - 1
	tmp = (a + b - 2^61 - 1) + (2^122 + 2^61 - 2^61*a - 2^61*b + a*b)
		= 2^122 - 1 + a*(1 - 2^61) + b*(1 - 2^61) + a*b
		= q^2 + 2q - a*q - b*q + a*b
		= 2q - (a - q)*(q - b) <= 2q

	Having tmp <= 2q = 2^62 - 2 always gives an output <= q (tmp = q gives q.)

UPSHOT: for 0 <= a <= 2^61 + sqrt(2^61),  0 <= b < 2^64, (b divisible by 8),
	        0 <= prodq8(a,b) <= 2q.

The sum (or difference) of any two such numbers cannot overflow into
the sign bit, i.e. if we are adding or subtracting two such products
we require no intermediate normalization step, except to add 2q to the
result of the subtract to put it into the same (nonnegative) range as
the result of the add. We thus have:

STANDARD COMPLEX MULTIPLY:
The standard complex multiply is of form

	xout = a0*b0 - a1*b1	(Real);		a0, a1, b0, b1 in [0, B]
	yout = a1*b0 + a0*b1	(Imag).

which we compute as

	xout = prodq8(a0, 8*b0) - prodq8(a1, 8*b1) + 2*q;	xout in [0, 4q]
	yout = prodq8(a1, 8*b0) + prodq8(a0, 8*b1);			yout in [0, 4q]

Assuming b0 and b1 have been premultiplied by 8, this costs 4 128-bit
multiplications and 8 IOPs for the four PRODQ8s, and a further 4 IOPs for
the adds and subtracts.

We choose to leave the two output digits unreduced, since they may be
put through one additional add/subtract cycle before needing reduction,
which may be beneficial to keep calls to qreduce to a minimum. In other
words, we provide positive-digit outputs both of which are in [0, 4q],
and leave it up to the user to reduce these or not.

On our two main target architectures, we thus have:

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  STANDARD MULTIPLY: Inputs: a0,a1,8*b0,8*b1. xout, yout in [0, 4q].	*
*																		*
*  MIPS:  4 128-bit multiplications (DMULTUs), 12 IOPS;					*
*  Alpha: 8 half-multiplies (4 MULQ, 4 UMULH), 12 IOPS.					*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

The scheme is quite multiply-heavy on the older Alpha 21064 and 21164, which
can start an integer half-multiply (MULQ or UMULH) only every 8th clock cycle,
only a bit better-balanced on the MIPS (where we need just 4 integer
multiplies, but these are not pipelined at all), but ideally suited to
the Alpha 21264 (where multiplies are fully pipelined.)


KARATSUBA MULTIPLY:
A complex product can be calculated using the Karatsuba trick to
reduce the number of 128-bit multiplies from 4 to 3:

	xout = a0*b0 - a1*b1
	yout = (a0 + a1)*(b0 + b1) - a0*b0  - a1*b1	(Imag2).

The computation of xout is the same as for the standard complex multiply.
For yout, the premultiplied-by-8 b0 and b1 terms maybe as large as 2^64-8,
so (b0 + b1) may overflow - to correct for this, we have two options:

(1) multiply the carry bit by 8 (since 2^64 mod q = 8) and add it back into
    the (mod-2^64) result. We don`t assume that there is a carry bit; rather
    we emulate it via the unsigned compare (b0 + b1) < b0 (which returns 1
    if there was a carry, 0 if not). Thus the overflow correction goes like

	tmp = b0 + b1;
	return( ((tmp < b0) << 3) + tmp )	(cost: 3 IOPS).

(2) Better is to instead premultiply b0 and b1 by 4 (or allocate
    a third register to hold 4*(b0 + b1)), then calculate

	t00 = prodq8(2*a0, 4*b0);			t00 in [0, 2q]
	t11 = prodq8(2*a1, 4*b1);			t11 in [0, 2q]
	t01 = prodq8(2*a0 + 2*a1, 4*b0 + 4*b1);		t01 in [0, 5q] <*** check! ***

Given a0, a1, 4*b0, 4*b1, this is 4 adds and 3 prodq8, for a total of
4 + 3*2 = 10 IOPs, in addition to 3 128-bit multiplies. Given a0, a1, b0, b1,
this costs 6 + 3*2 = 12 IOPs, in addition to 3 128-bit multiplies.

The real part of the result is t00 - t11, which is in [-2q, 2q]. We can add
2*q to bring this into [0, 4q]. t00 - t11 + 2*q costs 3 IOPs.

For the imaginary part, t01 - t00 - t11 appears to have a potential range of
5q + 2q + 2q = 9q > 2^64 if we use only the individual upper bounds
on the three variables.  But the three are not independent, since
t01 - t00 - t11 = (a0 + a1) * (b0 + b1) - a0*b0 - a1*b1 = a0*b1 + a1*b0,
which is in [0, 4q]. t01 - t00 - t11 costs 2 IOPs.

As with the standard complex multiply, we leave the outputs unreduced.

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*  KARATSUBA MULTIPLY: Inputs: a0,a1,4*b0,4*b1. xout, yout in [0, B].	*
*																		*
*  MIPS:  3 128-bit multiplications (DMULTUs), 15 IOPS;					*
*  Alpha: 6 half-multiplies (3 MULQ, 3 UMULH), 15 IOPS.					*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

Thus, Karatsuba seems preferable on all architectures except possibly
the Alpha 21264.

NOTE: If we can afford to fetch three constants from meory,
	we can use 8*b0, 8*b1, 4*(b0 + b1).  We save the add
	in 4*b0 + 4*b1 and can use 2*(a0 + a1) instead of doubling a0 and a1
	separately. We trade two adds for a load.

ESTIMATED CYCLES:

On MIPS, roughly assuming an 11-cycle latency and no pipelining of DMULTU,
we estimate 40 cycles (35 if successive CMUL_MODs overlap.)

On Alpha 21064 and 21164, assuming one multiply starts every eighth clock)
and MULQ takes 12 cycles to complete (i.e. we assume our final multiply is
a MULQ), with no overlap between complex modular multiplies, the minimum
number of cycles is obviously set by the multiply non-pipelineability and
is about 55 (50 if successive CMUL_MODs overlap.)

On Alpha 21264, assuming perfect pipelining and both integer units
completely busy, we could do the entire Karatsuba multiply in about 10
cycles, if successive CMUL_MODs overlap. 15-20 seems more realistic.
In either event, it`s a big gain over the previous Alpha architectures.

FLOATING-POINT VS. MODULAR:

Given that a well-pipelined sequence of floating-point complex multiplies
cost about 5-6 clocks per CMUL, the above modular timing estimates raise
the question: how can we expect the modular FGT to execute in anywhere
near the time of the floating-point FFT? Here is the idea: in each FFT
pass, there is only one phase where CMUL_MODs dominate, namely during
the twiddle multiply phase at the start of the pass. Here, the modular
algorithm will lag significantly behind the floating-point. However,
90-95% of the rest of the FFT pass consists of adds and subtracts -
multiplies by sqrt(-1) (primitive roots of order 2^2) can be effected
via adds and subtracts alone, and the only other multiplies (for radix-8)
are by the modular analog of sqrt(1/2), i.e. primitive roots of order 2^3,
which for q=2^p-1 have the simple form (-2)^[(p-1)/2]*(1-i), e.g. 2^30*(1-i)
for q = M61, and can thus be effected via a few simple shifts, adds and
ANDs, and don`t need the slower integer multiplies.

A radix-16 pass also needs 4 complex multiplies by roots of order 2^4,
and since these don`t seem to have any simple modular representation,
we`d have 15 complex twiddle multiplies and 4 more complex multiplies
in the finishing phase - this may tip the scales in favor of radix-8,
but note that there will also be savings due to 25% fewer FFT passes
using radix-16.

Anyway, in this large add/subtract-dominated phase, the floating adder
can only do one add per clock even with perfect pipelining, whereas
the two integer units can produce two results per cycle. We hope that
the modular FGT can thus "catch up" with the floating FFT at the finish.

Opcounts:
STANDARD MULTIPLY: 4 MUL_LOHI64, 15 IOPs;
KARATSUBA VARIANT: 3 MUL_LOHI64, 20 IOPs.

The CMUL_MODQ8 variant assumes the inputs are premultiplied by 8 and thus cuts IOPs to 11.
*/
#define KARATSUBA_CMUL	0	// Set = 1 to use Karatsuba variant

void cmul_modq(const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout, uint64*yout)
{
	uint64 t00,t01,t11;
	ASSERT((a0 <= bb && a1 <= bb && b0 <= q && b1 <= q), "ERROR: CMUL_MODQ input out of range!");

	// Bounds: b0,b1 in [0,q], so 4(b0+b1) in [0,8q]; prodq8lo/8 always in [0, q].
	t00 = prodq8(a0, b0<<3);		// a0    in [0, B]: prodq8hi(a0,8b0) in [0, q], t00 in [0, 2q]
	t11 = prodq8(a1, b1<<3);		// a1    in [0, B]: prodq8hi(a1,8b1) in [0, q], t11 in [0, 2q]
	*xout = t00 - t11 + q2;		// xout in [0, 4q]
	ASSERT(*xout <= q4, "ERROR: xout > 4q in CMUL_MODQ!");

#if !KARATSUBA_CMUL
	// Standard complex 4-multiply:
	t01   = prodq8(a0, b1<<3);			// a0    in [0, B]: prodq8hi(a0,8b1) in [0, q], t01 in [0, 2q]
	*yout = prodq8(a1, b0<<3) + t01;	// a1    in [0, B]: prodq8hi(a1,8b0) in [0, q], t10 in [0, 2q]; yout in [0, 4q]
	ASSERT(*yout <= q4, "ERROR: yout > 4q in CMUL_MODQ!");
#else
	// Karatsuba variant:
	t01 = prodq8((a0 + a1)<<1, (b0 + b1)<<2);	// prodq8hi( 2(a0+a1) , 4(b0+b1) ) in [0,4q], t01 in [0, 5q]
	*yout = qreduce(t01 - t00 - t11 + q4);	// t01 in [0, 5q] but t01-t00-t11 in [0,4q], so no overflow in t01-t00-t11+q4.
	ASSERT(t01 <= (q4 + q), "ERROR: t01 > 5q in Karatsuba-part of CMUL_MODQ!");
	// This version reduces both parts of the output:
  #if 0
	uint64 tmp = (q<<2) - t11;					// tmp in [2q, 4q]
	*xout = qreduce(t00 + tmp);			// t00 + tmp       in [2q, 6q]--->[0, B]
	*yout = qreduce(t01 - t00 + tmp);	// t01 - t00 + tmp in [2q, 6q]--->[0, B]
  #endif
#endif

	return;
}

// Variant of CMUL_MODQ which assumes the B-inputs are premultiplied by 8 and thus cuts IOPs to 11.
// This function is ideal for CMUL-by-a-constant, e.g. in FFT twiddle MULs.
// There is no Karatsuba option for this function, since (b0 + b1) overflows if b0,b1 are premultiplied by 8:
//
void cmul_modq8(const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout, uint64*yout)
{
	uint64 t00,t01,t11;
	ASSERT(( a0 <= bb && a1 <= bb ), "ERROR: CMUL_MODQ8 A-input out of range!");
	ASSERT((!(b0 & 7) && !(b1 & 7)), "ERROR: CMUL_MODQ8 B-input not divisible by 8!");

	// Bounds: b0,b1 in [0,q], so 4(b0+b1) in [0,8q]; prodq8lo/8 always in [0, q].
	t00 = prodq8(a0, b0);		// a0    in [0, B]: prodq8hi(a0,8b0) in [0, q], t00 in [0, 2q]
	t11 = prodq8(a1, b1);		// a1    in [0, B]: prodq8hi(a1,8b1) in [0, q], t11 in [0, 2q]
	*xout = t00 - t11 + q2;		// xout in [0, 4q]
	ASSERT(*xout <= q4, "ERROR: xout > 4q in CMUL_MODQ!");

	// Standard complex 4-multiply is only option here:
	t01   = prodq8(a0, b1);			// a0    in [0, B]: prodq8hi(a0,8b1) in [0, q], t01 in [0, 2q]
	*yout = prodq8(a1, b0) + t01;	// a1    in [0, B]: prodq8hi(a1,8b0) in [0, q], t10 in [0, 2q]; yout in [0, 4q]
	ASSERT(*yout <= q4, "ERROR: yout > 4q in CMUL_MODQ!");
	return;
}

/***************/

/* Complex-squaring:
Given input a0 + I*a1, computes the square xout + I*yout = (a0^2-a1^2) + I*(2*a0*a1) via

	xout = (a0+a1)*(a0-a1)
	yout = 2*a0*a1

Cost: 2 MUL_LOHI64, 11 IOPs.
*/
void csqr_modq(const uint64 a0, const uint64 a1, uint64*xout, uint64*yout)
{
	// This version reduces both parts of the output...
	*xout = prodq8((a0 + a1)<<1, (a0 - a1 + q)<<2);	// prodq8hi( 2(a0+a1) , 4(a0-a1+q) ) in [0,4q]; xout in [0,5q]
	ASSERT(*xout <= (q4+q), "ERROR: xout >= 5q in CSQR_MODQ!");

	*yout = prodq8(a0<<2, a1<<2);					// prodq8hi(     4*a0 , 4*a1       ) in [0,2q]; yout in [0,3q]
	ASSERT(*yout < (q4-q), "ERROR: yout > 3q in CSQR_MODQ!");
}

/***************/

/* Complex-plus-conjugate version of CMUL_MODQ. Given inputs A=a0 + I*a1 and B=b0 + I*b1,
computes
    xout1 + I*yout1 =      A *B = (a0*b0 - a1*b1) + I*(a0*b1 + a1*b0)
and
	xout2 + I*yout2 = c.c.(A)*B = (a0*b0 + a1*b1) + I*(a0*b1 - a1*b0).

Cost: 4 MUL_LOHI64, 16 IOPs, only a few IOPs more expensive than a single CMUL_MODQ.
*/
void cmul_conj_modq(const uint64 a0, const uint64 a1, const uint64 b0, const uint64 b1, uint64*xout1, uint64*yout1, uint64*xout2, uint64*yout2)
{
	uint64 t00,t11,t01 = b1<<3,t10 = b0<<3;	// Use t01,t10 as temps in first stages

	t00 = prodq8(a0, t10);	// a0    in [0, B]: prodq8hi(a0,8b0) in [0, q], t00 in [0, 2q]
	t11 = prodq8(a1, t01);	// a1    in [0, B]: prodq8hi(a1,8b1) in [0, q], t11 in [0, 2q]
	*xout1 = t00 - t11 + q2;	// xout1 in [0, 4q]
	*xout2 = t00 + t11;			// xout2 in [0, 4q]

	t01 = prodq8(a0, t01);	// a0    in [0, B]: prodq8hi(a0,8b1) in [0, q], t01 in [0, 2q]
	t10 = prodq8(a1, t10);	// a1    in [0, B]: prodq8hi(a1,8b0) in [0, q], t10 in [0, 2q]
	*yout1 = t01 + t10;			// yout1 in [0, 4q]
	*yout2 = t01 - t10 + q2;	// yout2 in [0, 4q]
}

/***************/

// Given an input (ord) returns a primitive root of q of that order via the root_re, root_im pointers.
// We use the fundamental primitive root (6+I) of order (q^2-1) = 2^62 * (2^60 - 1),
// which order must be divisible by (ord), obviously. We don't expect to need ord > 64-bits.
//
void prim_root_q(const uint64 ord, uint64*root_re, uint64*root_im)
{
	int i,zbits = trailz64(ord);
	uint64 r0,i0,rm,im,rtmp,itmp,pow;

	// Maximal order (q^2-1) = 2^62 * (2^60-1), allowing power-of-2 roots up to 2^62:
	ASSERT(zbits < 63, "PRIM_ROOT_Q: Maximal power-of-2 roots = 2^62!");

	// First raise result to the [(2^60-1)/(ord >> trailz(ord))]th power using LR binary powering:
	itmp = (1ull << 60) - 1;
	pow = itmp/(ord >> zbits);		// Odd component of the needed power; this should have 0 remainder for legal ord values
	ASSERT(itmp == pow*(ord >> zbits), "pow does not divide 2^60-1!");
	pow = pow << (leadz64(pow)+1);	// Left-justify pow and shift leftmost bit off.

	// 6 + I is a primitive root of full order q^2 - 1:
	r0 = rm = 6ull;		i0 = im = 1ull;
	while(pow) {
		// Squaring is here:
		csqr_modq(rm,im, &rtmp,&itmp);
		rm = qreduce(qreduce(rtmp)); im = qreduce(qreduce(itmp));
		// If leftmost bit = 1, multiply by (r0,i0):
		if((int64)pow < 0) {
			cmul_modq(rm,im,r0,i0, &rtmp,&itmp);
			rm = qreduce(qreduce(rtmp)); im = qreduce(qreduce(itmp));
		}
		// Promote the next-lower bit into the MSB slot. Again, this works because pow is odd:
		pow <<= 1;
	}

	// q^2 - 1 = 2^62*(2^60-1), so need to do [62-trailz(ord)] squarings...
	zbits = 62 - trailz64(ord);
	for(i = 0; i < zbits; i++) {
		cmul_modq(rm,im,rm,im, &rtmp,&itmp);
		rm = qreduce(qreduce(rtmp)); im = qreduce(qreduce(itmp));
	}

	// ...and normalize result one more time, to ensure outputs are in [0,q-1]:
	rm -= (-(uint64)(rm >= q)) & q;	im -= (-(uint64)(im >= q)) & q;
	*root_re = rm;	*root_im = im;
	return;
}

/***************/

// Raise input to the desired power using left-to-right binary powering.
//
void pow_modq(const uint64 power, const uint64 re_in, const uint64 im_in, uint64*re_out, uint64*im_out)
{
	int i = leadz64(power), zbits = 63-i;	// zbits = Number of bits to the right of leading ones bit.
	uint64 r0,i0,rm,im,rtmp,itmp, pow = power << (i+1);	// pow = Left-justify input power and shift leftmost bit off.
	// Special handling for 0th power:
	if(!power) {
		*re_out = 1ull;	*im_out = 0ull;
		return;
	}
	r0 = rm = re_in; i0 = im = im_in;
	for(i = 0; i < zbits; i++) {
		// Squaring is here:
		csqr_modq(rm,im, &rtmp,&itmp);
		rm = qreduce(qreduce(rtmp)); im = qreduce(qreduce(itmp));
		// If leftmost bit = 1, multiply by (r0,i0):
		if((int64)pow < 0) {
			cmul_modq(rm,im,r0,i0, &rtmp,&itmp);
			rm = qreduce(qreduce(rtmp)); im = qreduce(qreduce(itmp));
		}
		// Promote the next-lower bit into the MSB slot. Again, this works because pow is odd:
		pow <<= 1;
	}

	// ...and normalize result one more time, to ensure outputs are in [0,q-1]:
	rm -= (-(uint64)(rm >= q)) & q;	im -= (-(uint64)(im >= q)) & q;
	*re_out = rm;	*im_out = im;
	return;
}

#undef q
#undef q2
#undef q4
#undef bb

#endif	// USE_FGT61
