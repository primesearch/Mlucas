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

#define MI64_DEBUG	0

#define ARRAYS_OVERLAP(x, lenx, y, leny)	( (x <= y) && (x+len) > y ) || ( (x > y) && (y+len) > x )

#if MI64_DEBUG
	char string1[STR_MAX_LEN], string2[STR_MAX_LEN];
#endif

#include "mi64.h"

/*******************/

/* Set-X-equal-to-Y: */
void	mi64_set_eq(uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;
	DBG_ASSERT(HERE, len != 0, "mi64_set_eq: zero-length array!");
	if(x == y) return;
	for(i = 0; i < len; ++i)
	{
		x[i] = y[i];
	}
}

/* Set-X-equal-to-scalar-A: */
void	mi64_set_eq_scalar(uint64 x[], const uint64 a, uint32 len)
{
	uint32 i;

	DBG_ASSERT(HERE, len != 0, "mi64_set_eq: zero-length array!");
	x[0] = a;
	for(i = 1; i < len; ++i)
	{
		x[i] = 0ull;
	}
}

/* Set-X-equal-to-0: */
void	mi64_clear(uint64 x[], uint32 len)
{
	uint32 i;
	for(i = 0; i < len; ++i)
	{
		x[i] = 0ull;
	}
}

/* Set selected bit: */
void	mi64_set_bit(uint64 x[], uint32 bit)
{
	uint32 bit_in_word = bit&63, word = (bit >> 6);
	x[word] |= (0x0000000000000001ull << bit_in_word);
}

/* Return 1 if selected bit is set, 0 otherwise: */
int		mi64_test_bit(const uint64 x[], uint32 bit)
{
	uint32 bit_in_word = bit&63, word = (bit >> 6);
	return (int)(x[word] >> bit_in_word) & 0x1;
}

/*******************/

/*
Left-shift a base-2^64 int x[] by (nbits) bits, returning the result in y[].
Any off-shifted bits are lost. No array-bounds checking is done.

Allows in-place operation, i.e. x == y.
*/
void	mi64_shl(const uint64 x[], uint64 y[], uint32 nbits, uint32 len)
{
	int i;
	uint32 nwshift = (nbits >> 6), rembits = (nbits & 63);

	DBG_ASSERT(HERE, len != 0, "mi64_shl: zero-length array!");

	/* Special-casing for shift-into-Bolivian: */
	if(nwshift >= len)
	{
		nwshift = len;
		rembits = 0;
	}

	/* Take care of the whole-word part of the shift: */
	for(i = len-1; i >= (int)nwshift; i--)
	{
		y[i] = x[i-nwshift];
	}
	for(i = nwshift-1; i >= 0; i--)
	{
		y[i] = 0;
	}

	/* If nbits not an exact multiple of the wordlength, take care of remainder: */
	if(rembits)
	{
		/* Process all but the least-significant element, in reverse order: */
		for(i = len-1; i > 0; i--)
		{
			y[i] = (y[i] << rembits) + (y[i-1] >> (64-rembits));
		}
		/* Least-significant element gets zeros shifted in from the right: */
		y[0] <<= rembits;
	}
}

/*******************/

/*
Logical-right-shift a base-2^64 int x[] by (nbits) bits, returning the result in y[].
Any off-shifted bits are lost. No array-bounds checking is done.

Allows in-place operation, i.e. x == y.
*/
void	mi64_shrl(const uint64 x[], uint64 y[], uint32 nbits, uint32 len)
{
	int i;
	uint32 nwshift = (nbits >> 6), rembits = (nbits & 63);

	DBG_ASSERT(HERE, len != 0, "mi64_shrl: zero-length array!");

	/* Special-casing for shift-into-Bolivian: */
	if(nwshift >= len)
	{
		nwshift = len;
		rembits = 0;
	}

	/* Take care of the whole-word part of the shift: */
	for(i = 0; i < len-nwshift; i++)
	{
		y[i] = x[i+nwshift];
	}
	for(i = len-nwshift; i < len; i++)
	{
		y[i] = 0;
	}

	/* If nbits not an exact multiple of the wordlength, take care of remainder: */
	if(rembits)
	{
		/* Process all but the most-significant element, in reverse order: */
		for(i = 0; i < len-1; i++)
		{
			y[i] = (y[i] >> rembits) + (y[i+1] << (64-rembits));
		}
		/* Least-significant element gets zeros shifted in from the right: */
		y[len-1] >>= rembits;
	}
}

/*******************/
/* unsigned compare: these are all built from just two elemental functions: < and == */
uint32	mi64_cmpult(const uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;

	DBG_ASSERT(HERE, len != 0, "mi64_cmpult: zero-length array!");
	for(i=len-1; i!=0; i--)	/* Loop over all but the 0 elements while equality holds.... */
	{
		if(x[i] < y[i])
			return TRUE;
		if(x[i] > y[i])
			return FALSE;
	}
	return x[0] < y[0];
}

uint32	mi64_cmpeq(const uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;

	DBG_ASSERT(HERE, len != 0, "mi64_cmpeq: zero-length array!");
	for(i = 0; i < len; i++)
	{
		if(x[i] != y[i])
			return FALSE;
	}
	return TRUE;
}

uint32	mi64_cmplt_scalar(const uint64 x[], uint64 a, uint32 len)
{
	return ( (mi64_getlen(x, len) <= 1) && (x[0] < a) );
}

uint32	mi64_cmpgt_scalar(const uint64 x[], uint64 a, uint32 len)
{
	return ( (x[0] > a) || (mi64_getlen(x, len) > 1) );
}

uint32	mi64_cmpeq_scalar(const uint64 x[], uint64 a, uint32 len)
{
	uint32 i;

	DBG_ASSERT(HERE, len != 0, "mi64_cmpeq_scalar: zero-length array!");
	/* Structure things so we examine the low digits first */
	if(x[0] == a)
	{
		for(i = 1; i < len; i++)
		{
			if(x[i])
				return FALSE;
		}
		return TRUE;
	}
	else
		return FALSE;
}

/*******************/

/* Count trailing zeros of a 64-bit array-int (specify signed, but that's not crucial.)

	RETURNS: # of trailing zeros if x nonzero; 0 otherwise.
*/
uint32	mi64_trailz(const uint64 x[], uint32 len)
{
	int i;
	uint32 tz = 0, iszero = TRUE;

	DBG_ASSERT(HERE, len != 0, "mi64_trailz: zero-length array!");
	for(i=0; i<len; i++)
	{
		if(x[i]==0)
			tz += 64;
		else
		{
			iszero = FALSE;
			tz += trailz64(x[i]);
			break;
		}
	}

	if(iszero)
		return 0;
	else
		return tz;
}

/*******************/

/* Returns # of leading zeros of a 64-bit multiword int. */
uint32	mi64_leadz(const uint64 x[], uint32 len)
{
	int i;
	uint32 lz = 0, iszero = TRUE;

	DBG_ASSERT(HERE, len != 0, "mi64_leadz: zero-length array!");
	for(i=len-1; i>=0; i--)
	{
		if(x[i]==0)
			lz += 64;
		else
		{
			iszero = FALSE;
			lz += leadz64(x[i]);
			break;
		}
	}
	return lz;
}

/*******************/

/* Extract leading 64 significant bits, or 128-most-significant-bits-starting-at-a-specified-bit-position.
Note that despite the similarity of their names these functions do slightly different things - mi64_extract_lead64()
starts at a given word of an mi64, *finds* the leading nonzero bit in that word or [if none there] in the remaining [len - 1] words.
mi64_extract_lead128() extracts the leading 128 bits from the length-len vector X, starting at the *specified* MSB position of
X[len-1] passed in leading-zeros form in (nshift), filling in any missing low-order bits in lead_x[] with zero.
*/

/* Extract leading 64 significant bits, filling in any missing low-order bits in the result with zero.
Returns:
	- The leading 64 bits in *result;
	- Bitlength of x[] ( = 64*len - mi64_leadz(x) ) in function value. If this is < 64 (i.e. only the x[0] word, if any, is nonzero),
	  the absolute value of the difference between it and 64 equals the number of "borrowed" to fill in missing low-order bits.
*/
uint32 mi64_extract_lead64(const uint64 x[], uint32 len, uint64*result)
{
	uint32 i,nshift,nwshift,rembits;

	ASSERT(HERE, len != 0, "mi64_extract_lead64: zero-length array!");

	nshift = mi64_leadz(x, len);
	nwshift = (nshift >> 6), rembits = (nshift & 63);
	/* shift-word count may == len, but only if x[] = 0: */
	if(nwshift >= len)
	{
		DBG_ASSERT(HERE, nwshift == len, "mi64_extract_lead64: nwshift out of range!");
		DBG_ASSERT(HERE, mi64_iszero(x, len), "mi64_extract_lead64: expected zero-valued array!");
		*result = 0ull;
	}
	else
	{
		i = len-1-nwshift;
		if(rembits == 0)
		{
			*result = x[i];
		}
		else
		{
			*result  = (x[i] << rembits);
			if(i > 0)
			*result += (x[i-1] >> (64-rembits));
		}
	}
	return (len << 6) - nshift;
}

/* Convert mi64 to double, rounding lowest-order mantissa bit: */
double	mi64_cvt_double(const uint64 x[], uint32 len)
{
	uint64 itmp64, lead64, lead64_rnd;
	double retval;
	int pow2 = mi64_extract_lead64(x, len, &lead64);
	if(lead64 == 0ull) return 0.0;
	DBG_ASSERT(HERE,(lead64 >> 63) == 1ull, "mi64_cvt_double: lead64 lacks leftmost ones bit!");
	/*  round based on 1st neglected bit: */
	lead64_rnd = (lead64 >> 11) + ((lead64 >> 10) & 0x0000000000000001ull);
	/* exponent: */
	itmp64 = ((0x3FD + (uint64)pow2) << 52) + lead64_rnd;
	/* Add in mantissa, with hidden bit made explicit, hence the 0x3FD (rather than 0x3FE) initializer */
	itmp64 += lead64_rnd;
	ASSERT(HERE, itmp64 > lead64_rnd , "mi64_cvt_double: Exponent overflows IEEE64 field");
	retval = *(double *)&itmp64;
	/* GCC compiler bug: needed to insert the explicit range-check here, otherwise compiler 'optimizes' the (*(double *)&itmp64) to zero: */
	if(retval < 0.0)
	{
		sprintf(cbuf, "rng_isaac_rand_double_norm_pos: lead64 = %16llx, itmp64 = %16llx, retval = %lf not in [0,1]!\n", lead64, itmp64, retval);
		ASSERT(HERE, 0, cbuf);
	}
	return retval;
}

/*	Extract leading 128 bits from the length-len vector X, starting at the
	designated MSB position of X[len-1] passed in leading-zeros form in (nshift), which is required to be < 64.

	RETURNS: the desired 128 bits, left-justified, in the 2-array lead_x.
	If X has fewer than 128 bits below nshift, i.e. if (len*64 - nshift) < 128,
	missing low-order bits in lead_x[] are filled in with zero.

	Visually we have this:

    [                          x[len-1]                             ] [                          x[len-2]                             ] [                          x[len-3] (if nshift > 0)             ]
     321098765432109876543210987654321098765432109876543210987654321   321098765432109876543210987654321098765432109876543210987654321   321098765432109876543210987654321098765432109876543210987654321
     |-------- nshift --------|                                                                                                                                    |--------- (64-nshift) -------------|
                               |-------------------- 128 bits to be copied into lead_x[], right-padded with zeros if (len*64 - nshift) < 128 ---------------------|
*/
void mi64_extract_lead128(const uint64 x[], uint32 len, uint32 nshift, uint64 lead_x[])
{
	lead_x[0] = lead_x[1] = 0;

	DBG_ASSERT(HERE, len != 0, "mi64_extract_lead128: zero-length array!");
	DBG_ASSERT(HERE, nshift < 64, "mi64_extract_lead128: illegal nshift value!");

	/* Syntax reminder:
		MVBITS(from_integer,low_bit_of_from_integer,num_bits,to_integer,insert_bits_in_to_integer_starting_at_this_low_bit)
	*/
		mvbits64(x[len-1],0        ,64-nshift,&lead_x[1], nshift);	/* move nonzero bits of leading digit into high word (2) of 128-bit slot, left-justifying...*/

	if(len > 1)
	{
	  if(nshift)
	  {
		mvbits64(x[len-2],64-nshift,nshift   ,&lead_x[1], 0     );	/* if leading digit had any leftmost zero bits, fill low-order bits of high word with high bits of next digit...*/
	  }
		mvbits64(x[len-2],0        ,64-nshift,&lead_x[0], nshift);	/* move leading bits of (lenU-1)st digit into low word (1) of 128-bit slot, left-justifying...*/
	}

	if(len > 2)
	{
	  if(nshift)
	  {
		mvbits64(x[len-3],64-nshift,nshift   ,&lead_x[0], 0     );	/* if necessary, fill low-order bits of low word with high bits of (lenU-2)nd digit. */
	  }
	}
}

/*******************/

uint32	mi64_iszero(const uint64 x[], uint32 len)
{
	return mi64_getlen(x, len) == 0;
}

/*******************/
/* Given a vector int x[] and a dimension (len), returns the actual length of x[] based on the position of the most-significant word.
If all zero words, returns 0.
*/
uint32	mi64_getlen(const uint64 x[], uint32 len)
{
	int i;
	for(i = len-1; i >= 0; i--)
	{
		if(x[i] != 0)
			break;
	}
	return i+1;	/* Returns proper result (0) if len = 0, since i will be init'ed to -1 and the loop exit criterion immediately satisfied */
}

/*******************/
/* Clear high (nclear) words of a vector int x[start_word]. */
void	mi64_setlen(uint64 x[], uint32 oldlen, uint32 newlen)
{
	if(newlen <= oldlen) return;
	mi64_clear(&x[oldlen], (newlen - oldlen));
}

/*******************/

/*
	Unsigned add of two base-2^64 vector ints X + Y.
	Any or all of X, Y and Z can point to the same array object.
	Return any exit carry - up to user to determine what to do if this is nonzero.
*/
uint64	mi64_add(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 i;
	uint64 tmp, cy = 0;

	DBG_ASSERT(HERE, len != 0, "mi64_add: zero-length array!");
	for(i=0; i<len; i++)
	{
		tmp = x[i] + cy;
		cy  = (tmp < x[i]);
		tmp = tmp + y[i];
		cy += (tmp < y[i]);
		z[i] = tmp;
	}
	return cy;
}

/*
	Unsigned sub of two base-2^64 vector ints X - Y.
	Any or all of X, Y and Z can point to the same array object.
	Return any exit borrow - up to user to determine what to do if this is nonzero.
*/
uint64	mi64_sub(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 i;
	uint64 tmp, tmp2, bw = 0;

	DBG_ASSERT(HERE, len != 0, "mi64_sub: zero-length array!");
	for(i=0; i<len; i++)
	{
		tmp = x[i] - bw;
		bw  = (tmp > x[i]);
//bw  = ((uint64)tmp > (uint64)x[i]);
		DBG_ASSERT(HERE, bw == ((uint64)tmp > (uint64)x[i]), "mi64_sub: compiler using signed compare (tmp > x[i])!");
		/* Need an extra temp here due to asymmetry of subtract: */
		tmp2= tmp - y[i];
		bw += (tmp2 > tmp);
//bw += ((uint64)tmp2 > (uint64)tmp);
		DBG_ASSERT(HERE, (tmp2 > tmp) == ((uint64)tmp2 > (uint64)tmp), "mi64_sub: compiler using signed compare (tmp2 > tmp)!");
		z[i] = tmp2;
	}
	return bw;
}

/*******************/

/* arithmetic negation: */
void	mi64_nega(const uint64 x[], uint64 y[], uint32 len)
{
	uint32 i;
	DBG_ASSERT(HERE, len != 0, "mi64_add: zero-length array!");
	y[0] = -x[0];
	for(i=1; i<len; i++)
	{
		y[i] = ~x[i];
	}
}

void	mi64_negl(const uint64 x[], uint64 y[], uint32 len)
{
	uint32 i;
	DBG_ASSERT(HERE, len != 0, "mi64_add: zero-length array!");
	for(i=0; i<len; i++)
	{
		y[i] = ~x[i];
	}
}


/*******************/
/*
	Add of 64-bit scalar A to base-2^64 vector int X. Allows In-place addition.
	Return any exit carry - up to user to determine what to do if this is nonzero.
	(Typically one would store the return value in x[len] and then increment len.)
*/
uint64	mi64_add_scalar(const uint64 x[], uint64 a, uint64 y[], uint32 len)
{
	uint32 i;
	uint64 cy = a;

	DBG_ASSERT(HERE, len != 0, "mi64_add_scalar: zero-length array!");
	if(x == y)
	{
		/* In-place: Only need to proceed until carry peters out: */
		for(i=0; i<len && cy != 0; i++)
		{
			y[i] = x[i] + cy;
			cy = (y[i] < cy);
		}
	}
	else
	{
		for(i=0; i<len; i++)
		{
			y[i] = x[i] + cy;
			cy = (y[i] < cy);
		}
	}
	return cy;
}

/*
	Sub of 64-bit scalar A from base-2^64 vector int X. Allows In-place subtraction.
	Return any exit borrow - up to user to determine what to do if this is nonzero.
*/
uint64	mi64_sub_scalar(const uint64 x[], uint64 a, uint64 y[], uint32 len)
{
	uint32 i;
	uint64 bw = a, tmp;

	DBG_ASSERT(HERE, len != 0, "mi64_sub_scalar: zero-length array!");
	if(x == y)
	{
		/* In-place: Only need to proceed until borrow peters out: */
		for(i=0; i<len && bw != 0; i++)
		{
			tmp = x[i] - bw;
			/* Need an extra temp here due to asymmetry of subtract, since x[i] and y[i] may point to the same memloc: */
			bw = (tmp > x[i]);
			y[i] = tmp;
		}
	}
	else
	{
		for(i=0; i<len; i++)
		{
			tmp = x[i] - bw;
			bw = (tmp > x[i]);
			y[i] = tmp;
		}
	}
	return bw;
}

/*******************/

/*
	Unsigned multiply base-2^64 vector int x[] by scalar A, returning result in y[].
	Return any exit carry, rather than automatically storing it in a hypothetical
	[len+1]st array slot - up to user to determine what to do if this is nonzero.
	(Typically one would store the return value in y[len] and then increment len.)
	This has the advantage that we need assume nothing about the allocated length
	of the x[] array within the function.

	There is no code to check for special values of a, e.g. 0 or 1.

	Allows in-place, i.e. x == y.
*/
uint64	mi64_mul_scalar(const uint64 x[], uint64 a, uint64 y[], uint32 len)
{
	uint32 i;
	uint64 lo, hi, cy = 0;

	DBG_ASSERT(HERE, len != 0, "mi64_mul_scalar: zero-length array!");
	for(i=0; i<len; i++)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		y[i] = lo + cy;
		/*
		A*x[i] is at most (2^64 - 1)^2, i.e. hi <= 2^64 - 2. Since the
		carry out of the low 64 bits of the 128-bit add is at most 1,
		cy + hi <= 2^64 - 1, hence we don't need to check for a carry there,
		only when adding the result to x[i+1] on the next pass through the loop:
		*/
		cy = hi + (y[i] < lo);
	}
	return cy;
}

/*******************/

/*
	Unsigned multiply of base-2^64 vector ints X * Y, having respective lengths lenX, lenY.
	Result is returned in vector int Z. For simplicity, NEITHER X NOR Y MAY OVERLAP Z,
	though X and Y may point to the same (or to overlapping) memory addresses.
	An added restriction is that Z must be large enough to store the result of the multiply,
	i.e. has at least (lenX + lenY) allocated 64-bit integer elements.

	The routine first find the larger of the 2 input vectors (in terms of number of elements,
	i.e. lenA = max(lenX, lenY), lenB = min(lenX, lenY) and pointers A and B set to point to
	the corresponding input vector. It then performs a grammar-school multiply in (lenB) passes,
	during each of which the entire vector A is multiplied by the (pass)th 64-bit scalar element
	of vector B and the resulting partial product shifted and added into vector Z.

	There is no code to check for special values of the inputs, e.g. 0 or 1; However, if one of
	the inputs is specified to have zero length, lenZ is set to zero and an immediate return effected.

	RETURNS:
		- product in vector Z;
		- actual product length via pointer *lenZ;
*/
void	mi64_mul_vector(const uint64 x[], uint32 lenX, const uint64 y[], uint32 lenY, uint64 z[], uint32 *lenZ)
{
	uint32 i, j, lenA, lenB;
	const uint64 *A, *B;
	uint64 tmp64;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u;
	static uint32 dimU = 0;

	DBG_ASSERT(HERE, lenX != 0, "mi64_mul_vector: zero-length X-array!");
	DBG_ASSERT(HERE, lenY != 0, "mi64_mul_vector: zero-length Y-array!");
	DBG_ASSERT(HERE, x != z, "mi64_mul_vector: X and Z point to same array object!");
	DBG_ASSERT(HERE, y != z, "mi64_mul_vector: Y and Z point to same array object!");
	DBG_ASSERT(HERE, lenZ != 0x0, "mi64_mul_vector: Null lenZ pointer!");

	/* Init z[] = 0: */
	for(i = 0; i < lenX + lenY; i++)
	{
		z[i] = 0;
	}

	/* Find larger of the 2 inputs (in terms of actual length, not nominal length: */
	lenX = mi64_getlen(x, lenX);
	lenY = mi64_getlen(y, lenY);

	if(lenX >= lenY)
	{
		lenA = lenX;	lenB = lenY;
		A = x;			B = y;
	}
	else
	{
		lenA = lenY;	lenB = lenX;
		A = y;			B = x;
	}
	/* If the actual length of the smaller argument = 0, nothing to do, return 0: */
	if(lenB == 0)
	{
		DBG_WARN(HERE, "mi64_mul_vector: Zero-length input vector detected!", "", TRUE);
		*lenZ = 0;
		return;
	}
	else
		*lenZ = lenA;

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (lenA+1))
	{
		free((void *)u);
		u = (uint64 *)calloc((lenA+1), sizeof(uint64));
	}

	/* Loop over remaining (lenB-1) elements of B[], multiplying A by each, and
	using u[] as a scratch array to store B[i]*A[] prior to adding to z[]: */
	for(j = 0; j < lenB; j++)
	{
		u[lenA] = mi64_mul_scalar(A, B[j], u, lenA);

		/* Add j-word-left-shifted u[] to z[]: we don't need an explicit shift for this, we simply add
		u to z beginning at the (j)th slot of z. Remember, need to specify vector adds of length (lenA+1):
		*/
		tmp64 = mi64_add(&z[j], u, &z[j], lenA+1);	/* Any carryout would go into z[lenA+j+1] */
		/* Product of (1-word) x (N-word) should never exceed (N+1) words: */
		DBG_ASSERT(HERE, tmp64 == 0, "mi64_add(&z[j], u, &z[j], lenA+1) != 0");	(*lenZ)++;
	}

	/* Return actual length of result vector in the function argument, so that if one or
	more leading terms of the result is zero, caller can adjust vector length accordingly:
	*/
	*lenZ = mi64_getlen(z, *lenZ);
	DBG_ASSERT(HERE, *lenZ <= lenA + lenB, "mi64_mul_vector: *lenZ > (lenA + lenB)!");
}

/* Low and high-half mul allow in-place: */
void	mi64_mul_vector_lo_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 i, j;
	uint64 tmp64;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;

	DBG_ASSERT(HERE, len != 0, "`mi64_mul_vector_lo_half: zero-length X-array!");

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		free((void *)u);
		free((void *)v);
		u = (uint64 *)calloc((len+1), sizeof(uint64));
		v = (uint64 *)calloc((len*2), sizeof(uint64));
	}

	/* Loop over the elements of y[], multiplying x[] by each, and
	using u[] as a scratch array to store x[]*y[j] prior to adding to z[].

	For the high-half version, only want terms x[i]*y[j] with (i + j) >= len,
	plus the high halves of the (i + j) = len-1 terms for the requisite carryins.
	*/
	for(j = 0; j < len; j++)
	{
		if(y[j] == 0)
			continue;
		u[len] = mi64_mul_scalar(x, y[j], u, len);

		/* Add j-word-left-shifted u[] to v[]: we don't need an explicit shift for this, we simply add
		u to v beginning at the (j)th slot of v. Remember, need to specify vector adds of length (len+1):
		*/
		tmp64 = mi64_add(&v[j], u, &v[j], len+1);	/* Any carryout would go into v[len+j+1] */

		/* Product of (1-word) x (N-word) should never exceed (N+1) words: */
		DBG_ASSERT(HERE, tmp64 == 0, "mi64_add(&z[j], u, &v[j], len+1) != 0");
	}

	/* Copy v[0:len-1] into z[0:len-1]: */
	for(i = 0; i < len; i++)
	{
		z[i] = v[i];
	}
}

void	mi64_mul_vector_hi_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 i, j;
	uint64 tmp64;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;

	DBG_ASSERT(HERE, len != 0, "`mi64_mul_vector_hi_half: zero-length X-array!");

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		free((void *)u);
		free((void *)v);
		u = (uint64 *)calloc((len+1), sizeof(uint64));
		v = (uint64 *)calloc((len*2), sizeof(uint64));
	}

	/* Loop over the elements of y[], multiplying x[] by each, and
	using u[] as a scratch array to store x[]*y[j] prior to adding to z[].

	For the high-half version, only want terms x[i]*y[j] with (i + j) >= len,
	plus the high halves of the (i + j) = len-1 terms for the requisite carryins.
	*/
	for(j = 0; j < len; j++)
	{
		if(y[j] == 0)
			continue;
		u[len] = mi64_mul_scalar(x, y[j], u, len);

		/* Add j-word-left-shifted u[] to v[]: we don't need an explicit shift for this, we simply add
		u to v beginning at the (j)th slot of v. Remember, need to specify vector adds of length (len+1):
		*/
		tmp64 = mi64_add(&v[j], u, &v[j], len+1);	/* Any carryout would go into v[len+j+1] */

		/* Product of (1-word) x (N-word) should never exceed (N+1) words: */
		DBG_ASSERT(HERE, tmp64 == 0, "mi64_add(&z[j], u, &v[j], len+1) != 0");
	}

	/* Copy v[len:2*len-1] into z[0:len-1]: */
	for(i = 0; i < len; i++)
	{
		z[i] = v[i+len];
	}
}


/*************************************************************/
/****** STUFF RELATED TO FFT-BASED MULTIPLY BEGINS HERE ******/
/*************************************************************/

/*
!	Functions to convert a generic multiword int (i.e. not w.r.to a particular
!	modulus) from/to balanced-digit fixed-base floating-point form and uint64 form.
*/

/*
	EFFECTS: Converts pair of uint64 multiword ints x[],y[] ==> double (balanced-digit base-FFT_MUL_BASE
			representation) and stores the resulting floating values in packed-complex form in
			a[] (x[] terms in real/even-index terms of a[], y[] terms in imaginary/odd-index.)

	RETURNS: length of the result (i.e. number of significant complex floating-point words).

	ASSUMES: User has allocated enough room in a[] to store the complete result in
			balanced-digit form, and that there is no remaining output carry.
*/
uint32	mi64_cvt_uint64_double(const uint64 x[], const uint64 y[], uint32 len, double a[])
{
	uint32 i, j = 0, jpad, k;
	 int64 cyi, cyj, itmp, jtmp;
	uint64 curr_re64, curr_im64, bitsm1 = FFT_MUL_BITS-1, basem1 = FFT_MUL_BASE-1;

	DBG_ASSERT(HERE, len != 0, "mi64_cvt_uint64_double: zero-length array!");

	/* Only constant base 2^16 is supported for this conversion at present: */
	DBG_ASSERT(HERE, FFT_MUL_BITS == 16, "mi64_cvt_uint64_double: FFT_MUL_BITS != 16");

	/* Redo the quicker checks of those done in util.c::check_nbits_in_types() */
	DBG_ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_double_uint64: FFT_MUL_BASE not pure-integer!");
	DBG_ASSERT(HERE, FFT_MUL_BASE < TWO54FLOAT, "mi64_cvt_double_uint64: FFT_MUL_BASE >= maximum allowed value of 2^54!");

	/* As we extract each floating-point word, balance it and set
	resulting carry into next FP word: */
	cyi = cyj = 0;
	for(i=0; i<len; i++)
	{
		curr_re64 = x[i];
		curr_im64 = y[i];

		for(k=0; k<4; k++)
		{
			j = (i<<2) + k;
			j *= 2;	/* Since have two input arrays */
			jpad = j + ( (j >> DAT_BITS) << PAD_BITS );	/* Padded-array index */
			/* Extract the current 16-bit subfield: */
			itmp = cyi + ((curr_re64>>(k<<4)) & basem1);
			jtmp = cyj + ((curr_im64>>(k<<4)) & basem1);

			/* Is the current digit >= (base/2)? Note the extra != 0 check here which forces the carries to be at most 1 -
			that's needed for the extreme case tmp == base, in which the simple right-shift by bitsm1 gives cy = 2, which would cause
			(cy << FFT_MUL_BITS) = base*2, giving a "notmalized" result = -base, rather than the desired 0: */
			cyi = (int64)(((uint64)itmp>>bitsm1) != 0);	/* Cast to unsigned to ensure logical right-shift */
			cyj = (int64)(((uint64)jtmp>>bitsm1) != 0);	/* Cast to unsigned to ensure logical right-shift */
			/* If yes, balance it by subtracting the base: */
			/* RHS terms must be signed to prevent integer underflow-on-subtract: */
			a[jpad  ] = (double)(itmp - (cyi<<FFT_MUL_BITS));
			a[jpad+1] = (double)(jtmp - (cyj<<FFT_MUL_BITS));
		}
	}

	DBG_ASSERT(HERE, cyi <= 1,"mi64_cvt_double_uint64: Output Real carry out of range!");
	DBG_ASSERT(HERE, cyj <= 1,"mi64_cvt_double_uint64: Output Imag carry out of range!");

	/* It is desirable to not have the FP vector length exceed 4*len,
	so suppress any output carry by folding back into MS array element:
	*/
	if(cyi)
	{
		DBG_ASSERT(HERE, a[j  ] < 0,"mi64_cvt_double_uint64: MS array element >= 0!");
		a[j  ] += FFT_MUL_BASE;
	}

	if(cyj)
	{
		DBG_ASSERT(HERE, a[j+1] < 0,"mi64_cvt_double_uint64: MS array element >= 0!");
		a[j+1] += FFT_MUL_BASE;
	}

	return (4*len);
}

/*
	Takes a pair of multiprecision real vectors stored in packed-complex form in
	the double[] A-array, where each real vector has (n) significant words, and
	converts the real vectors to uint64[] form, storing uint64[] form of the
	vector stored in the real (even-index) terms of A[] in x[], imaginary in y[].

	ASSUMES:

		- Each real input vector correponds to a nonnegative multiword integer,
		i.e. even though the individual digits of Re(A[]) and Im(A[]) may be of either
		sign (i.e. due to balanced-digit representation), the most-significant digit
		of each real vector (assuming the vector is nonzero) MUST BE > 0.

		- Input array has been properly rounded (i.e. all elements have zero fractional
		part) and normalized, i.e. examination of MSW suffices to determine sign.

		- Floating-point base may be > 2^32 (e.g. for mixed FFT/FGT),
		but is < 2^54, i.e. properly-normalized balanced digits are
		in [-2^53, +2^53] and thus exactly representable via an IEEE double.

	RETURNS: length of the result vectors (i.e. number of significant uint64 words).
*/
uint32	mi64_cvt_double_uint64(const double a[], uint32 n, uint64 x[], uint64 y[])
{
	uint32 len;
	int curr_bits,i,j,nbits;
	int64 cyi, cyj, itmp, jtmp;
	uint64 curr_re64, curr_im64;

	DBG_ASSERT(HERE, n != 0, "mi64_cvt_double_uint64: zero-length array!");

	/* Redo the quicker checks of those done in util.c::check_nbits_in_types() */
	DBG_ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_double_uint64: FFT_MUL_BASE not pure-integer!");
	DBG_ASSERT(HERE, FFT_MUL_BASE < TWO54FLOAT, "mi64_cvt_double_uint64: FFT_MUL_BASE >= maximum allowed value of 2^54!");
	/*
	!...Make sure MSW of Re(A[]) and Im(A[]) in the
	balanced-representation form are both >= 0:
	*/
	/* Re(A[]) stored in even terms: */
	for(i=2*n-2; i >= 0; i-=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );
		if(a[j] != 0.0)
		{
			DBG_ASSERT(HERE, a[j] > 0.0, "mi64_cvt_double_uint64: MSW(Re(A[])) < 0!");
			break;
		}
	}
	/* Im(A[]) stored in odd terms: */
	for(i=2*n-1; i >= 1; i-=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );
		if(a[j] != 0.0)
		{
			DBG_ASSERT(HERE, a[j] > 0.0, "mi64_cvt_double_uint64: MSW(Im(A[])) < 0!");
			break;
		}
	}

	/*...Now form the terms of the uint64 representation:	*/

	len   = 0;		/* Full 64-bit words accumulated so far in the residue	*/
	nbits = 0;		/* Total bits accumulated so far in the x-array	*/
	curr_bits = 0;	/* Total bits currently stored in the 64-bit accumulator word */
	curr_re64 = 0;	/* Current value of the 64-bit accumulator word */
	curr_im64 = 0;

	cyi = cyj = 0;	/* Init carry */
	for(i = 0; i < 2*n; i+=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );

		DBG_ASSERT(HERE, curr_bits < 64,"mi64_cvt_double_uint64: curr_bits < 64");
		itmp = (uint64)1<<curr_bits;
		DBG_ASSERT(HERE, curr_re64 < itmp,"mi64_cvt_double_uint64: curr_wd64 !< (1<<curr_bits)");
		DBG_ASSERT(HERE, curr_im64 < itmp,"mi64_cvt_double_uint64: curr_wd64 !< (1<<curr_bits)");
		DBG_ASSERT(HERE, DNINT(a[j]) == a[j], "mi64_cvt_double_uint64: a[j] not pure-integer!");
		DBG_ASSERT(HERE, ABS(a[j]) < TWO54FLOAT, "mi64_cvt_double_uint64: |a[j]| >= maximum allowed value of 2^54!");

		itmp = (int64)a[j  ] + cyi;	/* current digit in int64 form, subtracting any borrow from previous digit.	*/
		if(itmp < 0)	/* If current digit < 0, add the base and set carry = -1	*/
		{
			itmp += FFT_MUL_BASE;
			cyi = -1;
		}
		else
		{
			cyi = 0;
		}
		DBG_ASSERT(HERE, itmp >= 0,"mi64_cvt_double_uint64: itmp < 0!");
		DBG_ASSERT(HERE, (curr_re64>>curr_bits) == 0,"mi64_cvt_double_uint64: (curr_re64>>curr_bits) != 0!");

		jtmp = (int64)a[j+1] + cyj;
		if(jtmp < 0)
		{
			jtmp += FFT_MUL_BASE;
			cyj = -1;
		}
		else
		{
			cyj = 0;
		}
		DBG_ASSERT(HERE, jtmp >= 0,"mi64_cvt_double_uint64: jtmp < 0!");
		DBG_ASSERT(HERE, (curr_im64>>curr_bits) == 0,"mi64_cvt_double_uint64: (curr_re64>>curr_bits) != 0!");

		/* Copy bits of the current residue word into the accumulator, starting
		at the (curr_bits)th bit. The resulting total number of accumulated bits
		will be curr_bits += FFT_MUL_BITS. If this is > 64, move the low 64 bits of
		the accumulator (which are in curr_wd64) into the x-array, increment (len),
		and move the overflow bits (which are in itmp >> (curr_bits - 64) ) into curr_wd64.
		*/
	#if 0
		EXAMPLE: curr_bits=40, FFT_MUL_BITS = 53: lower (64-40)=24 bits of itmp get copied to
		curr_wd64 starting at bit 40, accum now has (40+53) = 93 > 64 bits, so dump curr_wd64
		into x[len++], move upper (93-64) = 29 bits of itmp into curr_wd64, reset curr_bits=29.
	#endif

		curr_re64 += (itmp<<curr_bits);
		curr_im64 += (jtmp<<curr_bits);

		curr_bits += FFT_MUL_BITS;		/* Ex: 40+=53 => 93 */
		if(curr_bits >= 64)
		{
			x[len  ] = curr_re64;
			y[len++] = curr_im64;

			nbits += 64;
			curr_bits -= 64;	/* # of bits now left in curr_wd64 */	/* Ex: 93-=64 => 29 */

			/* Ex: Upper 29 bits of itmp via (unsigned) 24-bit rightt-shift: */
			curr_re64 = (uint64)itmp >> (FFT_MUL_BITS - curr_bits);
			curr_im64 = (uint64)jtmp >> (FFT_MUL_BITS - curr_bits);
		}
	}

	/* Dump any remaining partially-filled word: */
	if(curr_bits)
	{
		x[len  ] = curr_re64;
		y[len++] = curr_im64;
		nbits += curr_bits;
	}

	DBG_ASSERT(HERE, nbits == n*FFT_MUL_BITS,"mi64_cvt_double_uint64: nbits == n*FFT_MUL_BASE!");
	DBG_ASSERT(HERE, len > 0                ,"mi64_cvt_double_uint64: len has overflowed!");

	return len;
}

/***************/

/* returns 1 if the multiword unsigned integer p is a base-z Fermat pseudoprime, 0 otherwise. The base is assumed to fit in a 64-bit unsigned scalar. */
uint32 mi64_pprimeF(const uint64 p[], uint64 z, uint32 len)
{
/* The allocatable-array version of this gave bizarre stack faults under MSVC, so use fixed size,
especially as this routine is too slow to be useful for inputs larger than a couple hundred words anyway:
*/
/*	static uint64 *y, *n, *zvec, *ppad, *prod, *tmp;	*/
	uint64 y[1024], n[1024], zvec[1024], ppad[2048], prod[2048], tmp[2048];
	uint64 flag;
	uint32 curr_len, retval, plen;
#define MI64_PRP_DBG 0
#if MI64_PRP_DBG
	uint128 q128, r128;
	uint192 q192, r192;
	char str0[STR_MAX_LEN], str1[STR_MAX_LEN];
#endif

	DBG_ASSERT(HERE, len <= 1024, "mi64_pprimeF: Max 1024 words allowed at present!");

/*
	y    = (uint64 *)calloc((  len), sizeof(uint64));
	n    = (uint64 *)calloc((  len), sizeof(uint64));
	zvec = (uint64 *)calloc((  len), sizeof(uint64));
	ppad = (uint64 *)calloc((2*len), sizeof(uint64));
	prod = (uint64 *)calloc((2*len), sizeof(uint64));
	tmp  = (uint64 *)calloc((2*len), sizeof(uint64));
*/
	mi64_clear(y, len);
	mi64_clear(n, len);
	mi64_clear(zvec, len);
	mi64_clear(ppad, 2*len);
	mi64_clear(prod, 2*len);
	mi64_clear(tmp, 2*len);

	y[0] = 1ull;						/* y = 1 */
	zvec[0] = z;						/* vector version of z for powering */
	mi64_set_eq(ppad, p, len);			/* ppad is zero-padded */
	mi64_set_eq(n   , p, len);
	mi64_sub_scalar(n, 1ull, n, len);		/* n = p-1 */
	curr_len = mi64_getlen(n, len);

	while(!mi64_iszero(n, curr_len))
	{
		flag = n[0] & 1;
		mi64_shrl(n, n, 1, curr_len);	/* n >>= 1, also update curr_len */
		curr_len = mi64_getlen(n, curr_len);

	#if MI64_PRP_DBG
		if(len == 2)
		{
			q128[0] = y[0]; q128.d1 = y[1];
			r128[0] = zvec[0]; r128.d1 = zvec[1];
			fprintf(stderr, "mi64_pprimeF: flag = %1d, y = %s, z = %s\n", (uint32)flag, &str0[convert_uint128_base10_char(str0, q128)], &str1[convert_uint128_base10_char(str1, r128)]);
		}
		else if(len == 3)
		{
			q192[0] = y[0]; q192.d1 = y[1]; q192.d2 = y[2];
			r192[0] = zvec[0]; r192.d1 = zvec[1]; r192.d2 = zvec[2];
			fprintf(stderr, "mi64_pprimeF: flag = %1d, y = %s, z = %s\n", (uint32)flag, &str0[convert_uint192_base10_char(str0, q192)], &str1[convert_uint192_base10_char(str1, r192)]);
		}
	#endif

		if(flag)
		{
			mi64_mul_vector(y, len, zvec, len, prod, &plen);	/* y*z */
			mi64_div(prod, ppad, tmp, prod, 2*len);				/* Discard quotient here, overwrite prod with remainder */
			DBG_ASSERT(HERE, mi64_getlen(prod, 2*len) <= len, "mi64_pprimeF: (y*z)%p illegal length");
			mi64_set_eq(y, prod, len);	/* y = (y*z)%p */
		}
		mi64_mul_vector(zvec, len, zvec, len, prod, &plen);		/* z^2 */
		mi64_div(prod, ppad, tmp, prod, 2*len);
		DBG_ASSERT(HERE, mi64_getlen(prod, 2*len) <= len, "mi64_pprimeF: (z^2)%p illegal length");
		mi64_set_eq(zvec, prod, len);	/* z = (z^2)%p */
		if(mi64_iszero(zvec, len))
		{
			retval=0;
		}
	}

	retval = mi64_cmpeq_scalar(y, 1ull, len);
/*
	free((void *)y);
	free((void *)n);
	free((void *)zvec);
	free((void *)ppad);
	free((void *)prod);
	free((void *)tmp);
*/
	return retval;
}

/****************/

/*
Slow bit-at-a-time method to get quotient q = x/y and remainder r = x%y, where all vectors are of length [len] and nonoverlapping.
We formally expect y[] to be a const, but for practical reasons we leave the const qualifier off it so we can use its shifted values during execution.
*/
void mi64_div(const uint64 x[], uint64 y[], uint64 q[], uint64 r[], uint32 len)
{
	int i, nshift;
	uint32 lz_x, lz_y, xlen, ylen, max_len;
#define MI64_DIV_DBG 0
#if MI64_DIV_DBG
	char s0[1024],s1[1024];
#endif
	DBG_ASSERT(HERE, x && y && q && r, "mi64_div: At least one of X, Y, Q, R is null!");
	DBG_ASSERT(HERE, x != y, "mi64_div: X and Y arrays overlap!");
	DBG_ASSERT(HERE, r != y, "mi64_div: Y and Rem arrays overlap!");
	DBG_ASSERT(HERE, q != x && q != y && q != r, "mi64_div: Quotient array overlaps one of X, Y ,Rem!");

	/* Init Q = 0; don't do similarly for R since we allow X and R to point to same array: */
	mi64_clear(q, len);
	/* And now find the actual lengths of the divide operands and use those for the computation: */
	xlen = mi64_getlen(x, len);
	ylen = mi64_getlen(y, len);
	DBG_ASSERT(HERE, ylen != 0, "mi64_div: divide by 0!");

	/* Zero any extra high words R may have set on entry (only really needed if X and R not the same array): */
	if(r != x) mi64_set_eq(r, x, xlen);

	max_len = MAX(xlen, ylen);
	if((xlen < ylen) || mi64_cmpult(x, y, max_len))
	{
		return;		/* If x < y, return x */
	}
	lz_x = mi64_leadz(x, max_len);
	lz_y = mi64_leadz(y, max_len);
	nshift = lz_y - lz_x;
	DBG_ASSERT(HERE, nshift >= 0, "mi64_div: nshift < 0");

	/* Left-justify the modulus y to match x's leading bit: */
	mi64_shl(y, y, nshift, max_len);	ylen = max_len;
	for(i = nshift; i >= 0; --i)
	{
	#if MI64_DIV_DBG
		printf("I = %3d: r = %s, y = %s\n", i,&s0[convert_mi64_base10_char(s0, r, max_len)],&s1[convert_mi64_base10_char(s1, y, max_len)]);
	#endif
		if(mi64_cmpuge(r, y, max_len))
		{
			DBG_ASSERT(HERE, xlen == max_len,"mi64_div: xlen != max_len");
			mi64_sub(r, y, r, max_len);	/* r -= ytmp */
			DBG_ASSERT(HERE, mi64_cmpult(r, y, max_len),"mi64_div: r >= ytmp");
			xlen = mi64_getlen(r, max_len);
			mi64_set_bit(q, i);
		}
		if(i > 0)		/* Want y returned unchanged from its input value */
		{
			mi64_shrl(y, y, 1, ylen);
			ylen = mi64_getlen(y, ylen);
		}
		max_len = MAX(xlen, ylen);
	}
}

/* Fast is-divisible-by-32-bit-scalar using Montgomery modmul and right-to-left modding: */
int mi64_is_div_by_scalar32(const uint32 x[], uint32 q, uint32 len)
{
	int i,j,nshift;
	uint32 dlen,qinv,tmp,cy;

	DBG_ASSERT(HERE, q > 0, "mi64_is_div_by_scalar32: 0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz32(q);
	if(nshift)
	{
		if(trailz32(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint32)2;
	for(j = 0; j < 3; j++)
	{
		qinv = qinv*((uint32)2 - q*qinv);
	}
	cy = (uint32)0;
	dlen = len+len;	/* Since are processing a uint64 array cast to uint32[], double the #words parameter */
	for(i = 0; i < dlen; ++i)
	{
		tmp  = x[i] - cy;
		/* Add q if had a borrow - Since we low=half multiply tmp by qinv below and q*qinv == 1 (mod 2^32),
		   can simply add 1 to tmp*qinv result instead of adding q to the difference here:
		*/
		cy = (cy > x[i]); /* Comparing this rather than (tmp > x[i]) frees up tmp for the multiply */
		/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
		tmp *= qinv;
		tmp += cy;
	#ifdef MUL_LOHI32_SUBROUTINE
		cy = __MULH32(q,tmp);
	#else
		MULH32(q,tmp, cy);
	#endif
	}
	return (cy == 0);
}

uint32	mi64_is_div_by_scalar32_x4(const uint32 x[], uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 len)
{
	int i,j,nshift0,nshift1,nshift2,nshift3;
	uint32 retval=0,dlen = len+len, qinv0,qinv1,qinv2,qinv3,tmp0,tmp1,tmp2,tmp3,cy0,cy1,cy2,cy3;
	uint32 xcur,trailx;

	DBG_ASSERT(HERE, q0 && q1 && q2 && q3, "mi64_is_div_by_scalar32_x4: 0 modulus!");
	if(q0 + q1 + q2 + q3 == 4) return TRUE;
	if(len == 0) return TRUE;

	trailx = trailz32(x[0]);

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift0 = trailz32(q0);
	nshift1 = trailz32(q1);
	nshift2 = trailz32(q2);
	nshift3 = trailz32(q3);

	q0 >>= nshift0;
	q1 >>= nshift1;
	q2 >>= nshift2;
	q3 >>= nshift3;

	qinv0 = (q0+q0+q0) ^ (uint32)2;
	qinv1 = (q1+q1+q1) ^ (uint32)2;
	qinv2 = (q2+q2+q2) ^ (uint32)2;
	qinv3 = (q3+q3+q3) ^ (uint32)2;
	for(j = 0; j < 3; j++)
	{
		qinv0 = qinv0*((uint32)2 - q0*qinv0);
		qinv1 = qinv1*((uint32)2 - q1*qinv1);
		qinv2 = qinv2*((uint32)2 - q2*qinv2);
		qinv3 = qinv3*((uint32)2 - q3*qinv3);
	}
	cy0 = (uint32)0;
	cy1 = (uint32)0;
	cy2 = (uint32)0;
	cy3 = (uint32)0;
	for(i = 0; i < dlen; ++i)
	{
		xcur = x[i];

		tmp0 = xcur - cy0;
		tmp1 = xcur - cy1;
		tmp2 = xcur - cy2;
		tmp3 = xcur - cy3;
		/* Add q if had a borrow - Since we low=half multiply tmp by qinv below and q*qinv == 1 (mod 2^32),
		   can simply add 1 to tmp*qinv result instead of adding q to the difference here:
		*/
		cy0 = (cy0 > xcur);
		cy1 = (cy1 > xcur);
		cy2 = (cy2 > xcur);
		cy3 = (cy3 > xcur);
		/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
		tmp0 *= qinv0;
		tmp1 *= qinv1;
		tmp2 *= qinv2;
		tmp3 *= qinv3;

		tmp0 += cy0;
		tmp1 += cy1;
		tmp2 += cy2;
		tmp3 += cy3;
	#ifdef MUL_LOHI32_SUBROUTINE
		cy0 = __MULH32(q0,tmp0);
		cy1 = __MULH32(q1,tmp1);
		cy2 = __MULH32(q2,tmp2);
		cy3 = __MULH32(q3,tmp3);
	#else
		MULH32(q0,tmp0, cy0);
		MULH32(q1,tmp1, cy1);
		MULH32(q2,tmp2, cy2);
		MULH32(q3,tmp3, cy3);
	#endif
	}
	retval += ((cy0 == 0) && (nshift0 <= trailx));
	retval += ((cy1 == 0) && (nshift1 <= trailx)) << 1;
	retval += ((cy2 == 0) && (nshift2 <= trailx)) << 2;
	retval += ((cy3 == 0) && (nshift3 <= trailx)) << 3;
	return retval;
}

uint32	mi64_is_div_by_scalar32_x8(const uint32 x[], uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 q4, uint32 q5, uint32 q6, uint32 q7, uint32 len)
{
	int i,j,nshift0,nshift1,nshift2,nshift3,nshift4,nshift5,nshift6,nshift7;
	uint32 retval=0,dlen = len+len, qinv0,qinv1,qinv2,qinv3,qinv4,qinv5,qinv6,qinv7,tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7;
	uint32 xcur,trailx;

	DBG_ASSERT(HERE, q0 && q1 && q2 && q3 && q4 && q5 && q6 && q7, "mi64_is_div_by_scalar32_x8: 0 modulus!");
	if(q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7 == 8) return TRUE;
	if(len == 0) return TRUE;

	trailx = trailz32(x[0]);

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift0 = trailz32(q0);						nshift4 = trailz32(q4);
	nshift1 = trailz32(q1);						nshift5 = trailz32(q5);
	nshift2 = trailz32(q2);						nshift6 = trailz32(q6);
	nshift3 = trailz32(q3);						nshift7 = trailz32(q7);

	q0 >>= nshift0;								q4 >>= nshift4;
	q1 >>= nshift1;								q5 >>= nshift5;
	q2 >>= nshift2;								q6 >>= nshift6;
	q3 >>= nshift3;								q7 >>= nshift7;

	qinv0 = (q0+q0+q0) ^ (uint32)2;				qinv4 = (q4+q4+q4) ^ (uint32)2;
	qinv1 = (q1+q1+q1) ^ (uint32)2;				qinv5 = (q5+q5+q5) ^ (uint32)2;
	qinv2 = (q2+q2+q2) ^ (uint32)2;				qinv6 = (q6+q6+q6) ^ (uint32)2;
	qinv3 = (q3+q3+q3) ^ (uint32)2;				qinv7 = (q7+q7+q7) ^ (uint32)2;
	for(j = 0; j < 3; j++)
	{
		qinv0 = qinv0*((uint32)2 - q0*qinv0);	qinv4 = qinv4*((uint32)2 - q4*qinv4);
		qinv1 = qinv1*((uint32)2 - q1*qinv1);	qinv5 = qinv5*((uint32)2 - q5*qinv5);
		qinv2 = qinv2*((uint32)2 - q2*qinv2);	qinv6 = qinv6*((uint32)2 - q6*qinv6);
		qinv3 = qinv3*((uint32)2 - q3*qinv3);	qinv7 = qinv7*((uint32)2 - q7*qinv7);
	}
	cy0 = (uint32)0;							cy4 = (uint32)0;
	cy1 = (uint32)0;							cy5 = (uint32)0;
	cy2 = (uint32)0;							cy6 = (uint32)0;
	cy3 = (uint32)0;							cy7 = (uint32)0;
	for(i = 0; i < dlen; ++i)
	{
		xcur = x[i];

		tmp0 = xcur - cy0;						tmp4 = xcur - cy4;
		tmp1 = xcur - cy1;						tmp5 = xcur - cy5;
		tmp2 = xcur - cy2;						tmp6 = xcur - cy6;
		tmp3 = xcur - cy3;						tmp7 = xcur - cy7;

		cy0 = (cy0 > xcur);						cy4 = (cy4 > xcur);
		cy1 = (cy1 > xcur);						cy5 = (cy5 > xcur);
		cy2 = (cy2 > xcur);						cy6 = (cy6 > xcur);
		cy3 = (cy3 > xcur);						cy7 = (cy7 > xcur);

		tmp0 *= qinv0;							tmp4 *= qinv4;
		tmp1 *= qinv1;							tmp5 *= qinv5;
		tmp2 *= qinv2;							tmp6 *= qinv6;
		tmp3 *= qinv3;							tmp7 *= qinv7;

		tmp0 += cy0;							tmp4 += cy4;
		tmp1 += cy1;							tmp5 += cy5;
		tmp2 += cy2;							tmp6 += cy6;
		tmp3 += cy3;							tmp7 += cy7;
	#ifdef MUL_LOHI32_SUBROUTINE
		cy0 = __MULH32(q0,tmp0);				cy4 = __MULH32(q4,tmp4);
		cy1 = __MULH32(q1,tmp1);				cy5 = __MULH32(q5,tmp5);
		cy2 = __MULH32(q2,tmp2);				cy6 = __MULH32(q6,tmp6);
		cy3 = __MULH32(q3,tmp3);				cy7 = __MULH32(q7,tmp7);
	#else
		MULH32(q0,tmp0, cy0);					MULH32(q4,tmp4, cy4);
		MULH32(q1,tmp1, cy1);					MULH32(q5,tmp5, cy5);
		MULH32(q2,tmp2, cy2);					MULH32(q6,tmp6, cy6);
		MULH32(q3,tmp3, cy3);					MULH32(q7,tmp7, cy7);
	#endif
	}
	retval += ((cy0==0)&&(nshift0<=trailx));	retval += ((cy4==0)&&(nshift4<=trailx))<<4;
	retval += ((cy1==0)&&(nshift1<=trailx))<<1;	retval += ((cy5==0)&&(nshift5<=trailx))<<5;
	retval += ((cy2==0)&&(nshift2<=trailx))<<2;	retval += ((cy6==0)&&(nshift6<=trailx))<<6;
	retval += ((cy3==0)&&(nshift3<=trailx))<<3;	retval += ((cy7==0)&&(nshift7<=trailx))<<7;
	return retval;
}


/* Fast is-divisible-by-64-bit scalar using Montgomery modmul and right-to-left modding: */
int mi64_is_div_by_scalar64(const uint64 x[], uint64 q, uint32 len)
{
	int i,j,nshift;
	uint64 qinv,tmp,cy;

	DBG_ASSERT(HERE, q > 0, "mi64_is_div_by_scalar64: 0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz64(q);
	if(nshift)
	{
		if(trailz64(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint64)2;
	for(j = 0; j < 4; j++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}
	cy = (uint64)0;
	for(i = 0; i < len; ++i)
	{
		tmp  = x[i] - cy;
		/* Add q if had a borrow - Since we low=half multiply tmp by qinv below and q*qinv == 1 (mod 2^64),
		   can simply add 1 to tmp*qinv result instead of adding q to the difference here:
		*/
		cy = (cy > x[i]); /* Comparing this rather than (tmp > x[i]) frees up tmp for the multiply */
		/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
		tmp *= qinv;
		tmp += cy;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy = __MULH64(q,tmp);
	#else
		MULH64(q,tmp, cy);
	#endif
	}
	return (cy == 0);
}

/*
Divide-with-Remainder of x by y, where x is a multiword base 2^64 integer and y a 32-bit unsigned scalar.
(Optionally) Returns quotient (x - x%y)/y in q, 32-bit remainder in the function result. Allows division-in-place, i.e. x == q.
*/
uint32 mi64_div_y32(uint64 x[], uint32 y, uint64 q[], uint32 len)
{
	int i;
	uint64 cy, rem, xlomody, tsum;
	static uint32 ysave = 0;
	static uint64 two64divy, two64mody;

	if(y != ysave)
	{
		ysave = y;
		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
	}

	rem = 0;
	for(i = len-1; i >= 0; --i)
	{
		/* Current-term remainder (must calculate this before modifying (q[i], in case q and x point to same array) */
		xlomody = (x[i])%y;
		tsum = rem*two64mody + xlomody;

		/* Low digit of result: we must separately divide (x[i]) by y
		(making sure to add (x[i])%y to  cy*two64mody first, so as not to drop a digit)
		because x[i] may be as large as 2^64-1, and adding rem*two64mody
		prior to dividing risks unsigned integer overflow:
		*/
		if(q)
		{
			q[i] = rem*two64divy + tsum/y + (x[i])/y;
		}
		rem = tsum%y;
	}
	if(rem == 0)
	{
		DBG_ASSERT(HERE, mi64_is_div_by_scalar32(x, y, len), "Results of mi64_div_y32 and mi64_is_div_by_scalar32 differ!");
		return 0;
	}
	return (uint32)rem;
}

/*
Returns decimal character representation of a base-2^64 multiword unsigned int in char_buf,
and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
left-justified form) in the function result.
*/
int	convert_mi64_base10_char(char char_buf[], const uint64 x[], uint32 len)
{
	uint32 MAX_DIGITS;
	uint32 i, curr_len, n_dec_digits = 0;
	char c;
	const double log10_base = 19.26591972249479649367928926;	/* log10(2^64) */
	const double ln10 = log(10.0);
	uint64 *temp;

	/* Estimate # of decimal digits: */
	curr_len = mi64_getlen(x, len);	/* this checks that len > 0; need at least one digit, even if it = 0. curr_len guaranteed > 0. */
	curr_len = MAX(curr_len, 1);
	temp = (uint64 *)calloc(curr_len, sizeof(uint64));
	mi64_set_eq(temp, x, curr_len);
	MAX_DIGITS = ceil( (curr_len-1)*log10_base + log((double)x[curr_len-1])/ln10 );
	MAX_DIGITS = MAX(MAX_DIGITS, 1);

	DBG_ASSERT(HERE, MAX_DIGITS < STR_MAX_LEN, "convert_mi64_base10_char: Output string overflows buffer");
	char_buf[MAX_DIGITS-1]='0';
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i=0; i < MAX_DIGITS; i++)
	{
		/* Only print leading zero if q = 0, in which case we print a right-justifed single zero: */
		/* Since the x***_div_y32 routines return the mod *and* the divided input,
		   don't call the function until *after* performing the if() test:
		*/
		if(!mi64_iszero(temp, curr_len) || n_dec_digits == 0)
		{
			c = mi64_div_y32(temp, (uint32)10, temp, curr_len) + CHAROFFSET;
			curr_len = mi64_getlen(temp, curr_len);
			n_dec_digits++;
		}
		else
			c = ' ';

		char_buf[(MAX_DIGITS - 1) - i] = c;
	}
	free((void *)temp);
	return (int)MAX_DIGITS-n_dec_digits;
}

/****************/

/* This is essentially an mi64 version of the <stdlib.h> strtoul function: Takes an input string,
converts it - if numeric - to mi64 form, returning a pointer to the mi64 array allocated to
store the result in the function value, and the nominal length of the result in terms of 64-bit
words via argument #2.

If the input string is not a valid numeric or some other conversion problem is encountered,
both the
*/
//uint32 convert_base10_char_mi64(const char*char_buf, uint64 **out_array)
uint64 *convert_base10_char_mi64(const char*char_buf, uint32 *len)
{
	uint64 *mi64_vec;
	uint32 LEN_MAX;
	uint64 tmp = 0;
	uint32 i, imin, imax;
	char c;
	uint64 curr_digit;
	const double log10_base = 19.26591972249479649367928926;	/* log10(2^64) */

	/* Read the decimal digits from the string from left to right,
	skipping any leading whitespace, and stopping if either non-leading
	whitespace or '\0' is encountered:
	*/
	imax = strlen(char_buf);
	for(i=0; i < imax; i++)
	{
		c = char_buf[i];
		if(!isspace(c))
		{
			break;
		}
	}
	/* Estimate # of 64-bit words based on length of non-whitespace portion of the string: */
	*len = 1;
	LEN_MAX = (uint32)ceil( (imax-i)/log10_base );
	mi64_vec = (uint64 *)calloc(LEN_MAX+1, sizeof(uint64));	/* 01/09/2009: Add an extra zero-pad element here as workaround for bug in mi64_div called with differing-length operands */
	imin = i;
	for(i=imin; i < imax; i++)
	{
		c = char_buf[i];
		if(!isdigit(c))
		{
			free((void *)mi64_vec);	*len = 0;	return 0x0;
		}
		curr_digit = (uint64)(c - CHAROFFSET);
		DBG_ASSERT(HERE, curr_digit < 10,"util.c: curr_digit < 10");
		/* currsum *= 10, and check for overflow: */
		tmp = mi64_mul_scalar(mi64_vec, (uint64)10, mi64_vec, *len);
		if(tmp != 0)
		{
			if(*len == LEN_MAX)
			{
				fprintf(stderr, "ERROR: Mul-by-10 overflows in convert_base10_char_mi64: Offending input string = %s\n", char_buf);
				DBG_ASSERT(HERE, 0,"0");
			}
			mi64_vec[(*len)++] = tmp;
		}

		*len += mi64_add_scalar(mi64_vec, curr_digit, mi64_vec, *len);
		DBG_ASSERT(HERE, *len <= LEN_MAX,"len <= LEN_MAX");
	}
	*len = LEN_MAX;	/* Nominal length, so user knows how much memory was allocated */
	return mi64_vec;
}

/****************/

/*
Function to find 2^(-p) mod q, where p and q are both base-2^64 multiword unsigned integers.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^len. Result (optionally) returned in res[],
assumed to have allocated length at least as large as q[]:

The key 3-operation sequence here is as follows:

	SQR_LOHI(x,lo,hi);	// Input x has len words, lo/hi have len words each
	MULL(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have len words
	MULH(q,lo,lo);	// Inputs q & lo, and output (overwrites lo) have len words.
*/
uint32 mi64_twopmodq(const uint64 p[], uint32 len_p, uint64 q[], uint32 len, uint64*res)
{
	 int32 j, lenp_save = -1, lenq_save = -1;	/* Current-Bit index j needs to be signed because of the LR binary exponentiation. */
	uint32 idum, pbits;
	uint64 lead_chunk, lo64, cyout;
	static uint64 *qhalf = 0x0, *qinv = 0x0, *x = 0x0, *lo = 0x0, *hi = 0x0;
	static uint64 *psave = 0x0, *pshift = 0x0;
	static uint32 lenP, lenQ, qbits, log2_numbits, start_index, zshift, first_entry = TRUE;
#define MI64_POW_DBG	0
#if MI64_POW_DBG
	char s0[1024];
#endif

	lenP = mi64_getlen(p, len_p);
	DBG_ASSERT(HERE, lenP > 0, "mi64_twopmodq: 0 exponent");
	lenQ = mi64_getlen(q, len);
	DBG_ASSERT(HERE, lenQ > 0, "mi64_twopmodq: 0 modulus!");
	if(len_p != lenp_save)
	{
		lenp_save = len_p;
		free((void *)psave );
		free((void *)pshift);
		psave  = (uint64 *)calloc((len_p  ), sizeof(uint64));
		pshift = (uint64 *)calloc((len_p+1), sizeof(uint64));
	}
	if(len != lenq_save)
	{
		lenq_save = len;
		free((void *)psave );
		free((void *)pshift);
		free((void *)qhalf );
		free((void *)qinv  );
		free((void *)x     );
		free((void *)lo    );
		psave  = (uint64 *)calloc((len_p  ), sizeof(uint64));
		pshift = (uint64 *)calloc((len_p+1), sizeof(uint64));
		qhalf  = (uint64 *)calloc((lenQ   ), sizeof(uint64));
		qinv   = (uint64 *)calloc((lenQ   ), sizeof(uint64));
		x      = (uint64 *)calloc((lenQ   ), sizeof(uint64));
		lo     = (uint64 *)calloc((2*lenQ ), sizeof(uint64));
		hi     = lo + lenQ;	/* Pointer to high half of double-wide product */
	}

#if MI64_POW_DBG
	printf("mi64_twopmodq: len = %u, lenQ = %u\n", len, lenQ);
#endif
	qbits = lenQ << 6;
	mi64_shrl(q, qhalf, 1, lenQ);	/* (q >> 1) = (q-1)/2, since q odd. */

	if(first_entry || !mi64_cmpeq(p, psave, len_p))
	{
		first_entry = FALSE;
		mi64_set_eq(psave, p, len_p);
		/* pshift = p + len*64 */
		pshift[lenP] = mi64_add_scalar(p, lenP*64, pshift, lenP);
		DBG_ASSERT(HERE, !pshift[lenP], "mi64_twopmodq: pshift overflows!");

#if MI64_POW_DBG
	printf("mi64_twopmodq: lenP = %u, pshift = p+(64*len)\n", lenP);
#endif
		log2_numbits = ceil(log(1.0*qbits)/log(2.0));
	/*
	Find position of the leftmost ones bit in pshift, and subtract log2_numbits-1 or log2_numbits
	to account for the fact that we can do the powering for the leftmost log2_numbits-1 or log2_numbits
	bits (depending on whether the leftmost log2_numbits >= qbits or not) via a simple shift.
	*/
		/* Leftward bit at which to start the l-r binary powering, assuming
		the leftmost 7/8 bits have already been processed via a shift.

		Since 7 bits with leftmost bit = 1 is guaranteed
		to be in [64,127], the shift count here is in [0, 63].
		That means that zstart < 2^64. Together with the fact that
		squaring a power of two gives another power of two, we can
		simplify the modmul code sequence for the first iteration.
		Every little bit counts (literally in this case :), right?
		*/
	/* Extract leftmost log2_numbits bits of pshift (if >= qbits, use the leftmost log2_numbits-1) and subtract from qbits: */
		pbits = mi64_extract_lead64(pshift,len_p,&lo64);
		DBG_ASSERT(HERE, pbits >= log2_numbits, "mi64_twopmodq: leadz64!");
		if(pbits >= 64)
			lead_chunk = lo64>>(64-log2_numbits);
		else
			lead_chunk = lo64>>(pbits-log2_numbits);

		if(lead_chunk >= qbits)
		{
			lead_chunk >>= 1;
	#if MI64_POW_DBG
		printf("mi64_twopmodq: lead%u = %u\n", log2_numbits-1,lead_chunk);
	#endif
			start_index = pbits-(log2_numbits-1);	/* Use only the leftmost log2_numbits-1 bits */
		}
		else
		{
	#if MI64_POW_DBG
		printf("mi64_twopmodq: lead%u = %u\n", log2_numbits  ,lead_chunk);
	#endif
			start_index = pbits-log2_numbits;
		}

		zshift = (qbits-1) - lead_chunk;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		mi64_negl(pshift, pshift, lenP);	/* ~pshift[] */
	#if MI64_POW_DBG
		printf("mi64_twopmodq: pshift = %s\n", &s0[convert_mi64_base10_char(s0, pshift, lenP)]);
	#endif
	}

	/*
	Find modular inverse (mod 2^qbits) of q in preparation for modular multiply.
	q must be odd for Montgomery-style modmul to work.

	Init qinv = q. This formula returns the correct bottom 4 bits of qinv,
	and we double the number of correct bits on each of the subsequent iterations.
	*/
	DBG_ASSERT(HERE, (q[0] & (uint64)1) == 1, "mi64_twopmodq : q must be odd!");
	mi64_clear(qinv, lenQ);
	qinv[0] = (q[0] + q[0] + q[0]) ^ (uint64)2;

	/* Newton iteration involves repeated steps of form

		qinv = qinv*(2 - q*qinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as qinv_0 = 3*q ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 2; j < 6; j++)	/* At each step, have 2^j correct low-order bits in qinv */
	{
		lo64 = q[0]*qinv[0];
		qinv[0] = qinv[0]*((uint64)2 - lo64);
	}

	/* Now that have bottom 64 = 2^6 bits of qinv, do as many more Newton iterations as needed to get the full [qbits] of qinv: */
	for(j = 6; j < log2_numbits; j++)
	{
		mi64_mul_vector_lo_half(q, qinv, x, lenQ);
		mi64_nega              (x, x, lenQ);
		cyout = mi64_add_scalar(x, 2ull, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		mi64_mul_vector_lo_half(qinv, x, qinv, lenQ);
	}
	#if MI64_POW_DBG
		printf("mi64_twopmodq: q    = %s\n", &s0[convert_mi64_base10_char(s0, q   , lenQ)]);
		printf("mi64_twopmodq: qinv = %s\n", &s0[convert_mi64_base10_char(s0, qinv, lenQ)]);
		printf("start_index = %3d\n", start_index);
	#endif

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL192(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if MI64_POW_DBG
	printf("mi64_twopmodq: zshift = %d, lo = (qinv << zshift)...\n", zshift);
#endif
	mi64_shl(qinv, lo, zshift, lenQ);
#if MI64_POW_DBG
	printf("mi64_twopmodq: lo = %s\n", &s0[convert_mi64_base10_char(s0, lo, lenQ)]);
#endif
	mi64_mul_vector_hi_half(q, lo, lo, lenQ);
#if MI64_POW_DBG
	printf("mi64_twopmodq: q*lo/2^%u = %s\n", (lenQ<<6), &s0[convert_mi64_base10_char(s0, lo, lenQ)]);
#endif

	/* hi = 0 in this instance, which simplifies things. */
	cyout = mi64_sub(q, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");

	if(mi64_test_bit(pshift, j))
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(mi64_cmpugt(x, qhalf, lenQ))
		{
			cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			cyout = mi64_sub(x, q, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
		else
		{
			cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
	}
#if MI64_POW_DBG
	printf("mi64_twopmodq: x0 = %s\n", &s0[convert_mi64_base10_char(s0, x, lenQ)]);
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	#if MI64_POW_DBG
		printf("J = %3d: x = %s\n", j,&s0[convert_mi64_base10_char(s0, x, lenQ)]);
	#endif
		/*...x^2 mod q is returned in x. */
		mi64_mul_vector(x, lenQ, x, lenQ, lo, &idum);
		mi64_mul_vector_lo_half(lo, qinv, lo, lenQ);
		mi64_mul_vector_hi_half(q, lo, lo, lenQ);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(mi64_cmpult(hi, lo, lenQ))
		{
			cyout = mi64_sub(q, lo, lo, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			cyout = mi64_add(lo, hi, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
		else
		{
			cyout = mi64_sub(hi, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
#if MI64_POW_DBG
	if((j&0xf)==0)
	printf("");
#endif
		if(mi64_test_bit(pshift, j))
		{
			DBG_ASSERT(HERE, mi64_cmpult(x, q, lenQ), "mi64_twopmodq : x >= q");
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(mi64_cmpugt(x, qhalf, lenQ))
			{
				cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
				cyout = mi64_sub(x, q, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			}
			else
			{
				cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			}
	#if MI64_POW_DBG
		printf("2x = %s\n",&s0[convert_mi64_base10_char(s0, x, lenQ)]);
	#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");	/* In the case of interest, x = (q+1)/2, so x + x cannot overflow. */
	cyout = mi64_sub(x, q, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
	if(res != 0x0)
	{
		mi64_set_eq(res, x, lenQ);
	}
	return mi64_cmpeq_scalar(x, 1ull, lenQ);
}

