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

#define MI64_DEBUG	0

#include "mi64.h"
#include "align.h"
#include "qfloat.h"

#if MI64_DEBUG
	#include "Mdata.h"
	char s0[STR_MAX_LEN], s1[STR_MAX_LEN];
	char str_10k[10*1024];

//	#define YOU_WANT_IT_TO_BE_SLOW
#endif

#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

/*******************/

/* Set-X-equal-to-Y ... argument order may seem reversed, but leads to calls whose arg-order matches mnemonic "x = y": */
void	mi64_set_eq(uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	if(len && x != y) {
		for(i = 0; i < len; ++i)
		{
			x[i] = y[i];
		}
	}
}

/* Set-X-equal-to-scalar-A: */
void	mi64_set_eq_scalar(uint64 x[], const uint64 a, uint32 len)
{
	uint32 i;
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	if(len) {
		x[0] = a;
		for(i = 1; i < len; ++i)
		{
			x[i] = 0ull;
		}
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
Any off-shifted bits aside from the least-significant 64 are lost. No array-bounds checking is done.

Returns low 64 bits of off-shifted portion.
Allows in-place operation, i.e. x == y.
*/
uint64	mi64_shl(const uint64 x[], uint64 y[], uint32 nbits, uint32 len)
{
	int i;
	uint32 nwshift = (nbits >> 6), rembits = (nbits & 63), m64bits;
	uint64 lo64 = 0ull;

	DBG_ASSERT(HERE, len != 0, "mi64_shl: zero-length array!");
	/* Special-casing for shift-into-Bolivian (that's Mike-Tyson-ese for "oblivion"): */
	if(nwshift >= len)
	{
		// If nwshift == len, save low word in lo64...
		if(nwshift == len)
		{
			lo64 = x[0];
		}
		// ... because after this next line, we can only distinguish between (nwhsift < len) and (nwshift >= len):
		nwshift = len;
		rembits = 0;
	}

	/* Take care of the whole-word part of the shift: */
	for(i = len-1; i >= (int)nwshift; i--)
	{
		y[i] = x[i-nwshift];
	}
	// lo64 plays the part of "y[len]"; see above comment for why we deal with (nwshift >= len) above rather than here:
	if(nwshift && (nwshift < len)) {
		lo64 = x[len-nwshift];
	}
	for(i = nwshift-1; i >= 0; i--)
	{
		y[i] = 0ull;
	}

	/* If nbits not an exact multiple of the wordlength, take care of remainder: */
	if(rembits)
	{
		m64bits = (64-rembits);
		/* Process all but the least-significant element, in reverse order: */
		for(i = len-1; i > 0; i--)
		{
			y[i] = (y[i] << rembits) + (y[i-1] >> m64bits);
		}
		// lo64 plays the part of "y[len]":
		lo64 = (lo64 << rembits) + (y[len-1] >> m64bits);
		/* Least-significant element gets zeros shifted in from the right: */
		y[0] <<= rembits;
	}
	return lo64;
}

/*******************/

/*
Logical-right-shift a base-2^64 int x[] by (nbits) bits, returning the result in y[].
Any off-shifted bits aside from the most-significant 64 are lost. No array-bounds checking is done.

Returns high 64 bits of off-shifted portion.
Allows in-place operation, i.e. x == y.
*/
uint64	mi64_shrl(const uint64 x[], uint64 y[], uint32 nbits, uint32 len)
{
	int i;
	uint32 nwshift = (nbits >> 6), rembits = (nbits & 63), m64bits;
	uint64 hi64 = 0ull;

	DBG_ASSERT(HERE, len != 0, "mi64_shrl: zero-length array!");
	/* Special-casing for shift-into-Bolivian (that's Mike-Tyson-ese for "oblivion"): */
	if(nwshift >= len)
	{
		// If nwshift == len, save high word in hi64...
		if(nwshift == len)
		{
			hi64 = x[len-1];
		}
		// ... because after this next line, we can only distinguish between (nwhsift < len) and (nwshift >= len):
		nwshift = len;
		rembits = 0;
	}

	// hi64 plays the part of "y[-1]"; see above comment for why we deal with (nwshift >= len) above rather than here:
	// Must do this *before* the loop because of in-place possibility:
	if(nwshift && (nwshift < len)) {
		hi64 = x[nwshift-1];
	}
	/* Take care of the whole-word part of the shift: */
	for(i = 0; i < len-nwshift; i++)
	{
		y[i] = x[i+nwshift];
	}
	for(i = len-nwshift; i < len; i++)
	{
		y[i] = 0ull;
	}

	/* If nbits not an exact multiple of the wordlength, take care of remaining shift bits: */
	if(rembits)
	{
		m64bits = (64-rembits);
		hi64 = (hi64 >> rembits) + (y[0] << m64bits);
		/* Process all but the most-significant element, in reverse order: */
		for(i = 0; i < len-1; i++)
		{
			y[i] = (y[i] >> rembits) + (y[i+1] << m64bits);
		}
		/* Most-significant element gets zeros shifted in from the left: */
		y[len-1] >>= rembits;
	}
	return hi64;
}

/* Specialized 1-bit left and rightward shifts, also allow in-place: */
uint64	mi64_mul2(const uint64 x[], uint64 y[], uint32 len)
{
	uint32 i;
	uint64 cy = 0;
	if(len) {
		cy = (int64)x[len-1] < 0;
		for(i = len-1; i != 0; --i)
		{
			y[i] = (x[i] << 1) + ((int64)x[i-1] < 0);
		}
		y[0] = x[0] << 1;
	}
	return cy;
}

uint64	mi64_div2(const uint64 x[], uint64 y[], uint32 len)
{
	uint32 i;
	uint64 cy = 0;
	if(len) {
		cy = x[0] & 1ull;
		for(i = 0; i < len-1; ++i)
		{
			y[i] = (x[i] >> 1) + (x[i+1] << 63);
		}
		y[len-1] = x[len-1] >> 1;
	}
	return cy;
}


/*******************/
/* unsigned compare: these are all built from just two elemental functions: < and == */
uint32	mi64_cmpult(const uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;
	// Need hard-assert here due to zero-element default compare:
	ASSERT(HERE, len != 0, "mi64_cmpult: zero-length array!");
	for(i = len-1; i !=0 ; i--)	/* Loop over all but the 0 elements while equality holds.... */
	{
		if(x[i] < y[i]) {
			return TRUE;
		} else if(x[i] > y[i]) {
			return FALSE;
		}
	}
	return x[0] < y[0];
}

uint32	mi64_cmp_eq(const uint64 x[], const uint64 y[], uint32 len)
{
	uint32 i;
	// Allow for zero-length here with default return TRUE,
	// according to the convention that a zero-length mi64 object = 0:
	DBG_ASSERT(HERE, len != 0, "mi64_cmp_eq: zero-length array!");	// DBG_ allows us to catch zero-length cases in debug build & test
	for(i = 0; i < len; i++)
	{
		if(x[i] != y[i])
			return FALSE;
	}
	return TRUE;
}

uint32	mi64_cmplt_scalar(const uint64 x[], uint64 a, uint32 len)
{
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	return ( (mi64_getlen(x, len) <= 1) && (x[0] < a) );
}

uint32	mi64_cmpgt_scalar(const uint64 x[], uint64 a, uint32 len)
{
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	return ( (x[0] > a) || (mi64_getlen(x, len) > 1) );
}

uint32	mi64_cmp_eq_scalar(const uint64 x[], uint64 a, uint32 len)
{
	DBG_ASSERT(HERE, len != 0, "mi64_cmp_eq_scalar: zero-length array!");
	return ( (x[0] == a) && (mi64_getlen(x+1, len-1) == 0) );
}

/*******************/

/* Returns number of trailing zeros of a base-2^64 multiword int x if x != 0, 0 otherwise. */
uint32	mi64_trailz(const uint64 x[], uint32 len)
{
	uint32 i, tz = 0;
	DBG_ASSERT(HERE, len != 0, "mi64_trailz: zero-length array!");
	for(i = 0; i < len; i++, tz += 64)
	{
		if(x[i]) {
			return tz + trailz64(x[i]);
		}
	}
	return 0;
}

/*******************/

/* Returns number of leading zeros of a base-2^64 multiword int x, defaulting to len*64 if x == 0. */
uint32	mi64_leadz(const uint64 x[], uint32 len)
{
	int i;	// Loop index signed here
	uint32 lz = 0;
	DBG_ASSERT(HERE, len != 0, "mi64_leadz: zero-length array!");
	for(i = len-1; i >= 0; i--, lz += 64)
	{
		if(x[i]) {
			return lz + leadz64(x[i]);
		}
	}
	return lz;
}

/*******************/

/* Pair of functions to extract leading 64 significant bits, and 128-most-significant-bits-starting-at-a-specified-bit-position.
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
	nwshift = (nshift >> 6);
	rembits = (nshift & 63);
	/* shift-word count may == len, but only if x[] = 0: */
	if(nwshift >= len)
	{
		DBG_ASSERT(HERE, nwshift == len, "mi64_extract_lead64: nwshift out of range!");
		DBG_ASSERT(HERE, mi64_iszero(x, len), "mi64_extract_lead64: expected zero-valued array!");
		*result = 0ull;
	} else {
		i = len-1-nwshift;
		if(rembits) {
			*result  = (x[i] << rembits);
			if(i) {
				*result += (x[i-1] >> (64-rembits));
			}
		} else {
			*result = x[i];
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
	if(lead64 == 0ull) {
		return 0.0;
	}
	DBG_ASSERT(HERE,(lead64 >> 63) == 1ull, "mi64_cvt_double: lead64 lacks leftmost ones bit!");
	/*  round based on 1st neglected bit: */
	lead64_rnd = (lead64 >> 11) + ((lead64 >> 10) & 0x0000000000000001ull);
	/* exponent: */
	itmp64 = (((uint64)0x3FD + (uint64)pow2) << 52);
	/* Add in mantissa, with hidden bit made explicit, hence the 0x3FD (rather than 0x3FE) initializer */
	itmp64 += lead64_rnd;
	ASSERT(HERE, itmp64 > lead64_rnd , "mi64_cvt_double: Exponent overflows IEEE64 field");
	/* GCC bug: needed to add the explicit sign-check below, otherwise GCC 'optimizes' away the (*(double *)&itmp64): */
	retval = *(double *)&itmp64;
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
		if(x[i]) {
			return i+1;
		}
	}
	return 0;
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
#ifndef YES_ASM

	uint32 i;
	uint64 tmp, cy = 0;
	DBG_ASSERT(HERE, len != 0, "mi64_add: zero-length array!");

	for(i = 0; i < len; i++)
	{
		tmp = x[i] + cy;
		cy  = (tmp < x[i]);
		tmp = tmp + y[i];
		cy += (tmp < y[i]);
		z[i] = tmp;
	}

#elif 1//!defined(YES_ASM)	// 2-unrolled C loop gives a nice boost. (4-unrolled no better).

	uint64 tmp, cy, c2 = 0;
	uint32 i, odd = (len&1), len2 = len - odd;

	for(i = 0; i < len2; i++){
		tmp = x[i] + y[i];
		cy = (tmp < x[i]);
		z[i] = tmp + c2;
		cy += (tmp > z[i]);
		i++;

		tmp = x[i] + y[i];
		c2 = (tmp < x[i]);
		z[i] = tmp + cy;
		c2 += (tmp > z[i]);
	}

	if(odd) {
		tmp = x[i] + y[i];
		cy = (tmp < x[i]);
		z[i] = tmp + c2;
		cy += (tmp > z[i]);
	} else {
		cy = c2;
	}

#else	// defined(YES_ASM)

  #if MI64_DEBUG
	#define DEBUG_ADD_ASM	0
  #endif

  #if DEBUG_ADD_ASM
	uint32 i;
	uint64 tmp, c2 = 0, *u = 0x0;
	// Compute result using above 'nondestructive' way, to compare with ASM result later:
	u = (uint64 *)calloc(len+1, sizeof(uint64));
	for(i = 0; i < len; i++)
	{
		tmp = x[i] + c2;
		c2  = (tmp < x[i]);
		tmp = tmp + y[i];
		c2 += (tmp < y[i]);
		u[i] = tmp;
	}
  #endif

	uint64 cy = 0;
	/* x86_64 ASM implementation of the add/carry loop: */
	__asm__ volatile (\
		"movq	%[__x0],%%rax	\n\t"/* &x[0] */\
		"movq	%[__y0],%%rbx	\n\t"/* &y[0] */\
		"movq	%[__z0],%%rdx	\n\t"/* &z[0] */\
		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"xorq   %%rsi, %%rsi    \n\t"/* Index into the 3 arrays (really a uint64-array pointer offset) */\
		"movq	%%rcx,%%rdi	\n\t"/* Copy of array-length in rdi */\
		"shrq	$2,%%rcx	\n\t"/* How many iterations thru the 4-way loop */\
		"andq	$3,%%rdi	\n\t"/* Clear CF. Prepare for middle of loop jump */\
		"jz     1f		\n\t"/* Woohoo! Perfect multiple of 4 */\
		"incq   %%rcx		\n\t"/* Compensate for jumping into the middle of the loop */\
		"decq   %%rdi		\n\t"/* residue 1? */\
		"jz	4f		\n\t"/* Just do 1 addition in the first partial */\
		"decq   %%rdi		\n\t"/* residue 2? */\
		"jz	3f		\n\t"\
		"jmp	2f		\n\t"/* residue 3 */\
	"1:					\n\t"\
		"movq	(%%rax,%%rsi),%%r8 	\n\t"\
		"adcq	(%%rbx,%%rsi),%%r8	\n\t"\
		"movq	%%r8,(%%rdx,%%rsi)	\n\t"\
		"leaq   0x8(%%rsi), %%rsi	\n\t"\
	"2:					\n\t"\
		"movq	(%%rax,%%rsi),%%r8 	\n\t"\
		"adcq	(%%rbx,%%rsi),%%r8	\n\t"\
		"movq	%%r8,(%%rdx,%%rsi)	\n\t"\
		"leaq   0x8(%%rsi), %%rsi	\n\t"\
	"3:					\n\t"\
		"movq	(%%rax,%%rsi),%%r8 	\n\t"\
		"adcq	(%%rbx,%%rsi),%%r8	\n\t"\
		"movq	%%r8,(%%rdx,%%rsi)	\n\t"\
		"leaq   0x8(%%rsi), %%rsi	\n\t"\
	"4:					\n\t"\
		"movq	(%%rax,%%rsi),%%r8 	\n\t"\
		"adcq	(%%rbx,%%rsi),%%r8	\n\t"\
		"movq	%%r8,(%%rdx,%%rsi)	\n\t"\
		"leaq   0x8(%%rsi), %%rsi	\n\t"\

	"decq	%%rcx \n\t"\
	"jnz 1b 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"adcq	%%rcx,%[__cy]	\n\t"/* Carryout. RCX is guaranteed to be zero at this point */\
		: [__cy] "=m" (cy) /* outputs: cy */\
		: [__x0] "g" (x)	/* All inputs from memory/register here */\
		 ,[__y0] "g" (y)	\
		 ,[__z0] "g" (z)	\
		 ,[__len] "g" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8"	/* Clobbered registers */\
	);

  #if DEBUG_ADD_ASM
	if(!mi64_cmp_eq(u,z,len) || (cy != c2)) {
		for(i = 0; i < len; i++)
		{
			if(u[i] != z[i]) printf("i = %u Error: U = %20llu, Z = %20llu, Diff = %20lld\n",i,u[i],z[i],(int64)(u[i]-z[i]) );
		}
		if(cy != c2) printf("Carry Error: c2 = %20llu, cy = %20llu, Diff = %20lld\n",c2,cy,(int64)(c2-cy) );
		ASSERT(HERE, 0, "mi64_add ASM result incorrect!");
	}
	free((void *)u);
  #endif
  #undef DEBUG_ADD_ASM

#endif

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
	for(i = 0; i < len; i++)
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

/* Arithmetic negation uses that -x = ~x + 1: */
void	mi64_nega(const uint64 x[], uint64 y[], uint32 len)
{
	if(len) {
		mi64_negl(x,y,len);
		mi64_add_scalar(y,1,y,len);
	}
}

/* Logical negation: */
void	mi64_negl(const uint64 x[], uint64 y[], uint32 len)
{
	uint32 i;
	for(i = 0; i < len; i++)
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
		for(i = 0; i < len; i++)
		{
			y[i] = x[i] + cy;
			cy = (y[i] < cy);
			if(!cy) {
				break;
			}
		}
	}
	else
	{
		for(i = 0; i < len; i++)
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
		for(i = 0; i < len; i++)
		{
			tmp = x[i] - bw;
			/*  Since x[i] and y[i] point to the same memloc, need an extra temp here due to asymmetry of subtract: */
			bw = (tmp > x[i]);
			y[i] = tmp;	// This is really assigning to x[], but saying so much gives "error: assignment of read-only location"
			if(!bw) {
				break;
			}
		}
	}
	else
	{
		for(i = 0; i < len; i++)
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
	uint32 i = 0;
	uint64 lo, hi, cy = 0;

#define PIPELINED	1
#if PIPELINED
	uint32 lmod4 = (len&0x3), ihi = len - lmod4;
	for(; i < ihi; ++i)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		y[i] = lo + cy;
		cy = hi + (y[i] < lo);
		i++;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		y[i] = lo + cy;
		cy = hi + (y[i] < lo);
		i++;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		y[i] = lo + cy;
		cy = hi + (y[i] < lo);
		i++;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		y[i] = lo + cy;
		cy = hi + (y[i] < lo);
	}
	// Cleanup loop for remaining terms:
#endif
#undef PIPELINED
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	for(; i < len; i++)
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
	Unsigned multiply (vector * scalar) and add result to vector, a*X + Y = Z .
	Return any exit carry, rather than automatically storing it in a hypothetical
	[len+1]st array slot - up to user to determine what to do if this is nonzero.

	There is no code to check for special values of a, e.g. 0 or 1.

	Allows in-place, i.e. x == y.
*/
#if MI64_DEBUG
	#define MI64_MSAV2	0
#endif
uint64	mi64_mul_scalar_add_vec2(const uint64 x[], uint64 a, const uint64 y[], uint64 z[], uint32 len)
{
	uint64 cy = 0;

#if MI64_MSAV2
	uint64 *u = 0x0, *v = 0x0, c2;
	u = (uint64 *)calloc(len, sizeof(uint64));
	v = (uint64 *)calloc(len, sizeof(uint64));		memcpy(v,y,(len<<3));	// Save copy of x[]
	c2  = mi64_mul_scalar(x, a, u, len);
	c2 += mi64_add(u, y, u, len);
#endif

#if 1//ndef YES_ASM	// Toggle for x86_64 inline ASM

	uint32 i = 0;
	uint64 lo, hi, tmp;	// Oddly, the tmp-involving sequence below runs faster than the tmp-less one used in the ASM
	uint32 lmod2 = (len&0x1), ihi = len - lmod2;	// Unlike mi64_mul_scalar, 2x-unrolled works best here
	for(; i < ihi; ++i)			// "Oddly" enough, using a for-loop for cleanup rather than if(odd) is best
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		tmp = lo + cy;
		cy	= hi + (tmp < lo);
		z[i] = tmp + y[i];
		cy += (z[i] < tmp);
		i++;
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		tmp = lo + cy;
		cy	= hi + (tmp < lo);
		z[i] = tmp + y[i];
		cy += (z[i] < tmp);
	}
	// Cleanup loop for remaining terms:
	DBG_ASSERT(HERE, len != 0, "zero-length array!");
	for(; i < len; i++)
	{
	#ifdef MUL_LOHI64_SUBROUTINE
		MUL_LOHI64(a, x[i],&lo,&hi);
	#else
		MUL_LOHI64(a, x[i], lo, hi);
	#endif
		tmp = lo + cy;
		cy	= hi + (tmp < lo);
		z[i] = tmp + y[i];
		cy += (z[i] < tmp);
		/*
		A*x[i] is at most (2^64 - 1)^2, i.e. hi <= 2^64 - 2. Since the
		carry out of the low 64 bits of the 128-bit add is at most 1,
		cy + hi <= 2^64 - 1, hence we don't need to check for a carry there,
		only when adding the result to x[i+1] on the next pass through the loop:
		*/
	}

#else	// defined(YES_ASM)

	/* x86_64 ASM implementation of the a*X + Y = Z loop: */
	__asm__ volatile (\
		"movq	%[__a],%%rbx	\n\t"/* Constant scalar multiplier*/\
		"movq	%[__x0],%%r11	\n\t"/* &x[0] */\
		"movq	%[__y0],%%r12	\n\t"/* &y[0] */\
		"movq	%[__z0],%%r13	\n\t"/* &z[0] */\
		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"xorq	%%r9 , %%r9 	\n\t"/* Clear the '64-bit carry' */\
		"xorq	%%rsi, %%rsi	\n\t"/* Index into the 3 arrays (really a uint64-array pointer offset) */\
		"movq	%%rcx,%%rdi	\n\t"/* Copy of array-length in rdi */\
		"shrq	$2,%%rcx	\n\t"/* How many iterations thru the 4-way loop */\
		"andq	$3,%%rdi	\n\t"/* Clear CF. Prepare for middle of loop jump */\
		"jz	 1f		\n\t"/* Woohoo! Perfect multiple of 4 */\
		"incq	%%rcx		\n\t"/* Compensate for jumping into the middle of the loop */\
		"decq	%%rdi		\n\t"/* residue 1? */\
		"jz	4f		\n\t"/* Just do 1 addition in the first partial */\
		"decq	%%rdi		\n\t"/* residue 2? */\
		"jz	3f		\n\t"\
		"jmp	2f		\n\t"/* residue 3 */\
	"1:					\n\t"\
		"movq	(%%r11,%%rsi),%%rax \n\t"/* x[i] */\
		"mulq	%%rbx				\n\t"/* MUL_LOHI64(a, x[i]); lo:hi in rax:rdx; CF set according to (hi == 0) */\
		"addq	(%%r12,%%rsi),%%rax	\n\t"/* z[i] = lo + y[i], and set CF if overflow */\
		"adcq	$0   ,%%rdx 		\n\t"\
		"adcq	%%r9 ,%%rax			\n\t"/* lo += 64-bit carryin from pvs word, and set CF if overflow */\
		"adcq	$0   ,%%rdx 		\n\t"/* Carryout into next-higher word = hi + CF */\
		"movq	%%rdx,%%r9 			\n\t"/* Save carryout in r9 for next-higher word computation */\
		"movq	%%rax,(%%r13,%%rsi)	\n\t"/* Store z[i] */\
		"leaq	0x8(%%rsi), %%rsi	\n\t"\
	"2:					\n\t"\
		"movq	(%%r11,%%rsi),%%rax \n\t"\
		"mulq	%%rbx				\n\t"\
		"addq	(%%r12,%%rsi),%%rax	\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"adcq	%%r9 ,%%rax			\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"movq	%%rdx,%%r9 			\n\t"\
		"movq	%%rax,(%%r13,%%rsi)	\n\t"\
		"leaq	0x8(%%rsi), %%rsi	\n\t"\
	"3:					\n\t"\
		"movq	(%%r11,%%rsi),%%rax \n\t"\
		"mulq	%%rbx				\n\t"\
		"addq	(%%r12,%%rsi),%%rax	\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"adcq	%%r9 ,%%rax			\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"movq	%%rdx,%%r9 			\n\t"\
		"movq	%%rax,(%%r13,%%rsi)	\n\t"\
		"leaq	0x8(%%rsi), %%rsi	\n\t"\
	"4:					\n\t"\
		"movq	(%%r11,%%rsi),%%rax \n\t"\
		"mulq	%%rbx				\n\t"\
		"addq	(%%r12,%%rsi),%%rax	\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"adcq	%%r9 ,%%rax			\n\t"\
		"adcq	$0   ,%%rdx 		\n\t"\
		"movq	%%rdx,%%r9 			\n\t"\
		"movq	%%rax,(%%r13,%%rsi)	\n\t"\
		"leaq	0x8(%%rsi), %%rsi	\n\t"\

	"decq	%%rcx \n\t"\
	"jnz 1b 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"adcq	%%r9 ,%[__cy]	\n\t"/* Carryout. */\
		: [__cy] "=m" (cy) /* outputs: cy */\
		: [__x0] "g" (x)	/* All inputs from memory/register here */\
		 ,[__y0] "g" (y)	\
		 ,[__z0] "g" (z)	\
		 ,[__a] "g" (a)	\
		 ,[__len] "g" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r9","r11","r12","r13"	/* Clobbered registers */\
	);

#endif	// YES_ASM

#if MI64_MSAV2
	if(!mi64_cmp_eq(u,z,len) || (cy != c2)) {
		for(i = 0; i < len; i++)
		{
		//	if(u[i] != z[i])
			printf("i = %u Error: U = %20llu, Z = %20llu, Diff = %20lld\n",i,u[i],z[i],(int64)(u[i]-z[i]) );
		}
		if(cy != c2) printf("Carry Error: c2 = %20llu, cy = %20llu, Diff = %20lld\n",c2,cy,(int64)(c2-cy) );
		ASSERT(HERE, 0, "mi64_add ASM result incorrect!");
	}
	free((void *)u);
#endif
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
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u;
	static uint32 dimU = 0;

	ASSERT(HERE, x && y && z, "mi64_mul_vector: Null array x/y/z!");
	ASSERT(HERE, lenX != 0, "mi64_mul_vector: zero-length X-array!");
	ASSERT(HERE, lenY != 0, "mi64_mul_vector: zero-length Y-array!");
	ASSERT(HERE, x != z, "mi64_mul_vector: X and Z point to same array object!");
	ASSERT(HERE, y != z, "mi64_mul_vector: Y and Z point to same array object!");
	ASSERT(HERE, lenZ != 0x0, "mi64_mul_vector: Null lenZ pointer!");

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
		dimU = (lenA+1);
		free((void *)u);
		u = (uint64 *)calloc((lenA+1), sizeof(uint64));
	}

	/* Loop over remaining (lenB-1) elements of B[], multiplying A by each, and
	using u[] as a scratch array to store B[i]*A[] prior to adding to z[]: */
	for(j = 0; j < lenB; j++)
	{
		u[lenA] = mi64_mul_scalar(A, B[j], u, lenA);

		/* Add j-word-left-shifted u[] to z[]: */
		z[lenA+j] = u[lenA] + mi64_add(&z[j], u, &z[j], lenA);
	}
	*lenZ += lenB;

	/* Return actual length of result vector in the function argument, so that if one or
	more leading terms of the result is zero, caller can adjust vector length accordingly:
	*/
	*lenZ = mi64_getlen(z, *lenZ);
	DBG_ASSERT(HERE, *lenZ <= lenA + lenB, "mi64_mul_vector: *lenZ > (lenA + lenB)!");
}

#if MI64_DEBUG
	#define DEBUG_SQUARE	0
#endif

/* Squaring-specialized version of above. By way of example, consider a length-10 input vector and
examine full-double-width square in term-by-term fashion, powers of the base b = 2^64 in left col:

b^n:	Coefficient (unnormalized)
----	--------------------------
n = 0	x0^2
1		x0.x1.2
2		x0.x2.2 + x1^2
3		x0.x3.2 + x1.x2.2
4		x0.x4.2 + x1.x3.2 + x2^2
5		x0.x5.2 + x1.x4.2 + x2.x3.2
6		x0.x6.2 + x1.x5.2 + x2.x4.2 + x3^2
7		x0.x7.2 + x1.x6.2 + x2.x5.2 + x3.x4.2
8		x0.x8.2 + x1.x7.2 + x2.x6.2 + x3.x5.2 + x4^2
9		x0.x9.2 + x1.x8.2 + x2.x7.2 + x3.x6.2 + x4.x5.2
10		          x1.x9.2 + x2.x8.2 + x3.x7.2 + x4.x6.2 + x5^2
11		                    x2.x9.2 + x3.x8.2 + x4.x7.2 + x5.x6.2
12		                              x3.x9.2 + x4.x8.2 + x5.x7.2 + x6^2
13		                                        x4.x9.2 + x5.x8.2 + x6.x7.2
14		                                                  x5.x9.2 + x6.x8.2 + x7^2
15		                                                            x6.x9.2 + x7.x8.2
16		                                                                      x7.x9.2 + x8^2
17		                                                                                x8.x9.2
18		                                                                                   x9^2
19		                                                                                      0

We can thus first compute the "off-diagonal" terms:

b^n:	Coefficient (unnormalized)
----	--------------------------
n = 0
1		x0.x1
2		x0.x2
3		x0.x3 + x1.x2
4		x0.x4 + x1.x3
5		x0.x5 + x1.x4 + x2.x3
6		x0.x6 + x1.x5 + x2.x4
7		x0.x7 + x1.x6 + x2.x5 + x3.x4
8		x0.x8 + x1.x7 + x2.x6 + x3.x5
9		x0.x9 + x1.x8 + x2.x7 + x3.x6 + x4.x5
10		        x1.x9 + x2.x8 + x3.x7 + x4.x6
11		                x2.x9 + x3.x8 + x4.x7 + x5.x6
12		                        x3.x9 + x4.x8 + x5.x7
13		                                x4.x9 + x5.x8 + x6.x7
14		                                        x5.x9 + x6.x8
15		                                                x6.x9 + x7.x8
16		                                                        x7.x9
17		                                                                x8.x9
18		                                                                    0
19		                                                                    0

Which we can do as set of ever-shorter scalar-vector muls, with intermediate vector add-with-carry:

	memset(z+len, 0ull, (len<<3));	// Clear upper half of output vector
	// Initial scalar-vector mul can be used to populate low [len+1] slots of output vector:
	z[len] = mi64_mul_scalar(x+1, x[0], z+1, len-1);
	// Use loop to take care of the rest:
	z[len+1] = mi64_mul_scalar(x+2, x[1], u, len-2);	ASSERT(0 == mi64_add(z+ 3, u, z+ 3, len-2));
	z[len+2] = mi64_mul_scalar(x+3, x[2], u, len-3);	ASSERT(0 == mi64_add(z+ 5, u, z+ 5, len-3));
	z[len+3] = mi64_mul_scalar(x+4, x[3], u, len-4);	ASSERT(0 == mi64_add(z+ 7, u, z+ 7, len-4));
	z[len+4] = mi64_mul_scalar(x+5, x[4], u, len-5);	ASSERT(0 == mi64_add(z+11, u, z+11, len-5));
	z[len+5] = mi64_mul_scalar(x+6, x[5], u, len-6);	ASSERT(0 == mi64_add(z+13, u, z+13, len-6));
	z[len+6] = mi64_mul_scalar(x+7, x[6], u, len-7);	ASSERT(0 == mi64_add(z+15, u, z+15, len-7));
	z[len+7] = mi64_mul_scalar(x+8, x[7], u, len-8);	ASSERT(0 == mi64_add(z+17, u, z+17, len-8));
	z[len+8] = mi64_mul_scalar(x+9, x[8], u, len-9);	ASSERT(0 == mi64_add(z+19, u, z+19, len-9));

We then double the result, and compute the square terms and add them in. This can be done efficiently
by making a copy of the z-vector resulting from the above (the one needing double-and-add-to-square-terms),
initing a vector containing the (nonoverlapping) square terms, and then doing a full-length vector-add.
*/
void	mi64_sqr_vector(const uint64 x[], uint64 z[], uint32 len)
{
	uint32 i, j, len8 = (len<<3);
	uint64 sgn, cy;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;

	ASSERT(HERE, z != x, "Input and output arrays must be distinct!");
	DBG_ASSERT(HERE, len != 0, "zero-length X-array!");

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		dimU = (len+1);
		free((void *)u);
		free((void *)v);
		u = (uint64 *)calloc((len+1), sizeof(uint64));
		v = (uint64 *)calloc((len*2), sizeof(uint64));
	}
	z[0] = 0ull, memset(z+len, 0ull, len8);	// Clear z[0], z[len+1:2*len]
	memcpy(u, x, len8);			// Save copy of input vector in u

	// Initial scalar-vector mul can be used to populate low [len+1] slots of output vector:
	z[len] = mi64_mul_scalar(x+1, x[0], z+1, len-1);
	// Use loop to take care of the rest:
	for(i = 2, j = 1; i < len; j = i++)	// j stores i-1 throughout loop
	{
	#ifdef YOU_WANT_IT_TO_BE_SLOW
		z[len+j]  = mi64_mul_scalar(x+i, x[j], u, len-i);
		z[len+j] += mi64_add    (z+i+j, u, z+i+j, len-i);
	#else
		z[len+j]  = mi64_mul_scalar_add_vec2(x+i, x[j], z+i+j, z+i+j, len-i);
	#endif
	}
	// Init vector containing the (nonoverlapping) square terms at same time we do doubling of the cross-terms vector.
	// We avoid the need for add-with-carry by doing the doubling via left-shift of the current word 1 bit
	// and copying the (saved) high bit of the next-lower word into the thus-vacated low bit.
	cy = 0;
	for(i = j = 0; i < len; ++i, j +=2)
	{
		sgn = (int64)z[j  ] < 0;
		z[j  ] = (z[j  ] << 1) + cy;	cy = sgn;
	#ifdef MUL_LOHI64_SUBROUTINE
		SQR_LOHI64(x[i],v+j ,v+j+1 );
	#else
		SQR_LOHI64(x[i],v[j],v[j+1]);
	#endif
		sgn = (int64)z[j+1] < 0;
		z[j+1] = (z[j+1] << 1) + cy;	cy = sgn;
	}
	// then do a full-length vector-add to obtain the result:
	mi64_add(z,v,z,2*len);

#if DEBUG_SQUARE
	// Compare result to that using generic vector-mul:
	mi64_mul_vector(x, len, x, len, v, &i);
	if(!mi64_cmp_eq(v,z,2*len)) {
		for(i = 0; i < 2*len; ++i)
		{
			printf("i = %2u: v = %20llu, z = %20llu, error = %20lld\n",i,v[i],z[i],(int64)(z[i]-v[i]));
		}
		ASSERT(HERE, 0, "mi64_sqr_vector result incorrect!");
	}
#endif
}
/* Initial Debug:
																			error z[j]-v[j]
																			---------------
(gdb) p v[0]	=               515524		z[0]	=               514544	-980
(gdb) p v[1]	=                    0		z[1]	=                    1	+1
(gdb) p v[2]	=                    0		z[2]	=                    0
(gdb) p v[3]	=   695282096289087488		z[3]	=   695282096289087487	-1
(gdb) p v[4]	=  6561739352485612000		z[4]	=  6561739352485612000
(gdb) p v[5]	=  4280339908740614450		z[5]	=  4280339908740614450
(gdb) p v[6]	=  4655141245927292619		z[6]	=  4655141245927292618	-1
(gdb) p v[7]	= 16391347641773540250		z[7]	= 16391347641773540249	-1
(gdb) p v[8]	=  9660163612446427001		z[8]	=  9660163612446427001
(gdb) p v[9]	=  6118265205031034628		z[9]	=  6118265205031034627	-1
(gdb) p v[10]	=   949770555345385036		z[10]	=   949770555345385035	-1
(gdb) p v[11]	=  6888472502951621128		z[11]	=  6888472502951621128
(gdb) p v[12]	= 10263162476022843315		z[12]	= 10263162476022843315
(gdb) p v[13]	=  4557757262155789608		z[13]	=  4557757262155789608
(gdb) p v[14]	=  3367728947229070195		z[14]	=  3367728947229070195
(gdb) p v[15]	= 17118160634115880599		z[15]	= 17118160634115880599
(gdb) p v[16]	=  1896777173959327369		z[16]	=  1896777173959327369
(gdb) p v[17]	= 12756822023939931388		z[17]	= 12756822023939931387	-1
(gdb) p v[18]	= 14488044506087195387		z[18]	= 14488044506087195387
(gdb) p v[19]	=               129184		z[19]	=               129184
(gdb)
*/
/* Low and high-half mul allow in-place: */
void	mi64_mul_vector_lo_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 j;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;

	DBG_ASSERT(HERE, len != 0, "`mi64_mul_vector_lo_half: zero-length X-array!");

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		dimU = (len+1);
		free((void *)u);
		free((void *)v);
		u = (uint64 *)calloc((len+1), sizeof(uint64));
		v = (uint64 *)calloc((len*2), sizeof(uint64));
	}
	memset(v, 0ull, (len<<4));	// Accumulator v[] needs to be cleared each time

	/* Loop over the elements of y[], multiplying x[] by each, and
	using u[] as a scratch array to store x[]*y[j] prior to adding to z[].

	For the high-half version, only want terms x[i]*y[j] with (i + j) >= len,
	plus the high halves of the (i + j) = len-1 terms for the requisite carryins.
	*/
	for(j = 0; j < len; j++)
	{
		if(y[j] == 0)
			continue;
	#ifdef YOU_WANT_IT_TO_BE_SLOW
		// Only need y[j]*u[0:len-1-j], i.e. high term of order [len-1], and can discard carryouts here:
		mi64_mul_scalar(x, y[j], u, len-j);
		// Add j-word-left-shifted u[] to v[], retaining only terms with index < len in the sum:
		mi64_add(v+j, u, v+j, len-j);
	#else
		mi64_mul_scalar_add_vec2(x, y[j], v+j, v+j, len-j);
	#endif
	}
	/* Copy v[0:len-1] into z[0:len-1]: */
	memcpy(z,v,(len<<3));
}

void	mi64_mul_vector_hi_half	(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 j;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;

	DBG_ASSERT(HERE, len != 0, "`mi64_mul_vector_hi_half: zero-length X-array!");

	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		dimU = (len+1);
		free((void *)u);	free((void *)v);
		u = (uint64 *)calloc((len+1), sizeof(uint64));
		v = (uint64 *)calloc((len*2), sizeof(uint64));
	}
	memset(v, 0ull, (len<<4));	// Accumulator v[] needs to be cleared each time

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
		/* Add j-word-left-shifted u[] to v[]: */
		v[len+j] = u[len] + mi64_add(&v[j], u, &v[j], len);
	}
	/* Copy v[len:2*len-1] into z[0:len-1]: */
	memcpy(z,v+len,(len<<3));
}

// Fast version of above, which only computes hi-half product terms and the additional ones
// needed to approximate the exact carryin to the hi half by its upper 64 bits.
// This fails for certain well-known input types, e.g. when x and y have a significant fraction
// of their 64-bit 'digts' very close in value to 2^64-1:

#if MI64_DEBUG
	#define DEBUG_HI_MUL	0	// Set nonzero to enable debug of mi64_mul_vector_hi_fast
#endif

/* 4/22/2012: This rhombus-truncated fast version of the above routine needs further work to
implement an error-correction for the low-half carryout. Here are notes related to an example
illustrating this.

Inputs are in the conetxt of emulated 640-bit twos-complement integer arithmetic:

	x = 2^607-1 ==> x[0,9] = [2^64-1,...,2^64-1,2^31-1]
	y[0,9] = [0,0,0,0,0,0,2^64-2^60,2^64-1,2^64-1,2^64-1] ==> y =

y is output of mi64_mul_vector_lo_half(lo, qinv, lo) with lo = [0,0,0,2^30,0...,0] = 2^222, q = x = 2^607-1,
qinv[0,9] = [2^64-1,...,2^64-1,2^64-1-2^31] ==> qinv = 2^640-2^607-1, check Montgomery-inverse property:

q*qinv = [2^607-1]*[2^640-2^607-1]
		= 2^1247 - 2^640 - 2^1214 + 2^607 - 2^607 + 1
		= 2^1247 - 2^1214 - 2^640 + 1 == 1 (mod 2^640) .

Now examine full-double-width multiply in term-by-term fashion, powers of the base b = 2^64 in left col:

b^n:	Coefficient (unnormalized)
----	--------------------------
n = 0	x0.y0
1		x0.y1 + x1.y0
2		x0.y2 + x1.y1 + x2.y0
3		x0.y3 + x1.y2 + x2.y1 + x3.y0
4		x0.y4 + x1.y3 + x2.y2 + x3.y1 + x4.y0
5		x0.y5 + x1.y4 + x2.y3 + x3.y2 + x4.y1 + x5.y0
6		x0.y6 + x1.y5 + x2.y4 + x3.y3 + x4.y2 + x5.y1 + x6.y0
7		x0.y7 + x1.y6 + x2.y5 + x3.y4 + x4.y3 + x5.y2 + x6.y1 + x7.y0
8		x0.y8 + x1.y7 + x2.y6 + x3.y5 + x4.y4 + x5.y3 + x6.y2 + x7.y1 + x8.y0
9		x0.y9 + x1.y8 + x2.y7 + x3.y6 + x4.y5 + x5.y4 + x6.y3 + x7.y2 + x8.y1 + x9.y0
10		        x1.y9 + x2.y8 + x3.y7 + x4.y6 + x5.y5 + x6.y4 + x7.y3 + x8.y2 + x9.y1
11		                x2.y9 + x3.y8 + x4.y7 + x5.y6 + x6.y5 + x7.y4 + x8.y3 + x9.y2
12		                        x3.y9 + x4.y8 + x5.y7 + x6.y6 + x7.y5 + x8.y4 + x9.y3
13		                                x4.y9 + x5.y8 + x6.y7 + x7.y6 + x8.y5 + x9.y4
14		                                        x5.y9 + x6.y8 + x7.y7 + x8.y6 + x9.y5
15		                                                x6.y9 + x7.y8 + x8.y7 + x9.y6
16		                                                        x7.y9 + x8.y8 + x9.y7
17		                                                                x8.y9 + x9.y8
18		                                                                        x9.y9
19		                                                                          0

Does use of only the b^9 = x0.y9 + x1.y8 + x2.y7 + x3.y6 + x4.y5 + x5.y4 + x6.y3 + x7.y2 + x8.y1 + x9.y0 term
yield the correct carry into the b^10 term here?

bc
b=2^64;
x0=x1=x2=x3=x4=x5=x6=x7=x8=b-1; x9=2^31-1
y0=y1=y2=y3=y4=y5=0; y6=b-2^60; y7=y8=y9=b-1;

																								bn (normalized)			cy
b0-5 = 0																						0						0
b6 = cy+x0*y6+x1*y5+x2*y4+x3*y3+x4*y2+x5*y1+x6*y0                   ;m=b6%b ;m;cy=(b6 -m)/b;cy	1152921504606846976		17293822569102704639
b7 = cy+x0*y7+x1*y6+x2*y5+x3*y4+x4*y3+x5*y2+x6*y1+x7*y0             ;m=b7%b ;m;cy=(b7 -m)/b;cy	0						35740566642812256254
b8 = cy+x0*y8+x1*y7+x2*y6+x3*y5+x4*y4+x5*y3+x6*y2+x7*y1+x8*y0       ;m=b8%b ;m;cy=(b8 -m)/b;cy	0						54187310716521807869
b9 = cy+x0*y9+x1*y8+x2*y7+x3*y6+x4*y5+x5*y4+x6*y3+x7*y2+x8*y1+x9*y0 ;m=b9%b ;m;cy=(b9 -m)/b;cy	0						72634054790231359484
b10= cy+      x1*y9+x2*y8+x3*y7+x4*y6+x5*y5+x6*y4+x7*y3+x8*y2+x9*y1 ;m=b10%b;m;cy=(b10-m)/b;cy	18446744073709551615	72634054790231359484
b11= cy+            x2*y9+x3*y8+x4*y7+x5*y6+x6*y5+x7*y4+x8*y3+x9*y2 ;m=b11%b;m;cy=(b11-m)/b;cy
b12= cy+                  x3*y9+x4*y8+x5*y7+x6*y6+x7*y5+x8*y4+x9*y3 ;m=b12%b;m;cy=(b12-m)/b;cy
b13= cy+                        x4*y9+x5*y8+x6*y7+x7*y6+x8*y5+x9*y4 ;m=b13%b;m;cy=(b13-m)/b;cy
b14= cy+                              x5*y9+x6*y8+x7*y7+x8*y6+x9*y5 ;m=b14%b;m;cy=(b14-m)/b;cy
b15= cy+                                    x6*y9+x7*y8+x8*y7+x9*y6 ;m=b15%b;m;cy=(b15-m)/b;cy
b16= cy+                                          x7*y9+x8*y8+x9*y7 ;m=b16%b;m;cy=(b16-m)/b;cy
b17= cy+                                                x8*y9+x9*y8 ;m=b17%b;m;cy=(b17-m)/b;cy
b18= cy+                                                      x9*y9 ;m=b18%b;m;cy=(b18-m)/b;cy
b19= cy+

Using low-order-truncated-rhombus, get

b9  = x0*y9 + x1*y8 + x2*y7 + x3*y6 + x4*y5 + x5*y4 + x6*y3 + x7*y2 + x8*y1 + x9*y0 = 1152921504606846979 + b*72634054790231359481	<*** cy is 3 too low! ***
thus carry into b10 = 72634054790231359481 .

With the 3-too-low approx-carry, get

b10 = cy + x1*y9 + x2*y8 + x3*y7 + x4*y6 + x5*y5 + x6*y4 + x7*y3 + x8*y2 + x9*y1
	= 18446744073709551612 + b*72634054790231359484 ,
which has remainder also 3-too-low (18446744073709551612 vs the correct 18446744073709551615 = 2^64-1) but correct carryout.

UPSHOT: Could spend days coming up with error-coorection scheme, but no point, since for mm(p)-TF-related modmuls we already
have an exact fast O(n) alternative to compute the mi64_mul_vector_hi_half(q, lo) result.
*/
#if 0
#error mi64_mul_vector_hi_trunc not yet debugged!
void	mi64_mul_vector_hi_trunc(const uint64 x[], const uint64 y[], uint64 z[], uint32 len)
{
	uint32 j, lm1 = len-1;
	static uint64 *u, *v;	// Scratch arrays for storing intermediate scalar*vector products
	static uint32 dimU = 0;
	DBG_ASSERT(HERE, len != 0, "zero-length X-array!");

	/* Do scratch arrays need allocating or reallocating? */
	if(dimU < (len+1))
	{
		dimU = (len+1);
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
		// Only need y[j]*u[len-1-j:len-1], i.e. low term of order [len-1]
		u[lm1] = mi64_mul_scalar(&x[lm1-j], y[j], &u[lm1-j], j);
		// Add j-word-left-shifted u[] to v[], retaining only terms with index < len in the sum:
		v[lm1+j] = u[lm1] + mi64_add(&v[lm1], &u[lm1-j], &v[lm1], j);
		// Full-length: v[j, j+1, ..., j+len-1] gets added to u[0,1,...,len-1]
		// Half-length: v[ len-1, ..., j+len-1] gets added to u[len-1-j,...,len-1]
	}
#if DEBUG_HI_MUL
	// Test code for fast version of this function - re-use low half of v[] for output::
	mi64_mul_vector_hi_half(x,y,v,len);
	if(!mi64_cmp_eq(v,v+len,len)) {
		printf("mi64_mul_vector_hi_trunc result incorrect!");
	}
#endif

	/* Copy v[len:2*len-1] into z[0:len-1]: */
	memcpy(z,v+len,(len<<3));
}
#endif


/* Fast O(n) variant of the more-general O(n^2) mi64_mul_vector_hi_half() to compute UMULH(q, Y)
for the special case of MM(p)-TF-related modmul, i.e. with modulus q = 2.k.M(p) + 1, where M(p)
is a Mersenne prime. For such special-form factor candidates q we define

	Z = 2.k.Y, needing just O(N) work, where N = number of computer words needed to store Y,

And then use that

	UMULH(q,Y) = ((Z << p) - (2k-1).Y) >> B . [*]

where B is the emulated bits per vector-integer 'word' as defined by the associated MULL implementation;
that is, 2^B is the emulated base of the 2s-comp arithmetic being done.

Note that the result can NOT be computed like so:

			   = (Z >> (B-p)) - (2k-1).(Y >> B) ,

because this neglects the possible borrow from the bits of the difference which get shifted off.

The (Z << p) term in [*] needs O(n) hardware shifts and adds to compute and (2k-1).Y needs no
explicit shifts but rather O(n) hardware multiplies to evaluate the scalar-vector product,
followed by simply retaining the single-word carry out of that, discarding the low B bits.

We replace the literal vector-q-input in the function arglist with the single-word integers k and p.
*/
// emulated bits per 'word' as implemented by MULL = (len<<6):
void	mi64_mul_vector_hi_qmmp(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 bits)
{
	// If bits not a multiple of 64, len2 may = 2*len-1 rather than 2*len:
	uint32 i = (bits+63), len = (i >> 6), len2 = ((i+bits) >> 6), len8 = (len << 3);
	uint64 k2 = k+k, bw;
	/* Scratch array for storing intermediate scalar*vector products: */
	static uint64 *u, *v;
	static uint32 dimU = 0;
//====need to finish 200-bit support! =======================
	ASSERT(HERE, z != y, "Input and output arrays must be distinct!");
	ASSERT(HERE, p < bits, "shift parameters out of range!");
	DBG_ASSERT(HERE, len != 0, "zero-length X-array!");
	/* Does scratch array need allocating or reallocating? */
	if(dimU < (len+1))
	{
		dimU = (len+1);
		// U needs extra pad slot (beyond the one needed for carryout of Y*k2) to ensure mi64_shl can never grab an uninited high word
		free((void *)u);	u = (uint64 *)calloc(len+2, sizeof(uint64));
		free((void *)v);	v = (uint64 *)calloc(len*2, sizeof(uint64));
	} else {
		for(i = len+1; i < len2; i++) {
			u[i] = 0ull;	// With proper padding of U don't need any zeroing of V prior to V = (U << p) step below
		}
	}
	DBG_ASSERT(HERE, (k != 0) && ((k2>>1) == k), "2*k overflows!");	// Make sure 2*k did not overflow
	u[len] = mi64_mul_scalar(y,k2,u,len);	// u[] stores Z = 2.k.Y
	mi64_shl(u,v,p,len2);			// v[] stores (Z << p), store result in V
	u[len] -= mi64_sub(u,y,u,len);	// (2k-1).Y = Z-Y, store result in U
	bw = mi64_sub(v,u,v,len+1);
	DBG_ASSERT(HERE, !bw, "Unexpected borrow!");

	/* Right-shift by B bits to get UMULH(q,Y) = ((Z << p) - (2k-1).Y) >> B: */
	mi64_shrl(v,v,bits,len2);
	memcpy(z,v,len8);
//	memcpy(z,v+len,len8);

#if MI64_DEBUG
	#define DEBUG_HI_MUL	0	// Set nonzero to compare versus mi64_mul_vector_hi_fast
#endif
#if DEBUG_HI_MUL
//rror Check that you want DEBUG_HI_MUL set in mi64_mul_vector_hi_qmmp!
	// Compute q, store in u[]:
	memset(u, 0ull, (dimU<<3));
	u[0] = 1;
	mi64_shl(u, u, p, len);			// 2^p
	mi64_sub_scalar(u, 1, u, len);	// M(p) = 2^p-1
	ASSERT(HERE, 0 == mi64_mul_scalar(u, k2, u, len), "2.k.M(p) overflows!");	// 2.k.M(p)
	mi64_add_scalar(u, 1ull, u, len);	// q = 2.k.M(p) + 1
	// Test code for fast version of this function - re-use v[] for output::
//	mi64_mul_vector_hi_half(u,y,v,len);
	mi64_mul_vector_hi_fast(y,p,k,v,len);
	if(!mi64_cmp_eq(v,z,len)) {
		ASSERT(HERE, 0, "mi64_mul_vector_hi_qmmp/fast results differ!");
	}
#endif
}

/* Version #2 of the above, which further streamlines things. again start with the following 2 multiplicands:

	1. q = 2.k.M(p) + 1, where M(p) is a Mersenne prime;
	2. multiword multiplicand Y < q,

and take B as the emulated bits per vector-integer 'word' as defined by the associated MULL implementation; that
is, 2^B is the emulated base of the 2s-comp arithmetic being done; in terms of the hardware-integer wordsize W = 2^b,
B = n*b, i.e. the multiplicands are n-word vectors in terms of the hardware-integer wordsize.
In terms of relative sizes, y < q < B and k < W, i.e. k is a 1-word scalar multiplier. We further assume k < W/2,
thus 2*k fits into a single word, as well.

Define

	Z = 2.k.Y, which needs just O(n) work to compute via mi64_mul_scalar().

Again we use that

	UMULH(q,Y) = ((Z << p) - (2k-1).Y) >> B . [*]

The above version of this function computes the full-length scalar-vector product Z = 2.k.Y, then does a full-length
vector-vector subtract to obtain Z - Y = (2k-1).Y, followed by a double-wide vector-left-shift to get (Z << p), another
full-length vector-vector subtract to get (Z << p) - (2k-1).Y, followed by a vector-subtract-scalar to propagate the
resulting borrow into the high (n) words of the length-2n (Z << p), followed by a length-n vector-copy (in lieu of
a double-wide vector-right-shift of the result by B bits) in order to obtain the result [*].

This is O(n) work but still involves much wasted computation, specifically related to the fact that half the
words of the penultimate double-wide vector are discarded: These are only needed for the possible
subtraction-borrow they provide to the retained bits.

	Ex.: q = 2.k.M(127) + 1 with k = 7143819210136784550, i.e.
	q = 2430915709680614116949754105299803650411408301848040235701 ;
	y =  915005412744957807408012591600653057424688130286064771258 = y0 + 2^64*y1 + 2^128*y2,
	with y0 = 2294959606785646778; y1 = 10167084567166165345; y2 = 2688959234133783535 .
[Note that
2^192 = 6277101735386680763835789423207666416102355444464034512896 .
]
	Exact result:

		UMULH_192(q,y) = 354351598245602020483095922210514413558224553895064094733 = u0 + 2^64*u1 + 2^128*u2,
		with u0 = 141525868296128525, u1 = 4269430960237156763, u2 = 1041345754856384950 .

	b = 192 bits
	p = 127 bits
	(b-p) = 65
	Compute (2k-1).y, store in z:
	z' = (2k-1).y = 2^b*2082691509712769900 + [low 192 bits], compare that retained coefficient to direct high-64 bits of 128-bit product (2k-1).y2:
	umulh64(2k-1, y2) = 2082691509712769899, one too low, because this neglects the possibility of a carryin resulting from
	 mull64(2k-1, y2) + umulh64(2k-1, y1).

But we need the full-length product (2k-1).y *anyway* to get Z = 2.k.y (by adding another copy of y), so do as follows:

1. compute z' = (2k-1).y via vector-scalar mul, the carryout word cw = ((2k-1).Y >> B);

2. compute low n words of z = z' + y via vector-vector add, make sure there is no carryout of that (or if there is, add it to a 2nd copy of cw, say, cz);

3. compute low n words of z >> (b-p), then separately shift in cw from the left, via (2^b*cz) >> (b-p) = (cz << p).
[*** idea: define special version of mi64_shift for this ***]
4. subtract scalar cw from resulting vector to effect ... - (2k-1).(Y >> B) step in [*].
*/
void	mi64_mul_vector_hi_fast(const uint64 y[], const uint64 p, const uint64 k, uint64 z[], uint32 len)
{
	uint32 i, bits;
	uint64 k2m1 = k-1+k, tmp,bw0,bw1,bw,cw,cy,cz;
	uint64 *zptr;
	ASSERT(HERE, z != y, "Input and output arrays must be distinct!");
	ASSERT(HERE, (k != 0) && ((k2m1>>1) == k-1), "2*k-1 overflows!");
	DBG_ASSERT(HERE, len != 0, "zero-length X-array!");

// 1. compute z' = (2k-1).y via vector-scalar mul, the carryout word cw = ((2k-1).Y >> B);
	cw = mi64_mul_scalar(y,k2m1,z,len);	// z' = (2k-1).y
	bw0 = z[len-1];
//if(k==900) printf("Mi64: bw0 = %20llu, cw = %20llu, z` = %s\n", bw0,cw,&s0[convert_mi64_base10_char(s0, z, len)]);
// 2. compute low n words of z = z' + y via vector-vector add, any carryout of that gets added to a 2nd copy of cw, cz;
	cz = cw + mi64_add(y,z,z, len);	// z = z' + y
//if(k==900) printf("Mi64: cz = %20llu, z = %s\n", cz,&s0[convert_mi64_base10_char(s0, z, len)]);

// 3. compute low n words of z >> (b-p), then separately shift in cz from the left, via (2^b*cz) >> (b-p) = (cz << p).
	ASSERT(HERE, (len<<6) > p, "shift parameters out of range!");
	bw1 = mi64_shrl(z,z,(len<<6)-p,len);	// low n words of z >> (b-p); high 64 bits of off-shifted portion saved in bw1
//if(k==900) printf("Mi64: bw1 = %20llu, z>> = %s\n", bw1,&s0[convert_mi64_base10_char(s0, z, len)]);

/* Check for borrow-on-subtract of to-be-off-shifted sections: have a borrow if
	z' (result from above mul_scalar, not including the carryout word cw) >	((z << p) % 2^b) (off-shifted portion of z = z' + y above, left-justified to fill a b-bit field)
In order to compute this borrow with high probability without having to store the entire b-bit vectors in this comparison,
we save just the high word of z', bw0 := z'[len-1], and the high word of the off-shifted portion of z (call it bw1), and compare.
The only tricky case is is these word-length proxies for the full-length vectors compare equal, in which case we should abort
and tell the user to call the slow exact version of this function, currently inactivated inside the above #if 0.
*/
	bw = (bw0 > bw1);
	// Add (cw << (b-p)) to result.
	bits = (p&63);
	i = (p >> 6);	// Low word into which cw will get added
	zptr = z+i;
	// If (b-p) == 0 (mod 64) all of cz goes into z[i], with i = (b-p)/64;
	if(bits == 0) {
		ASSERT(HERE, 0 == mi64_add_scalar(zptr,cz,zptr,len-i), "unexpected carryout of ( + cw)!");
	// Otherwise cz gets split between z[i] and z[i+1]:
	} else {
		// low 64-(p%64) bits of cz = (cz << bits) go into z[i]:
		tmp = (cz << bits);
		*zptr += tmp; cy = (*zptr++ < tmp);
		// high (p%64) bits of cw = (cw >> bits) go into z[i+1]
		ASSERT(HERE, 0 == mi64_add_scalar(zptr,(cz >> (64-bits)) + cy,zptr,len-i-1), "unexpected carryout of ( + cw).hi!");
	}

// 4. subtract scalar (bw + cw) from resulting vector to effect ... - (2k-1).Y step in [*].
	ASSERT(HERE, 0 == mi64_sub_scalar(z,(bw + cw),z,len), "unexpected carryout of (... - cw) !");
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
	DBG_ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "mi64_cvt_uint64_double: FFT_MUL_BASE not pure-integer!");
	DBG_ASSERT(HERE, FFT_MUL_BASE < TWO54FLOAT, "mi64_cvt_uint64_double: FFT_MUL_BASE >= maximum allowed value of 2^54!");

	/* As we extract each floating-point word, balance it and set
	resulting carry into next FP word: */
	cyi = cyj = 0;
	for(i = 0; i < len; i++)
	{
		curr_re64 = x[i];
		curr_im64 = y[i];

		for(k = 0; k < 4; k++)
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

	DBG_ASSERT(HERE, cyi <= 1,"mi64_cvt_uint64_double: Output Real carry out of range!");
	DBG_ASSERT(HERE, cyj <= 1,"mi64_cvt_uint64_double: Output Imag carry out of range!");

	/* It is desirable to not have the FP vector length exceed 4*len,
	so suppress any output carry by folding back into MS array element:
	*/
	if(cyi)
	{
		DBG_ASSERT(HERE, a[j  ] < 0,"mi64_cvt_uint64_double: MS array element >= 0!");
		a[j  ] += FFT_MUL_BASE;
	}

	if(cyj)
	{
		DBG_ASSERT(HERE, a[j+1] < 0,"mi64_cvt_uint64_double: MS array element >= 0!");
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
	int64 cy_re, cy_im, itmp, jtmp;
	uint64 curr_re64, curr_im64;

	ASSERT(HERE, n != 0, "zero-length array!");

	/* Redo the quicker checks of those done in util.c::check_nbits_in_types() */
	ASSERT(HERE, DNINT(FFT_MUL_BASE) == FFT_MUL_BASE, "FFT_MUL_BASE not pure-integer!");
	ASSERT(HERE, FFT_MUL_BASE < TWO54FLOAT, "FFT_MUL_BASE >= maximum allowed value of 2^54!");
	/*
	!...Make sure MSW of Re(A[]) and Im(A[]) in the
	balanced-representation form are both >= 0:
	*/
	/* Re(A[]) stored in even terms: */
	for(i = 2*n-2; i >= 0; i-=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );
		if(a[j] != 0.0)
		{
			ASSERT(HERE, a[j] > 0.0, "MSW(Re(A[])) < 0!");
			break;
		}
	}
	/* Im(A[]) stored in odd terms: */
	for(i = 2*n-1; i >= 1; i-=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );
		if(a[j] != 0.0)
		{
			ASSERT(HERE, a[j] > 0.0, "MSW(Im(A[])) < 0!");
			break;
		}
	}

	/*...Now form the terms of the uint64 representation:	*/

	len   = 0;		/* Full 64-bit words accumulated so far in the residue	*/
	nbits = 0;		/* Total bits accumulated so far in the x-array	*/
	curr_bits = 0;	/* Total bits currently stored in the 64-bit accumulator word */
	curr_re64 = 0;	/* Current value of the 64-bit accumulator word */
	curr_im64 = 0;

	cy_re = cy_im = 0;	/* Init carry */
	for(i = 0; i < 2*n; i+=2)
	{
		j = i + ( (i >> DAT_BITS) << PAD_BITS );

		ASSERT(HERE, curr_bits < 64,"curr_bits < 64");
		itmp = (uint64)1<<curr_bits;
		ASSERT(HERE, curr_re64 < itmp,"curr_wd64 !< (1<<curr_bits)");
		ASSERT(HERE, curr_im64 < itmp,"curr_wd64 !< (1<<curr_bits)");
		ASSERT(HERE, DNINT(a[j]) == a[j], "a[j] not pure-integer!");
		ASSERT(HERE, ABS(a[j]) < TWO54FLOAT, "|a[j]| >= maximum allowed value of 2^54!");

		itmp = (int64)a[j  ] + cy_re;	/* current digit in int64 form, subtracting any borrow from previous digit.	*/
		if(itmp < 0)	/* If current digit < 0, add the base and set carry = -1	*/
		{
			itmp += FFT_MUL_BASE;
			cy_re = -1;
		}
		else
		{
			cy_re = 0;
		}
		ASSERT(HERE, itmp >= 0,"itmp < 0!");
		ASSERT(HERE, (curr_re64>>curr_bits) == 0,"(curr_re64>>curr_bits) != 0!");

		jtmp = (int64)a[j+1] + cy_im;
		if(jtmp < 0)
		{
			jtmp += FFT_MUL_BASE;
			cy_im = -1;
		}
		else
		{
			cy_im = 0;
		}
		ASSERT(HERE, jtmp >= 0,"jtmp < 0!");
		ASSERT(HERE, (curr_im64>>curr_bits) == 0,"(curr_re64>>curr_bits) != 0!");

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

	ASSERT(HERE, nbits == n*FFT_MUL_BITS,"nbits == n*FFT_MUL_BASE!");
	ASSERT(HERE, len > 0                ,"len has overflowed!");

	return len;
}

/***************/

/* returns 1 if the multiword unsigned integer p is a base-z Fermat pseudoprime, 0 otherwise. The base is assumed to fit in a 64-bit unsigned scalar. */
#if MI64_DEBUG
	#define MI64_PRP_DBG 0
#endif
uint32 mi64_pprimeF(const uint64 p[], uint64 z, uint32 len)
{
/* The allocatable-array version of this gave bizarre stack faults under MSVC, so use fixed size,
especially as this routine is too slow to be useful for inputs larger than a couple hundred words anyway:
*/
/*	static uint64 *y, *n, *zvec, *ppad, *prod, *tmp;	*/
	uint64 y[1024], n[1024], zvec[1024], ppad[2048], prod[2048], tmp[2048];
	uint64 flag;
	uint32 len2 = 2*len, curr_len, retval, plen;
#if MI64_PRP_DBG
	uint128 q128, r128;
	uint192 q192, r192;
#endif
//if(p[0] == 4961431981940124743ull) {
//	printf("here!\n");
//}
	ASSERT(HERE, len <= 1024, "mi64_pprimeF: Max 1024 words allowed at present!");

/*
	y    = (uint64 *)calloc((  len), sizeof(uint64));
	n    = (uint64 *)calloc((  len), sizeof(uint64));
	zvec = (uint64 *)calloc((  len), sizeof(uint64));
	ppad = (uint64 *)calloc((len2), sizeof(uint64));
	prod = (uint64 *)calloc((len2), sizeof(uint64));
	tmp  = (uint64 *)calloc((len2), sizeof(uint64));
*/
	mi64_clear(y, len);
	mi64_clear(n, len);
	mi64_clear(zvec, len);
	mi64_clear(ppad, len2);
	mi64_clear(prod, len2);
	mi64_clear(tmp, len2);

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
			fprintf(stderr, "mi64_pprimeF: flag = %1d, y = %s, z = %s\n", (uint32)flag, &s0[convert_uint128_base10_char(s0, q128)], &s1[convert_uint128_base10_char(s1, r128)]);
		}
		else if(len == 3)
		{
			q192[0] = y[0]; q192.d1 = y[1]; q192.d2 = y[2];
			r192[0] = zvec[0]; r192.d1 = zvec[1]; r192.d2 = zvec[2];
			fprintf(stderr, "mi64_pprimeF: flag = %1d, y = %s, z = %s\n", (uint32)flag, &s0[convert_uint192_base10_char(s0, q192)], &s1[convert_uint192_base10_char(s1, r192)]);
		}
	#endif

		if(flag)
		{
			mi64_mul_vector(y, len, zvec, len, prod, &plen);	/* y*z */
			mi64_div(prod, ppad, len2, len2, tmp, prod);		/* Discard quotient here, overwrite prod with remainder */
			DBG_ASSERT(HERE, mi64_getlen(prod, len2) <= len, "mi64_pprimeF: (y*z)%p illegal length");
			mi64_set_eq(y, prod, len);	/* y = (y*z)%p */
		}
		mi64_mul_vector(zvec, len, zvec, len, prod, &plen);		/* z^2 */
		mi64_div(prod, ppad, len2, len2, tmp, prod);
		DBG_ASSERT(HERE, mi64_getlen(prod, len2) <= len, "mi64_pprimeF: (z^2)%p illegal length");
		mi64_set_eq(zvec, prod, len);	/* z = (z^2)%p */
		if(mi64_iszero(zvec, len))
		{
			retval=0;
		}
	}

	retval = mi64_cmp_eq_scalar(y, 1ull, len);
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

/* Fast division using Montgomery-modmul. First implemented: May 2012 */
void mi64_div(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[])
{
	uint32 xlen, ylen, max_len;
	// Only the quotient array is optional:
	ASSERT(HERE, lenX && lenY, "illegal 0 dimension!");
	ASSERT(HERE, x && y && r, "At least one of X, Y, R is null!");
	ASSERT(HERE, x != y, "X and Y arrays overlap!");
	ASSERT(HERE, r != y, "Y and Rem arrays overlap!");
	ASSERT(HERE, q != x && q != y && q != r, "Quotient array overlaps one of X, Y ,Rem!");

	/* Init Q = 0; don't do similarly for R since we allow X and R to point to same array: */
	if(q && (q != x)) {
		mi64_clear(q, lenX);
	}
	/* And now find the actual lengths of the divide operands and use those for the computation: */
	xlen = mi64_getlen(x, lenX);
	ylen = mi64_getlen(y, lenY);
	DBG_ASSERT(HERE, ylen != 0, "divide by 0!");

	/* If x < y then q = 0, r = x */
	max_len = MAX(xlen, ylen);
	if((xlen < ylen) || mi64_cmpult(x, y, max_len))
	{
		/* If no modular reduction needed and X and R not the same array, set R = X: */
		if(r != x) {
			mi64_set_eq(r, x, lenX);
		}
		return;
	}
	// If single-word divisor, use specialized single-word-divisor version:
	if(ylen == 1)
	{
		r[0] = mi64_div_by_scalar64(x, y[0], xlen, q);
		if(r == x) mi64_clear(r+1, lenY-1);	// In single-word special case, only clear high words of remainder array if rem-in-place (x == r)
	} else {
		mi64_div_mont(x, y, xlen, ylen, q, r);
		mi64_clear(r+ylen, lenX-ylen);	// This must wait until after divide completes in case R, X point to same memory
	}
	return;
}

/* Fast div-with-remainder with arbitrary-length divisor using Montgomery modmul, my right-to-left mod
algorithm, and my true-remainder postprocessing step of May 2012.

Returns x % y in r (required), and x / y in q (optional).
*/
#if MI64_DEBUG
	#define MI64_DIV_MONT	0
#endif
void mi64_div_mont(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[])
{
#if MI64_DIV_MONT
	uint32 dbg = 0;//STREQ(&s0[convert_mi64_base10_char(s0, x, lenX)], "184027369222821018387581692039");
	uint64 *qref = 0x0,*rref = 0x0;
#endif
	int i,j;	// i and j must be signed
	uint32 lenW,lenD, ybits, log2_numbits, n,p,nshift,nws=0,nbs=0;
	uint64 bw,mask,lo64,itmp64;
	double fquo;
	// pointers to local-storage:
	static int lens = 0;	// # of 64-bit ints allocated for current scratch space
	static uint64 *yinv = 0x0,*cy = 0x0, *tmp = 0x0, *itmp = 0x0, *lo = 0x0, *hi = 0x0, *w = 0x0, *rem_save = 0x0;
	static uint64 *scratch = 0x0;	// "base pointer" for local storage shared by all of the above subarrays
	static uint64 *v = 0x0;		// This one gets alloc'ed separately from the above ones, declare separately as reminder

	ASSERT(HERE, lenX && lenY, "illegal 0 dimension!");
	ASSERT(HERE, (lenY > 1) || (y[0] > 0), "Divide by zero!");
	ASSERT(HERE, (x && y) && (x != y), "Bad x or y array!");
	ASSERT(HERE, (r != 0) && (q != r), "Bad remainder array!");	// q may be 0x0, but must not overlap r

#if MI64_DIV_MONT
	if(dbg)	// Compute result using slow binary-div algo, use that as reference
	{
		qref = (uint64 *)calloc((lenX), sizeof(uint64));	ASSERT(HERE, qref != 0x0, "alloc fail!");
		rref = (uint64 *)calloc((lenX), sizeof(uint64));	ASSERT(HERE, rref != 0x0, "alloc fail!");
		lo   = (uint64 *)calloc((lenX), sizeof(uint64));	ASSERT(HERE,   lo != 0x0, "alloc fail!");
		hi   = (uint64 *)calloc((lenX), sizeof(uint64));	ASSERT(HERE,   hi != 0x0, "alloc fail!");
		mi64_set_eq(lo,x,lenX);
		mi64_set_eq(hi,y,lenY);
		mi64_div_binary(lo,hi,lenX,lenY,qref,rref);
		free((void *)lo);
		free((void *)hi);
	}
#endif

	/* Modulus y must be odd for Montgomery-style modmul to work, so first shift off any low 0s.
	Unlike the binary-result is-divisible-by routines, we must restore the shift to the
	quotient and remainder at the end. */
	lenD = mi64_getlen(y, lenY);
#if MI64_DIV_MONT
	if(dbg) {
		if(lenX <= 50)printf("x = %s\n", &s0[convert_mi64_base10_char(s0, x, lenX)]);
	//	printf("x = %s\n", &str_10k[__convert_mi64_base10_char(str_10k, 10<<10, x, lenX)]);
		printf("y = %s\n", &s0[convert_mi64_base10_char(s0, y, lenD)]);	// Leave length-check off this so if y too large for print we assert right here
	}
#endif
	if(lenD > lens) {
		free((void *)scratch);	scratch = yinv = cy = tmp = itmp = lo = hi = w = rem_save = 0x0;
		/* (re)Allocate the needed auxiliary storage: */
		scratch = (uint64 *)calloc((6*lenD), sizeof(uint64));	ASSERT(HERE, scratch != 0x0, "alloc fail!");
		yinv    = scratch;	// These ptrs just point to various disjoint length-lenD sections of the shared local-storage chunk
		cy      = yinv    + lenD;
		tmp     = cy      + lenD;
		itmp    = tmp     + lenD;
		lo      = itmp    + lenD;
		hi      = lo      + lenD;
		w       = hi      + lenD;	// This gets reset on each call and may end up == y, this is reminder what its copy-mode value should be
		rem_save= w       + lenD;
	}

	// Special case: x and y have same number of machine words, x > y (i.e. x does indeed need modding).
	// In such cases call binary-div to do modular reduction:
	if(lenD == lenX && !mi64_cmpult(x, y, lenX)) {
		// Here first attempt fast-mod based on FP approximation of the
		// quotient x/y, which only defaults to bin-div if x/y > 2^53
		i = mi64_extract_lead64(x, lenX, &lo64);
		j = mi64_extract_lead64(y, lenX, &itmp64);
		bw = 0x1ull << (i-j);
		fquo  = (double)bw;	// Power-of-2 scaling needed on ratio of left-justified lead64 fields
		fquo *= (double)lo64/(double)itmp64;
		if(fquo < TWO54FLOAT) {
			itmp64 = (uint64)fquo;
			// Since x,v,r may all point to same memory, need local-storage to hold y*fquo - use yinv to point to that:
			yinv = (uint64 *)calloc((lenD), sizeof(uint64));	ASSERT(HERE, yinv != 0x0, "alloc fail!");
			mi64_mul_scalar(y,itmp64,yinv,lenX);
			bw = mi64_sub(x,yinv,r,lenX);
		#if MI64_DIV_MONT
			if(dbg)printf("fquo*x = %s\n", &s0[convert_mi64_base10_char(s0, yinv, lenD)]);
			if(dbg)printf("     r = %s\n", &s0[convert_mi64_base10_char(s0, r   , lenD)]);
		#endif
			ASSERT(HERE, !bw, "Unexpected borrow!");
			if(!mi64_cmpult(r, y, lenX)) {
			#if MI64_DIV_MONT
				printf("Remainder exceeds modulus:");
				printf("x = %s\n", &s0[convert_mi64_base10_char(s0, x, lenX)]);
				printf("y = %s\n", &s0[convert_mi64_base10_char(s0, y, lenX)]);
				printf("r = %s\n", &s0[convert_mi64_base10_char(s0, r, lenX)]);
			#endif
				// Hack: subtract (y-r) and see if that satisfies:
				mi64_sub(r,y,r,lenX);
				ASSERT(HERE, mi64_cmpult(r, y, lenX), "Remainder should be < modulus!");
			}
			// At this point are done with x, so set low word of quotient array and clear rest:
			if(q) {
				mi64_clear(q, lenX);	q[0] = itmp64;
			}
			free((void *)yinv);
		} else {
			mi64_div_binary(x,y,lenX,lenY,q,r);
		}
	#if MI64_DIV_MONT
		if(dbg)	// Compute result using slow binary-div algo, use that as reference
		{
			ASSERT(HERE, mi64_cmp_eq(qref,q,lenX), "bzzt!");
			ASSERT(HERE, mi64_cmp_eq(rref,r,lenY), "bzzt!");
			free((void *)qref);
			free((void *)rref);
		}
	#endif
		return;
	}

	nshift = mi64_trailz(y,lenD);
	if(nshift)
	{
		w = hi + lenD;
		// If need to right-justify x and y, make a copy of x and work with that.
		if(q) {		// If q-array supplied, use that for (x >> nshift):
			v = q;
		} else {	// Otherwise allocate scratch space
			v = (uint64 *)calloc((lenX), sizeof(uint64));	ASSERT(HERE, v != 0x0, "alloc fail!");
		}
		nws = (nshift+63)>>6;	// Any partial word counts as 1 here
		nbs = nshift&63;
		// Save the bottom nshift bits of x (we work with copy of x saved in v) prior to right-shifting:
		rem_save = (uint64 *)calloc(nws, sizeof(uint64));
		mi64_set_eq(rem_save, x, nws);
		mask = ((uint64)-1 >> (64 - nshift));
		rem_save[nws-1] &= mask;
		mi64_shrl(x,v,nshift,lenX);
		mi64_shrl(y,w,nshift,lenD);
	#if MI64_DIV_MONT
		if(dbg && (lenX <= 50))printf("x >> %u = %s\n",nshift, &s0[convert_mi64_base10_char(s0, v, lenD)]);
		if(dbg)printf("y >> %u = %s\n",nshift, &s0[convert_mi64_base10_char(s0, w, lenD)]);
	#endif
	} else {
		v = x;	// Otherwise just point v at x, since from here on v is treated as read-only (even if not declared so)
		w = y;	// Do similarly for w and y.
		// GCC: "warning: assignment discards qualifiers from pointer target type"; we wish C had a C++ - style <const_cast>.
	}
	// We don't bother trimming leading zeros in x, only the divisor y:
	lenD = mi64_getlen(w, lenD);

	/*
	Find modular inverse (mod 2^nbits) of w in preparation for modular multiply.
	w must be odd for Montgomery-style modmul to work.

	Init yinv = 3*w ^ 2. This formula returns the correct bottom 4 bits of yinv,
	and we double the number of correct bits on each of the subsequent iterations.
	*/
	ASSERT(HERE, (w[0] & (uint64)1) == 1, "modulus must be odd!");
	ybits = lenD << 6;
	log2_numbits = ceil(log(1.0*ybits)/log(2.0));
	DBG_ASSERT(HERE, (w[0] & (uint64)1) == 1, "w must be odd!");
	mi64_clear(yinv, lenD);
	yinv[0] = (w[0] + w[0] + w[0]) ^ (uint64)2;

	/* Newton iteration involves repeated steps of form

		yinv = yinv*(2 - w*yinv);

	Number of significant bits at the bottom doubles on each iteration, starting from 4 for the initial seed
	defined as yinv_0 = 3*w ^ 2. The doubling continues until we reach the bitwidth set by the MULL operation.
	*/
	for(j = 2; j < 6; j++)	/* At each step, have 2^j correct low-order bits in yinv */
	{
		lo64 = w[0]*yinv[0];
		yinv[0] = yinv[0]*((uint64)2 - lo64);
	}

	/* Now that have bottom 64 = 2^6 bits of yinv, do as many more Newton iterations as needed to get the full [ybits] of yinv.
	Number of good bits doubes each iteration, and we only need to compute the new good bits, via

		yinv.hi = MULL(-yinv.lo, MULL(w.hi, yinv.lo) + UMULH(w.lo, yinv.lo))	//========== do this optimization later===========

	where lo = converged bits, and hi = bits converged on current iteration
	*/
	i = 1;	// I stores number of converged 64-bit words
	for(j = 6; j < log2_numbits; j++, i <<= 1)
	{
		mi64_mul_vector_lo_half(w, yinv,tmp, lenD);
		mi64_nega              (tmp,tmp, lenD);
		bw = mi64_add_scalar(tmp, 2ull,tmp, lenD);	DBG_ASSERT(HERE, !bw, "");
		mi64_mul_vector_lo_half(yinv,tmp, yinv, lenD);
	}
	// Check the computed inverse:
	mi64_mul_vector_lo_half(w, yinv, tmp, lenD);
	ASSERT(HERE, mi64_cmp_eq_scalar(tmp, 1ull, lenD), "Bad Montmul inverse!");
#if MI64_DIV_MONT
	if(dbg)printf("yinv = %s\n", &s0[convert_mi64_base10_char(s0, yinv, lenD)]);
#endif

	// Process the x-array data in lenD-sized chunks. If last chunk < lenD-sized it gets special handling.
	// Number of lenD-sized chunks to do (incl. 0-padding of x if its number of machine words not an exact multiple of lenD):
	lenW = (lenX+lenD-1)/lenD;
	for(i = 0; i < lenX-lenD+1; i += lenD)
	{
	#if MI64_DIV_MONT
	//	if(dbg)printf("i = %u; v = %s\n", i,&s0[convert_mi64_base10_char(s0, v+i, lenD)]);
	#endif
		/* Add w if had a borrow - Since we low-half multiply tmp by yinv below and w*yinv == 1 (mod 2^64),
		   can simply add 1 to tmp*yinv result instead of adding w to the difference here: */
		bw = mi64_sub(v+i,hi,tmp,lenD);

		// Compute expected value of low-half of MUL_LOHI for sanity-checking
		if(bw) {
			mi64_add(tmp,w,itmp,lenD);	// itmp = tmp + w
		} else {
			mi64_set_eq(itmp,tmp,lenD);	// itmp = tmp
		}
	#if MI64_DIV_MONT
	//	if(dbg)printf("v-cy = %s, bw = %llu\n", &s0[convert_mi64_base10_char(s0, tmp, lenD)], bw);
	#endif

		// Now do the Montgomery mod: cy = umulh( w, mull(tmp, yinv) );
		mi64_mul_vector_lo_half(tmp,yinv,tmp, lenD);	// tmp = tmp*yinv + bw;
	#if MI64_DIV_MONT
	//	if(dbg)printf("MULL = %s\n", &s0[convert_mi64_base10_char(s0, tmp, lenD)]);
	#endif
		tmp[0] += bw;	ASSERT(HERE, tmp[0] >= bw, "bw!");	// Since bw = 0 or 1, check that bw=1 does not propagate is simply check (sum >= bw)

		// Do double-wide product. Fast-divisibility test needs just high half (stored in hi); low half (lo) useful to extract true-mod
		mi64_mul_vector(tmp,lenD,w,lenD,lo, (uint32*)&j);	// lo:hi = MUL_LOHI(q, tmp)

	#if MI64_DIV_MONT
	//	if(dbg)printf("  lo = %s\n", &s0[convert_mi64_base10_char(s0,   lo, lenD)]);
	//	if(dbg)printf("  hi = %s\n", &s0[convert_mi64_base10_char(s0,   hi, lenD)]);
	#endif

		if(!mi64_cmp_eq(lo,itmp,lenD))
		{
		#if MI64_DIV_MONT
		//	if(dbg)printf("itmp = %s\n", &s0[convert_mi64_base10_char(s0, itmp, lenD)]);
		#endif
			ASSERT(HERE, 0, "Low-half product check mismatch!");
		}
	}

	// Last term gets special handling:
	// Zero-pad the final x-array section (by copying into cy and zeroing remaining higher-order terms of cy) if lenX != 0 (mod lenD):
	j = lenX-i;
	if(j) {
		// Set cy = {x[i],x[i+1],...,x[lenX-1],0,...,0}
		mi64_set_eq(cy,v+i,j);
	#if MI64_DIV_MONT
	//	if(dbg)printf("i+ = %u; v = %s\n", i,&s0[convert_mi64_base10_char(s0, cy, j)]);	// use 'i+' here to indicate this is the post-loop code
	#endif
		for(i = j; i < lenD; i++) {
			cy[i] = 0ull;
		}
		// only need the low half of the final product, which we can obtain sans explicit MUL:
		bw = mi64_sub(cy,hi,cy,lenD);
		if(bw) {
			mi64_add(cy,w,cy,lenD);	// cy += w
		}
	} else {
		mi64_set_eq(cy,lo,lenD);
	}
#if MI64_DIV_MONT
//	if(dbg)printf("MR = %s\n", &s0[convert_mi64_base10_char(s0, cy, lenD)]);	// MR = "Montgomery remainder"
#endif

//----------------------------------

	// Prepare to transform back out of "Montgomery space" ... first compute B^2 mod q.
	// B^2 overflows our double-wide scratch [lo]-array field, so compute B^2/2 mod q...
	mi64_clear(lo,2*lenD);	lo[2*lenD-1] = 0x8000000000000000ull;
	mi64_div_binary(lo,w,2*lenD,lenD,0,tmp);	// B^2/2 mod q returned in tmp
	// ...and mod-double the result:
	itmp64 = mi64_mul2(tmp,tmp,lenD);
	if(itmp64 || mi64_cmpugt(tmp,w,lenD)) {
		mi64_sub(tmp,w,tmp,lenD);
	}
#if MI64_DIV_MONT
if(dbg) {
	printf("B^2 mod q = %s\n", &s0[convert_mi64_base10_char(s0, tmp, lenD)]);
}
#endif
	/* tmp holds B^2 mod q - Now compute sequence of powers needed to obtain B^len mod q via Montgomery-muls.
	See the 64-bit-modulus-specialized version of this routine in mi64_div_by_scalar64() for details.
	*/
	p = lenW;

	// If p == 2, tmp already contains the needed power
	if(p == 3) {
		MONT_SQR_N(tmp,lo,hi,w,yinv, tmp,lenD);
	}
	else if(p > 3)
	{
		/*
		We always start with p = 2 and M-square that to get p = 3:
		*/
		MONT_SQR_N(tmp,lo,hi,w,yinv,itmp,lenD);

		n = 0;		// Init the bitstring
		for(j = 0; p > 5; j++)		// j counts number of bits processed
		{
			BIT_SETC(n,j,IS_EVEN(p));
			// Each M-mul includes another inverse power B^(-1) with the product, so we add 1 to the
			// current power p after each halving step here to account for that:
			p = (p >> 1) + 1;
		}
		ASSERT(HERE, j <= 32, "Need 64-bit bitstring!");
		/*
		Now do the needed powering. We always start with p = 2 and M-square that to get p = 3:
		*/
		MONT_SQR_N(tmp,lo,hi,w,yinv,itmp,lenD);	// tmp has p = 2, itmp has p = 3

		/* Starting off we have the following 2 choices:
			A. Mp(2,3) -> p=4;
			B. Mp(3,3) -> p=5,
		*/
		if(p == 4) {
			MONT_MUL_N(itmp,tmp,lo,hi,w,yinv,tmp,lenD);
		} else if(p == 5) {
			MONT_SQR_N(itmp,lo,hi,w,yinv,tmp,lenD);
		} else {
			ASSERT(HERE, 0,"Bad starting value for power p!");
		}
		for(i = j-1; i >= 0; i--)
		{
			if(BIT_TEST(n,i)) {
				mi64_set_eq(itmp,tmp,lenD);	// 'itmp = tmp', but copy *data*, not the pointers :)
				MONT_UNITY_MUL_N(itmp,w,yinv,itmp,lenD);	// Reduce power of B by 1 in one of the 2 multiplicands...
				MONT_MUL_N(itmp,tmp,lo,hi,w,yinv,tmp,lenD);	// ...and multiply `em.
			} else {
				MONT_SQR_N(tmp,lo,hi,w,yinv,tmp,lenD);
			}
		}
	}

	/*
	Now multiply the Montgomery residue from the mod-loop by the mod-power of the base.
	Properly scaled, right-justified remainder output in cy ... need the r-j remainder
	if doing a further quotient computation anyway, and since function allows x/r-pointers
	to refer to same memloc, keep remainder in local-array cy for now, defer copy-to-r until last:
	*/
	MONT_MUL_N(cy,tmp,lo,hi,w,yinv,cy,lenD);
#if MI64_DIV_MONT
	if(dbg && lenW > 2)printf("B^%u mod q = %s\n", lenW,&s0[convert_mi64_base10_char(s0, tmp, lenD)]);
	if(dbg)printf("Remainder = %s\n", &s0[convert_mi64_base10_char(s0,cy, lenD)]);
#endif

//----------------------------------

	// If q-array supplied, compute quotient and return it in that:
	if(q) {
	#if MI64_DIV_MONT
		if(dbg)printf("Computing quotient...\n");
	#endif
		// If even modulus, right-justified copy of input array already in v.
		// Now can use a simple loop and a sequence of word-size MULLs to obtain quotient.
		// Fusing the functionality of mi64_sub_scalar and the quotient extraction is fastest here:
		bw = 0;
		mi64_set_eq(hi,cy,lenD);	// Copy remainder (stored in cy) into hi-array
		for(i = 0; i < lenX-lenD+1; i += lenD)
		{
		#if MI64_DIV_MONT
			if(dbg)printf("i = %u; v = %s\n", i,&s0[convert_mi64_base10_char(s0, v+i, lenD)]);
		#endif
			bw = mi64_sub(v+i,hi,tmp,lenD);	// tmp = x[i] - (bw+cy);

			// Compute expected value of low-half of MUL_LOHI for sanity-checking
			mi64_set_eq(itmp,tmp,lenD);	// itmp = tmp

			// Now do the Montgomery mod: cy = umulh( y, mull(tmp, yinv) );
			mi64_mul_vector_lo_half(tmp,yinv,tmp, lenD);	// tmp = tmp*yinv + bw;
		#if MI64_DIV_MONT
			if(dbg)printf("tmp*yinv = %s, bw = %llu\n", &s0[convert_mi64_base10_char(s0, tmp, lenD)], bw);
		#endif
			// Do double-wide product. Fast-divisibility test needs just high half (stored in hi); low half (lo) useful to extract true-mod
			mi64_mul_vector(tmp,lenD,w,lenD,lo, (uint32*)&j);	// lo:hi = MUL_LOHI(q, tmp); cy is in hi half
			hi[0] += bw;	ASSERT(HERE, hi[0] >= bw, "bw!");// (cy + bw); Since bw = 0 or 1, check that bw=1 does not propagate is simply check (sum >= bw)

		#if MI64_DIV_MONT
			if(dbg)printf("  lo = %s\n", &s0[convert_mi64_base10_char(s0,   lo, lenD)]);
			if(dbg)printf("  hi = %s\n", &s0[convert_mi64_base10_char(s0,   hi, lenD)]);
		#endif

			if(!mi64_cmp_eq(lo,itmp,lenD))
			{
			#if MI64_DIV_MONT
				printf("itmp = %s\n", &s0[convert_mi64_base10_char(s0, itmp, lenD)]);
			#endif
				ASSERT(HERE, 0, "Low-half product check mismatch!");
			}
			mi64_set_eq(q+i,tmp,lenD);	// Equivalent to the y[i] = tmp step of the scalar routine
		}
		// Any words remaining due to lenD not exactly dividing lenX are guaranteed to end up = 0, use that as a sanity check:
		j = lenX-i;
		if(j) {
			// Check cy = {v[i],v[i+1],...,v[lenX-1],0,...,0}
			ASSERT(HERE, mi64_cmp_eq(hi,v+i,j), "cy check!");
			mi64_clear(q+i,j);	// Do after above check since v may == q
			for(i = j; i < lenD; i++) {
				ASSERT(HERE, hi[i] == 0ull, "cy check!");
			}
		}
	#if MI64_DIV_MONT
		if(dbg) {
		//	if(lenX <= 50)printf("q = %s\n", &s0[convert_mi64_base10_char(s0, q, lenX)]);
			printf("q = %s\n", &str_10k[__convert_mi64_base10_char(str_10k, 10<<10, q, lenX)]);
		}
	#endif
	}

	// Copy remainder from local cy-array to output r-array (which may == x).
	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift)
	{
		// rem = (rem << nshift) + rem_save:
		mi64_shl(cy,r,nshift,lenD);
		// No carryout here since we are filling in the just-vacated low bits of rem with rem_save:
		mi64_add(r,rem_save,r,nws);
	} else {
		mi64_set_eq(r,cy,lenD);
	}

//----------------------------------

#if MI64_DIV_MONT
if(dbg)	// Compute result using slow binary-div algo, use that as reference
{
	if(!mi64_cmp_eq(rref,r,lenY)) {
		printf("rref = %s\n", &s0[convert_mi64_base10_char(s0, rref, lenD)]);
		printf("rewm = %s\n", &s0[convert_mi64_base10_char(s0, r   , lenD)]);
		printf("bzzt!\n");
	}
	if(!mi64_cmp_eq(qref,q,lenX)) {
		printf("qref = %s\n", &s0[convert_mi64_base10_char(s0, qref, lenX)]);
		printf("qewm = %s\n", &s0[convert_mi64_base10_char(s0, q   , lenX)]);
		printf("bzzt!\n");
	}

	free((void *)qref);
	free((void *)rref);
}
#endif
	return;
}

/*
Slow bit-at-a-time method to get quotient q = x/y and remainder r = x%y.
[optional] Output Q-array, if supplied, must have dimension at least as large as that of X-array, even if actual  quotient is smaller.
[required] Output R-array **** must have dimension at least as large as that of Y-array ****, even if actual remainder is smaller.
*/
#if MI64_DEBUG
	#define MI64_DIV_DBG	0
#endif
void mi64_div_binary(const uint64 x[], const uint64 y[], uint32 lenX, uint32 lenY, uint64 q[], uint64 r[])
{
#if MI64_DIV_DBG
	uint32 dbg = 0;//(lenX== 6 && lenY==3) && STREQ(&s0[convert_mi64_base10_char(s0, y, lenY)], "53625112691923843508117942311516428173021903300344567");
#endif
	int i, nshift;
	uint32 lz_x, lz_y, xlen, ylen, max_len;
	// pointers to local storage:
	static uint64 *xloc = 0x0, *yloc = 0x0;
	static uint64 *scratch = 0x0;	// "base pointer" for local storage shared by all of the above subarrays
	static int lens = 0;	// # of 64-bit ints allocated for current scratch space

#if MI64_DIV_DBG
if(dbg) {
	printf("mi64_div_binary: x = %s, y = %s\n",&s0[convert_mi64_base10_char(s0, x, lenX)],&s1[convert_mi64_base10_char(s1, y, lenY)]);
}
#endif
	ASSERT(HERE, lenX && lenY, "illegal 0 dimension!");
	ASSERT(HERE, x && y && r, "At least one of X, Y, R is null!");
	ASSERT(HERE, x != y, "X and Y arrays overlap!");
	ASSERT(HERE, r != y, "Y and Rem arrays overlap!");
	ASSERT(HERE, q != x && q != y && q != r, "Quotient array overlaps one of X, Y ,Rem!");

	/* Init Q = 0; don't do similarly for R since we allow X and R to point to same array: */
	if(q) {
		mi64_clear(q, lenX);
	}

	/* And now find the actual lengths of the divide operands and use those for the computation: */
	xlen = mi64_getlen(x, lenX);
	ylen = mi64_getlen(y, lenY);
	ASSERT(HERE, ylen != 0, "divide by 0!");
	/* Allocate the needed auxiliary storage: */
	if(lens < (xlen << 1)) {
		lens = (xlen << 1);	// Alloc yloc same as x to allow for left-justification of y-copy
		scratch = (uint64 *)realloc(scratch, lens*sizeof(uint64));	ASSERT(HERE, scratch != 0x0, "alloc fail!");
	}

	max_len = MAX(xlen, ylen);
	if((xlen < ylen) || ((xlen == ylen) && mi64_cmpult(x, y, xlen)))
	{
		return;		/* If x < y, return x */
	}

	xloc = scratch       ;	mi64_set_eq(xloc, x, xlen);
	yloc = scratch + xlen;	mi64_set_eq(yloc, y, ylen);	mi64_clear(yloc+ylen,xlen-ylen);

	lz_x = mi64_leadz(xloc, max_len);
	lz_y = mi64_leadz(yloc, max_len);
	nshift = lz_y - lz_x;
	DBG_ASSERT(HERE, nshift >= 0, "nshift < 0");

	/* Left-justify the modulus (copy) y to match x's leading bit: */
	mi64_shl(yloc, yloc, nshift, max_len);	ylen = max_len;
	for(i = nshift; i >= 0; --i)
	{
	#if MI64_DIV_DBG
	//	if(dbg)printf("I = %3d: r = %s, yshift = %s\n", i,&s0[convert_mi64_base10_char(s0, xloc, max_len)],&s1[convert_mi64_base10_char(s1, yloc, max_len)]);
	#endif
		if(mi64_cmpuge(xloc, yloc, max_len))
		{
			DBG_ASSERT(HERE, xlen == max_len,"xlen != max_len");
			mi64_sub(xloc, yloc, xloc, max_len);	/* r -= yshift */
			DBG_ASSERT(HERE, mi64_cmpult(xloc, yloc, max_len),"r >= yshift");
			xlen = mi64_getlen(xloc, max_len);
			if(q) {
				mi64_set_bit(q, i);
			}
		}
		if(i > 0)
		{
			mi64_shrl(yloc, yloc, 1, ylen);
			ylen = mi64_getlen(yloc, ylen);
		}
		max_len = MAX(xlen, ylen);
	}
	// Remainder in xloc - do some sanity checks piror to copying into r[]:
	xlen = mi64_getlen(xloc, lenX);
	ASSERT(HERE, xlen <= ylen && mi64_cmpugt(y,xloc,ylen), "Remainder should be < modulus!");
	mi64_set_eq(r, xloc, ylen);	mi64_clear(r+ylen,lenY-ylen);

	/* Final value of yloc is unchanged from its (unshifted) starting value == y */
	ASSERT(HERE, mi64_cmp_eq(yloc,y,ylen), "Final value of y-copy differs from original!");
#if MI64_DIV_DBG
if(dbg) {
	if(q)printf("mi64_div_binary: quotient  = %s\n",&s0[convert_mi64_base10_char(s0, q, lenX)]);
	printf("mi64_div_binary: remainder = %s\n",&s0[convert_mi64_base10_char(s0, r, lenX)]);
	printf("\n");
}
#endif
}

/* Fast is-divisible-by-32-bit-scalar using Montgomery modmul and right-to-left modding: */
int mi64_is_div_by_scalar32(const uint32 x[], uint32 q, uint32 len)
{
	uint32 i,j,nshift,dlen,qinv,tmp,cy;

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

/* Same as above, but assumes q and its modular inverse have been precomputed: */
int		mi64_is_div_by_scalar32p(const uint32 x[], uint32 q, uint32 qinv, uint32 len)
{
	uint32 i,dlen,tmp,cy;

	DBG_ASSERT(HERE, qinv == qinv*((uint32)2 - q*qinv), "mi64_is_div_by_scalar32p: bad qinv!");
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

/* Same as above, but assumes q and its modular inverse have been precomputed: */
int		mi64_is_div_by_scalar32p_x8(
	const uint32 a[],
	const uint32 b[],
	const uint32 c[],
	const uint32 d[],
	const uint32 e[],
	const uint32 f[],
	const uint32 g[],
	const uint32 h[],
	uint32 q, uint32 qinv, uint32 len, uint32*sum)
{
	int retval=0;
	uint32 tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7;
	cy0 = cy1 = cy2 = cy3 = cy4 = cy5 = cy6 = cy7 = (uint32)0;

	DBG_ASSERT(HERE, qinv == qinv*((uint32)2 - q*qinv), "mi64_is_div_by_scalar32p: bad qinv!");

	tmp0 = a[0] * qinv;
	tmp1 = b[0] * qinv;
	tmp2 = c[0] * qinv;
	tmp3 = d[0] * qinv;
	tmp4 = e[0] * qinv;
	tmp5 = f[0] * qinv;
	tmp6 = g[0] * qinv;
	tmp7 = h[0] * qinv;

  #ifdef MUL_LOHI32_SUBROUTINE
	cy0 = __MULH32(q,tmp0);
	cy1 = __MULH32(q,tmp1);
	cy2 = __MULH32(q,tmp2);
	cy3 = __MULH32(q,tmp3);
	cy4 = __MULH32(q,tmp4);
	cy5 = __MULH32(q,tmp5);
	cy6 = __MULH32(q,tmp6);
	cy7 = __MULH32(q,tmp7);
  #else
	MULH32(q,tmp0, cy0);
	MULH32(q,tmp1, cy1);
	MULH32(q,tmp2, cy2);
	MULH32(q,tmp3, cy3);
	MULH32(q,tmp4, cy4);
	MULH32(q,tmp5, cy5);
	MULH32(q,tmp6, cy6);
	MULH32(q,tmp7, cy7);
  #endif

	tmp0 = a[1] - cy0;
	tmp1 = b[1] - cy1;
	tmp2 = c[1] - cy2;
	tmp3 = d[1] - cy3;
	tmp4 = e[1] - cy4;
	tmp5 = f[1] - cy5;
	tmp6 = g[1] - cy6;
	tmp7 = h[1] - cy7;
	/* Add q if had a borrow - Since we low=half multiply tmp by qinv below and q*qinv == 1 (mod 2^32),
	   can simply add 1 to tmp*qinv result instead of adding q to the difference here:
	*/
	cy0 = (cy0 > a[1]);
	cy1 = (cy1 > b[1]);
	cy2 = (cy2 > c[1]);
	cy3 = (cy3 > d[1]);
	cy4 = (cy4 > e[1]);
	cy5 = (cy5 > f[1]);
	cy6 = (cy6 > g[1]);
	cy7 = (cy7 > h[1]);
	/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
	tmp0 *= qinv;
	tmp1 *= qinv;
	tmp2 *= qinv;
	tmp3 *= qinv;
	tmp4 *= qinv;
	tmp5 *= qinv;
	tmp6 *= qinv;
	tmp7 *= qinv;

	tmp0 += cy0;
	tmp1 += cy1;
	tmp2 += cy2;
	tmp3 += cy3;
	tmp4 += cy4;
	tmp5 += cy5;
	tmp6 += cy6;
	tmp7 += cy7;

  #ifdef MUL_LOHI32_SUBROUTINE
	cy0 = __MULH32(q,tmp0);
	cy1 = __MULH32(q,tmp1);
	cy2 = __MULH32(q,tmp2);
	cy3 = __MULH32(q,tmp3);
	cy4 = __MULH32(q,tmp4);
	cy5 = __MULH32(q,tmp5);
	cy6 = __MULH32(q,tmp6);
	cy7 = __MULH32(q,tmp7);
  #else
	MULH32(q,tmp0, cy0);
	MULH32(q,tmp1, cy1);
	MULH32(q,tmp2, cy2);
	MULH32(q,tmp3, cy3);
	MULH32(q,tmp4, cy4);
	MULH32(q,tmp5, cy5);
	MULH32(q,tmp6, cy6);
	MULH32(q,tmp7, cy7);
  #endif

	retval += (cy0 == 0);
	retval += (cy1 == 0) << 1;
	retval += (cy2 == 0) << 2;
	retval += (cy3 == 0) << 3;
	retval += (cy4 == 0) << 4;
	retval += (cy5 == 0) << 5;
	retval += (cy6 == 0) << 6;
	retval += (cy7 == 0) << 7;
	*sum = cy0+cy1+cy2+cy3+cy4+cy5+cy6+cy7;
	return retval;
}

#if(defined(USE_SSE2) && defined(COMPILER_TYPE_MSVC))
/* Same as above, but assumes q and its modular inverse have been precomputed: */
int		mi64_is_div_by_scalar32p_x8_SSE2(
	const uint32 a[],
	const uint32 b[],
	const uint32 c[],
	const uint32 d[],
	const uint32 e[],
	const uint32 f[],
	const uint32 g[],
	const uint32 h[],
	uint32 q, uint32 qinv, uint32 len, uint32*sum)
{
	int retval=0;
	static uint32*cy;
#define OLD_ALGO	0
#if OLD_ALGO
	static uint64 *sm_arr,*sm_ptr;
	static uint64 *ablo,*cdlo,*eflo,*ghlo,*abhi,*cdhi;
	uint64 itmp64;
	static int first_entry=TRUE;
	if(first_entry)
	{
		first_entry=FALSE;
		sm_arr = ALLOC_UINT64(20);
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
		ablo = sm_ptr + 0x00;			/* Each of these needs 2 64-bit slots... */
		cdlo = sm_ptr + 0x02;
		eflo = sm_ptr + 0x04;
		ghlo = sm_ptr + 0x06;
		abhi = sm_ptr + 0x08;
		cdhi = sm_ptr + 0x0a;
		cy = (uint32*)(sm_ptr + 0xc);	/* This needs 4x64-bits, so allocate 2*6+4 = 16 64-bit slots, plus 4 slots to allow for the alignment */
	}
	DBG_ASSERT(HERE, qinv == qinv*((uint32)2 - q*qinv), "mi64_is_div_by_scalar32p: bad qinv!");
	*ablo++ = a[0];							*ablo-- = b[0];
	*cdlo++ = c[0];							*cdlo-- = d[0];
	*eflo++ = e[0];							*eflo-- = f[0];
	*ghlo++ = g[0];							*ghlo-- = h[0];
	 itmp64 = a[1]+(((uint64)c[1]) << 32);
	*abhi++ = itmp64;
	 itmp64 = b[1]+(((uint64)d[1]) << 32);
	*abhi-- = itmp64;
	 itmp64 = e[1]+(((uint64)g[1]) << 32);
	*cdhi++ = itmp64;
	 itmp64 = f[1]+(((uint64)h[1]) << 32);
	*cdhi-- = itmp64;
	__asm	mov	eax, ablo
	__asm	lea	ebx, q
	__asm	lea	ecx, qinv
	__asm	movaps	xmm0,[eax     ]	/* ablo: d3210 = [  0|blo|  0|alo] */
	__asm	movaps	xmm1,[eax+0x10]	/* cdlo: d3210 = [  0|dlo|  0|clo] */
	__asm	movaps	xmm2,[eax+0x20]	/* eflo: d3210 = [  0|flo|  0|elo] */
	__asm	movaps	xmm3,[eax+0x30]	/* ghlo: d3210 = [  0|hlo|  0|glo] */
#else
	DBG_ASSERT(HERE, qinv == qinv*((uint32)2 - q*qinv), "mi64_is_div_by_scalar32p: bad qinv!");
	DBG_ASSERT(HERE, ((uint32)&a[0] & 0x3f) == 0, "A-array not 64-byte aligned!");
	__asm	mov	eax, a	// Assumes inputs a,b,c,d,... are 64-bit separated and &a[0} is 64-byte aligned
	__asm	lea	ebx, q
	__asm	lea	ecx, qinv
	__asm	movaps	xmm0,[eax     ]	/* ab: d3210 = [bhi|blo|ahi|alo] */
	__asm	movaps	xmm1,[eax+0x10]	/* cd: d3210 = [dhi|dlo|chi|clo] */
	__asm	movaps	xmm2,[eax+0x20]	/* ef: d3210 = [fhi|flo|ehi|elo] */
	__asm	movaps	xmm3,[eax+0x30]	/* gh: d3210 = [hhi|hlo|ghi|glo] */
	__asm	movaps	xmm6,xmm0	// Circularly-permute [4,6,7] -> [6,7,4] here so the 2 packed outputs end up in xmm6,7
	__asm	movaps	xmm5,xmm1
	__asm	movaps	xmm7,xmm2
	__asm	movaps	xmm4,xmm3
	__asm	psrlq	xmm6, 32		/* d3210 = [  0|bhi|  0|ahi] */
	__asm	psrlq	xmm5, 32		/* d3210 = [  0|dhi|  0|chi] */
	__asm	psrlq	xmm7, 32		/* d3210 = [  0|fhi|  0|ehi] */
	__asm	psrlq	xmm4, 32		/* d3210 = [  0|hhi|  0|ghi] */
	__asm	psllq	xmm5, 32		/* d3210 = [dhi|  0|chi|  0] */
	__asm	psllq	xmm4, 32		/* d3210 = [hhi|  0|ghi|  0] */
	__asm	paddd	xmm6,xmm5		/* d3210 = [dhi|bhi|chi|ahi], xmm5 FREE */
	__asm	paddd	xmm7,xmm4		/* d3210 = [hhi|fhi|ghi|ehi], xmm4 FREE */
#endif
	__asm	movd	xmm4,[ebx]
	__asm	movd	xmm5,[ecx]
	__asm	pshufd	xmm4,xmm4,0x44	/* Broadcast q    to slots 0,2 of xmm4 */
	__asm	pshufd	xmm5,xmm5,0x44	/* Broadcast qinv to slots 0,2 of xmm5 */
	/* (a-h)[0]*qinv; Alas SSE2 has no 32-bit low-half packed MUL, so use 32x32->64 -bit and discard high halves */
	__asm	pmuludq	xmm0,xmm5
	__asm	pmuludq	xmm1,xmm5
	__asm	pmuludq	xmm2,xmm5
	__asm	pmuludq	xmm3,xmm5
	/* cy[0-7] = MULH32(tmp[0-7]*q) - high halves of above MULQs automatically get overwritten: */
	__asm	pmuludq	xmm0,xmm4
	__asm	pmuludq	xmm1,xmm4
	__asm	pmuludq	xmm2,xmm4
	__asm	pmuludq	xmm3,xmm4
	__asm	psrlq	xmm0, 32		/* d3210 = [  0|cy1|  0|cy0] */
	__asm	psrlq	xmm1, 32		/* d3210 = [  0|cy3|  0|cy2] */
	__asm	psrlq	xmm2, 32		/* d3210 = [  0|cy5|  0|cy4] */
	__asm	psrlq	xmm3, 32		/* d3210 = [  0|cy7|  0|cy6] */
	__asm	psllq	xmm1, 32		/* d3210 = [cy3|  0|cy2|  0] */
	__asm	psllq	xmm3, 32		/* d3210 = [cy7|  0|cy6|  0] */
	__asm	paddd	xmm0,xmm1		/* d3210 = [cy3|cy1|cy2|cy0], xmm1 FREE */
	__asm	paddd	xmm2,xmm3		/* d3210 = [cy7|cy5|cy6|cy4], xmm3 FREE */
#if OLD_ALGO
	/* 		tmp[0-7] = (a-h)[1] - cy[0-7]; */
	__asm	mov	edx, abhi
	__asm	movaps	xmm6,[edx     ]	/* acbd[1]: d3210 = [dhi|bhi|chi|ahi] */
	__asm	movaps	xmm7,[edx+0x10]	/* egfh[1]: d3210 = [hhi|fhi|ghi|ehi] */
#endif
	__asm	movaps	xmm3,xmm6		/* Copy of acbd[1] */
	__asm	movaps	xmm1,xmm7		/* Copy of efgh[1] */
	__asm	psubd	xmm6,xmm0		/* acbd[1] - cy0213, xmm0 FREE */
	__asm	psubd	xmm7,xmm2		/* egfh[1] - cy4657, xmm2 FREE */
	__asm	movaps	xmm2,xmm6		/* Copy of acbd[1] - cy0213 */
	__asm	movaps	xmm0,xmm7		/* Copy of efgh[1] - cy4657 */
	/* Had a borrow? Frickin' SSE2 only gives us signed packed-integer compares,
	so need to emulate unsigned (x > y) via signed (x ^ 0x80000000) < (y ^ 0x80000000): */
	__asm	pcmpeqd	xmm4,xmm4		/* All 1s  - will need to restore q to this register later */
	__asm	pslld	xmm4, 31		/* 4-way 0x80000000 */
	__asm	pxor	xmm6,xmm4		/* (acbd[1]-cy0213) ^ 0x80000000 */
	__asm	pxor	xmm7,xmm4		/* (egfh[1]-cy4657) ^ 0x80000000 */
	__asm	pxor	xmm3,xmm4		/* (acbd[1]) ^ 0x80000000 */
	__asm	pxor	xmm1,xmm4		/* (egfh[1]) ^ 0x80000000 */
	__asm	pcmpgtd	xmm6,xmm3		/* cy0213 = (acbd[1]-cy0213) > abcd[1], xmm3 FREE */
	__asm	pcmpgtd	xmm7,xmm1		/* cy4657 = (egfh[1]-cy4657) > efgh[1], xmm1 FREE */
	__asm	pshufd	xmm3,xmm2,0x31	/* xmm2 = [----|tmp1|----|tmp0], xmm3 = [----|tmp3|----|tmp2], don't care what's in ---- slots */
	__asm	pshufd	xmm1,xmm0,0x31	/* xmm0 = [----|tmp5|----|tmp4], xmm1 = [----|tmp7|----|tmp6], don't care what's in ---- slots */
	__asm	movd	xmm4,[ebx]		/* Restore q to xmm4 */
	__asm	pshufd	xmm4,xmm4,0x44	/* Broadcast q    to slots 0,2 of xmm4 */
	/* tmp[0-7]*qinv; Alas SSE2 has no 32-bit low-half packed MUL, so use 32x32->64 -bit and discard high halves */
	__asm	pmuludq	xmm3,xmm5
	__asm	pmuludq	xmm1,xmm5
	__asm	pmuludq	xmm2,xmm5
	__asm	pmuludq	xmm0,xmm5
	/* Add carries 01/45, scatter carries 23/67 into slots of 01/45, add those...Since SSE2 compare result is ~()ed, add really means sub: */
	__asm	psubd	xmm2,xmm6		/* xmm6 = [----|tmp1|----|tmp0], don't care what's in ---- slots */
	__asm	psubd	xmm0,xmm7		/* xmm7 = [----|tmp5|----|tmp4], don't care what's in ---- slots */
	__asm	pshufd	xmm6,xmm6,0x31
	__asm	pshufd	xmm7,xmm7,0x31
	__asm	psubd	xmm3,xmm6		/* xmm3 = [----|tmp3|----|tmp2], don't care what's in ---- slots */
	__asm	psubd	xmm1,xmm7		/* xmm1 = [----|tmp7|----|tmp6], don't care what's in ---- slots */
	/* cy[0-7] = MULH32(tmp[0-7]*q) - high halves of above MULQs automatically get overwritten: */
	__asm	pmuludq	xmm2,xmm4
	__asm	pmuludq	xmm0,xmm4
	__asm	pmuludq	xmm3,xmm4
	__asm	pmuludq	xmm1,xmm4
	__asm	psrlq	xmm2, 32		/* d3210 = [  0|cy1|  0|cy0] */
	__asm	psrlq	xmm0, 32		/* d3210 = [  0|cy5|  0|cy4] */
	__asm	psrlq	xmm3, 32		/* d3210 = [  0|cy3|  0|cy2] */
	__asm	psrlq	xmm1, 32		/* d3210 = [  0|cy7|  0|cy6] */
	__asm	pshufd	xmm2,xmm2,0x58	/* [  0|  0|cy1|cy0] */
	__asm	pshufd	xmm0,xmm0,0x58	/* [  0|  0|cy5|cy4] */
	__asm	pshufd	xmm3,xmm3,0x85	/* [cy3|cy2|  0|  0] */
	__asm	pshufd	xmm1,xmm1,0x85	/* [cy7|cy6|  0|  0] */
	__asm	paddd	xmm2,xmm3		/* d3210 = [cy3|cy1|cy2|cy0] */
	__asm	paddd	xmm0,xmm1		/* d3210 = [cy7|cy5|cy6|cy4] */
//__asm	movaps	[eax+0x60],xmm2
//__asm	movaps	[eax+0x70],xmm0
	__asm	pcmpgtd	xmm7,xmm7		/* All 0s */
	__asm	pcmpeqd	xmm2,xmm7		/* retval[0-3] */
	__asm	pcmpeqd	xmm0,xmm7		/* retval[4-7] */
	__asm	movmskps eax,xmm2		/* retval[0-3] */
	__asm	movmskps ebx,xmm0		/* retval[4-7] */
	__asm	shl		 ebx, 4		/* retval[4-7] << 4 */
	__asm	add		 eax,ebx	/* retval[0-7] */
	__asm	mov	retval,  eax
//  	*sum = cy[0]+cy[1]+cy[2]+cy[3]+cy[4]+cy[5]+cy[6]+cy[7];

	return retval;
}
#endif

uint32	mi64_is_div_by_scalar32_x4(const uint32 x[], uint32 q0, uint32 q1, uint32 q2, uint32 q3, uint32 len)
{
	uint32 i,j,nshift0,nshift1,nshift2,nshift3;
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
	uint32 i,j,nshift0,nshift1,nshift2,nshift3,nshift4,nshift5,nshift6,nshift7;
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

/*
Data: Integer n >= 0, unsigned 64-bit integers q; qinv, q odd and q * qinv == 1 (mod R), with R = 2^64 here.
Returns: The (n)th modular power of the twos-complement radix, R^n (mod q).
*/
#if MI64_DEBUG
	#define RADIX_POWER64	1
#endif

uint64 radix_power64(const uint64 q, const uint64 qinv, uint32 n)
{
#if MI64_DEBUG
	int dbg = q == 16357897499336320049ull;//87722769297534671ull;
#endif
	int i,j,bmap,p;	// i and j must be signed
	uint64 x,iquo,lo,hi,rem64,itmp64;
	double fquo,fqinv;
#ifdef RADIX_POWER64
	const uint64 two96[2] = {0ull, (1ull << 32)};
	uint64 rsqr,t128,ntry = 0, nbad = 0, nerr = 0;
	const double ILG2 = 1.0/log(2.0);
	const struct qfloat qtwo96 = {0x45F0000000000000ull, 0x0000000000000000ull};	// exp-field = 0x3ff + 96
	struct qfloat qx;
#endif
	ASSERT(HERE, (uint64)TWO48FLOAT == 0x0001000000000000ull, "Util.c const-init check fails!");
	DBG_ASSERT(HERE, (q & 1) == 0, "modulus must be odd!");
	if(n == 0) return 1ull;
	if(n == 1) return 0ull;
	// The minimum nontrivial power of B we need is B^2, so first obtain that.

#ifndef RADIX_POWER64
	x = q;
	// Floating-point (approximate) inverse:
	fqinv = 1.0/q;
#else
// test code for radix_power64 R^2 (mod q) computation:
x = 1;
for(;;)
{
	// If adding x/64 to x will overflow 64 bits, break:
	if((x >> 6) > (1ull << 58)) {
		break;
	}
	if(x < 1000) {
		x += 2;
	} else {
		x += (x >> 6);
		x += (1ull - (x&1));	// x must be odd
	}
	// Floating-point (approximate) inverse:
	fqinv = 1.0/x;
#endif
	if(q < (1ull << 48)) {	// q < 2^48
	// Method A: Use mod-doublings to get 2^72 (mod q), followed by 3 MONT_SQR64:
	#if 0
		itmp64 = 0x8000000000000000ull % q;	// 2^63 % q
		// 5 mod-doublings yield 2^68 (mod q)
		itmp64 = itmp64 + itmp64 - q;	itmp64 += (-((int64)itmp64 < 0) & q);
		itmp64 = itmp64 + itmp64 - q;	itmp64 += (-((int64)itmp64 < 0) & q);
		itmp64 = itmp64 + itmp64 - q;	itmp64 += (-((int64)itmp64 < 0) & q);
		itmp64 = itmp64 + itmp64 - q;	itmp64 += (-((int64)itmp64 < 0) & q);
		itmp64 = itmp64 + itmp64 - q;	itmp64 += (-((int64)itmp64 < 0) & q);
		MONT_SQR64(itmp64,q,qinv,itmp64);	// 2^(2*68-64) == 2^72 (mod q)
		MONT_SQR64(itmp64,q,qinv,itmp64);	// 2^80 (mod q)
		MONT_SQR64(itmp64,q,qinv,itmp64);	// 2^96 (mod q)
	#else
	// Method B uses fast floating-point mod:
		fquo = TWO48FLOAT*fqinv;
		iquo = (uint64)fquo;
		itmp64 = TWO48FLOAT - q*iquo;	// 2^48 % q, which may be as large as 48 bits
		SQR_LOHI64(itmp64, lo,hi);
		fquo = (hi*TWO64FLOAT + lo)*fqinv;
		iquo = (uint64)fquo;	// Need both of these for check below
		itmp64 = lo - q*iquo;	// Only need low 64 bits of difference
		// Floating-point method here fails in a tiny fraction of cases - see notes in "else" block for details:
		if(itmp64 > q)
		{
			// This check allows us to differentiate between incorrect upward-rounded and (rarer) downward-rounded cases:
			if(DNINT(fquo) == (double)iquo) {	// Incorrect   upward-rounded
				itmp64 += q;	//	printf("A...");
			} else {							// Incorrect downward-rounded
				itmp64 -= q;	//	printf("B...");
			}
		}
	#endif
	#ifdef RADIX_POWER64
	// Compute ref-value using stdlib mod:
		if(x < 0x10000ull) {	// If x < 2^16, 2^48/x > 2^32, hence cannot compute (2^48/x)^2 using 64-bit int; instead compute 2^32/x and modsqr-twice
			rsqr = (1ull << 32) % x;	rem64 = rsqr;
			rsqr = (rsqr*rsqr) % x;		ASSERT(HERE, rsqr < (1ull << 32), "gah!");
			rsqr = (rem64*rsqr) % x;	// 2^96 % x
			t128 = (rsqr*rsqr) % x;		// 2^128 % x
		} else {	// Compute 2^128 (mod x) via 3-step modmul: (2^48/x)^2 * (2^32/x) (mod x)
			rsqr = (1ull << 48) % x;
			rsqr = (rsqr*rsqr) % x;
			t128 = (rsqr << 32) % x;	// 2^128 % x
		}
		if(itmp64 != rsqr) {
			printf("radix_power64: <48 fails for x = %llu, float-result = %20llu, true remainder = %20llu\n",x, itmp64, rsqr);
		}
	#endif
	}
	else	// q in [2^48, 2^64)
	{
	#ifdef RADIX_POWER64
		++ntry;
	#endif
		fquo = TWO48FLOAT*TWO48FLOAT*fqinv;
		rem64  = (uint64)fquo;
		// Bottom 64 bits of 2^96 - q*(2^96/q)
		itmp64 = -(rem64 * x);
		// Floating-point method here fails in a tiny fraction of cases (which on x86 is highly dependent on build mode, likely due to in-register
		// precision effects: in my test, ~35000 of 10^9 fail for debug-build, only 19 for opt-build), but these failures easily spotted via this check.
		// If quotient has (exact) fractional part = 0.99...., may end up with fquo = ceiling(2^96/q) rather than floor; thus want x - (2^64 - itmp64):
		if(itmp64 > x)
		{
		#ifdef RADIX_POWER64
			++nbad;
		#endif
			// This check allows us to differentiate between incorrect upward-rounded and (rarer) downward-rounded cases:
			if(DNINT(fquo) == (double)rem64) {	// Incorrect   upward-rounded, e.g. fquo = 1084809392143.0001, exact = 1084809392142.999...
			//	printf("radix_power64: x = %20llu < itmp64 = (int64)%20lld, fquo = %20.4f, (double)rem64 = %20.4f\n",x,(int64)itmp64, fquo, (double)rem64);
				itmp64 += x;
			} else {							// Incorrect downward-rounded, e.g. fquo = 7344640876302.9990, exact = 7344640876303.0000002...
			//	printf("radix_power64: x = %20llu < itmp64 = (int64)%20lld, fquo = %20.4f *** Bad Downward ***\n",x,(int64)itmp64, fquo);
				itmp64 -= x;
			}
		#ifdef RADIX_POWER64
			mi64_div_binary(two96, &x, 2,1, 0x0, &rsqr);
			if(itmp64 != rsqr) {
				++nerr;
				printf("radix_power64: x = %20llu fails, float-result = %20llu, true remainder = %20llu\n",x, itmp64, rsqr);
				qx = qfdiv(qtwo96, i64_to_q(x));
				printf("    xbits = %5.2f, fquo = %25.8f, qquo = %25.8f, q_hi = ...%8X, q_lo = %16llX\n",log((double)x)*ILG2,fquo,qfdbl(qx),(uint32)qx.hi,qx.lo);
			//	ASSERT(HERE, 0, "Bad float-result!");
			}
		#endif
		}
	#ifdef RADIX_POWER64
		if(ntry == 1000000000) break;
	#endif
	}
#ifdef RADIX_POWER64
}
if(nerr == 0) {
	printf("radix_power64: >48-bit: %llu incorrectly-rounded of %llu tried, all corrected successfully.\n", nbad, ntry);
} else {
	printf("radix_power64: >48-bit: %llu incorrectly-rounded of %llu tried, %llu remain uncorrected.\n", nbad, ntry, nerr);
}
exit(0);
#endif

	// Now that have B^(3/2), do a Mont-square to get B^2 % q:
	MONT_SQR64(itmp64,q,qinv,rem64);

#if MI64_DEBUG
	if(dbg)printf("B^2 mod q = %20llu\n",rem64);
#endif

	/* rem64 holds B^2 mod q - Now compute sequence of powers needed to obtain B^len mod q via Montgomery-muls: */
	p = n;	// Add 1 since used high-MUL version of the scaled-remainder algo ( = Algorithm A in the paper)
	// If p == 2, rem64 already contains the needed power
	if(p == 3) {
		MONT_SQR64(rem64,q,qinv, rem64);
	}
	else if(p > 3)
	{
		// We always start with p = 2 and M-square that to get p = 3:
		MONT_SQR64(rem64,q,qinv,itmp64);
		bmap = 0;		// Init bitstring
		for(j = 0; p > 5; j++)		// j counts #bits processed
		{
			BIT_SETC(bmap,j,IS_EVEN(p));
			p = (p >> 1) + 1;
		}
		// Now do the needed powering. We always start with p = 2 and M-square that to get p = 3:
		MONT_SQR64(rem64,q,qinv,itmp64);	// rem64 has p = 2, itmp64 has p = 3
		if(p == 4) {
			MONT_MUL64(itmp64,rem64,q,qinv,rem64);
		} else if(p == 5) {
			MONT_SQR64(itmp64,q,qinv,rem64);
		} else {
			ASSERT(HERE, 0,"Bad starting value for power p!");
		}
		for(i = j-1; i >= 0; i--)
		{
			if(BIT_TEST(bmap,i)) {
				itmp64 = rem64;
				MONT_UNITY_MUL64(itmp64,q,qinv,itmp64);	// Reduce power of B by 1 in one of the 2 multiplicands...
				MONT_MUL64(itmp64,rem64,q,qinv,rem64);	// ...and multiply `em.
			} else {
				MONT_SQR64(rem64,q,qinv,rem64);
			}
		}
	}
#if MI64_DEBUG
	if(dbg && p > 2)printf("B^%u mod q = %20llu\n",n,rem64);
#endif
	return rem64;
}

/* Fast is-divisible-by-64-bit scalar using Montgomery modmul and right-to-left modding: */
int mi64_is_div_by_scalar64(const uint64 x[], uint64 q, uint32 len)
{
#ifndef YES_ASM
	uint64 tmp;
#endif
	uint32 i,nshift;
	uint64 qinv,cy;

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
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

#ifndef YES_ASM

	cy = (uint64)0;
	for(i = 0; i < len; ++i)
	{
		tmp  = x[i] - cy;
		/* Add q if had a borrow - Since we low=half multiply tmp by qinv below and q*qinv == 1 (mod 2^64),
		   can simply add 1 to tmp*qinv result instead of adding q to the difference here:
		*/
		cy = (cy > x[i]); /* Comparing this rather than (tmp > x[i]) frees up tmp for the multiply */
		/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
		tmp = tmp*qinv + cy;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy = __MULH64(q,tmp);
	#else
		MULH64(q,tmp, cy);
	#endif
	}

#elif 1

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Will use RDI/SIL to save CF */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
	"loop_start:					\n\t"\
		"movq	(%%r10),%%rax	\n\t"/* load x[i] */\
		"addq   $0x8, %%r10	\n\t"/* Increment array pointer */\
		"subq	%%rdx,%%rax		\n\t"/* tmp  = x[i] - cy */\
		"setc	%%dil			\n\t"/* save CF */\
		"imulq	%%rbx,%%rax 	\n\t"/* tmp *= qinv */\
		"addq	%%rdi,%%rax		\n\t"/* tmp = tmp*qinv + CF */\
		"mulq	%%rsi		\n\t"/* cy = __MULH64(q,tmp) */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop_start 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdx,%[__cy]	\n\t"\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__cy] "m" (cy)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10"	/* Clobbered registers */\
		);
#else	// This variant of my original ASM due to Robert Holmes:

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
	"loop_start:					\n\t"\
		"movq	(%%r10),%%rax	\n\t"/* load x[i] */\
		"addq   $0x8, %%r10	\n\t"/* Increment array pointer */\
		"subq	%%rdx,%%rax		\n\t"/* tmp  = x[i] - cy */\
		"sbbq	%%rdi,%%rdi		\n\t"/* save -CF */\
		"imulq	%%rbx,%%rax 	\n\t"/* tmp *= qinv */\
		"subq	%%rdi,%%rax		\n\t"/* tmp = tmp*qinv - (-CF) */\
		"mulq	%%rsi		\n\t"/* cy = __MULH64(q,tmp) */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop_start 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdx,%[__cy]	\n\t"\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__cy] "m" (cy)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10"	/* Clobbered registers */\
		);
#endif
	return (cy == 0);
}

// 4 trial divisors at a time:
int mi64_is_div_by_scalar64_x4(const uint64 x[], uint64 q0, uint64 q1, uint64 q2, uint64 q3, uint32 len)
{
	int retval = 0;
#ifndef YES_ASM
	uint64 tmp0,tmp1,tmp2,tmp3;
#endif
#if MI64_DEBUG
	int dbg = 0;
#endif
	uint32 i,trailx;
	uint32 nshift0,nshift1,nshift2,nshift3;
	uint64 qinv0,qinv1,qinv2,qinv3,cy0,cy1,cy2,cy3;

	ASSERT(HERE, (len == 0), "0 length!");
	trailx = trailz64(x[0]);
	ASSERT(HERE, trailx < 64, "0 low word!");

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift0 = trailz64(q0);
	nshift1 = trailz64(q1);
	nshift2 = trailz64(q2);
	nshift3 = trailz64(q3);

	q0 >>= nshift0;
	q1 >>= nshift1;
	q2 >>= nshift2;
	q3 >>= nshift3;
	ASSERT(HERE, q1 > 1 && q1 > 1 && q2 > 1 && q3 > 1 , "modulus must be > 1!");
	ASSERT(HERE, q0 & 1 && q1 & 1 && q2 & 1 && q3 & 1 , "even modulus!");

	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}

#ifndef YES_ASM
	cy0 = cy1 = cy2 = cy3 = (uint64)0;
	for(i = 0; i < len; ++i)
	{
		tmp0 = x[i] - cy0;			tmp1 = x[i] - cy1;			tmp2 = x[i] - cy2;			tmp3 = x[i] - cy3;
		cy0 = (cy0 > x[i]);			cy1 = (cy1 > x[i]);			cy2 = (cy2 > x[i]);			cy3 = (cy3 > x[i]);
		tmp0 = tmp0*qinv0 + cy0;	tmp1 = tmp1*qinv1 + cy1;	tmp2 = tmp2*qinv2 + cy2;	tmp3 = tmp3*qinv3 + cy3;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q0,tmp0);	cy1 = __MULH64(q1,tmp1);	cy2 = __MULH64(q2,tmp2);	cy3 = __MULH64(q3,tmp3);
	#else
		MULH64(q0,tmp0, cy0);		MULH64(q1,tmp1, cy1);		MULH64(q2,tmp2, cy2);		MULH64(q3,tmp3, cy3);
	#endif
	}

#else

	// 4-way version leaves registers RSI and RBX unused:
	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Init cy0 = 0 */\
		"xorq	%%r8 ,%%r8 		\n\t"/* Init cy1 = 0 */\
		"xorq	%%r9 ,%%r9 		\n\t"/* Init cy2 = 0 */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy3 = 0 */\

		"movslq	%[__len], %%r14	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1      , %%r14	\n\t"/* len/2 */\
		"leaq	(%%r10,%%r14,8),%%r15	\n\t"/* x+len/2 */\
		"shrq	$1      , %%r14	\n\t"/* len/4 */\
		"movq	%%r14,%%rcx	\n\t"/* Copy len/4 into loop counter*/\
	"loop4x:		\n\t"\
		"movq	(%%r10),%%rax	\n\t	leaq	(%%r10,%%r14,8),%%r11	\n\t"/* load x0,&x1 */\
		"addq	$0x8 ,%%r10		\n\t	movq	(%%r11),%%r11	\n\t"/* Increment x0-ptr, load x1 */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - cy0, save -CF */\
		"subq	%%r8 ,%%r11		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp1 = x1 - cy1, save -CF */\
		"movq	(%%r15),%%r12	\n\t	leaq	(%%r15,%%r14,8),%%r13	\n\t"/* load x2,&x3 */
		"addq	$0x8 ,%%r15		\n\t	movq	(%%r13),%%r13	\n\t"/* Increment x2-ptr, load x3 */\
		"subq	%%r9 ,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp2 = x2 - cy2, save -CF */\
		"subq	%%rdx,%%r13		\n\t	sbbq	%%rdx,%%rdx		\n\t"/* tmp3 = x3 - cy3, save -CF */\
		"imulq	%[__qinv0],%%rax \n\t	imulq	%[__qinv2],%%r11 	\n\t"/* tmp0,1 *= qinv */\
		"imulq	%[__qinv1],%%r12 \n\t	imulq	%[__qinv3],%%r13 	\n\t"/* tmp2,3 *= qinv */\
		"subq	%%rdi,%%rax		\n\t	subq	%%r8 ,%%r11		\n\t"/* tmp0,1 = tmp0,1 * qinv + CF */\
		"subq	%%r9 ,%%r12		\n\t	subq	%%rdx,%%r13		\n\t"/* tmp2,3 = tmp2,3 * qinv + CF */\
		"								mulq	%[__q0]	\n\t	movq	%%rdx,%%rdi	\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r11,%%rax		\n\t	mulq	%[__q1]	\n\t	movq	%%rdx,%%r8 	\n\t"/* cy1 = MULH64(q,tmp1), then move cy1 out of rdx in prep for cy2 computation */\
		"movq	%%r12,%%rax		\n\t	mulq	%[__q2]	\n\t	movq	%%rdx,%%r9 	\n\t"/* cy2 = MULH64(q,tmp2), then move cy2 out of rdx in prep for cy3 computation */\
		"movq	%%r13,%%rax		\n\t	mulq	%[__q3]	\n\t"/* load x3 into rax; cy3 = MULH64(q,tmp3), leave result in rdx */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop4x 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdi,%[__cy0]	\n\t	movq	%%r8 ,%[__cy1]	\n\t	movq	%%r9 ,%[__cy2]	\n\t	movq	%%rdx,%[__cy3]	\n\t"\
	:	/* outputs: none */\
	: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
	 ,[__q1] "m" (q1)	\
	 ,[__q2] "m" (q2)	\
	 ,[__q3] "m" (q3)	\
	 ,[__qinv0] "m" (qinv0)	\
	 ,[__qinv1] "m" (qinv1)	\
	 ,[__qinv2] "m" (qinv2)	\
	 ,[__qinv3] "m" (qinv3)	\
	 ,[__x] "m" (x)	\
	 ,[__cy0] "m" (cy0)	\
	 ,[__cy1] "m" (cy1)	\
	 ,[__cy2] "m" (cy2)	\
	 ,[__cy3] "m" (cy3)	\
	 ,[__len] "m" (len)	\
	: "cc","memory","rax","rcx","rdx","rdi","r8","r9","r10","r11","r12","r13","r14","r15"	/* Clobbered registers */\
	);
#endif

#if MI64_DEBUG
	if(dbg)printf("4-way carryouts: cy0-3 = %20llu, %20llu, %20llu, %20llu\n",cy0,cy1,cy2,cy3);
#endif
	retval += ((cy0 == 0) && (nshift0 <= trailx));
	retval += ((cy1 == 0) && (nshift1 <= trailx)) << 1;
	retval += ((cy2 == 0) && (nshift2 <= trailx)) << 2;
	retval += ((cy3 == 0) && (nshift3 <= trailx)) << 3;
	return retval;
}


// 2-way loop splitting:
int mi64_is_div_by_scalar64_u2(const uint64 x[], uint64 q, uint32 len)
{
#if MI64_DEBUG
	int dbg = q == 16357897499336320049ull;//87722769297534671ull;
#endif
#ifndef YES_ASM
	uint64 tmp0,tmp1,bw0,bw1;
#endif
	uint32 i,len2 = (len>>1),nshift;
	uint64 qinv,cy0,cy1,rpow;

	DBG_ASSERT(HERE, q > 0, "mi64_is_div_by_scalar64: 0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;
	ASSERT(HERE, (len&1) == 0, "odd length!");
	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz64(q);
ASSERT(HERE, !nshift, "2-way folded DIV requires odd q!");
	if(nshift)
	{
		if(trailz64(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

#ifndef YES_ASM
	cy0 = cy1 = (uint64)0;
	for(i = 0; i < len2; ++i)
	{
		tmp0  = x[i] - cy0;				tmp1 = x[i+len2] - cy1;
		cy0 = (cy0 > x[i]);				cy1 = (cy1 > x[i+len2]);
		tmp0 = tmp0*qinv + cy0;			tmp1 = tmp1*qinv + cy1;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);			cy1 = __MULH64(q,tmp1);
	#else
		MULH64(q,tmp0, cy0);				MULH64(q,tmp1, cy1);
	#endif
	}

#else

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Init cy0 = 0 */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy1 = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1, %%rcx	\n\t"/* len/2 */\
		"leaq	(%%r10,%%rcx,8),%%r11	\n\t"/* x+len/2 */\
	"loop2a:		\n\t"\
		"movq	(%%r10),%%rax	\n\t	movq	(%%r11),%%r12	\n\t"/* load x0,x1 */\
		"addq	$0x8 ,%%r10		\n\t	addq	$0x8 ,%%r11		\n\t"/* Increment array pointers */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - cy0, save -CF */\
		"subq	%%rdx,%%r12		\n\t	sbbq	%%rdx,%%rdx		\n\t"/* tmp1 = x1 - cy1, save -CF */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r12 	\n\t"/* tmp *= qinv */\
		"subq	%%rdi,%%rax		\n\t	subq	%%rdx,%%r12		\n\t"/* tmp = tmp*qinv + CF */\
		"mulq	%%rsi			\n\t	movq	%%rdx,%%rdi		\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r12,%%rax		\n\t"/* load x1 into rax */\
		"mulq	%%rsi			\n\t"/* cy1 = MULH64(q,tmp1), leave result in rdx */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop2a 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdi,%[__cy0]	\n\t	movq	%%rdx,%[__cy1]	"\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__cy0] "m" (cy0)	\
		 ,[__cy1] "m" (cy1)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12"		/* Clobbered registers */\
		);
#endif

#if MI64_DEBUG
	if(dbg)printf("Half-length carryouts: cy0 = %20llu, cy1 = %20llu\n",cy0,cy1);
#endif
	// Compute radix-power; add 1 since used high-MUL version of the scaled-remainder algo ( = Algorithm A in the paper)
	rpow = radix_power64(q,qinv,len2+1);

	// Multiply the Montgomery residue from the mod-loop by the mod-power of the base:
	MONT_MUL64(cy1,rpow,q,qinv,cy1);	// cy1*B^p (mod q)
#if MI64_DEBUG
	if(dbg) {
		printf("s1     mod q) = %20llu\n",cy0);
		printf("s2*B^p mod q) = %20llu\n",cy1);
	}
#endif
	// Sum the scaled partial remainders:
	cy0 += cy1;
	if(cy0 < cy1 || cy0 >= q) cy0 -= q;
	// Negation (mod q) needed for Algo A scaled remainder
	if(cy0) cy0 = q-cy0 ;
#if MI64_DEBUG
	if(dbg)printf("(s1 + s2*B^p) mod q = %20llu, q = %20llu\n",cy0,q);
#endif
	// One more modmul of sum by same power of the base gives true remainder - may as well, since we already have B^p handy:
	MONT_MUL64(cy0,rpow,q,qinv,cy0);
#if MI64_DEBUG
	if(dbg) {
		printf("True mod x mod q = %20llu\n",cy0);
		exit(0);
	}
#endif

	return (cy0 == 0);
}

// 4-way loop splitting:
int mi64_is_div_by_scalar64_u4(const uint64 x[], uint64 q, uint32 len)
{
#ifndef YES_ASM
	uint32 i0,i1,i2,i3;
	uint64 tmp0,tmp1,tmp2,tmp3;
#endif
#if MI64_DEBUG
	int dbg = 0;
#endif
	uint32 i,len4 = (len>>2),nshift;
	uint64 qinv,cy0,cy1,cy2,cy3,rpow;

	DBG_ASSERT(HERE, q > 0, "mi64_is_div_by_scalar64: 0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;
	ASSERT(HERE, (len&1) == 0, "odd length!");
	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz64(q);
ASSERT(HERE, !nshift, "2-way folded DIV requires odd q!");
	if(nshift)
	{
		if(trailz64(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

#ifndef YES_ASM
	cy0 = cy1 = cy2 = cy3 = (uint64)0;
	for(i0 = 0, i1 = len4, i2 = i1+i1, i3 = i1+i2; i0 < len4; ++i0, ++i1, ++i2, ++i3)
	{
		tmp0 = x[i0] - cy0;			tmp1 = x[i1] - cy1;			tmp2 = x[i2] - cy2;			tmp3 = x[i3] - cy3;
		cy0 = (cy0 > x[i0]);		cy1 = (cy1 > x[i1]);		cy2 = (cy2 > x[i2]);		cy3 = (cy3 > x[i3]);
		tmp0 = tmp0*qinv + cy0;		tmp1 = tmp1*qinv + cy1;		tmp2 = tmp2*qinv + cy2;		tmp3 = tmp3*qinv + cy3;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);		cy1 = __MULH64(q,tmp1);		cy2 = __MULH64(q,tmp2);		cy3 = __MULH64(q,tmp3);
	#else
		MULH64(q,tmp0, cy0);		MULH64(q,tmp1, cy1);		MULH64(q,tmp2, cy2);		MULH64(q,tmp3, cy3);
	#endif
	}

#else

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Init cy0 = 0 */\
		"xorq	%%r8 ,%%r8 		\n\t"/* Init cy1 = 0 */\
		"xorq	%%r9 ,%%r9 		\n\t"/* Init cy2 = 0 */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy3 = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%r14	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1      , %%r14	\n\t"/* len/2 */\
		"leaq	(%%r10,%%r14,8),%%r15	\n\t"/* x+len/2 */\
		"shrq	$1      , %%r14	\n\t"/* len/4 */\
		"movq	%%r14,%%rcx	\n\t"/* Copy len/4 into loop counter*/\
	"loop4u:		\n\t"\
		"movq	(%%r10),%%rax	\n\t	leaq	(%%r10,%%r14,8),%%r11	\n\t"/* load x0,&x1 */\
		"addq	$0x8 ,%%r10		\n\t	movq	(%%r11),%%r11	\n\t"/* Increment x0-ptr, load x1 */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - cy0, save -CF */\
		"subq	%%r8 ,%%r11		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp1 = x1 - cy1, save -CF */\
		"movq	(%%r15),%%r12	\n\t	leaq	(%%r15,%%r14,8),%%r13	\n\t"/* load x2,&x3 */
		"addq	$0x8 ,%%r15		\n\t	movq	(%%r13),%%r13	\n\t"/* Increment x2-ptr, load x3 */\
		"subq	%%r9 ,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp2 = x2 - cy2, save -CF */\
		"subq	%%rdx,%%r13		\n\t	sbbq	%%rdx,%%rdx		\n\t"/* tmp3 = x3 - cy3, save -CF */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r11 	\n\t"/* tmp0,1 *= qinv */\
		"imulq	%%rbx,%%r12 	\n\t	imulq	%%rbx,%%r13 	\n\t"/* tmp2,3 *= qinv */\
		"subq	%%rdi,%%rax		\n\t	subq	%%r8 ,%%r11		\n\t"/* tmp0,1 = tmp0,1 * qinv + CF */\
		"subq	%%r9 ,%%r12		\n\t	subq	%%rdx,%%r13		\n\t"/* tmp2,3 = tmp2,3 * qinv + CF */\
		"								mulq	%%rsi	\n\t	movq	%%rdx,%%rdi	\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r11,%%rax		\n\t	mulq	%%rsi	\n\t	movq	%%rdx,%%r8 	\n\t"/* cy1 = MULH64(q,tmp1), then move cy1 out of rdx in prep for cy2 computation */\
		"movq	%%r12,%%rax		\n\t	mulq	%%rsi	\n\t	movq	%%rdx,%%r9 	\n\t"/* cy2 = MULH64(q,tmp2), then move cy2 out of rdx in prep for cy3 computation */\
		"movq	%%r13,%%rax		\n\t	mulq	%%rsi	\n\t"/* load x3 into rax; cy3 = MULH64(q,tmp3), leave result in rdx */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop4u 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdi,%[__cy0]	\n\t	movq	%%r8 ,%[__cy1]	\n\t	movq	%%r9 ,%[__cy2]	\n\t	movq	%%rdx,%[__cy3]	\n\t"\
	:	/* outputs: none */\
	: [__q] "m" (q)	/* All inputs from memory addresses here */\
	 ,[__qinv] "m" (qinv)	\
	 ,[__x] "m" (x)	\
	 ,[__cy0] "m" (cy0)	\
	 ,[__cy1] "m" (cy1)	\
	 ,[__cy2] "m" (cy2)	\
	 ,[__cy3] "m" (cy3)	\
	 ,[__len] "m" (len)	\
	: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"	/* Clobbered registers */\
	);
#endif

#if MI64_DEBUG
	if(dbg)printf("Half-length carryouts: cy0-3 = %20llu, %20llu, %20llu, %20llu\n",cy0,cy1,cy2,cy3);
#endif
	// Compute radix-power; add 1 since used high-MUL version of the scaled-remainder algo ( = Algorithm A in the paper)
	rpow = radix_power64(q,qinv,len4+1);

	// Build up the sum of the scaled partial remainders, two terms at a time:
	MONT_MUL64(cy3,rpow,q,qinv,cy3);	cy2 += cy3;	if(cy2 < cy3 || cy2 >= q) cy2 -= q;	//               cy2 + cy3*B^p           (mod q)
	MONT_MUL64(cy2,rpow,q,qinv,cy2);	cy1 += cy2;	if(cy1 < cy2 || cy1 >= q) cy1 -= q;	//        cy1 + (cy2 + cy3*B^p)*B^p      (mod q)
	MONT_MUL64(cy1,rpow,q,qinv,cy1);	cy0 += cy1;	if(cy0 < cy1 || cy0 >= q) cy0 -= q;	// cy0 + (cy1 + (cy2 + cy3*B^p)*B^p)*B^p (mod q)
	// Negation (mod q) needed for Algo A scaled remainder
	if(cy0) cy0 = q-cy0 ;
#if MI64_DEBUG
	if(dbg) printf("(sum0-3) mod q = %20llu, q = %20llu\n",cy0,q);
#endif
	// One more modmul of sum by same power of the base gives true remainder:
	MONT_MUL64(cy0,rpow,q,qinv,cy0);
#if MI64_DEBUG
	if(dbg) {
		printf("True mod x mod q = %20llu\n",cy0);
		exit(0);
	}
#endif

	return (cy0 == 0);
}

/* Fast div-with-remainder with divisor 64-bit scalar using Montgomery modmul, my right-to-left mod
algorithm (as used in the strictly binary is-divisble-by test above), and my true-remainder postprocessing
step of May 2012.

Returns x % q in function result, and quotient x / q in y-vector, if one is provided. (Permits in-place, i.e. y == x).
*/
#if MI64_DEBUG
	#define MI64_DIV_MONT64	1
#endif
uint64 mi64_div_by_scalar64(const uint64 x[], uint64 q, uint32 len, uint64 y[])
{
#if MI64_DIV_MONT64//MI64_DEBUG
	int dbg = q == 16357897499336320049ull;
#endif
	uint32 i,nshift;
	uint64 qinv,tmp = 0,bw,cy,lo,rem64,rem_save = 0,itmp64,mask;
	double fquo,fqinv;

	ASSERT(HERE, x != 0, "Null input array!");
	DBG_ASSERT(HERE, q > 0, "0 modulus!");
	if((q == 1) || (len == 0)) return TRUE;

	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s.
	Unlike the binary-result is-divisible-by routines, we must restore the shift to the
	quotient and remainder at the end.

	Ex 1: 195 % 44 = 19, 195/44 = 4.
	After right-justifying divisor to get q = (44 >> 2) = 11, we must also right-shift 195 the same amount and save any off-shifted bits for later:
	(195 >> 2) = (192 + 3)/4 = 48, offshift = 3.
	We then compute 48/11 = 4 which is the true divisor, and re-add 3 to the resulting shift-restored remainder to get rem = (4 << 2) + 3 = 19.

	Ex 2: x=53:
	q =  5 -> x%q = 3
	q = 10 -> x%q = 3, above algo gives shifted rem = [53/2]%5 = 26%5 = 1, then restore as 2*1 + 53%2 = 2+1 = 3, ok
	q = 20 -> x%q =13, above algo gives shifted rem = [53/4]%5 = 13%5 = 3, then restore as 4*3 + 53%4 =12+1 =13, ok
	q = 40 -> x%q =13, above algo gives shifted rem = [53/8]%5 =  6%5 = 1, then restore as 8*1 + 53%8 = 8+5 =13, ok
	*/
	nshift = trailz64(q);
	if(nshift)
	{
		mask = ((uint64)-1 >> (64 - nshift));	// Save the bits which would be off-shifted if we actually right-shifted x[]
		rem_save = x[0] & mask;		// (Which we don`t do since x is read-only; thus we are forced into accounting tricks :)
		q >>= nshift;
	}
	ASSERT(HERE, (q & (uint64)1) == 1, "q must be odd!");

	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}
#if MI64_DIV_MONT64
	if(dbg) {
		printf("mi64_div_by_scalar64: q = %20llu, qinv = %20llu\n",q,qinv);
	}
#endif

	cy = (uint64)0;
	if(!nshift)
	{
		for(i = 0; i < len; ++i)
		{
			tmp  = x[i] - cy;	// Need to add q if had a borrow - Since we do MULL(tmp,qinv) below and MULL(q*qinv) = 1,
								// can simply add 1 to tmp*qinv result instead of adding q to the difference here.
			cy = (cy > x[i]);	// Comparing this rather than (tmp > x[i]) frees up tmp for the MULL-by-qinv
			/* Now do the Montgomery mod: cy = umulh( q, mulq(tmp, qinv) ); */
		#if MI64_DIV_MONT64
			bw = cy;	// Save a copy of the borrow flag for debug-printing
			itmp64 = tmp + ((-cy)&q);	// Expected value of low-half of MUL_LOHI
	//		if(dbg)printf("i = %4u, tmp*qinv = %20llu\n",i,tmp*qinv);
		#endif
			tmp = tmp*qinv + cy;
			// Do double-wide product. Fast-divisibility test needs just high half (stored in cy); low half (tmp) needed to extract true-mod
		#ifdef MUL_LOHI64_SUBROUTINE
			MUL_LOHI64(q, tmp, &tmp,&cy);
		#else
			MUL_LOHI64(q, tmp,  tmp, cy);
		#endif
		#if MI64_DIV_MONT64
			if(dbg)printf("i = %4u, lo = %20llu, hi = %20llu, bw = %1u\n",i,tmp,cy,(uint32)bw);
			ASSERT(HERE, itmp64 == tmp, "Low-half product check mismatch!");
		#endif
		}
	} else if(!y) {	// Even modulus, no quotient computation
		// Inline right-shift of x-vector with modding, store x >> nshift in y:
		for(i = 0; i < len-1; ++i)
		{
			itmp64 = (x[i] >> nshift) + (x[i+1] << (64-nshift));
			tmp  = itmp64 - cy;
			cy = (cy > itmp64);
		#if MI64_DIV_MONT64
			bw = cy;	// Save a copy of the borrow flag for debug-printing
			itmp64 = tmp + ((-cy)&q);	// Expected value of low-half of MUL_LOHI
		#endif
			tmp = tmp*qinv + cy;
		#ifdef MUL_LOHI64_SUBROUTINE
			MUL_LOHI64(q, tmp, &tmp,&cy);
		#else
			MUL_LOHI64(q, tmp,  tmp, cy);
		#endif
		#if MI64_DIV_MONT64
			if(dbg)printf("i = %4u, lo = %20llu, hi = %20llu, bw = %1u\n",i,tmp,cy,(uint32)bw);
			ASSERT(HERE, itmp64 == tmp, "Low-half product check mismatch!");
		#endif
		}
		// Last element has no shift-in from next-higher term, so can compute just the low-half output term, sans explicit MULs:
		itmp64 = (x[i] >> nshift);
		tmp  = itmp64 - cy;
		cy = (cy > itmp64);
		tmp = tmp + ((-cy)&q);
	#if MI64_DIV_MONT64
		if(dbg)printf("i = %4u, lo_out = %20llu\n",i,tmp);
	#endif
	} else {	// Even modulus, with quotient computation
		for(i = 0; i < len-1; ++i)
		{
			y[i] = (x[i] >> nshift) + (x[i+1] << (64-nshift));
			tmp  = y[i] - cy;
			cy = (cy > y[i]);
		#if MI64_DIV_MONT64
			bw = cy;	// Save a copy of the borrow flag for debug-printing
			itmp64 = tmp + ((-cy)&q);	// Expected value of low-half of MUL_LOHI
		#endif
			tmp = tmp*qinv + cy;
		#ifdef MUL_LOHI64_SUBROUTINE
			MUL_LOHI64(q, tmp, &tmp,&cy);
		#else
			MUL_LOHI64(q, tmp,  tmp, cy);
		#endif
		#if MI64_DIV_MONT64
			if(dbg)printf("i = %4u, lo = %20llu, hi = %20llu, bw = %1u\n",i,tmp,cy,(uint32)bw);
			ASSERT(HERE, itmp64 == tmp, "Low-half product check mismatch!");
		#endif
		}
		y[i] = (x[i] >> nshift);
		tmp  = y[i] - cy;
		cy = (cy > y[i]);
		tmp = tmp + ((-cy)&q);
	#if MI64_DIV_MONT64
		if(dbg)printf("i = %4u, lo_out = %20llu, bw = %1u\n",i,tmp,(uint32)cy);
	#endif
	}

	// Prepare to transform back out of "Montgomery space" ... first compute B^2 mod q using successive FP approximation:
	// Compute how many 2^48s there are in the quotient 2^96/q:
	fqinv = 1.0/q;

	// If len = 1, simply return x[0] % q = tmp % q:
	if(len == 1) {
		fquo = tmp*fqinv;
		itmp64 = (uint64)fquo;
		rem64 = tmp - q*itmp64;
		// May need a 2nd pass to clean up any ROE in 1st iteration, and
		// must account for ROE which leads to a borrow in the above subtraction to get rem64:
		if(rem64 > tmp) {	// Had a borrow
			fquo = -rem64*fqinv + 1;	// Add one to FP quotient to effectround-toward-zero
			itmp64 -= (uint64)fquo;
			rem64 = rem64 + q*(uint64)fquo;
		} else {
			fquo = rem64*fqinv;
			itmp64 += (uint64)fquo;
			rem64 = rem64 - q*(uint64)fquo;
		}
		DBG_ASSERT(HERE, rem64 == tmp%q, "Bad floating-mod mod!");
		if(y) {
			y[0] = itmp64;
		}
		return rem64;
	}

	// Compute radix-power; no add-1 here since use scaled-remainder Algorithm B:
	rem64 = radix_power64(q,qinv,len);

	// And multiply the Montgomery residue from the mod-loop by the mod-power of the base:
	MONT_MUL64(tmp,rem64,q,qinv,rem64);

#if MI64_DIV_MONT64
	if(dbg)printf("True mod x mod q = %20llu\n",rem64);
#endif

	// If we applied an initial right-justify shift to the modulus, restore the shift to the
	// current (partial) remainder and re-add the off-shifted part of the true remainder.
	if(nshift)
	{
		rem64 = (rem64 << nshift) + rem_save;
	}

	// If y-array supplied, compute quotient and return it in that:
	if(y) {
	#if MI64_DIV_MONT64
		if(dbg)printf("Computing quotient...\n");
	#endif
		// Now can use a simple loop and a sequence of word-size MULLs to obtain quotient.
		// Fusing the functionality of mi64_sub_scalar and the quotient extraction is fastest here:
		if(!nshift) {
			// If odd modulus, have not yet copied input array to y...
			bw = 0;	cy = rem64;
			for(i = 0; i < len; ++i)
			{
#if MI64_DEBUG
	if(dbg && i%(len>>2) == 0)printf("mi64_div_by_scalar64: bw = %1llu, cy%1u = %20llu\n",bw,i/(len>>2),cy);	// Use to debug loop-folded implemntation
#endif
				tmp = x[i] - bw - cy;
				/*  Since may be working in-place, need an extra temp here due to asymmetry of subtract: */
				bw = (tmp > x[i]);
			#if MI64_DIV_MONT64
				itmp64 = tmp;	// Expected value of low-half of MUL_LOHI; here the borrow gets subtracted from the next-higher word so no mod-q
			#endif
				tmp *= qinv;
			#ifdef MUL_LOHI64_SUBROUTINE
				MUL_LOHI64(q, tmp, &lo,&cy);
			#else
				MUL_LOHI64(q, tmp,  lo, cy);
			#endif
			#if MI64_DIV_MONT64
				if(dbg)printf("i = %4u, quot[i] = %20llu, lo1 = %20llu, lo2 = %20llu, hi = %20llu, bw = %1u\n",i,tmp,itmp64,lo,cy,(uint32)bw);
				ASSERT(HERE, itmp64 == lo, "Low-half product check mismatch!");
			#endif
				y[i] = tmp;
			}
		} else {
			// If even modulus, right-justified copy of input array already in y.
			bw = 0;	cy = rem64>>nshift;
			for(i = 0; i < len; ++i)
			{
				tmp = y[i] - bw - cy;
				/*  Since may be working in-place, need an extra temp here due to asymmetry of subtract: */
				bw = (tmp > y[i]);
				tmp *= qinv;
			#ifdef MUL_LOHI64_SUBROUTINE
				MUL_LOHI64(q, tmp, &lo,&cy);
			#else
				MUL_LOHI64(q, tmp,  lo, cy);
			#endif
			#if MI64_DIV_MONT64
				if(dbg)printf("i = %4u, quot[i] = %20llu\n",i,tmp);
			#endif
				y[i] = tmp;
			}
		}
		ASSERT(HERE, bw == 0 && cy == 0, "bw/cy check!");
	}
#if MI64_DIV_MONT64
	if(dbg) {
		return rem64;
	}
#endif
	return rem64;
}

// 2-way loop splitting:
uint64 mi64_div_by_scalar64_u2(const uint64 x[], uint64 q, uint32 len, uint64 y[])
{
#ifndef YES_ASM
	uint64 tmp0,tmp1,bw0,bw1;
#endif
#if MI64_DEBUG
	int dbg = q == 16357897499336320049ull;//87722769297534671ull;
#endif
	int i,len2 = (len>>1),nshift;
	uint64 qinv,cy0,cy1,rpow;

	DBG_ASSERT(HERE, q > 0, "0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;
	ASSERT(HERE, (len&1) == 0, "odd length!");
	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz64(q);
ASSERT(HERE, !nshift, "2-way folded DIV requires odd q!");
	if(nshift)
	{
		if(trailz64(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

/* 07/30/2012: Compare full-DIV timings for the standard test, 10000 length-34420 DIvs:
Pure C:  3.585 sec --> 21 cycles
GCC ASM: 2.953 sec --> 17 cycles.
Bizarrely, switch just ONE of the 2 steps (either one) from C to ASM gives a negligible speedup; only doing both helps a lot.
[Possible code/loop alignment issue?]
See similar behavior for 4-way-split version of the algorithm.
*/
#ifndef YES_ASM
	cy0 = cy1 = (uint64)0;
	for(i = 0; i < len2; ++i)
	{
		tmp0  = x[i] - cy0;				tmp1 = x[i+len2] - cy1;
		cy0 = (cy0 > x[i]);				cy1 = (cy1 > x[i+len2]);
		tmp0 = tmp0*qinv + cy0;			tmp1 = tmp1*qinv + cy1;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);			cy1 = __MULH64(q,tmp1);
	#else
		MULH64(q,tmp0, cy0);				MULH64(q,tmp1, cy1);
	#endif
	}

#else

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Init cy0 = 0 */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy1 = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1, %%rcx	\n\t"/* len/2 */\
		"leaq	(%%r10,%%rcx,8),%%r11	\n\t"/* x+len/2 */\
	"loop2b:		\n\t"\
		"movq	(%%r10),%%rax	\n\t	movq	(%%r11),%%r12	\n\t"/* load x0,x1 */\
		"addq	$0x8 ,%%r10		\n\t	addq	$0x8 ,%%r11		\n\t"/* Increment array pointers */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - cy0, save -CF */\
		"subq	%%rdx,%%r12		\n\t	sbbq	%%rdx,%%rdx		\n\t"/* tmp1 = x1 - cy1, save -CF */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r12 	\n\t"/* tmp *= qinv */\
		"subq	%%rdi,%%rax		\n\t	subq	%%rdx,%%r12		\n\t"/* tmp = tmp*qinv + CF */\
		"mulq	%%rsi			\n\t	movq	%%rdx,%%rdi		\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r12,%%rax		\n\t"/* load x1 into rax */\
		"mulq	%%rsi			\n\t"/* cy1 = MULH64(q,tmp1), leave result in rdx */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop2b 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdi,%[__cy0]	\n\t	movq	%%rdx,%[__cy1]	"\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__cy0] "m" (cy0)	\
		 ,[__cy1] "m" (cy1)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","r11","r12"		/* Clobbered registers */\
		);
#endif

#if MI64_DEBUG
	if(dbg)printf("Half-length carryouts: cy0 = %20llu, cy1 = %20llu\n",cy0,cy1);
#endif
	// Compute radix-power; add-1 here since use scaled-remainder Algorithm A:
	rpow = radix_power64(q,qinv,len2+1);

	/*** Here is loop-split quotient computation: ***/

	// Negation (mod q) needed for Algo A true remainder computation:
	if(cy0) cy0 = q-cy0;
	if(cy1) cy1 = q-cy1;
	MONT_MUL64(cy1,rpow,q,qinv,cy1);	// cy1*B^p (mod q)
	// Sum the scaled partial remainders:
	cy0 += cy1;
	if(cy0 < cy1 || cy0 >= q) cy0 -= q;
	// One more modmul of sum by same power of the base gives full remainder in cy0, where we need it:
	MONT_MUL64(cy0,rpow,q,qinv,cy0);
#if MI64_DEBUG
	if(dbg) {
		printf("cy0 = %20llu\n",cy0);
		printf("cy1 = %20llu\n",cy1);
	}
#endif
	// Since re-use cy0-3 in quotient loop, save copy of full-vector remainder (cy0) to use for return value:
	rpow = cy0;

#ifndef YES_ASM
	bw0 = bw1 = (uint64)0;
	for(i = 0; i < len2; ++i)
	{
		tmp0  = x[i] - bw0 - cy0;		tmp1 = x[i+len2] - bw1 - cy1;
		/*  Since may be working in-place, need an extra temp here due to asymmetry of subtract: */
		bw0 = (tmp0 > x[i]);			bw1 = (tmp1 > x[i+len2]);
		tmp0 *= qinv;					tmp1 *= qinv;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);			cy1 = __MULH64(q,tmp1);
	#else
		MULH64(q,tmp0, cy0);			MULH64(q,tmp1, cy1);
	#endif
	#if MI64_DEBUG
		if(dbg)printf("quot[%2u] = %20llu, quot[%2u] = %20llu, bw0,1 = %1u,%1u, cy0,1 = %20llu,%20llu\n",i,tmp0,i+len2,tmp1,(uint32)bw0,(uint32)bw1,cy0,cy1);
	#endif
		// Write quotient word(s):
		y[i] = tmp0;					y[i+len2] = tmp1;
	}

#elif 0

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"movq	%[__y],%%r13	\n\t"\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\
		"xorq	%%r8 ,%%r8 		\n\t	xorq	%%r9 ,%%r9 		\n\t"/* Init bw0,1 = 0 */\
		"movq	%[__cy0],%%rdi	\n\t	movq	%[__cy1],%%rdx	\n\t"/* Load cy0, cy1 */\
		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1, %%rcx	\n\t"/* len/2 */\
		"leaq	(%%r10,%%rcx,8),%%r11	\n\t"/* x+len/2 */\
		"leaq	(%%r13,%%rcx,8),%%r14	\n\t"/* y+len/2 */\
	"loop2c:		\n\t"
		"subq	%%r8 ,%%rdi		\n\t	subq	%%r9 ,%%rdx		\n\t"/* [rdi,r12] = (cy + bw)[0,1] via cy -(-bw), no carryout */\
		"movq	(%%r10),%%rax	\n\t	movq	(%%r11),%%r12	\n\t"/* load x0,x1 */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp0 = x0 - (bw0 + cy0), save -CF */\
		"subq	%%rdx,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp1 = x1 - (bw1 + cy1), save -CF */\
		"addq	$0x8 ,%%r10		\n\t	addq	$0x8 ,%%r11		\n\t"/* Increment x-array pointers */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r12 	\n\t"/* tmp *= qinv */\
		"movq	%%rax,(%%r13)	\n\t	movq	%%r12,(%%r14)	\n\t"/* store tmp0,1 */\
		"addq	$0x8 ,%%r13		\n\t	addq	$0x8 ,%%r14		\n\t"/* Increment y-array pointers */\
		"mulq	%%rsi			\n\t	movq	%%rdx,%%rdi		\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r12,%%rax		\n\t	mulq	%%rsi			\n\t"/* load tmp1 into rax, then cy1 = MULH64(q,tmp1) */\
	"subq	$1,%%rcx \n\t"\
	"jnz loop2c 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \
		"movq	%%rdi,%[__cy0]	\n\t	movq	%%rdx,%[__cy1]	"\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__y] "m" (y)	\
		 ,[__cy0] "m" (cy0)	\
		 ,[__cy1] "m" (cy1)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14"		/* Clobbered registers */\
		);

#else	// Low-register version; use similar for 4-folded quotient loop:

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"movq	%[__y],%%r13	\n\t"\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\
		"xorq	%%r8 ,%%r8 		\n\t	xorq	%%r9 ,%%r9 		\n\t"/* Init bw0,1 = 0 */\
		"subq	%[__cy0],%%r8	\n\t	movq	%[__cy1],%%rdx	\n\t"/* Load -(bw0 + cy0), +(bw1 + cy1) */\
		"movslq	%[__len], %%rcx	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1, %%rcx	\n\t"/* len/2 */\
		"leaq	(%%r10,%%rcx,8),%%r11	\n\t"/* x+len/2 */\
		"leaq	(%%r13,%%rcx,8),%%r14	\n\t"/* y+len/2 */\
	"loop2c:		\n\t"
		"negq	%%r8			\n\t	subq	%%r9 ,%%rdx		\n\t"/* r8 = +(bw0 + cy0), r12 = +(bw1 + cy1) via cy1 -(-bw1), no carryout possible */\
		"movq	(%%r10),%%rax	\n\t	movq	(%%r11),%%r12	\n\t"/* load x0,x1 */\
		"subq	%%r8 ,%%rax		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp0 = x0 - (bw0 + cy0), save -CF */\
		"subq	%%rdx,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp1 = x1 - (bw1 + cy1), save -CF */\
		"addq	$0x8 ,%%r10		\n\t	addq	$0x8 ,%%r11		\n\t"/* Increment x-array pointers */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r12 	\n\t"/* tmp *= qinv */\
		"movq	%%rax,(%%r13)	\n\t	movq	%%r12,(%%r14)	\n\t"/* store tmp0,1 */\
		"addq	$0x8 ,%%r13		\n\t	addq	$0x8 ,%%r14		\n\t"/* Increment y-array pointers */\
		"mulq	%%rsi			\n\t	subq	%%rdx,%%r8		\n\t"/* cy0 = MULH64(q,tmp0), move cy0 out of rdx into -(bw0 + cy0) in prep for cy1 computation */\
		"movq	%%r12,%%rax		\n\t	mulq	%%rsi			\n\t"/* load tmp1 into rax, then cy1 = MULH64(q,tmp1) */\
	"subq	$1,%%rcx \n\t"\
	"jnz loop2c 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \
		"						\n\t	movq	%%rdx,%[__cy1]	"/* Only useful carryout is cy1, check-equal-to-zero */\
		:	/* outputs: none */\
		: [__q] "m" (q)	/* All inputs from memory addresses here */\
		 ,[__qinv] "m" (qinv)	\
		 ,[__x] "m" (x)	\
		 ,[__y] "m" (y)	\
		 ,[__cy0] "m" (cy0)	\
		 ,[__cy1] "m" (cy1)	\
		 ,[__len] "m" (len)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","r9","r10","r11","r12","r13","r14"		/* Clobbered registers */\
		);
#endif
//	ASSERT(HERE, bw0 == 0 && bw1 == 0, "bw check!");	** in fact do not require this **
	ASSERT(HERE, cy1 == 0, "cy check!");	// all but the uppermost carryout are generally nonzero
	return rpow;
}

// 4-way loop splitting:
uint64 mi64_div_by_scalar64_u4(const uint64 x[], uint64 q, uint32 len, uint64 y[])
{
#ifndef YES_ASM
	uint32 i0,i1,i2,i3;
	uint64 tmp0,tmp1,tmp2,tmp3,bw0,bw1,bw2,bw3;
#endif
#if MI64_DEBUG
	int dbg = q == 16357897499336320049ull;
#endif
	int i,len4 = (len>>2),nshift;
	uint64 qinv,cy0,cy1,cy2,cy3,rpow;
	uint64*xy_ptr_diff = (uint64*)((uint64)y - (uint64)x);	// 2-step cast to avoid GCC "initialization makes pointer from integer without a cast" warning
						// The inner (uint64) casts are needed for the compiler to emit correctly functioning code.
	DBG_ASSERT(HERE, q > 0, "mi64_is_div_by_scalar64: 0 modulus!");
	if(q == 1) return TRUE;
	if(len == 0) return TRUE;
	ASSERT(HERE, (len&1) == 0, "odd length!");
	/* q must be odd for Montgomery-style modmul to work, so first shift off any low 0s: */
	nshift = trailz64(q);
ASSERT(HERE, !nshift, "2-way folded DIV requires odd q!");
	if(nshift)
	{
		if(trailz64(x[0]) < nshift) return FALSE;
		q >>= nshift;
	}

	qinv = (q+q+q) ^ (uint64)2;
	for(i = 0; i < 4; i++)
	{
		qinv = qinv*((uint64)2 - q*qinv);
	}

/* 07/30/2012: Compare full-DIV timings for the standard test, 10000 length-34420 DIvs:
Pure C       : 2.651 sec --> 15.4 cycles
Rem-ASM/Div-C: 2.763 sec --> 16.0 cycles	<*** slower! Bizarre... ***
Both ASM     : 2.166 sec --> 12.5 cycles.	<*** Ignore incorrct output of -bw-cy version of loop for now ***
*/
#ifndef YES_ASM
	cy0 = cy1 = cy2 = cy3 = (uint64)0;
	for(i0 = 0, i1 = len4, i2 = i1+i1, i3 = i1+i2; i0 < len4; ++i0, ++i1, ++i2, ++i3)
	{
		tmp0 = x[i0] - cy0;			tmp1 = x[i1] - cy1;			tmp2 = x[i2] - cy2;			tmp3 = x[i3] - cy3;
		cy0 = (cy0 > x[i0]);		cy1 = (cy1 > x[i1]);		cy2 = (cy2 > x[i2]);		cy3 = (cy3 > x[i3]);
		tmp0 = tmp0*qinv + cy0;		tmp1 = tmp1*qinv + cy1;		tmp2 = tmp2*qinv + cy2;		tmp3 = tmp3*qinv + cy3;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);		cy1 = __MULH64(q,tmp1);		cy2 = __MULH64(q,tmp2);		cy3 = __MULH64(q,tmp3);
	#else
		MULH64(q,tmp0, cy0);		MULH64(q,tmp1, cy1);		MULH64(q,tmp2, cy2);		MULH64(q,tmp3, cy3);
	#endif
	}

#else

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t"/* Init cy0 = 0 */\
		"xorq	%%r8 ,%%r8 		\n\t"/* Init cy1 = 0 */\
		"xorq	%%r9 ,%%r9 		\n\t"/* Init cy2 = 0 */\
		"xorq	%%rdx,%%rdx		\n\t"/* Init cy3 = 0 */\
		"movq	%[__q],%%rsi	\n\t"\
		"movq	%[__qinv],%%rbx	\n\t"\

		"movslq	%[__len], %%r14	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1      , %%r14	\n\t"/* len/2 */\
		"leaq	(%%r10,%%r14,8),%%r15	\n\t"/* x+len/2 */\
		"shrq	$1      , %%r14	\n\t"/* len/4 */\
		"movq	%%r14,%%rcx	\n\t"/* Copy len/4 into loop counter*/\
	"loop4b:		\n\t"\
		"movq	(%%r10),%%rax	\n\t	leaq	(%%r10,%%r14,8),%%r11	\n\t"/* load x0,&x1 */\
		"addq	$0x8 ,%%r10		\n\t	movq	(%%r11),%%r11	\n\t"/* Increment x0-ptr, load x1 */\
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - cy0, save -CF */\
		"subq	%%r8 ,%%r11		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp1 = x1 - cy1, save -CF */\
		"movq	(%%r15),%%r12	\n\t	leaq	(%%r15,%%r14,8),%%r13	\n\t"/* load x2,&x3 */
		"addq	$0x8 ,%%r15		\n\t	movq	(%%r13),%%r13	\n\t"/* Increment x2-ptr, load x3 */\
		"subq	%%r9 ,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp2 = x2 - cy2, save -CF */\
		"subq	%%rdx,%%r13		\n\t	sbbq	%%rdx,%%rdx		\n\t"/* tmp3 = x3 - cy3, save -CF */\
		"imulq	%%rbx,%%rax 	\n\t	imulq	%%rbx,%%r11 	\n\t"/* tmp0,1 *= qinv */\
		"imulq	%%rbx,%%r12 	\n\t	imulq	%%rbx,%%r13 	\n\t"/* tmp2,3 *= qinv */\
		"subq	%%rdi,%%rax		\n\t	subq	%%r8 ,%%r11		\n\t"/* tmp0,1 = tmp0,1 * qinv + CF */\
		"subq	%%r9 ,%%r12		\n\t	subq	%%rdx,%%r13		\n\t"/* tmp2,3 = tmp2,3 * qinv + CF */\
		"								mulq	%%rsi	\n\t	movq	%%rdx,%%rdi	\n\t"/* cy0 = MULH64(q,tmp0), then move cy0 out of rdx in prep for cy1 computation */\
		"movq	%%r11,%%rax		\n\t	mulq	%%rsi	\n\t	movq	%%rdx,%%r8 	\n\t"/* cy1 = MULH64(q,tmp1), then move cy1 out of rdx in prep for cy2 computation */\
		"movq	%%r12,%%rax		\n\t	mulq	%%rsi	\n\t	movq	%%rdx,%%r9 	\n\t"/* cy2 = MULH64(q,tmp2), then move cy2 out of rdx in prep for cy3 computation */\
		"movq	%%r13,%%rax		\n\t	mulq	%%rsi	\n\t"/* load x3 into rax; cy3 = MULH64(q,tmp3), leave result in rdx */\

	"subq	$1,%%rcx \n\t"\
	"jnz loop4b 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \

		"movq	%%rdi,%[__cy0]	\n\t	movq	%%r8 ,%[__cy1]	\n\t	movq	%%r9 ,%[__cy2]	\n\t	movq	%%rdx,%[__cy3]	\n\t"\
	:	/* outputs: none */\
	: [__q] "m" (q)	/* All inputs from memory addresses here */\
	 ,[__qinv] "m" (qinv)	\
	 ,[__x] "m" (x)	\
	 ,[__cy0] "m" (cy0)	\
	 ,[__cy1] "m" (cy1)	\
	 ,[__cy2] "m" (cy2)	\
	 ,[__cy3] "m" (cy3)	\
	 ,[__len] "m" (len)	\
	: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"	/* Clobbered registers */\
	);

#endif

#if MI64_DEBUG
	if(dbg)printf("Half-length carryouts: cy0-3 = %20llu, %20llu, %20llu, %20llu\n",cy0,cy1,cy2,cy3);
#endif
	// Compute radix-power; add-1 here since use scaled-remainder Algorithm A:
	rpow = radix_power64(q,qinv,len4+1);

	/*** Here is loop-split quotient computation: ***/

	// Negation (mod q) needed for Algo A true remainder computation:
	if(cy3) cy3 = q-cy3;
	if(cy2) cy2 = q-cy2;
	if(cy1) cy1 = q-cy1;
	if(cy0) cy0 = q-cy0;
	// Build up the sum of the scaled partial remainders, two terms at a time:
	MONT_MUL64(cy3,rpow,q,qinv,cy3);	cy2 += cy3;	if(cy2 < cy3 || cy2 >= q) cy2 -= q;	//               cy2 + cy3*B^p           (mod q)
	MONT_MUL64(cy2,rpow,q,qinv,cy2);	cy1 += cy2;	if(cy1 < cy2 || cy1 >= q) cy1 -= q;	//        cy1 + (cy2 + cy3*B^p)*B^p      (mod q)
	MONT_MUL64(cy1,rpow,q,qinv,cy1);	cy0 += cy1;	if(cy0 < cy1 || cy0 >= q) cy0 -= q;	// cy0 + (cy1 + (cy2 + cy3*B^p)*B^p)*B^p (mod q)
	// One more modmul of sum by same power of the base gives full remainder in cy0, where we need it:
	MONT_MUL64(cy0,rpow,q,qinv,cy0);
#if MI64_DEBUG
	if(dbg) {
		printf("cy0 = %20llu\n",cy0);
		printf("cy1 = %20llu\n",cy1);
		printf("cy2 = %20llu\n",cy2);
		printf("cy3 = %20llu\n",cy3);
	}
#endif
	// Since re-use cy0-3 in quotient loop, save copy of full-vector remainder (cy0) to use for return value:
	rpow = cy0;

#ifndef YES_ASM

	bw0 = bw1 = bw2 = bw3 = (uint64)0;
	for(i0 = 0, i1 = len4, i2 = i1+i1, i3 = i1+i2; i0 < len4; ++i0, ++i1, ++i2, ++i3)
	{
		tmp0 = x[i0] - bw0 - cy0;	tmp1 = x[i1] - bw1 - cy1;	tmp2 = x[i2] - bw2 - cy2;	tmp3 = x[i3] - bw3 - cy3;
		/*  Since may be working in-place, need an extra temp here due to asymmetry of subtract: */
		bw0 = (tmp0 > x[i0]);		bw1 = (tmp1 > x[i1]);		bw2 = (tmp2 > x[i2]);		bw3 = (tmp3 > x[i3]);
		tmp0 = tmp0*qinv;			tmp1 = tmp1*qinv;			tmp2 = tmp2*qinv;			tmp3 = tmp3*qinv;
	#ifdef MUL_LOHI64_SUBROUTINE
		cy0 = __MULH64(q,tmp0);		cy1 = __MULH64(q,tmp1);		cy2 = __MULH64(q,tmp2);		cy3 = __MULH64(q,tmp3);
	#else
		MULH64(q,tmp0, cy0);		MULH64(q,tmp1, cy1);		MULH64(q,tmp2, cy2);		MULH64(q,tmp3, cy3);
	#endif
	#if MI64_DEBUG
		if(dbg)printf("quot[%2u,%2u,%2u,%2u] = %20llu,%20llu,%20llu,%20llu, bw0-3 = %1u,%1u,%1u,%1u, cy0-3 = %20llu,%20llu,%20llu,%20llu\n",i0,i1,i2,i3,tmp0,tmp1,tmp2,tmp3,(uint32)bw0,(uint32)bw1,(uint32)bw2,(uint32)bw3,cy0,cy1,cy2,cy3);
	#endif
		// Write quotient words:
		y[i0] = tmp0;				y[i1] = tmp1;				y[i2] = tmp2;				y[i3] = tmp3;
	}

#else

	__asm__ volatile (\
		"movq	%[__x],%%r10	\n\t"\
		"xorq	%%rdi,%%rdi		\n\t	subq	%[__cy0],%%rdi	\n\t"/* Init bw0 = 0, -(bw0 + cy0) */\
		"xorq	%%r8 ,%%r8 		\n\t	subq	%[__cy1],%%r8 	\n\t"/* Init bw1 = 0, -(bw1 + cy1) */\
		"xorq	%%r9 ,%%r9 		\n\t	subq	%[__cy2],%%r9 	\n\t"/* Init bw2 = 0, -(bw2 + cy2) */\
		"xorq	%%rbx,%%rbx		\n\t	movq	%[__cy3],%%rdx	\n\t"/* Init bw3 = 0, +(bw3 + cy3) */\
	/* Quotient loop needs a register for output-array ptr; timings of rem-loop show replacing %%rsi with %[__q] in MULQs is timing-neutral: */\
	/*	"movq	%[__q],%%rsi	\n\t"\	*/\
	/*	"movq	%[__qinv],%%rbx	\n\t"\	... similar timing-neutrality holds for [__q] */\
		"movq	%[__xy_ptr_diff],%%rsi	\n\t"/* Load (y-x) pointer-diff */\

		"movslq	%[__len],%%r14	\n\t"/* ASM loop structured as for(j = len; j != 0; --j){...} */\
		"shrq	$1      ,%%r14	\n\t"/* len/2 */\
		"leaq	(%%r10,%%r14,8),%%r15	\n\t"/* x+len/2 */\
		"shrq	$1      ,%%r14	\n\t"/* len/4 */\
		"movq	%%r14,%%rcx	\n\t"/* Copy len/4 into loop counter*/\
		"shlq	$3,%%r14	\n\t"/* r14 holds pointer offset corr. to len/4 array elements */\
	"loop4c:		\n\t"
		"leaq	(%%r10,%%r14),%%r11	\n\t	movq	(%%r10),%%rax	\n\t	movq	(%%r11),%%r11	\n\t"/* Load &x1,*x0,*x1 */\
		"leaq	(%%r15,%%r14),%%r13	\n\t	movq	(%%r15),%%r12	\n\t	movq	(%%r13),%%r13	\n\t"/* Load &x3,*x2,*x3 */
		"negq	%%rdi			\n\t	negq	%%r8			\n\t"/* rdi = +(bw0 + cy0), r8  = +(bw1 + cy1) */\
		"negq	%%r9			\n\t	subq	%%rbx,%%rdx		\n\t"/* r9  = +(bw2 + cy2), rbx = +(bw3 + cy3) via cy3 -(-bw3), no carryout possible */
		"subq	%%rdi,%%rax		\n\t	sbbq	%%rdi,%%rdi		\n\t"/* tmp0 = x0 - (bw0 + cy0), save -CF */\
		"subq	%%r8 ,%%r11		\n\t	sbbq	%%r8 ,%%r8 		\n\t"/* tmp1 = x1 - (bw1 + cy1), save -CF */\
		"subq	%%r9 ,%%r12		\n\t	sbbq	%%r9 ,%%r9 		\n\t"/* tmp2 = x2 - (bw2 + cy2), save -CF */\
		"subq	%%rdx,%%r13		\n\t	sbbq	%%rbx,%%rbx		\n\t"/* tmp3 = x3 - (bw3 + cy3), save -CF */\
		"imulq	%[__qinv],%%rax \n\t	imulq	%[__qinv],%%r11 \n\t"/* tmp0,1 *= qinv */\
		"imulq	%[__qinv],%%r12 \n\t	imulq	%[__qinv],%%r13 \n\t"/* tmp2,3 *= qinv */\
		"addq	%%rsi,%%r10		\n\t	addq	%%rsi,%%r15		\n\t"/* Add (y-x) pointer-diff to x0,2 ptrs to get y0,2 */\
		"movq	%%rax,(%%r10)	\n\t	movq	%%r12,(%%r15)	\n\t"/* store tmp0,2 */\
		"addq	%%r14,%%r10		\n\t	addq	%%r14,%%r15		\n\t"/* Add len/4 pointer-diff to y0,2 ptrs to get y1,3 */\
		"movq	%%r11,(%%r10)	\n\t	movq	%%r13,(%%r15)	\n\t"/* store tmp1,3 */\
		"subq	%%r14,%%r10		\n\t	subq	%%r14,%%r15		\n\t"/* Sub len/4 pointer-diff from y1,3 ptrs to get y0,2 */\
		"subq	%%rsi,%%r10		\n\t	subq	%%rsi,%%r15		\n\t"/* Sub (y-x) pointer-diff from y0,2 ptrs to get x0,2 */\
		"								mulq	%[__q]	\n\t	subq	%%rdx,%%rdi	\n\t"/* cy0 = MULH64(q,tmp0), move cy0 from rdx -> -bw0-cy0 */\
		"movq	%%r11,%%rax		\n\t	mulq	%[__q]	\n\t	subq	%%rdx,%%r8 	\n\t"/* cy1 = MULH64(q,tmp1), move cy1 from rdx -> -bw1-cy1 */\
		"movq	%%r12,%%rax		\n\t	mulq	%[__q]	\n\t	subq	%%rdx,%%r9 	\n\t"/* cy2 = MULH64(q,tmp2), move cy2 from rdx -> -bw2-cy2 */\
		"movq	%%r13,%%rax		\n\t	mulq	%[__q]	\n\t"/* load x3 into rax; cy3 = MULH64(q,tmp3), leave result in rdx */\
		"addq	$0x8 ,%%r10		\n\t	addq	$0x8 ,%%r15		\n\t"/* Increment x0,2-pointers */\
	"subq	$1,%%rcx \n\t"\
	"jnz loop4c 	\n\t"/* loop1 end; continue is via jump-back if rcx != 0 */ \
		"movq	%%rdx,%[__cy3]	\n\t"\

	:	/* outputs: none */\
	: [__q] "m" (q)	/* All inputs from memory addresses here */\
	 ,[__qinv] "m" (qinv)	\
	 ,[__x] "m" (x)	\
	 ,[__xy_ptr_diff] "m" (xy_ptr_diff)	\
	 ,[__cy0] "m" (cy0)	\
	 ,[__cy1] "m" (cy1)	\
	 ,[__cy2] "m" (cy2)	\
	 ,[__cy3] "m" (cy3)	\
	 ,[__len] "m" (len)	\
	: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"	/* Clobbered registers */\
	);
#endif
//	ASSERT(HERE, bw0 == 0 && bw1 == 0 && bw2 == 0 && bw3 == 0, "bw check!");
	ASSERT(HERE, cy3 == 0, "cy check!");	// all but the uppermost carryout are generally nonzero
	return rpow;
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
	uint32 n_alloc_chars = STR_MAX_LEN;
	return __convert_mi64_base10_char(char_buf, n_alloc_chars, x, len);
}

/* This is an arbitrary-string-length core routine which can be called directly, requires caller to supply allocated-length of input string: */
int	__convert_mi64_base10_char(char char_buf[], uint32 n_alloc_chars, const uint64 x[], uint32 len)
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
	ASSERT(HERE, MAX_DIGITS < n_alloc_chars, "Output string overflows buffer");
	char_buf[MAX_DIGITS-1]= '0';	// Init least-significant digit = 0, in case input = 0
	char_buf[MAX_DIGITS  ]='\0';

	/* Write the decimal digits into the string from right to left.
	This avoids the need to reverse the digits after calculating them.
	*/
	for(i = 0; i < MAX_DIGITS; i++)
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
	for(i = 0; i < imax; i++)
	{
		c = char_buf[i];
		if(!isspace(c)) {
			break;
		}
	}
	/* Estimate # of 64-bit words based on length of non-whitespace portion of the string: */
	*len = 1;
	LEN_MAX = (uint32)ceil( (imax-i)/log10_base );
	mi64_vec = (uint64 *)calloc(LEN_MAX+1, sizeof(uint64));	/* 01/09/2009: Add an extra zero-pad element here as workaround for bug in mi64_div called with differing-length operands */
	imin = i;
	for(i = imin; i < imax; i++)
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

#if MI64_DEBUG
	#define MI64_POW_DBG	0	// Set nonzero to enable debug-print in the mi64 modpow functions below
#endif
/*
Function to find 2^(-p) mod q, where p and q are both base-2^64 multiword unsigned integers.
Uses a Montgomery-style modmul with a power-of-2 modulus = 2^len. Result (optionally) returned in res[],
assumed to have allocated length at least as large as q[]:

The key 3-operation sequence here is as follows:

	SQR_LOHI(x,lo,hi);	// Input x has len words, lo/hi have len words each
	MULL(lo,qinv,lo);	// Inputs lo & qinv, and output (overwrites lo) have len words
	MULH(q,lo,lo);	// Inputs q & lo, and output (overwrites lo) have len words.
*/
uint32 mi64_twopmodq(const uint64 p[], uint32 len_p, const uint64 k, uint64 q[], uint32 len, uint64*res)
{
#if MI64_POW_DBG
	uint32 dbg = (k==4296452645);
#endif
	 int32 j, lenp_save = -1, lenq_save = -1;	/* Current-Bit index j needs to be signed because of the LR binary exponentiation. */
	uint32 idum, pbits;
	uint64 lead_chunk, lo64, cyout;
	static uint64 *qhalf = 0x0, *qinv = 0x0, *x = 0x0, *lo = 0x0, *hi = 0x0;
	static uint64 *psave = 0x0, *pshift = 0x0;
	static uint32 lenP, lenQ, qbits, log2_numbits, start_index, zshift, first_entry = TRUE;

	lenP = mi64_getlen(p, len_p);
	DBG_ASSERT(HERE, lenP > 0, "0 exponent");
	lenQ = mi64_getlen(q, len);
	DBG_ASSERT(HERE, lenQ > 0, "0 modulus!");
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
	if(dbg) { printf("mi64_twopmodq: k = %llu, len = %u, lenQ = %u\n",k,len,lenQ); }
#endif
	qbits = lenQ << 6;
	mi64_shrl(q, qhalf, 1, lenQ);	/* (q >> 1) = (q-1)/2, since q odd. */

	if(first_entry || !mi64_cmp_eq(p, psave, len_p))
	{
		first_entry = FALSE;
		mi64_set_eq(psave, p, len_p);
		/* pshift = p + len*64 */
		pshift[lenP] = mi64_add_scalar(p, lenP*64, pshift, lenP);
		DBG_ASSERT(HERE, !pshift[lenP], "pshift overflows!");

#if MI64_POW_DBG
	if(dbg) { printf("Init: k = %llu, lenP = %u, lenQ = %u\n",k,lenP,lenQ); }
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
		ASSERT(HERE, pbits >= log2_numbits, "leadz64!");
	//	if(pbits >= 64)
			lead_chunk = lo64>>(64-log2_numbits);
	//	else
	//		lead_chunk = lo64>>(pbits-log2_numbits);	**** lead_chunk now normalized to have >= 64 bits even if arg < 2^64 ****

		if(lead_chunk >= qbits)
		{
			lead_chunk >>= 1;
	#if MI64_POW_DBG
		if(dbg) { printf("lead%u = %llu\n", log2_numbits-1,lead_chunk); }
	#endif
			start_index = pbits-(log2_numbits-1);	/* Use only the leftmost log2_numbits-1 bits */
		}
		else
		{
	#if MI64_POW_DBG
		if(dbg) { printf("lead%u = %llu\n", log2_numbits  ,lead_chunk); }
	#endif
			start_index = pbits-log2_numbits;
		}

		zshift = (qbits-1) - lead_chunk;
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		mi64_negl(pshift, pshift, lenP);	/* ~pshift[] */
	#if MI64_POW_DBG
		if(dbg) { printf("pshift = %s\n", &s0[convert_mi64_base10_char(s0, pshift, lenP)]); }
	#endif
	}

	/*
	Find modular inverse (mod 2^qbits) of q in preparation for modular multiply.
	q must be odd for Montgomery-style modmul to work.

	Init qinv = q. This formula returns the correct bottom 4 bits of qinv,
	and we double the number of correct bits on each of the subsequent iterations.
	*/
	ASSERT(HERE, (q[0] & (uint64)1) == 1, "q must be odd!");
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
	// Check the computed inverse:
	mi64_mul_vector_lo_half(q, qinv, x, lenQ);
	ASSERT(HERE, mi64_cmp_eq_scalar(x, 1ull, lenQ), "Bad Montmul inverse!");
#if MI64_POW_DBG
	if(dbg) {
		printf("q    = %s\n", &s0[convert_mi64_base10_char(s0, q   , lenQ)]);
		printf("qinv = %s\n", &s0[convert_mi64_base10_char(s0, qinv, lenQ)]);
		printf("start_index = %3d\n", start_index);
	}
#endif

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if MI64_POW_DBG
	if(dbg) { printf("zshift = %d, lo = (qinv << zshift)...\n", zshift); }
#endif
	mi64_shl(qinv, lo, zshift, lenQ);
#if MI64_POW_DBG
	if(dbg) { printf("lo = %s\n", &s0[convert_mi64_base10_char(s0, lo, lenQ)]); }
#endif
	mi64_mul_vector_hi_half(q, lo, lo, lenQ);
#if MI64_POW_DBG
	if(dbg) { printf("q*lo/2^%u = %s\n", (lenQ<<6), &s0[convert_mi64_base10_char(s0, lo, lenQ)]); }
#endif

	/* hi = 0 in this instance, which simplifies things. */
	cyout = mi64_sub(q, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");

	if(mi64_test_bit(pshift, j))
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(mi64_cmpugt(x, qhalf, lenQ))
		{
			cyout = mi64_add(x, x, x, lenQ);
			cyout = mi64_sub(x, q, x, lenQ);
		}
		else
		{
			cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
	}
#if MI64_POW_DBG
	if(dbg) { printf("x0 = %s\n", &s0[convert_mi64_base10_char(s0, x, lenQ)] ); }
#endif

	for(j = start_index-2; j >= 0; j--)
	{
	#if MI64_POW_DBG
		if(dbg) { printf("J = %d:\n",j); }
	#endif
		/*...x^2 mod q is returned in x. */
		mi64_mul_vector(x, lenQ, x, lenQ, lo, &idum);
		mi64_mul_vector_lo_half(lo, qinv, lo, lenQ);
		mi64_mul_vector_hi_half(q, lo, lo, lenQ);

	#if MI64_POW_DBG
		if(dbg) { printf("lo = %s\n", &s0[convert_mi64_base10_char(s0, lo, lenQ)] ); }
	#endif

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(mi64_cmpult(hi, lo, lenQ))
		{
			cyout = mi64_sub(q, lo, lo, lenQ);
			cyout = mi64_add(lo, hi, x, lenQ);
		}
		else
		{
			cyout = mi64_sub(hi, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
	#if MI64_POW_DBG
		if(dbg) { printf("x = %s\n",&s0[convert_mi64_base10_char(s0, x, lenQ)]); }
	#endif

		if(mi64_test_bit(pshift, j))
		{
			DBG_ASSERT(HERE, mi64_cmpult(x, q, lenQ), "x >= q");
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(mi64_cmpugt(x, qhalf, lenQ))
			{
				cyout = mi64_add(x, x, x, lenQ);
				cyout = mi64_sub(x, q, x, lenQ);
			}
			else
			{
				cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			}
		#if MI64_POW_DBG
			if(dbg) { printf("2x...\n"); }
		#endif
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	cyout = mi64_add(x, x, x, lenQ);	/* In the case of interest, x = (q+1)/2, so x + x cannot overflow. */
	cyout = mi64_sub(x, q, x, lenQ);
	if(res != 0x0)
	{
		mi64_set_eq(res, x, lenQ);
	}
#if MI64_POW_DBG
	if(dbg) { printf("xout = %s\n",&s0[convert_mi64_base10_char(s0, x, lenQ)]); }
#endif
	return mi64_cmp_eq_scalar(x, 1ull, lenQ);
}


/****************************************************************************************************/
/* Specialized version of mi64_twopmodq for moduli q = 2.k.M(p) + 1, where M(p) is a Mersenne prime */
/****************************************************************************************************/
uint32 mi64_twopmodq_qmmp(const uint64 p, const uint64 k, uint64*res)//, uint32 half_mul_words)	<*** 05/02/2012: Compute required array length needed for modulus on the fly, adjust as needed
{
#if MI64_POW_DBG
	uint32 dbg = (k==4296452645);
#endif
	 int32 j;
	uint32 lenP, lenQ, pbits;
	uint64 k2 = k+k, lead_chunk, lo64, cyout;
	static uint64 *q = 0x0, *qhalf = 0x0, *qinv = 0x0, *x = 0x0, *lo = 0x0, *hi = 0x0;
	static uint64 psave = 0, *pshift = 0x0;
	static uint32 lenQ_save = 0, qbits, log2_numbits, start_index, zshift;
	static uint32  first_entry = TRUE;

	// Quick computation of number of uint64 needed to hold current q:
	ASSERT(HERE, (k != 0) && ((k2>>1) == k), "2*k overflows!");	// Make sure 2*k does not overflow
	j = (p+1)&63;	// p+1 mod 64, needed since q = 2*k*MMp+1 ~= k*MM(p+1)
	lenP = ((p+1) + 63)>>6;	// #64-bit words needed
	lo64 = k;		// Copy of k
	if( ((lo64 << j) >> j) != k ) {	// 2*k*MMp crosses a word boundary, need one more word than required by MMp to store it
		lenQ = lenP + 1;
	} else {
		lenQ = lenP;
	}
#if MI64_POW_DBG
	if(dbg) { printf("mi64_twopmodq_qmmp: k = %llu, lenP = %u, lenQ = %u\n",k,lenP,lenQ); }
#endif

	if(first_entry || (p != psave) || (lenQ != lenQ_save))
	{
		first_entry = FALSE;
		psave = p;
		free((void *)pshift);
		pshift = (uint64 *)calloc((lenP+1), sizeof(uint64));
		pshift[0] = 1;
		mi64_shl(pshift, pshift, p, lenP);	// 2^p
		mi64_sub_scalar(pshift, 1, pshift, lenP);	// M(p) = 2^p-1
		/* pshift = p + len*64: */
		pshift[lenP] = mi64_add_scalar(pshift, lenP*64, pshift, lenP);
		ASSERT(HERE, !pshift[lenP], "pshift overflows!");
	#if MI64_POW_DBG
		if(dbg) { printf("mi64_twopmodq_qmmpInit: k = %llu, lenP = %u, lenQ = %u\n",k,lenP,lenQ); }
	#endif
		lenQ_save = lenQ;
		free((void *)q    );
		free((void *)qhalf);
		free((void *)qinv );
		free((void *)x    );
		free((void *)lo   );
		q      = (uint64 *)calloc((lenQ), sizeof(uint64));
		qhalf  = (uint64 *)calloc((lenQ), sizeof(uint64));
		qinv   = (uint64 *)calloc((lenQ), sizeof(uint64));
		x      = (uint64 *)calloc((lenQ), sizeof(uint64));
		lo   = (uint64 *)calloc((2*lenQ), sizeof(uint64));
		hi   = lo + lenQ;	/* Pointer to high half of double-wide product */

		qbits = lenQ << 6;
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
		pbits = mi64_extract_lead64(pshift,lenP,&lo64);
		ASSERT(HERE, pbits >= log2_numbits, "leadz64!");
	//	if(pbits >= 64)
			lead_chunk = lo64>>(64-log2_numbits);
	//	else
	//		lead_chunk = lo64>>(pbits-log2_numbits);	**** lead_chunk now normalized to have >= 64 bits even if arg < 2^64 ****

		if(lead_chunk >= qbits)
		{
			lead_chunk >>= 1;
	#if MI64_POW_DBG
		if(dbg) { printf("lead%u = %llu\n", log2_numbits-1,lead_chunk); }
	#endif
			start_index = pbits-(log2_numbits-1);	/* Use only the leftmost log2_numbits-1 bits */
		}
		else
		{
	#if MI64_POW_DBG
		if(dbg) { printf("lead%u = %llu\n", log2_numbits  ,lead_chunk); }
	#endif
			start_index = pbits-log2_numbits;
		}

		zshift = (qbits-1) - lead_chunk;	// For MMp this = 2^(log2_numbits-1) = 100...000.
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		mi64_negl(pshift, pshift, lenP);	/* ~pshift[] */
	#if MI64_POW_DBG
		if(dbg) { printf("pshift = %s\n", &s0[convert_mi64_base10_char(s0, pshift, lenP)]); }
	#endif
	}

	/* compute q */
	memset(q, 0ull, (lenQ<<3));
#if 1	// Oddly, the "slower" way here gives a significantly faster binary:
	q[0] = 1; mi64_shl(q, q, p, lenQ);
	mi64_sub_scalar(q, 1, q, lenQ);	// M(p) = 2^p-1
	cyout = mi64_mul_scalar(q, k2, q, lenQ);	ASSERT(HERE, !cyout, "2.k.M(p) overflows!");	// 2.k.M(p)
	ASSERT(HERE, 0 != q[lenQ-1], "Excessive word size allocated for q!");
	mi64_add_scalar(q, 1ull, q, lenQ);	// q = 2.k.M(p) + 1
	mi64_shrl(q, qhalf, 1, lenQ);	/* (q >> 1) = (q-1)/2, since q odd. */
#else
	// 2^p; cheaper than the "long way": q[0] = 1; mi64_shl(q, q, p, lenQ);
	j = p>>6;	// p/64; the set-bit in 2^p goes into the (j)th word of q[]
	q[j] = ( 1ull << (p-(j<<6)) );
	mi64_sub_scalar(q, 1, q, lenQ);	// M(p) = 2^p-1
	cyout = mi64_mul_scalar(q, k2, q, lenQ);	ASSERT(HERE, !cyout, "2.k.M(p) overflows!");	// 2.k.M(p)
	ASSERT(HERE, 0 != q[lenQ-1], "Excessive word size allocated for q!");
	mi64_add_scalar(q, 1ull, q, lenQ);	// q = 2.k.M(p) + 1
	mi64_div2(q, qhalf, lenQ);	/* (q >> 1) = (q-1)/2, since q odd. */
#endif
	/*
	Find modular inverse (mod 2^qbits) of q in preparation for modular multiply.
	q must be odd for Montgomery-style modmul to work.

	Init qinv = q. This formula returns the correct bottom 4 bits of qinv,
	and we double the number of correct bits on each of the subsequent iterations.
	*/
	DBG_ASSERT(HERE, (q[0] & (uint64)1) == 1, "q must be odd!");
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
		mi64_add_scalar(x, 2ull, x, lenQ);
		mi64_mul_vector_lo_half(qinv, x, qinv, lenQ);
	}
	// Check the computed inverse:
	mi64_mul_vector_lo_half(q, qinv, x, lenQ);
	ASSERT(HERE, mi64_cmp_eq_scalar(x, 1ull, lenQ), "Bad Montmul inverse!");
#if MI64_POW_DBG
	if(dbg) {
		printf("q    = %s\n", &s0[convert_mi64_base10_char(s0, q   , lenQ)]);
		printf("qinv = %s\n", &s0[convert_mi64_base10_char(s0, qinv, lenQ)]);
		printf("start_index = %3d\n", start_index);
	}
#endif

	/* Since zstart is a power of two < 2^192, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
#if MI64_POW_DBG
	if(dbg) { printf("zshift = %d, lo = (qinv << zshift)...\n", zshift); }
#endif
	mi64_shl(qinv, lo, zshift, lenQ);
#if MI64_POW_DBG
	if(dbg) { printf("lo = %s\n", &s0[convert_mi64_base10_char(s0, lo, lenQ)]); }
#endif
	mi64_mul_vector_hi_half(q, lo, lo, lenQ);
#if MI64_POW_DBG
	if(dbg) { printf("q*lo/2^%u = %s\n", (lenQ<<6), &s0[convert_mi64_base10_char(s0, lo, lenQ)]); }
#endif

	/* hi = 0 in this instance, which simplifies things. */
	cyout = mi64_sub(q, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");

	// mi64_test_bit(pshift, j) always true for this portion of MMp powering
	DBG_ASSERT(HERE, mi64_test_bit(pshift, j), "pshift bit = 0 for pre-loop step!");
	DBG_ASSERT(HERE, mi64_cmpult(x, q, lenQ), "x >= q");
	/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
	if(mi64_cmpugt(x, qhalf, lenQ))
	{
		cyout = mi64_add(x, x, x, lenQ);
		cyout = mi64_sub(x, q, x, lenQ);
	}
	else
	{
		cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
	}

#if MI64_POW_DBG
	if(dbg) { printf("x0 = %s\n", &s0[convert_mi64_base10_char(s0, x, lenQ)] ); }
#endif

	for(j = start_index-2; j >= log2_numbits; j--)
	{
	#if MI64_POW_DBG
		if(dbg) { printf("J = %d:\n",j); }
	#endif
		/*...x^2 mod q is returned in x. */
		mi64_sqr_vector(x, lo, lenQ);
	#if MI64_POW_DBG
		if(dbg) {
			printf("lo = %s\n",&s0[convert_mi64_base10_char(s0, lo, lenQ)]);
			printf("hi = %s\n",&s0[convert_mi64_base10_char(s0, hi, lenQ)]);
		}
	#endif
		mi64_mul_vector_lo_half(lo, qinv, x, lenQ);
	#if MI64_POW_DBG
		if(dbg) { printf(" x = %s\n",&s0[convert_mi64_base10_char(s0,  x, lenQ)]); }
	#endif
		mi64_mul_vector_hi_qmmp(x, p, k, lo, (lenQ<<6));
	#if MI64_POW_DBG
		if(dbg) { printf("lo = %s\n",&s0[convert_mi64_base10_char(s0, lo, lenQ)]); }
	#endif

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(mi64_cmpult(hi, lo, lenQ))
		{
			cyout = mi64_sub(q, lo, lo, lenQ);
			cyout = mi64_add(lo, hi, x, lenQ);
		}
		else
		{
			cyout = mi64_sub(hi, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}

	#if MI64_POW_DBG
		if(dbg) { printf("x = %s\n",&s0[convert_mi64_base10_char(s0, x, lenQ)]); }
	#endif
		// mi64_test_bit(pshift, j) always true for this portion of MMp powering
		DBG_ASSERT(HERE, mi64_test_bit(pshift, j), "pshift bit = 0!");
	#if MI64_POW_DBG
		if(!mi64_cmpult(x, q, lenQ)) {
			printf("x < q test failed for k = %llu, j = %u!\n",k,j);
		}
		if(dbg) { printf("2x...\n"); }
	#else
		DBG_ASSERT(HERE, mi64_cmpult(x, q, lenQ), "x >= q");
	#endif

		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		if(mi64_cmpugt(x, qhalf, lenQ))
		{
			cyout = mi64_add(x, x, x, lenQ);
			cyout = mi64_sub(x, q, x, lenQ);
		}
		else
		{
			cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}
	}
	for(; j >= 0; j--)
	{
		/*...x^2 mod q is returned in x. */
		mi64_sqr_vector(x, lo, lenQ);
		mi64_mul_vector_lo_half(lo, qinv, x, lenQ);
		mi64_mul_vector_hi_qmmp(x, p, k, lo,(lenQ<<6));

		if(mi64_cmpult(hi, lo, lenQ))
		{
			cyout = mi64_sub(q, lo, lo, lenQ);
			cyout = mi64_add(lo, hi, x, lenQ);
		}
		else
		{
			cyout = mi64_sub(hi, lo, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
		}

		if((pshift[0] >> j) & (uint64)1)
		{
			DBG_ASSERT(HERE, mi64_cmpult(x, q, lenQ), "x >= q");
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(mi64_cmpugt(x, qhalf, lenQ))
			{
				cyout = mi64_add(x, x, x, lenQ);
				cyout = mi64_sub(x, q, x, lenQ);
			}
			else
			{
				cyout = mi64_add(x, x, x, lenQ);	DBG_ASSERT(HERE, cyout == 0ull, "");
			}
		}
	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	mi64_add(x, x, x, lenQ);	/* In the case of interest, x = (q+1)/2, so x + x cannot overflow. */
	mi64_sub(x, q, x, lenQ);
	if(res != 0x0)
	{
		mi64_set_eq(res, x, lenQ);
	}
#if MI64_POW_DBG
	if(dbg) { printf("mi64_twopmodq_qmmp: xout = %s\n", &s0[convert_mi64_base10_char(s0, x, lenQ)] ); }
#endif
	return mi64_cmp_eq_scalar(x, 1ull, lenQ);
}

