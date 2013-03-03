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

#include "factor.h"
#include "twopmodq80.h"

#define STREAMLINE_3WORD_MODMUL

#ifdef USE_GPU

	#include "gpu_iface.h"

	#undef FAC_DEBUG
	#define FAC_DEBUG	1	// Set nonzero and populate the relevant 'dbg =' lines in the code below with desired known-factor target k-values to enable debug printing.
	#if FAC_DEBUG

	/*
	Divide-with-Remainder of x by y, where x is a 128-bit unsigned (vector) integer and y a 32-bit unsigned scalar.
	Returns (x - x%y)/y in x, 32-bit remainder in the function result.

	If you only want the remainder, not to perform the divide, call x128_mod_y32 instead.
	*/
	__device__
	uint32 X128_DIV_Y32(uint128 *x, uint32 y)
	{
		uint64 cy, rem, xlomody, tsum;
		uint64 two64divy, two64mody;

		two64divy = 0x8000000000000000ull/y;
		two64divy = (two64divy + two64divy);

		two64mody = 0x8000000000000000ull%y;
		/* To save a second (expensive) standard library mod call,
		double and subtract y, then re-add y if the result underflows: */
		two64mody = (two64mody + two64mody) - y;
		cy = (two64mody >> 63);
		two64mody += (-cy) & y;
		two64divy += (cy == 0);
/*printf("INIT: two64divy, two64mody = %20llu %20llu\n\n", two64divy, two64mody); */

		/* Divide high digit by y, storing remainder in cy: */
		cy = (x->d1)%y;
		(x->d1) /= y;

		/* Remainder (must calculate this before modifying (x->d0), obviously): */
		xlomody = (x->d0)%y;
		tsum = cy*two64mody + xlomody;
		rem = tsum%y;

		/* Low digit of result: we must separately divide (x->d0) by y
		(making sure to add (x->d0)%y to  cy*two64mody first, so as not to drop a digit)
		because x->d0 may be as large as 2^64-1, and adding cy*two64mody
		prior to dividing risks unsigned integer overflow:
		*/
		(x->d0) = cy*two64divy + tsum/y + (x->d0)/y;
	/*printf("%20llu %20llu %2llu %2llu\n", x->d1, x->d0, cy, rem); */
		return (uint32)rem;
	}

	/*
	Returns decimal character representation of a base-2^64 2-word unsigned int in char_buf,
	and the position of the leftmost nonzero digit (e.g. if the caller wants to print in
	left-justified form) in the function result.
	*/
	__device__
	int	CONVERT_UINT128_BASE10_CHAR(char char_buf[], uint128 q128)
	{
		const int CHAROFFSET = '0';
		uint32 i, n_dec_digits = 0;
		char c;
		/* 2^128 has 39 decimal digits: */
		uint32 MAX_DIGITS = 39;

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
			if((q128.d0 || q128.d1) || n_dec_digits == 0)
			{
				c = X128_DIV_Y32(&q128, (uint32)10) + CHAROFFSET;
				n_dec_digits++;
			}
			else
				c = ' ';

			char_buf[(MAX_DIGITS - 1) - i] = c;
		}

		return (int)MAX_DIGITS-n_dec_digits;
	}

	__device__
	int	CONVERT_UINT96_BASE10_CHAR(char char_buf[], uint96 q96)
	{
		uint128 q128;
		q128.d0 = q96.d0;
		q128.d1 = (uint64)q96.d1;
		return CONVERT_UINT128_BASE10_CHAR(char_buf, q128);
	}
	#endif

#endif

#ifdef USE_GPU
	#include "gpu_iface.h"

	#define USE_KERN_INTRINSICS	1	// Use CUDA kernel IMUL intrinsics, rather than generic 32-bit C macros
	#if USE_KERN_INTRINSICS

		// Easier to just undef and redef these wide-int-mul macros rather than overhaul header-file infrastructure:
		#undef MUL_LOHI32
		#undef SQR_LOHI32
		#undef	__MULL32
		#undef	__MULH32
		#undef MULL32
		#undef MULH32

		#undef MUL_LOHI64
		#undef MUL64x32
		#undef SQR_LOHI64
		#undef	__MULL64
		#undef	__MULH64
		#undef MULL64
		#undef MULH64

		#undef leadz32
		#undef leadz64
		#undef trailz32
		#undef trailz64
		#undef popcnt32
		#undef popcnt64

		#define MUL_LOHI32(_x,_y,_lo,_hi) {uint32 _t = (uint32)(_x)*(uint32)(_y); _hi = __umulhi((uint32)(_x), (uint32)(_y)); _lo = _t;}
		#define SQR_LOHI32(_x,   _lo,_hi)	MUL_LOHI32(_x,_x,_lo,_hi)
		#define	__MULL32(	 _x,_y      )	         (uint32)(_x)*(uint32)(_y)
		#define	__MULH32(	 _x,_y      )	__umulhi((uint32)(_x),(uint32)(_y))
		#define MULL32(	 _x,_y, _lo     )	_lo =          (uint32)(_x)*(uint32)(_y)
		#define MULH32(  _x,_y,      _hi)	_hi = __umulhi((uint32)(_x),(uint32)(_y))

		#define MUL_LOHI64(_x,_y,_lo,_hi) {uint64 _t = (_x)*(_y); _hi = __umul64hi((_x), (_y));	_lo = _t;}
		#define MUL64x32(  _x,_y,_lo,_hi)	MUL_LOHI64(_x,(uint64)_y,_lo,_hi)	/* Need 32-bit-optimized version of t his */
		#define SQR_LOHI64(_x,   _lo,_hi)	MUL_LOHI64(_x,_x,_lo,_hi)
		#define	__MULL64(	 _x,_y      )	           (_x)*(_y)
		#define	__MULH64(	 _x,_y      )	__umul64hi((_x),(_y))
		#define MULL64(	 _x,_y, _lo     )	_lo =            (_x)*(_y)
		#define MULH64(  _x,_y,      _hi)	_hi = __umul64hi((_x),(_y))

		#define leadz32(_x)	( 32 - __ffs((_x)) )
		#define leadz64(_x)	( 64 - __ffsll((_x)) )
		#define trailz32(_x)	( 32 - __clz((_x)) )
		#define trailz64(_x)	( 64 - __clzll((_x)) )
		#define popcnt32(_x)	__popc((_x))
		#define popcnt64(_x)	__popcll((_x))

	#else	// Generic 32-bit C macros for wide IMUL

		#undef MUL_LOHI64
		#undef SQR_LOHI64
		#undef MULH64

		/* Generic 128-bit product macros: represent the inputs as
		x = a + b*2^32, y = c + d*2^32, and then do 4 MULs and a bunch of
		adds-with-carry to get x*y = b*d*2^64 + (a*d + b*c)*2^32 + a*c .
		Actual calling arguments are assumed to be 64-bit ints - user must
		make sure this is true. Result written into lo, hi.
		*/
		#define MUL_LOHI64(_x,_y,_lo,_hi)\
		{\
			uint64 _a,_b,_c,_d,_ac,_ad,_bc;				\
			/*_a = (_x) & (uint64)0x00000000ffffffff;*/		\
			_a = ((_x)<<32)>>32;								\
			_b =  (_x)>>32;									\
			/*_c = (_y) & (uint64)0x00000000ffffffff;*/		\
			_c = ((_y)<<32)>>32;								\
			_d =  (_y)>>32;									\
			/* Calculate 4 subproducts in order in which they are first used */\
			_ac = _a*_c;										\
			_bc = _b*_c;										\
			_hi = _b*_d;										\
			_ad = _a*_d;										\
			_lo  =  _ac;		/* use _lo to store copy of _ac */\
			_ac +=         (_bc<<32);	_hi += (_ac < _lo);	\
			_lo  =  _ac + (_ad<<32);	_hi += (_lo < _ac);	\
							_hi += (_bc>>32) + (_ad>>32);	\
		}

		/* Generic 128-bit squaring algorithm: represent the input as
		x = a + b*2^32, and then do 3 MULs and a bunch of add-with-carries
		to get x^2 = b^2*2^64 + a*b*2^33 + a^2 . Actual calling arguments
		are assumed to be 64-bit ints - user must make sure this is true.
		Result written into lo, hi.
		*/
		#define SQR_LOHI64(_x,_lo,_hi)\
		{\
			uint64 _a,_b,_aa,_ab;			\
			/*_a = (_x) & (uint64)0x00000000ffffffff;*/\
			_a = ((_x)<<32)>>32;				\
			_b =  (_x)>>32;					\
			_aa = _a*_a;						\
			_ab = _a*_b;						\
			_hi  = _b*_b;					\
			_lo  = _aa + (_ab<<33);			\
			_hi += (_ab>>31) + (_lo < _aa);	\
		}

		/* Generic 128-bit multiply algorithm, returning only the upper 64 bits (high part) of the result.
		Actual calling arguments are assumed to be 64-bit ints - user must make sure this is true.
		Result written into hi.
		*/
		#define MULH64(_x, _y,_hi)\
		{\
			uint64 _tt;					\
			MUL_LOHI64((_x), (_y), _tt, _hi);	\
		}

		__device__ uint32 leadz32(uint32 i)
		{
			uint32 lz, k, shift, ones_mask = 0xFFFFFFFF;
			if(i == 0) return 32;
			if((int32)i < 0) return 0;

			lz    =  0;
			shift = 16;
			k     = 16;
			while(k > 0)
			{
				if( (ones_mask >> shift) < i )
				{
					k >>= 1;
					shift -= k;
				}
				else
				{
					lz += k;
					k >>= 1;
					shift += k;
				}
			}

			DBG_ASSERT(HERE, ( (i << lz) >> lz == i ),"ERROR A in leadz32");
			DBG_ASSERT(HERE, ( i >> (32-lz)    == 0 ),"ERROR B in leadz32");

			return lz;
		}

	#endif	// USE_KERN_INTRINSICS

#endif

/***************/

/* Since GPU code tests many candidates per batch, adopt the following return-value convention:
			|	= 0		No factors found in current batch
	retval	|	> 0		 1 factor in current batch, whose k-value is returned
			|	< 0		n = (-retval) (>1) factors in current batch, user should re-sieve interval in 'slow' CPU-based mode to find the factor k's
*/

// GPU Vector-Kernel definition:
#ifdef USE_GPU
__global__
#endif
void GPU_TF78(uint64*checksum1, uint64*checksum2, uint64 p, const uint64 kvec[], uint32 n, int64 *retval)
{
	// Inlined-const byte-array storing popcount(i) for i = 0-15:
    //                                0000,0001,0010,0011,0100,0101,0110,0111,1000,1001,1010,1011,1100,1101,1110,1111;
	const uint8 nfac_in_4batch[16] = {   0,   1,   1,   2,   1,   2,   2,   3,   1,   2,   2,   3,   2,   3,   3,   4};
	 int64 nfac = 0;
	uint64 res;
	const uint64 *k0,*k1,*k2,*k3;
	int i;
#if defined(USE_GPU) && FAC_DEBUG

	i = blockDim.x * blockIdx.x + threadIdx.x;
//	if(i < 10) {
	printf("GPU block %d[dim %d], thread %d ==> seq-thread %d ... \n", blockIdx.x, blockDim.x, threadIdx.x, i);
	k0 = kvec;
	k1 = k0+1;
	k2 = k0+2;
	k3 = k0+3;
	printf("k0-3 = %20llu, %20llu, %20llu, %20llu\n", *k0, *k1, *k2, *k3);
//	}
#endif
	for(i = 0; i < n; i += 4) {
		// Each twopmodq*_q4 call processes 4 k's:
		k0 = kvec+i;
		k1 = k0+1;
		k2 = k0+2;
		k3 = k0+3;
		twopmodq78_q4_GPU(checksum1, checksum2, p, *k0, *k1, *k2, *k3, &res);
		if(res) {
			// If nonzero return value [res], proceed as follows:
			// - If nfac == 0 and res has precisely one bit set, copy the the associated k into nfac;
			if(!nfac && nfac_in_4batch[res] == 1) {
				if(res >> 3) {
					nfac = *k3;
				} else if(res >> 2) {
					nfac = *k2;
				} else if(res >> 1) {
					nfac = *k1;
				} else {
					nfac = *k0;
				}
			}
			// - Otherwise, test nfac and proceed as follows:
			//		o If nfac == 0 (i.e. we found > 1 factors in the current quartet and these are the first
			//				found on this invocation of the function), set nfac = -[popcnt(res)].
			else if(nfac == 0) {
				nfac = -nfac_in_4batch[res];
			}
			//		o If nfac > 0 (i.e. we previously found a total of 1 factor, whose k is stored in nfac),
			//				set nfac = -1 - [popcnt(res)] = -[total number of factors found so far];
			else if(nfac > 0) {
				nfac = -nfac_in_4batch[res] - 1;
			}
			//		o If nfac < 0 (i.e. we previously found n > 1 factors, with the total number n = -nfac),
			//				decrement nfac = nfac - [popcnt(res)] = -[total number of factors found so far].
			else if((int64)nfac < 0) {
				nfac -= nfac_in_4batch[res];
			}
		}
	}
	*retval = nfac;
}

/*** 4-trial-factor version ***/
#ifdef USE_GPU
__device__
#endif
void twopmodq78_q4_GPU(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 *retval)
{
#ifdef USE_GPU
	// These globals, normally defined in utll.c and types.c, need device-code-specific redefinitions
	const double TWO64FLOAT = (double)4.0*0x80000000*0x80000000;
//	const double TWO64FLINV = 1.0/TWO64FLOAT;
	const double TWO26FLOAT = (double)0x04000000;
	const double TWO26FLINV = 1.0/TWO26FLOAT;
	const uint96  ONE96  = {(uint64)1, (uint32)0};
#endif
	 int32 j;	/* This needs to be signed because of the LR binary exponentiation. */
	uint32 q32_0, qinv32_0, tmp32_0
		 , q32_1, qinv32_1, tmp32_1
		 , q32_2, qinv32_2, tmp32_2
		 , q32_3, qinv32_3, tmp32_3;
	uint64 lead7, tmp0, tmp1, tmp2, tmp3, r;
	uint96 q0, qinv0, qhalf0, x0
		 , q1, qinv1, qhalf1, x1
		 , q2, qinv2, qhalf2, x2
		 , q3, qinv3, qhalf3, x3;
	uint32 pshift;	/* Require this to be 32-bit here for 32-bit ASM support */
	uint32 start_index, zshift;
#if FAC_DEBUG
	char char_buf[1024];
	int dbg = (k0 == 292936666300) + ((k1 == 292936666300)<<1) + ((k2 == 292936666300)<<2) + ((k3 == 292936666300)<<3);
	dbg = dbg && (k0 != k1);	// Disable debug during self-test mode.
	if(dbg) {
		printf("Here it is: Target k is in k%1d\n", trailz32(dbg));
	}
#endif

#if !defined(STREAMLINE_3WORD_MODMUL) || FAC_DEBUG
	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2;
	double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2;
	double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2;
#else
	double fq0,fq1,fq2, fqinv0,fqinv1,fqinv2, fx0,fx1,fx2, flo0,flo1,flo2, fhi0;
	double gq0,gq1,gq2, gqinv0,gqinv1,gqinv2, gx0,gx1,gx2, glo0,glo1,glo2, ghi0;
	double hq0,hq1,hq2, hqinv0,hqinv1,hqinv2, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0;
	double iq0,iq1,iq2, iqinv0,iqinv1,iqinv2, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0;
#endif
#ifdef STREAMLINE_3WORD_MODMUL
	double fqhalf0,fqhalf1,fqhalf2,fqhalf3;
	double flohi52,glohi52,hlohi52,ilohi52;
	double fq_or_nil_lo26[2], fq_or_nil_hi52[2];
	double gq_or_nil_lo26[2], gq_or_nil_hi52[2];
	double hq_or_nil_lo26[2], hq_or_nil_hi52[2];
	double iq_or_nil_lo26[2], iq_or_nil_hi52[2];
	int fidx,gidx,hidx,iidx;
#endif

	DBG_ASSERT(HERE, (p >> 63) == 0, "twopmodq78_q4 : p must be < 2^63!");
	q0.d0 = q1.d0 = q2.d0 = q3.d0 = p+p;

	MUL_LOHI64(q0.d0, k0, q0.d0, q0.d1);
	MUL_LOHI64(q1.d0, k1, q1.d0, q1.d1);
	MUL_LOHI64(q2.d0, k2, q2.d0, q2.d1);
	MUL_LOHI64(q3.d0, k3, q3.d0, q3.d1);

	q0.d0 += 1;	/* Since 2*p*k even, no need to check for overflow here */
	q1.d0 += 1;
	q2.d0 += 1;
	q3.d0 += 1;
	DBG_ASSERT(HERE, (q0.d1 >> 14) == 0, "twopmodq78_q4 : (q0.d1 >> 14) != 0");
	DBG_ASSERT(HERE, (q1.d1 >> 14) == 0, "twopmodq78_q4 : (q1.d1 >> 14) != 0");
	DBG_ASSERT(HERE, (q2.d1 >> 14) == 0, "twopmodq78_q4 : (q2.d1 >> 14) != 0");
	DBG_ASSERT(HERE, (q3.d1 >> 14) == 0, "twopmodq78_q4 : (q3.d1 >> 14) != 0");
	*checksum1 += q0.d0 + q1.d0 + q2.d0 + q3.d0;

	q32_0 = (uint32)q0.d0;
	q32_1 = (uint32)q1.d0;
	q32_2 = (uint32)q2.d0;
	q32_3 = (uint32)q3.d0;

	/* Convert q to floating form: */
	CVT_UINT78_3WORD_DOUBLE(q0, fq0,fq1,fq2);
	CVT_UINT78_3WORD_DOUBLE(q1, gq0,gq1,gq2);
	CVT_UINT78_3WORD_DOUBLE(q2, hq0,hq1,hq2);
	CVT_UINT78_3WORD_DOUBLE(q3, iq0,iq1,iq2);

#ifdef STREAMLINE_3WORD_MODMUL
	fq_or_nil_lo26[0] = 0.0; fq_or_nil_lo26[1] = fq0; fq_or_nil_hi52[0] = 0.0; fq_or_nil_hi52[1] = fq1 + TWO26FLOAT*fq2;
	gq_or_nil_lo26[0] = 0.0; gq_or_nil_lo26[1] = gq0; gq_or_nil_hi52[0] = 0.0; gq_or_nil_hi52[1] = gq1 + TWO26FLOAT*gq2;
	hq_or_nil_lo26[0] = 0.0; hq_or_nil_lo26[1] = hq0; hq_or_nil_hi52[0] = 0.0; hq_or_nil_hi52[1] = hq1 + TWO26FLOAT*hq2;
	iq_or_nil_lo26[0] = 0.0; iq_or_nil_lo26[1] = iq0; iq_or_nil_hi52[0] = 0.0; iq_or_nil_hi52[1] = iq1 + TWO26FLOAT*iq2;
#endif

	RSHIFT_FAST96(q0, 1, qhalf0);	/* = (q-1)/2, since q odd. */
	RSHIFT_FAST96(q1, 1, qhalf1);
	RSHIFT_FAST96(q2, 1, qhalf2);
	RSHIFT_FAST96(q3, 1, qhalf3);

#ifdef STREAMLINE_3WORD_MODMUL
	fqhalf0 = qhalf0.d1*TWO64FLOAT + qhalf0.d0;
	fqhalf1 = qhalf1.d1*TWO64FLOAT + qhalf1.d0;
	fqhalf2 = qhalf2.d1*TWO64FLOAT + qhalf2.d0;
	fqhalf3 = qhalf3.d1*TWO64FLOAT + qhalf3.d0;
#endif

	DBG_ASSERT(HERE, ((p + 78) >> 32) == 0, "twopmodq78_q2 : (p+78) exceeds 32 bits!");
	pshift = (uint32)p + 78;
	j = leadz32(pshift);
	lead7 = ((pshift<<j) >> 25);	// Leading 7 bits, in [2^6, 2^7-1] = [64, 127]
	/* If leading 7 yield a shift count > MULL-set 'word' size,ß Use only the leftmost 6 bits */
	if(lead7 > 77)
	{
		lead7 >>= 1;			// Input in [78,127], Result in [39,64]
		start_index =  32-j-6;
	} else {
		start_index =  32-j-7;
	}
	zshift = 77 - lead7;	// lead7 now in [39,77], yielding zshift in [0,38]

	/* Streamlined 1st-iteration which uses that the starting value of x = (1 << zshift) to
	replace sqr_lohi(x, lo,hi) / mull(lo,qinv,lo) with (qinv << (zshift*2)),
	followed by a masking-off off the high bits beyond those left by the MULL boundary being used:
	*/
	zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
	DBG_ASSERT(HERE, zshift < 77, "Shift count exceeds limit!");	// 2*zshift in [0,76]
	/* Still need to do mulh(lo,q,x) once we have x in hand */

	pshift = ~pshift;

	/* This formula for inverse-inits has low 4 bits correct: */
	qinv32_0 = (q32_0 + q32_0 + q32_0) ^ (uint32)2;
	qinv32_1 = (q32_1 + q32_1 + q32_1) ^ (uint32)2;
	qinv32_2 = (q32_2 + q32_2 + q32_2) ^ (uint32)2;
	qinv32_3 = (q32_3 + q32_3 + q32_3) ^ (uint32)2;
	/* Do 3 NR-iterations using 32-bit arithmetic to get 32-bit approximation: */
	for(j = 0; j < 3; j++)
	{
		tmp32_0 = q32_0*qinv32_0;
		tmp32_1 = q32_1*qinv32_1;
		tmp32_2 = q32_2*qinv32_2;
		tmp32_3 = q32_3*qinv32_3;
		qinv32_0 = qinv32_0*((uint32)2 - tmp32_0);
		qinv32_1 = qinv32_1*((uint32)2 - tmp32_1);
		qinv32_2 = qinv32_2*((uint32)2 - tmp32_2);
		qinv32_3 = qinv32_3*((uint32)2 - tmp32_3);
	}
	/* Copy 32-bit result into 64-bit field: */
	qinv0.d0 = (uint64)qinv32_0;
	qinv1.d0 = (uint64)qinv32_1;
	qinv2.d0 = (uint64)qinv32_2;
	qinv3.d0 = (uint64)qinv32_3;
	/* And do 1 further iteration yielding a 64-bit result: */
	tmp0 = q0.d0*qinv0.d0;
	tmp1 = q1.d0*qinv1.d0;
	tmp2 = q2.d0*qinv2.d0;
	tmp3 = q3.d0*qinv3.d0;
	qinv0.d0 = qinv0.d0*((uint64)2 - tmp0);
	qinv1.d0 = qinv1.d0*((uint64)2 - tmp1);
	qinv2.d0 = qinv2.d0*((uint64)2 - tmp2);
	qinv3.d0 = qinv3.d0*((uint64)2 - tmp3);

	/* One more specialzed iteration to get desired 78-bit mod-inverse: */
	MULH64(q0.d0, qinv0.d0, tmp0);
	MULH64(q1.d0, qinv1.d0, tmp1);
	MULH64(q2.d0, qinv2.d0, tmp2);
	MULH64(q3.d0, qinv3.d0, tmp3);
	qinv0.d1 = -(uint32)qinv0.d0*(q0.d1*(uint32)qinv0.d0 + (uint32)tmp0);
	qinv1.d1 = -(uint32)qinv1.d0*(q1.d1*(uint32)qinv1.d0 + (uint32)tmp1);
	qinv2.d1 = -(uint32)qinv2.d0*(q2.d1*(uint32)qinv2.d0 + (uint32)tmp2);
	qinv3.d1 = -(uint32)qinv3.d0*(q3.d1*(uint32)qinv3.d0 + (uint32)tmp3);

	/* 78 bits: */
	qinv0.d1 &= 0x00003fff;	/* Only want the lower 14 bits here */
	qinv1.d1 &= 0x00003fff;
	qinv2.d1 &= 0x00003fff;
	qinv3.d1 &= 0x00003fff;

	/* Convert qinv to floating form: */
	CVT_UINT78_3WORD_DOUBLE(qinv0, fqinv0,fqinv1,fqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv1, gqinv0,gqinv1,gqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv2, hqinv0,hqinv1,hqinv2);
	CVT_UINT78_3WORD_DOUBLE(qinv3, iqinv0,iqinv1,iqinv2);

	/* Since zstart is a power of two < 2^78, use a streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL78(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	LSHIFT96(qinv0, zshift, x0);	x0.d1 &= 0x00003fff;	/* Only want the lower 14 bits here */
	LSHIFT96(qinv1, zshift, x1);	x1.d1 &= 0x00003fff;
	LSHIFT96(qinv2, zshift, x2);	x2.d1 &= 0x00003fff;
	LSHIFT96(qinv3, zshift, x3);	x3.d1 &= 0x00003fff;

	CVT_UINT78_3WORD_DOUBLE(x0, flo0,flo1,flo2);
	CVT_UINT78_3WORD_DOUBLE(x1, glo0,glo1,glo2);
	CVT_UINT78_3WORD_DOUBLE(x2, hlo0,hlo1,hlo2);
	CVT_UINT78_3WORD_DOUBLE(x3, ilo0,ilo1,ilo2);

	/* MULH96(q,lo,lo); */
	MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
		  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flohi52
		, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glohi52
		, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlohi52
		, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilohi52
	);

	/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
	/* hi = 0 in this instance so simply unconditionally compute (q-l), which simplifies things. */
	fx1 = fq_or_nil_hi52[1] - flohi52;
	gx1 = gq_or_nil_hi52[1] - glohi52;
	hx1 = hq_or_nil_hi52[1] - hlohi52;
	ix1 = iq_or_nil_hi52[1] - ilohi52;

	fx0 = fq0 - flo0;
	gx0 = gq0 - glo0;
	hx0 = hq0 - hlo0;
	ix0 = iq0 - ilo0;

	/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
		x = x + x - ((-(x > qhalf)) & q);
	In FP version replace integer and with array lookup.
	*/
	if((pshift >> j) & (uint64)1)
	{
		/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
		/* Fill FP approximation to x with as many SBs as possible prior to comparing */
		/* TRUE here means will need to subtract q from 2x: */
		fidx = (fqhalf0 < fx1*TWO26FLOAT + fx0);
		gidx = (fqhalf1 < gx1*TWO26FLOAT + gx0);
		hidx = (fqhalf2 < hx1*TWO26FLOAT + hx0);
		iidx = (fqhalf3 < ix1*TWO26FLOAT + ix0);

		fx1 *= 2;	fx0 += fx0;
		gx1 *= 2;	gx0 += gx0;
		hx1 *= 2;	hx0 += hx0;
		ix1 *= 2;	ix0 += ix0;

		fx1 -= fq_or_nil_hi52[fidx];
		gx1 -= gq_or_nil_hi52[gidx];
		hx1 -= hq_or_nil_hi52[hidx];
		ix1 -= iq_or_nil_hi52[iidx];

		fx0 -= fq_or_nil_lo26[fidx];
		gx0 -= gq_or_nil_lo26[gidx];
		hx0 -= hq_or_nil_lo26[hidx];
		ix0 -= iq_or_nil_lo26[iidx];
	}

	/* Normalize the result - use currently-unused x2 coeff as carry: */
	/* Digit 0: */
	fx2 = DNINT(fx0*TWO26FLINV);
	gx2 = DNINT(gx0*TWO26FLINV);
	hx2 = DNINT(hx0*TWO26FLINV);
	ix2 = DNINT(ix0*TWO26FLINV);

	fx0 -= fx2*TWO26FLOAT;
	gx0 -= gx2*TWO26FLOAT;
	hx0 -= hx2*TWO26FLOAT;
	ix0 -= ix2*TWO26FLOAT;

	/* Digit 1: */
	fx1 += fx2;
	gx1 += gx2;
	hx1 += hx2;
	ix1 += ix2;

	fx2 = DNINT(fx1*TWO26FLINV);
	gx2 = DNINT(gx1*TWO26FLINV);
	hx2 = DNINT(hx1*TWO26FLINV);
	ix2 = DNINT(ix1*TWO26FLINV);

	fx1 -= fx2*TWO26FLOAT;
	gx1 -= gx2*TWO26FLOAT;
	hx1 -= hx2*TWO26FLOAT;
	ix1 -= ix2*TWO26FLOAT;
	/* Digit 2 already in x2 term. */

#if FAC_DEBUG
	if(dbg)
	{
		printf("q0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, q0)]);
		printf("q1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, q1)]);
		printf("q2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, q2)]);
		printf("q3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, q3)]);
		printf("\n");
		printf("qinv0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, qinv0)]);
		printf("qinv1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, qinv1)]);
		printf("qinv2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, qinv2)]);
		printf("qinv3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, qinv3)]);
		printf("\n");
		printf("Initial ix0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
		printf("Initial fx0 = %20.5f, %20.5f, %20.5f\n", fx0,fx1,fx2);
		printf("\n");
	}
#endif

	// Don't use the streamlined-first-modmul code in this case because it uses too much int64-based stuff:
	for(j = start_index-2; j >= 0; j--)
	{
	#if FAC_DEBUG
		if(dbg) { printf("j = %2d:\n", j); }
	#endif
		/*...x^2 mod q is returned in x. */

	#ifdef STREAMLINE_3WORD_MODMUL	// Streamlined version of modmul sequence, leaving high 52 bits in one double whenever feasible:

		/* SQR_LOHI96(x,lo,hi); */
		SQR_LOHI78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fx1
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,gx1
			, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hx1
			, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ix1
		);
#if FAC_DEBUG
	if(dbg)
	{
		CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x0);	printf("Multiply output #1.lo0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
		CVT78_3WORD_DOUBLE_UINT96(glo0,glo1,glo2, x1);	printf("Multiply output #1.lo1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
		CVT78_3WORD_DOUBLE_UINT96(hlo0,hlo1,hlo2, x2);	printf("Multiply output #1.lo2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
		CVT78_3WORD_DOUBLE_UINT96(ilo0,ilo1,ilo2, x3);	printf("Multiply output #1.lo3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
		fhi1 = fx1;	fhi2 = 0.;	NORMALIZE78_3WORD_DOUBLE(fhi0,fhi1,fhi2);
		ghi1 = gx1;	ghi2 = 0.;	NORMALIZE78_3WORD_DOUBLE(ghi0,ghi1,ghi2);
		hhi1 = hx1;	hhi2 = 0.;	NORMALIZE78_3WORD_DOUBLE(hhi0,hhi1,hhi2);
		ihi1 = ix1;	ihi2 = 0.;	NORMALIZE78_3WORD_DOUBLE(ihi0,ihi1,ihi2);
		CVT78_3WORD_DOUBLE_UINT96(fhi0,fhi1,fhi2, x0);	printf("Multiply output #1.hi0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
		CVT78_3WORD_DOUBLE_UINT96(ghi0,ghi1,ghi2, x1);	printf("Multiply output #1.hi1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
		CVT78_3WORD_DOUBLE_UINT96(hhi0,hhi1,hhi2, x2);	printf("Multiply output #1.hi2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
		CVT78_3WORD_DOUBLE_UINT96(ihi0,ihi1,ihi2, x3);	printf("Multiply output #1.hi3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
	}
#endif
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE_q4(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
			, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
			, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
		);
#if FAC_DEBUG
	if(dbg)
	{
		CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x0);	printf("Multiply output #2.0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
		CVT78_3WORD_DOUBLE_UINT96(glo0,glo1,glo2, x1);	printf("Multiply output #2.1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
		CVT78_3WORD_DOUBLE_UINT96(hlo0,hlo1,hlo2, x2);	printf("Multiply output #2.2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
		CVT78_3WORD_DOUBLE_UINT96(ilo0,ilo1,ilo2, x3);	printf("Multiply output #2.3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
	}
#endif
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE_LEAVE_HIGH52_UNNORMALIZED_q4(
			  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flohi52
			, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glohi52
			, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlohi52
			, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilohi52
		);

		/* If h < l, then calculate q+(h-l) < q; otherwise calculate h-l. */
	/**** to-do: simply delay upper-52-bit-splitting in MULH78_3WORD_DOUBLE_q4 until after doing 52-bit compare step:
		fx1 = fhi1 + TWO26FLOAT*fhi2;	flohi52 = flo1 + TWO26FLOAT*flo2;
		gx1 = ghi1 + TWO26FLOAT*ghi2;	glohi52 = glo1 + TWO26FLOAT*glo2;
		hx1 = hhi1 + TWO26FLOAT*hhi2;	hlohi52 = hlo1 + TWO26FLOAT*hlo2;
		ix1 = ihi1 + TWO26FLOAT*ihi2;	ilohi52 = ilo1 + TWO26FLOAT*ilo2;
	****/
		fidx = (fx1 < flohi52);
		gidx = (gx1 < glohi52);
		hidx = (hx1 < hlohi52);
		iidx = (ix1 < ilohi52);
 #if FAC_DEBUG
	if(dbg)
	{
		flo1 = flohi52;	flo2 = 0.;	NORMALIZE78_3WORD_DOUBLE(flo0,flo1,flo2);
		glo1 = glohi52;	glo2 = 0.;	NORMALIZE78_3WORD_DOUBLE(glo0,glo1,glo2);
		hlo1 = hlohi52;	hlo2 = 0.;	NORMALIZE78_3WORD_DOUBLE(hlo0,hlo1,hlo2);
		ilo1 = ilohi52;	ilo2 = 0.;	NORMALIZE78_3WORD_DOUBLE(ilo0,ilo1,ilo2);
		CVT78_3WORD_DOUBLE_UINT96(flo0,flo1,flo2, x0);	printf("Multiply output #3.0 = %s, fidx = %u\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)], fidx);
		CVT78_3WORD_DOUBLE_UINT96(glo0,glo1,glo2, x1);	printf("Multiply output #3.1 = %s, gidx = %u\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)], gidx);
		CVT78_3WORD_DOUBLE_UINT96(hlo0,hlo1,hlo2, x2);	printf("Multiply output #3.2 = %s, hidx = %u\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)], hidx);
		CVT78_3WORD_DOUBLE_UINT96(ilo0,ilo1,ilo2, x3);	printf("Multiply output #3.3 = %s, iidx = %u\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)], iidx);
	}
#endif

		fx1 -= flohi52;
		gx1 -= glohi52;
		hx1 -= hlohi52;
		ix1 -= ilohi52;

		fhi0 -= flo0;
		ghi0 -= glo0;
		hhi0 -= hlo0;
		ihi0 -= ilo0;

		fx1 += fq_or_nil_hi52[fidx];
		gx1 += gq_or_nil_hi52[gidx];
		hx1 += hq_or_nil_hi52[hidx];
		ix1 += iq_or_nil_hi52[iidx];

		fx0 = fhi0 + fq_or_nil_lo26[fidx];
		gx0 = ghi0 + gq_or_nil_lo26[gidx];
		hx0 = hhi0 + hq_or_nil_lo26[hidx];
		ix0 = ihi0 + iq_or_nil_lo26[iidx];
		// Normalize result and store in x-vector...use hi0 terms as carries:
		fhi0 = DNINT(fx0*TWO26FLINV);
		ghi0 = DNINT(gx0*TWO26FLINV);
		hhi0 = DNINT(hx0*TWO26FLINV);
		ihi0 = DNINT(ix0*TWO26FLINV);

		fx0 -= fhi0*TWO26FLOAT;
		gx0 -= ghi0*TWO26FLOAT;
		hx0 -= hhi0*TWO26FLOAT;
		ix0 -= ihi0*TWO26FLOAT;

		fx1 += fhi0;
		gx1 += ghi0;
		hx1 += hhi0;
		ix1 += ihi0;

		/* Branchless version of the double-and-subtract-q-if-(2x >= q) sequence:
			x = x + x - ((-(x > qhalf)) & q);
		In FP version replace integer and with array lookup.
		*/
		if((pshift >> j) & (uint64)1)
		{
			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			/* Fill FP approximation to x with as many SBs as possible prior to comparing */
			/* TRUE here means will need to subtract q from 2x: */
			fidx = (fqhalf0 < fx1*TWO26FLOAT + fx0);
			gidx = (fqhalf1 < gx1*TWO26FLOAT + gx0);
			hidx = (fqhalf2 < hx1*TWO26FLOAT + hx0);
			iidx = (fqhalf3 < ix1*TWO26FLOAT + ix0);

			fx1 *= 2;	fx0 += fx0;
			gx1 *= 2;	gx0 += gx0;
			hx1 *= 2;	hx0 += hx0;
			ix1 *= 2;	ix0 += ix0;

			fx1 -= fq_or_nil_hi52[fidx];
			gx1 -= gq_or_nil_hi52[gidx];
			hx1 -= hq_or_nil_hi52[hidx];
			ix1 -= iq_or_nil_hi52[iidx];

			fx0 -= fq_or_nil_lo26[fidx];
			gx0 -= gq_or_nil_lo26[gidx];
			hx0 -= hq_or_nil_lo26[hidx];
			ix0 -= iq_or_nil_lo26[iidx];
		}

		/* Normalize the result - use currently-unused x2 coeff as carry: */
		/* Digit 0: */
		fx2 = DNINT(fx0*TWO26FLINV);
		gx2 = DNINT(gx0*TWO26FLINV);
		hx2 = DNINT(hx0*TWO26FLINV);
		ix2 = DNINT(ix0*TWO26FLINV);

		fx0 -= fx2*TWO26FLOAT;
		gx0 -= gx2*TWO26FLOAT;
		hx0 -= hx2*TWO26FLOAT;
		ix0 -= ix2*TWO26FLOAT;

		/* Digit 1: */
		fx1 += fx2;
		gx1 += gx2;
		hx1 += hx2;
		ix1 += ix2;

		fx2 = DNINT(fx1*TWO26FLINV);
		gx2 = DNINT(gx1*TWO26FLINV);
		hx2 = DNINT(hx1*TWO26FLINV);
		ix2 = DNINT(ix1*TWO26FLINV);

		fx1 -= fx2*TWO26FLOAT;
		gx1 -= gx2*TWO26FLOAT;
		hx1 -= hx2*TWO26FLOAT;
		ix1 -= ix2*TWO26FLOAT;
		/* Digit 2 already in x2 term. */

	  #if FAC_DEBUG
		if(dbg)
		{
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
			CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
			CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
			CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);
			if((pshift >> j) & (uint64)1)
			{
				printf("x0 = %s, *2 .\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
				printf("x1 = %s, *2 .\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
				printf("x2 = %s, *2 .\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
				printf("x3 = %s, *2 .\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
			} else {
				printf("x0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
				printf("x1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
				printf("x2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
				printf("x3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
			}
			printf("\n");
		}
	  #endif

	#else	// Basic impl of 3-word-double-based 78-bit-integer modmul:

		/* SQR_LOHI96(x,lo,hi); */
		SQR_LOHI78_3WORD_DOUBLE_q4(
			  fx0,fx1,fx2, flo0,flo1,flo2, fhi0,fhi1,fhi2
			, gx0,gx1,gx2, glo0,glo1,glo2, ghi0,ghi1,ghi2
			, hx0,hx1,hx2, hlo0,hlo1,hlo2, hhi0,hhi1,hhi2
			, ix0,ix1,ix2, ilo0,ilo1,ilo2, ihi0,ihi1,ihi2
		);
		/* MULL96(lo,qinv,lo); */
		MULL78_3WORD_DOUBLE_q4(
			  flo0,flo1,flo2, fqinv0,fqinv1,fqinv2, flo0,flo1,flo2
			, glo0,glo1,glo2, gqinv0,gqinv1,gqinv2, glo0,glo1,glo2
			, hlo0,hlo1,hlo2, hqinv0,hqinv1,hqinv2, hlo0,hlo1,hlo2
			, ilo0,ilo1,ilo2, iqinv0,iqinv1,iqinv2, ilo0,ilo1,ilo2
		);
		/* MULH96(q,lo,lo); */
		MULH78_3WORD_DOUBLE_q4(
			  fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2
			, gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2
			, hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2
			, iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2
		);

		/* If h < l, then calculate q-l+h < q; otherwise calculate h-l. */
		if(CMPLT78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2))
		{
			SUB78_3WORD_DOUBLE(fq0,fq1,fq2, flo0,flo1,flo2, flo0,flo1,flo2);
			ADD78_3WORD_DOUBLE(flo0,flo1,flo2, fhi0,fhi1,fhi2, fx0,fx1,fx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(fhi0,fhi1,fhi2, flo0,flo1,flo2, fx0,fx1,fx2);
		}
		NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

		if(CMPLT78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2))
		{
			SUB78_3WORD_DOUBLE(gq0,gq1,gq2, glo0,glo1,glo2, glo0,glo1,glo2);
			ADD78_3WORD_DOUBLE(glo0,glo1,glo2, ghi0,ghi1,ghi2, gx0,gx1,gx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(ghi0,ghi1,ghi2, glo0,glo1,glo2, gx0,gx1,gx2);
		}
		NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

		if(CMPLT78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2))
		{
			SUB78_3WORD_DOUBLE(hq0,hq1,hq2, hlo0,hlo1,hlo2, hlo0,hlo1,hlo2);
			ADD78_3WORD_DOUBLE(hlo0,hlo1,hlo2, hhi0,hhi1,hhi2, hx0,hx1,hx2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(hhi0,hhi1,hhi2, hlo0,hlo1,hlo2, hx0,hx1,hx2);
		}
		NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

		if(CMPLT78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2))
		{
			SUB78_3WORD_DOUBLE(iq0,iq1,iq2, ilo0,ilo1,ilo2, ilo0,ilo1,ilo2);
			ADD78_3WORD_DOUBLE(ilo0,ilo1,ilo2, ihi0,ihi1,ihi2, ix0,ix1,ix2);
		}
		else
		{
			SUB78_3WORD_DOUBLE(ihi0,ihi1,ihi2, ilo0,ilo1,ilo2, ix0,ix1,ix2);
		}
		NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);

		if((pshift >> j) & (uint64)1)
		{
			/* Convert fx to uint96 for purpose of this particular comparison: */
			CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
			CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
			CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
			CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x0, qhalf0))
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
				SUB78_3WORD_DOUBLE(fx0,fx1,fx2, fq0,fq1,fq2, fx0,fx1,fx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(fx0,fx1,fx2, fx0,fx1,fx2, fx0,fx1,fx2);
			}
			NORMALIZE78_3WORD_DOUBLE(fx0,fx1,fx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x1, qhalf1))
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
				SUB78_3WORD_DOUBLE(gx0,gx1,gx2, gq0,gq1,gq2, gx0,gx1,gx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(gx0,gx1,gx2, gx0,gx1,gx2, gx0,gx1,gx2);
			}
			NORMALIZE78_3WORD_DOUBLE(gx0,gx1,gx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x2, qhalf2))
			{
				ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
				SUB78_3WORD_DOUBLE(hx0,hx1,hx2, hq0,hq1,hq2, hx0,hx1,hx2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(hx0,hx1,hx2, hx0,hx1,hx2, hx0,hx1,hx2);
			}
			NORMALIZE78_3WORD_DOUBLE(hx0,hx1,hx2);

			/* Combines overflow-on-add and need-to-subtract-q-from-sum checks */
			if(CMPUGT96(x3, qhalf3))
			{
				ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
				SUB78_3WORD_DOUBLE(ix0,ix1,ix2, iq0,iq1,iq2, ix0,ix1,ix2);
			}
			else
			{
				ADD78_3WORD_DOUBLE(ix0,ix1,ix2, ix0,ix1,ix2, ix0,ix1,ix2);
			}
			NORMALIZE78_3WORD_DOUBLE(ix0,ix1,ix2);
		}

	#endif

	}

	/*...Double and return.	These are specialized for the case
	where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/
	CVT78_3WORD_DOUBLE_UINT96(fx0,fx1,fx2, x0);
	CVT78_3WORD_DOUBLE_UINT96(gx0,gx1,gx2, x1);
	CVT78_3WORD_DOUBLE_UINT96(hx0,hx1,hx2, x2);
	CVT78_3WORD_DOUBLE_UINT96(ix0,ix1,ix2, x3);

	ADD96(x0,x0,x0);	/* In the case of interest, x = (q+1)/2 < 2^95, so x + x cannot overflow. */
	ADD96(x1,x1,x1);
	ADD96(x2,x2,x2);
	ADD96(x3,x3,x3);

	SUB96(x0,q0,x0);
	SUB96(x1,q1,x1);
	SUB96(x2,q2,x2);
	SUB96(x3,q3,x3);

	tmp0 = CMPEQ96(x0, ONE96);	*checksum2 += x0.d0;
	tmp1 = CMPEQ96(x1, ONE96);	*checksum2 += x1.d0;
	tmp2 = CMPEQ96(x2, ONE96);	*checksum2 += x2.d0;
	tmp3 = CMPEQ96(x3, ONE96);	*checksum2 += x3.d0;
	r = tmp0;
	r += tmp1 << 1;
	r += tmp2 << 2;
	r += tmp3 << 3;
#if FAC_DEBUG
	if(dbg)
	{
		printf("xout0 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x0)]);
		printf("xout1 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x1)]);
		printf("xout2 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x2)]);
		printf("xout3 = %s\n", &char_buf[CONVERT_UINT96_BASE10_CHAR(char_buf, x3)]);
		printf("!");
	}
#endif

	*retval = r;
}

