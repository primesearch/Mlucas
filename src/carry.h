/*******************************************************************************
*                                                                              *
*   (C) 1997-2019 by Ernst W. Mayer.                                           *
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
#ifndef carry_h_included
#define carry_h_included

#include "Mdata.h"	// For e.g. decls of externs FFT_MUL_BASE

/*** Floating-point-based carry macros: ***/

#if EWM_DEBUG
	#define check_nint(temp, x)\
	{\
		if(temp != NINT(x))\
		{\
			printf("WARNING: rnd(x) != nint(x); %20.5f %20.5f %20.5f\n", x, temp, NINT(x));\
		}\
	}
#else
	#define check_nint(temp, x)	/* */
#endif

/*************************************************************/
/*************** GENERIC-FFT-MUL CARRY MACROS ****************/
/*************************************************************/
#if 0
	/*
	In the generic-FFT-mul (unmodded) carry scheme, real & imaginary parts
	are carried separately and there are no DWT weights:
	*/						// Idx-arg is debug-use only:vvv
	#define genfftmul_carry_norm_pow2_errcheck(x,y,cx,cy,idx)\
	{\
		/* Multiply the current transform output by any scale factor: */\
			rt = x*scale + cx;\
			it = y*scale + cy;\
		/* Real part: */\
			temp = DNINT(rt);\
			frac = fabs(rt - temp);\
			if(frac > maxerr) {\
	printf("Block %2u, full = %d, scale = %15.13f, j = %u: x = %20.5f, cx = %5.0f, maxerr = %7.5f\n",idx,full_pass,scale,j,x,cx,maxerr);\
				maxerr=frac;\
			}\
			cx   = DNINT(temp*FFT_MUL_BASE_INV);\
			x    = temp - cx*FFT_MUL_BASE;\
		/* Imag part: */\
			temp = DNINT(it);\
			frac = fabs(it - temp);\
			if(frac > maxerr) {\
	printf("Block %2u, full = %d, scale = %15.13f, j = %u: y = %20.5f, cy = %5.0f, maxerr = %7.5f\n",idx,full_pass,scale,j,y,cy,maxerr);\
				maxerr=frac;\
			}\
			cy   = DNINT(temp*FFT_MUL_BASE_INV);\
			y    = temp - cy*FFT_MUL_BASE;\
	if(j == jhi-2)printf("j = %4u ,jpad = %4u: Block %2u x,y = %10.5f, %10.5f, carryouts = %10.5f, %10.5f\n",j,j1, idx, x,y, cx,cy);\
	/*if(full_pass && (x != 0. || y != 0.))printf("Block %2u, j = %4u: u += %6d*bpow; v += %6d*bpow; bpow *= b;\n",idx,j,(int)x,(int)y);*/\
	}
#else
	/* V2 computes x-y (= a*u - b*v) on the fly and stores in the normalized x-output, setting y = 0:
	*/						// Idx-arg is debug-use only:vvv
	#define genfftmul_carry_norm_pow2_errcheck(x,y,cx,cy,idx)\
	{\
		/* Multiply the current transform output by any scale factor: */\
			rt = (x - y)*scale + cx;\
			temp = DNINT(rt);\
			frac = fabs(rt - temp);\
			if(frac > maxerr) {\
	printf("Block %2u, full = %d, scale = %15.13f, j = %u: x = %20.5f, cx = %5.0f, maxerr = %7.5f\n",idx,full_pass,scale,j,x,cx,maxerr);\
				maxerr=frac;\
			}\
			cx   = DNINT(temp*FFT_MUL_BASE_INV);\
			x    = temp - cx*FFT_MUL_BASE;\
			y    = 0.0;\
	}
#endif

/*************************************************************/
/**************** FERMAT-MOD CARRY MACROS ********************/
/*************************************************************/
/*
In the Fermat-mod negacyclic-DWT carry scheme, real & imaginary parts
are carried separately due to the right-angle transform:
*/
#define fermat_carry_norm_pow2_errcheck(x,y,cx,cy,idx_offset,NRTM1,NRT_BITS,prp_mult)\
{\
	/* Multiply the current transform output by any scale factor: */\
		x *= scale;\
		y *= scale;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		temp=rn0[k1].re;		wt_im=rn0[k1].im;\
		rt  =rn1[k2].re;		it   =rn1[k2].im;\
		wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im;\
		it = y*wt_re - x*wt_im;\
		temp = DNINT(rt);\
	/*** Feb 2018: Apply any PRP_BASE-multiplier to this rounded result, but only AFTER computing frac-part on next line. Further,
	must delay adding the carryins until the same spot, since these, being from the already-completed previous term in the carry
	chain, already have PRP_BASE-multiplier included in their computation - i.e. need to avoid double-mpying carries by PRP_BASE: ***/\
		frac = fabs(rt - temp);\
		temp = temp*prp_mult + cx;\
	if(frac > maxerr) maxerr=frac;\
		cx = DNINT(temp*baseinv[0]);\
		rt = temp-cx*base[0];\
		temp = DNINT(it);\
		frac = fabs(it - temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy = DNINT(temp*baseinv[0]);\
		it = temp-cy*base[0];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
}

/*
Non-power-of-2 runlengths combine the acyclic sincos weights
with the Mersenne-mod-style IBDWT roots-of-2 weights:
*/
#define fermat_carry_norm_errcheck(x,y,cx,cy,ii,bjmodn,idx_offset,NRTM1,NRT_BITS,prp_mult)\
{\
	/* For Fermat-mod case, combine inverse weight (same for real and imaginary */\
	/* parts of the output) with inverse-FFT scale factor: */\
		wt    =       wt0[ii];\
		wtinv = scale*wt1[ii];\
		ii += SW_DIV_N - nwt;\
		ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		temp=rn0[k1].re;		wt_im=rn0[k1].im;\
		rt  =rn1[k2].re;		it   =rn1[k2].im;\
		wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
	i = (bjmodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */\
		x *= wtinv;\
		y *= wtinv;\
	/* Inverse Negacyclic weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im;\
		it = y*wt_re - x*wt_im;\
		temp = DNINT(rt);\
	/*** Feb 2018: Apply any PRP_BASE-multiplier to this rounded result, but only AFTER computing frac-part on next line. Further,
	must delay adding the carryins until the same spot, since these, being from the already-completed previous term in the carry
	chain, already have PRP_BASE-multiplier included in their computation - i.e. need to avoid double-mpying carries by PRP_BASE: ***/\
		frac = fabs(rt - temp);\
		temp = temp*prp_mult + cx;\
	if(frac > maxerr) maxerr=frac;\
	bjmodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( ((int)i-1) & n);		/*       add 0 if a bigword,   N if a smallword */\
		cx = DNINT(temp*baseinv[i]);\
		rt = temp-cx*base[i];\
		temp = DNINT(it);\
		frac = fabs(it - temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy = DNINT(temp*baseinv[i]);\
		it = temp-cy*base[i];\
	/* Forward Negacyclic weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
	/* Forward IBDWT weight: */\
		x *= wt;\
		y *= wt;\
}

	/*
	Non-power-of-2 runlengths combine the acyclic sincos weights
	with the Mersenne-mod-style IBDWT roots-of-2 weights. Compare this sleeker version - which takes advantage of
	the short-cyclic nature of the Fermat-mod IBDWT weights and bases - with the original non-B one in carry.h:
	*/
	#define fermat_carry_norm_errcheckB(x,y,cx,cy,i,idx_offset,NRTM1,NRT_BITS,prp_mult)\
	{\
		/* For Fermat-mod case, combine inverse weight (same for real and imaginary */\
		/* parts of the output) with inverse-FFT scale factor: */\
			wt    = wt_arr   [i];\
			wtinv = wtinv_arr[i];\
			x *= wtinv;\
			y *= wtinv;\
		/* Get the needed Nth root of -1: */\
			l = ((j + idx_offset) >> 1);\
			k1=(l & NRTM1);\
			k2=(l >> NRT_BITS);\
			temp=rn0[k1].re;		wt_im=rn0[k1].im;\
			rt  =rn1[k2].re;		it   =rn1[k2].im;\
			wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
		/* Inverse root is (wt_re, -wt_im): */\
			rt = x*wt_re + y*wt_im;\
			it = y*wt_re - x*wt_im;\
			bs    = bs_arr   [i];\
			bsinv = bsinv_arr[i];\
			temp = DNINT(rt);\
			frac = fabs(rt - temp);\
			temp = temp*prp_mult + cx;\
		if(frac > maxerr) maxerr=frac;\
			cx = DNINT(temp*bsinv);\
			rt = temp-cx*bs;\
			temp = DNINT(it);\
			frac = fabs(it - temp);\
			temp = temp*prp_mult + cy;\
		if(frac > maxerr) maxerr=frac;\
			cy = DNINT(temp*bsinv);\
			it = temp-cy*bs;\
		/* Forward root is (wt_re, +wt_im): */\
			x = rt*wt_re - it*wt_im;\
			y = rt*wt_im + it*wt_re;\
		/* Forward IBDWT weight: */\
			x *= wt;\
			y *= wt;\
	}

/*************************************************************/
/**************** MERSENNE-MOD CARRY MACROS ******************/
/*************************************************************/
#define cmplx_carry_norm_pow2_errcheck0(x,y,cy,bjmodn,set,prp_mult)\
{\
		/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
		wtA = wt1[col+(set)];\
		wtB = wt1[co2-(set)];\
		wtC = wt1[co3-(set)];\
		wt_re = wtl*wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_pow2_errcheck */\
		wi_re = wtn*wtB;\
/* If it's the 0-word the usual i-computation gives I == 0 here, but want I == 1; force-bigword-ness for the
0-word only by XORing the i-result with ((j+set)==0) = 1. In all other cases (j+set) != 0, thus XOR = no-op: */\
		i  = ((uint32)(sw - bjmodn) >> 31) ^ ((j+set)==0);/* Don't want this for set 0. */\
		m  = ((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2 = 1 + ((uint32)(bjmodn - sinwt) >> 31);\
/*printf("Set %u, Re: bjmodn,sinwt  , m,m2 = %u,%u, %u,%u\n",set,bjmodn,sinwt  , m,m2);*/\
		wt_re *= one_half[m ];\
		wi_re *= one_half[m2];\
		x *= wi_re;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt_re;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt_im = wtlp1*wtA;\
		wi_im = wtnm1*wtC;\
		i  = ((uint32)(sw - bjmodn) >> 31);\
		m  = ((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2 = 1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
/*printf("Set %u, Im: bjmodn,sinwtm1, m,m2 = %u,%u, %u,%u\n",set,bjmodn,sinwtm1, m,m2);*/\
		wt_im *= one_half[m ];\
		wi_im *= one_half[m2];\
		y *= wi_im;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt_im;\
	bjmodn = (bjmodn + bw) & nm1;\
/*if(j >= 0 && j < 8)printf("A: Set %u: wt_re,inv = %12.10f, %15.10e, wt_im,inv = %12.10f, %15.10e\n",set,wt_re,wi_re,wt_im,wi_im);*/\
}

// Lower-accuracy cheap-carry macro which used chained-multipliers to compute current block's weights
// from those of the previous one. *Requires* cmplx_carry_norm_pow2_errcheck0 used on first block, as in HIACC case:
#define cmplx_carry_fast_pow2_errcheck(x,y,cy,bjmodn,set,prp_mult)\
{\
		i = (wt_re >= inv_mult[1]);\
		wi_re *= inv_mult[i];\
		wt_re *= wts_mult[i];\
		x *= wi_re;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt_re;\
	bjmodn = (bjmodn + bw) & nm1;\
		i = (wt_im >= inv_mult[1]);\
		wi_im *= inv_mult[i];\
		wt_im *= wts_mult[i];\
		y *= wi_im;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt_im;\
	bjmodn = (bjmodn + bw) & nm1;\
/*if(j >= 0 && j < 8)printf("B: Set %u: wt_re,inv = %12.10f, %15.10e, wt_im,inv = %12.10f, %15.10e\n",set,wt_re,wi_re,wt_im,wi_im);*/\
}

#define cmplx_carry_norm_pow2_errcheck(x,y,cy,bjmodn,set,prp_mult)\
{\
		/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
		wtA = wt1[col+(set)];\
		wtB = wt1[co2-(set)];\
		wtC = wt1[co3-(set)];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x *= wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y *= wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
}

// This shorter ...wtsinit() macro is actually for use by the SIMD analogs of this fast-carry macro.
// Need a special version in avx mode since there n_minus_sil,n_minus_silp1,sinwt,sinwtm1 are typed
// as (struct uint32x4 *) pointers rather than uint32 int-data:
#ifdef USE_AVX

	#define cmplx_carry_fast_pow2_wtsinit(__n_minus_sil,__n_minus_silp1,__sinwt,__sinwtm1, __wt_re,__wi_re,__wt_im,__wi_im, __bjmodn,set)\
	{\
		double __wtA,__wtB,__wtC;\
		uint32 __m,__m2;\
			/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
			__wtA = wt1[col+(set)];\
			__wtB = wt1[co2-(set)];\
			__wtC = wt1[co3-(set)];\
			__wt_re = wtl*__wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_pow2_errcheck */\
			__wi_re = wtn*__wtB;\
			__m  = ((uint32)(__n_minus_sil-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - __sinwt) >> 31);\
/*printf("Set %u, Re: bjmodn,sinwt  , m,m2 = %u,%u, %u,%u\n",(set),__bjmodn,__sinwt  ,__m,__m2);*/\
			__wt_re *= one_half[__m ];\
			__wi_re *= one_half[__m2];\
		__bjmodn = (__bjmodn + bw) & nm1;\
			__wt_im = wtlp1*__wtA;\
			__wi_im = wtnm1*__wtC;\
			__m  = ((uint32)(__n_minus_silp1-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - __sinwtm1) >> 31);\
/*printf("Set %u, Im: bjmodn,sinwtm1, m,m2 = %u,%u, %u,%u\n",(set),__bjmodn,__sinwtm1,__m,__m2);*/\
			__wt_im *= one_half[__m ];\
			__wi_im *= one_half[__m2];\
		__bjmodn = (__bjmodn + bw) & nm1;\
	}

#elif defined(USE_SSE2)

	#define cmplx_carry_fast_pow2_wtsinit(__wt_re,__wi_re,__wt_im,__wi_im,__bjmodn,set)\
	{\
		double __wtA,__wtB,__wtC;\
		uint32 __m,__m2;\
			/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
			__wtA = wt1[col+(set)];\
			__wtB = wt1[co2-(set)];\
			__wtC = wt1[co3-(set)];\
			__wt_re = wtl*__wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_pow2_errcheck */\
			__wi_re = wtn*__wtB;\
			__m  = ((uint32)(n_minus_sil-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - sinwt) >> 31);\
			__wt_re *= one_half[__m ];\
			__wi_re *= one_half[__m2];\
		__bjmodn = (__bjmodn + bw) & nm1;\
			__wt_im = wtlp1*__wtA;\
			__wi_im = wtnm1*__wtC;\
			__m  = ((uint32)(n_minus_silp1-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - sinwtm1) >> 31);\
			__wt_im *= one_half[__m ];\
			__wi_im *= one_half[__m2];\
		__bjmodn = (__bjmodn + bw) & nm1;\
	}

#endif	// USE_AVX ?

/*
Non-power-of-2 runlengths - here can't do a fast AND-based += bw (mod n),
so replace with a branchless -= sw (+= n if result < 0):
*/
#define cmplx_carry_norm_errcheck0(x,y,cy,bjmodn,set,prp_mult)\
{\
		/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
		wtA = wt1[col+(set)];\
		wtB = wt1[co2-(set)];\
		wtC = wt1[co3-(set)];\
		wt_re = wtl*wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_errcheck */\
		wi_re = wtn*wtB;\
/* If it's the 0-word the usual i-computation gives I == 0 here, but want I == 1; force-bigword-ness for the
0-word only by XORing the i-result with ((j+set)==0) = 1. In all other cases (j+set) != 0, thus XOR = no-op: */\
		i  = ((uint32)(sw - bjmodn) >> 31) ^ ((j+set)==0);/* Don't want this for set 0. */\
		m  = ((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2 = 1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt_re *= one_half[m ];\
		wi_re *= one_half[m2];\
		x *= wi_re;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt_re;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt_im = wtlp1*wtA;\
		wi_im = wtnm1*wtC;\
		i  = ((uint32)(sw - bjmodn) >> 31);\
		m  = ((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2 = 1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt_im *= one_half[m ];\
		wi_im *= one_half[m2];\
		y *= wi_im;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt_im;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

// Lower-accuracy cheap-carry macro which used chained-multipliers to compute current block's weights
// from those of the previous one. *Requires* cmplx_carry_norm_errcheck0 used on first block, as in HIACC case:
#define cmplx_carry_fast_errcheck(x,y,cy,bjmodn,set,prp_mult)\
{\
		i = (wt_re >= inv_mult[1]);\
		wi_re *= inv_mult[i];\
		wt_re *= wts_mult[i];\
		x *= wi_re;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt_re;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		i = (wt_im >= inv_mult[1]);\
		wi_im *= inv_mult[i];\
		wt_im *= wts_mult[i];\
		y *= wi_im;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt_im;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define cmplx_carry_norm_errcheck(x,y,cy,bjmodn,set,prp_mult)\
{\
		/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
		wtA = wt1[col+(set)];\
		wtB = wt1[co2-(set)];\
		wtC = wt1[co3-(set)];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x *= wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y *= wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		temp = temp*prp_mult + cy;\
	if(frac > maxerr) maxerr=frac;\
		i =((uint32)(sw - bjmodn) >> 31);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

// This shorter ...wtsinit() macro is actually for use by the SIMD analogs of this fast-carry macro:
// Need a special version in avx mode since there n_minus_sil,n_minus_silp1,sinwt,sinwtm1 are typed
// as (struct uint32x4 *) pointers rather than uint32 int-data:
#ifdef USE_AVX

	#define cmplx_carry_fast_wtsinit(__n_minus_sil,__n_minus_silp1,__sinwt,__sinwtm1, __wt_re,__wi_re,__wt_im,__wi_im, __bjmodn,set)\
	{\
		double __wtA,__wtB,__wtC;\
		uint32 __m,__m2;\
			/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
			__wtA = wt1[col+(set)];\
			__wtB = wt1[co2-(set)];\
			__wtC = wt1[co3-(set)];\
			__wt_re = wtl*__wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_pow2_errcheck */\
			__wi_re = wtn*__wtB;\
			__m  = ((uint32)(__n_minus_sil-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - __sinwt) >> 31);\
/*printf("Set %u, Re: bjmodn,sinwt  , m,m2 = %u,%u, %u,%u\n",(set),__bjmodn,__sinwt  ,__m,__m2);*/\
			__wt_re *= one_half[__m ];\
			__wi_re *= one_half[__m2];\
		__bjmodn -= sw;\
		__bjmodn += ( (-(int)((uint32)__bjmodn >> 31)) & n);\
			__wt_im = wtlp1*__wtA;\
			__wi_im = wtnm1*__wtC;\
			__m  = ((uint32)(__n_minus_silp1-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - __sinwtm1) >> 31);\
/*printf("Set %u, Im: bjmodn,sinwtm1, m,m2 = %u,%u, %u,%u\n",(set),__bjmodn,__sinwtm1,__m,__m2);*/\
			__wt_im *= one_half[__m ];\
			__wi_im *= one_half[__m2];\
		__bjmodn -= sw;\
		__bjmodn += ( (-(int)((uint32)__bjmodn >> 31)) & n);\
	}

#elif defined(USE_SSE2)

	#define cmplx_carry_fast_wtsinit(__wt_re,__wi_re,__wt_im,__wi_im,__bjmodn,set)\
	{\
		double __wtA,__wtB,__wtC;\
		uint32 __m,__m2;\
			/* Must parenthesize (set) to allow arg to be an expression, since e.g. (co2 - a+b) != (co2 - (a+b)) */\
			__wtA = wt1[col+(set)];\
			__wtB = wt1[co2-(set)];\
			__wtC = wt1[co3-(set)];\
			__wt_re = wtl*__wtA;	/* Use separate re/im wts in the 0-version of this macro for compatibility with cmplx_carry_fast_errcheck */\
			__wi_re = wtn*__wtB;\
			__m  = ((uint32)(n_minus_sil-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - sinwt) >> 31);\
			__wt_re *= one_half[__m ];\
			__wi_re *= one_half[__m2];\
		__bjmodn -= sw;\
		__bjmodn += ( (-(int)((uint32)__bjmodn >> 31)) & n);\
			__wt_im = wtlp1*__wtA;\
			__wi_im = wtnm1*__wtC;\
			__m  = ((uint32)(n_minus_silp1-__bjmodn) >> 31);\
			__m2 = 1 + ((uint32)(__bjmodn - sinwtm1) >> 31);\
			__wt_im *= one_half[__m ];\
			__wi_im *= one_half[__m2];\
		__bjmodn -= sw;\
		__bjmodn += ( (-(int)((uint32)__bjmodn >> 31)) & n);\
	}

#endif	// USE_AVX ?

/******************************************************************************************/
/*** SIMD implementation of the key carry macros is in linked header files: ***************/
/******************************************************************************************/
#ifdef USE_SSE2

	#if defined(COMPILER_TYPE_MSVC)	/* MSVC-style inline ASM: */

		#error MSVC SIMD builds not supported!

	#else	/* GCC-style inline ASM: */

		#include "carry_gcc64.h"

	#endif	/* MSVC or GCC */

#endif	/* USE_SSE2 */

#endif	/* carry_h_included */
