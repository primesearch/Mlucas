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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef carry_h_included
#define carry_h_included

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
/**************** FERMAT-MOD CARRY MACROS ********************/
/*************************************************************/
/*
In the Fermat-mod negacyclic-DWT carry scheme, real & imaginary parts
are carried separately due to the right-angle transform:
*/
#define fermat_carry_norm_pow2_errcheck(x,y,cx,cy,idx_offset,NRTM1,NRT_BITS)\
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
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp = DNINT(rt);\
		frac = fabs(rt - temp);\
		if(frac > maxerr) maxerr=frac;\
		cx = DNINT(temp*baseinv[0]);\
		rt = temp-cx*base[0];\
		temp = DNINT(it);\
		frac = fabs(it - temp);\
		if(frac > maxerr) maxerr=frac;\
		cy = DNINT(temp*baseinv[0]);\
		it = temp-cy*base[0];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
}

#define fermat_carry_norm_pow2_nocheck(x,y,cx,cy,idx_offset,NRTM1,NRT_BITS)\
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
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp = DNINT(rt);\
		cx = DNINT(temp*baseinv[0]);\
		rt = temp-cx*base[0];\
		temp = DNINT(it);\
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
#define fermat_carry_norm_errcheck(x,y,cx,cy,ii,bjmodn,idx_offset,NRTM1,NRT_BITS)\
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
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp = DNINT(rt);\
		frac = fabs(rt - temp);\
		if(frac > maxerr) maxerr=frac;\
	bjmodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( ((int)i-1) & n);		/*       add 0 if a bigword,   N if a smallword */\
		cx = DNINT(temp*baseinv[i]);\
		rt = temp-cx*base[i];\
		temp = DNINT(it);\
		frac = fabs(it - temp);\
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

#define fermat_carry_norm_nocheck(x,y,cx,cy,ii,bjmodn,idx_offset,NRTM1,NRT_BITS)\
{\
	/* For Fermat-mod case, combine inverse weight (same for real and imaginary */\
	/* parts of the output) with inverse-FFT scale factor: */\
		wt    =       wt0[ii];\
		wtinv = scale*wt1[ii];\
		ii += SW_DIV_N - nwt;\
		ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
		x *= wtinv;\
		y *= wtinv;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		temp=rn0[k1].re;		wt_im=rn0[k1].im;\
		rt  =rn1[k2].re;		it   =rn1[k2].im;\
		wt_re =temp*rt-wt_im*it;wt_im =temp*it+wt_im*rt;\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp = DNINT(rt);\
	i = (bjmodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */\
	bjmodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( ((int)i-1) & n);		/*       add 0 if a bigword,   N if a smallword */\
		cx = DNINT(temp*baseinv[i]);\
		rt = temp-cx*base[i];\
		temp = DNINT(it);\
		cy = DNINT(temp*baseinv[i]);\
		it = temp-cy*base[i];\
	/* Forward weight is (wt_re, +wt_im): */\
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
		#define fermat_carry_norm_errcheckB(x,y,cx,cy,i,idx_offset,NRTM1,NRT_BITS)\
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
				rt = x*wt_re + y*wt_im + cx;\
				it = y*wt_re - x*wt_im + cy;\
				bs    = bs_arr   [i];\
				bsinv = bsinv_arr[i];\
				temp = DNINT(rt);\
				frac = fabs(rt - temp);\
			if(frac > maxerr) maxerr=frac;\
				cx = DNINT(temp*bsinv);\
				rt = temp-cx*bs;\
				temp = DNINT(it);\
				frac = fabs(it - temp);\
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

#define cmplx_carry_norm_pow2_errcheck0(x,y,cy,bjmodn)\
{\
		wtA=wt1[col];\
		wtB=wt1[co2];\
		wtC=wt1[co3];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
/*		i =((uint32)(sw - bjmodn) >> 31);Don't want this for set 0. */\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
}

#define cmplx_carry_norm_pow2_errcheck(x,y,cy,bjmodn,set)\
{\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];\
		wtC=wt1[co3-set];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
}

#define cmplx_carry_norm_pow2_nocheck0(x,y,cy,bjmodn)\
{\
		wtA=wt1[col];\
		wtB=wt1[co2];\
		wtC=wt1[co3];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
/*		i =((uint32)(sw - bjmodn) >> 31);Don't want this for set 0. */\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
}

#define cmplx_carry_norm_pow2_nocheck(x,y,cy,bjmodn,set)\
{\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];\
		wtC=wt1[co3-set];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn = (bjmodn + bw) & nm1;\
}


/*
Non-power-of-2 runlengths:
*/
#define cmplx_carry_norm_errcheck0(x,y,cy,bjmodn)\
{\
		wtA=wt1[col];\
		wtB=wt1[co2];\
		wtC=wt1[co3];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
/*		i =((uint32)(sw - bjmodn) >> 31);Don't want this for set 0. */\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
/*if(j<5) printf("x_in[%d] = %15.5f, wt,inv = %15.14f,%15.14f, base,inv = %7d, %20.15e, temp,frac,cy_out = %15.14f,%15.14f,%15.14f, a[%d]* = %15.5f\n",j,x,wt,wtinv/scale,(int)base[i],baseinv[i],temp,frac,cy,j,(temp-cy*base[i]));*/\
		x = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
/*if(j<5) printf("y_in[%d] = %15.5f, wt,inv = %15.14f,%15.14f, base,inv = %7d, %20.15e, temp,frac,cy_out = %15.14f,%15.14f,%15.14f, a[%d]* = %15.5f\n",j+1,y,wt,wtinv/scale,(int)base[i],baseinv[i],temp,frac,cy,j+1,(temp-cy*base[i]));*/\
		y = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define cmplx_carry_norm_errcheck(x,y,cy,bjmodn,set)\
{\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];\
		wtC=wt1[co3-set];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

/* A version with a slightly different sequencing, attempting to achieve better pipelining: */
#define cmplx_carry_norm_errcheckB(x,y,cy,bjmodn,set)\
{\
		x *= wtn;\
		y *= wtnm1;\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];	x *= wtB;\
		wtC=wt1[co3-set];	y *= wtC;\
		wt   =wtl*wtA;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt =     wt*one_half[m];\
		x  = cy + x*one_half[m2];\
		temp = DNINT(x);				check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt =     wt*one_half[m];\
		y  = cy + y*one_half[m2];\
		temp = DNINT(y);				check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define cmplx_carry_norm_nocheck0(x,y,cy,bjmodn)\
{\
		wtA=wt1[col];\
		wtB=wt1[co2];\
		wtC=wt1[co3];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
/*		i =((uint32)(sw - bjmodn) >> 31);Don't want this for set 0. */\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
	bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
	bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define cmplx_carry_norm_nocheck(x,y,cy,bjmodn,set)\
{\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];\
		wtC=wt1[co3-set];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
		temp = DNINT(x);				check_nint(temp, x);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
		temp = DNINT(y);				check_nint(temp, y);\
		cy   = DNINT(temp*baseinv[i]);	check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}


/*** Integer-based carry macros (not currently used, as they're slow): ***/

#define Icmplx_carry_norm_errcheck0(x,y,cy,bjmodn)\
{\
const uint64 himask= ((uint64)1)<<63,\
    	     two52 = ((uint64)1)<<52,\
    	     ones  = ((uint64)((int64)-1)),\
		 mmask = ones>>12;\
const uint32 swbits = p/n;        /* Number of bits in a smallword. */\
uint64 ix,ifrac,ifracmax = 0;\
 int64 sign,mant,word,topbit;\
uint32 dexp,bits;\
 int32 shift;\
\
		wtA=wt1[col];\
		wtB=wt1[co2];\
		wtC=wt1[co3];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
/*if(j1==0)printf("x0,cy,wtinv = %20.15e  %20.15e  %20.15e\n",x,cy,wtinv);*/\
		x = cy+ x*wtinv;\
bits  = swbits + i;\
ix    = *((uint64 *)&x);\
sign  = ix & himask;\
mant  = ix & mmask;\
dexp  = (ix-sign)>>52;\
shift = 1074 - dexp;\
/*if(j1==0)printf("0xmant,shift,bits = %20llX  %10d  %10u\n",mant,shift,bits);*/\
if(shift<0)printf("WARN: j1 = %10d  %20.15e gives negative shift count = %10d\n",j1,x,shift);\
if(shift < 52)\
{\
sign>>=63;  /* Signed arithmetic left-shift here, i.e. get -1 if float had sign bit set. */\
ifrac = mant << (63-shift);\
if(ifrac > ifracmax) ifracmax=ifrac;\
mant += ((uint64)1)<<shift;\
mant  = (mant+two52)>>(shift+1);\
/*if(j1==0)printf("A: 0xmant = %20llX\n",mant);*/\
mant -= (mant & sign)<<1;\
/*if(j1==0)printf("B: 0xmant = %20llX\n",mant);*/\
word  = mant & (~(ones << bits));\
/*if(j1==0)printf("C: 0xword = %20llX\n",word);*/\
topbit= word >> (bits - 1);\
/*if(j1==0)printf("D: 0xtbit = %20llX\n",topbit);*/\
word -= topbit << bits;\
/*if(j1==0)printf("E: 0xword = %20llX\n",word);*/\
x     = wt*(double)word;\
cy    = (double)( (mant >> bits) + topbit );\
/*if(j1==0)printf("%20.4f  %20.4f\n",x,cy);*/\
}\
else\
{\
  x = 0; cy = 0;\
}\
		/*\
		temp = (x + RND_A) - RND_B;\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		x = (temp-cy*base[i])*wt;\
		*/\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
/*if(j1==0)printf("y0,cy,wtinv = %20.15e  %20.15e  %20.15e\n",y,cy,wtinv);*/\
		y = cy+ y*wtinv;\
bits  = swbits + i;\
ix    = *((uint64 *)&y);\
sign  = ix & himask;\
mant  = ix & mmask;\
dexp  = (ix-sign)>>52;\
shift = 1074 - dexp;\
if(shift<0)printf("WARN: j1 = %10d  %20.15e gives negative shift count = %10d\n",j1,y,shift);\
if(shift < 52)\
{\
sign>>=63;  /* Signed arithmetic left-shift here, i.e. get -1 if float had sign bit set. */\
ifrac = mant << (63-shift);\
if(ifrac > ifracmax) ifracmax=ifrac;\
mant += ((uint64)1)<<shift;\
mant  = (mant+two52)>>(shift+1);\
mant -= (mant & sign)<<1;\
word  = mant & (~(ones << bits));\
topbit= word >> (bits - 1);\
word -= topbit << bits;\
y     = wt*(double)word;\
cy    = (double)( (mant >> bits) + topbit );\
/*if(j1==0)printf("%20.4f  %20.4f\n",y,cy);*/\
}\
else\
{\
  y = 0; cy = 0;\
}\
		/*\
		temp = (y + RND_A) - RND_B;\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		y = (temp-cy*base[i])*wt;\
		*/\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define Icmplx_carry_norm_errcheck(x,y,cy,bjmodn,set)\
{\
		wtA=wt1[col+set];\
		wtB=wt1[co2-set];\
		wtC=wt1[co3-set];\
		wt   =wtl*wtA;\
		wtinv=wtn*wtB;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_sil-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwt) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		x = cy+ x*wtinv;\
bits  = swbits + i;\
ix    = *((uint64 *)&x);\
sign  = ix & himask;\
mant  = ix & mmask;\
dexp  = (ix-sign)>>52;\
shift = 1074 - dexp;\
if(shift<0)printf("WARN: j1 = %10d  %20.15e gives negative shift count = %10d\n",j1,x,shift);\
if(shift < 52)\
{\
sign>>=63;  /* Signed arithmetic left-shift here, i.e. get -1 if float had sign bit set. */\
ifrac = mant << (63-shift);\
if(ifrac > ifracmax) ifracmax=ifrac;\
mant += ((uint64)1)<<shift;\
mant  = (mant+two52)>>(shift+1);\
mant -= (mant & sign)<<1;\
word  = mant & (~(ones << bits));\
topbit= word >> (bits - 1);\
word -= topbit << bits;\
x     = wt*(double)word;\
cy    = (double)( (mant >> bits) + topbit );\
}\
else\
{\
  x = 0; cy = 0;\
}\
		/*\
		temp = (x + RND_A) - RND_B;\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		x = (temp-cy*base[i])*wt;\
		*/\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		y = cy+ y*wtinv;\
bits  = swbits + i;\
ix    = *((uint64 *)&y);\
sign  = ix & himask;\
mant  = ix & mmask;\
dexp  = (ix-sign)>>52;\
shift = 1074 - dexp;\
if(shift<0)printf("WARN: j1 = %10d  %20.15e gives negative shift count = %10d\n",j1,y,shift);\
if(shift < 52)\
{\
sign>>=63;  /* Signed arithmetic left-shift here, i.e. get -1 if float had sign bit set. */\
ifrac = mant << (63-shift);\
if(ifrac > ifracmax) ifracmax=ifrac;\
mant += ((uint64)1)<<shift;\
mant  = (mant+two52)>>(shift+1);\
mant -= (mant & sign)<<1;\
word  = mant & (~(ones << bits));\
topbit= word >> (bits - 1);\
word -= topbit << bits;\
y     = wt*(double)word;\
cy    = (double)( (mant >> bits) + topbit );\
}\
else\
{\
  y = 0; cy = 0;\
}\
		/*\
		temp = (y + RND_A) - RND_B;\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		y = (temp-cy*base[i])*wt;\
		*/\
	bjmodn -= sw;\
	bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}


	/******************************************************************************************/
	/*************** SSE2 implementation of the key carry macros ******************************/
	/******************************************************************************************/

#ifdef USE_SSE2

	#if defined(COMPILER_TYPE_MSVC)	/* MSVC-style inline ASM: */

		/*
		In the Fermat-mod negacyclic-DWT carry scheme, real & imaginary parts
		are carried separately due to the right-angle transform:
		*/
		#define SSE2_fermat_carry_norm_pow2_errcheck(__data,__cy,__idx_offset,__idx_incr,__NRTM1,__NRT_BITS)\
		{\
			__asm	mov		esi, __idx_offset	/* esi stores [j + idx_offset], idx_offset starts = 0, gets incremented by idx_incr each macro invocation */\
			__asm	mov		ecx, __NRT_BITS	\
			__asm	mov		eax, esi	/* j + idx_offset */	\
			__asm	shr		eax, 1	/* l = ((j + idx_offset) >> 1) */	\
			__asm	mov		ebx, eax	\
			__asm	and		eax, __NRTM1		/* k1 = (l & __NRTM1) */	\
			__asm	shr		ebx, cl	/* k2=(l >> __NRT_BITS) */	\
			__asm	shl		eax, 4	/* 16 bytes for array-of-complex */	\
			__asm	shl		ebx, 4	/* 16 bytes for array-of-complex */	\
			__asm	add		eax, add1	/* rn0[k1] */	\
			__asm	add		ebx, add2	/* rn1[k2] */	\
			__asm	movaps	xmm0,[eax]	/* [c0,s0] */	\
			__asm	movaps	xmm1,[ebx]	/* [x0,y0] */	\
			__asm	mov		eax, esi	\
			__asm	movaps	xmm2,xmm1	/* [x0,y0] copy */	\
			__asm	shufpd	xmm2,xmm2,1	/* [y0,x0] (swap re <--> im) */	\
			__asm	mulpd	xmm1,xmm0	/* [c0.x0,s0.y0] */	\
			__asm	mulpd	xmm2,xmm0	/* [c0.y0,s0.x0] 1,2 used */	\
			/* Get next root for interleaving with the first: */	\
			__asm	add		eax, 2	\
			__asm	shr		eax, 1	/* l = ((j + idx_offset) >> 1) */	\
			__asm	mov		ebx, eax	\
			__asm	and		eax, __NRTM1		/* k1 = (l & __NRTM1) */	\
			__asm	shr		ebx, cl	/* k2=(l >> __NRT_BITS) */	\
			__asm	shl		eax, 4	/* 16 bytes for array-of-complex */	\
			__asm	shl		ebx, 4	/* 16 bytes for array-of-complex */	\
			__asm	add		eax, add1	/* rn0[k1] */	\
			__asm	add		ebx, add2	/* rn1[k2] */	\
			__asm	movaps	xmm0,[eax]	/* [c1,s1] */	\
			__asm	movaps	xmm3,[ebx]	/* [x1,y1] 0-3 used*/	\
			__asm	mov		eax, esi	\
			__asm	movaps	xmm4,xmm3	/* [x1,y1] copy */	\
			__asm	shufpd	xmm4,xmm4,1	/* [y1,x1] (swap re <--> im) */	\
			__asm	mulpd	xmm3,xmm0	/* [c1.x1,s1.y1] */	\
			__asm	mulpd	xmm4,xmm0	/* [c1.y1,s1.x1] 1-4 used */	\
			__asm	movaps	xmm0,xmm1	/* xmm0 <- copy [c0.x0,s0.y0] */	\
			__asm	unpcklpd	xmm0,xmm3	/* [c0.x0,c1.x1] */	\
			__asm	unpckhpd	xmm1,xmm3	/* [s0.y0,s1.y1], 0-2,4 used */	\
			__asm	subpd	xmm0,xmm1	/* XMM0 = [wt_r0,wt_r1] 0,2,4 used */	\
			__asm	movaps	xmm1,xmm2	/* xmm1 <- copy [c0.y0,s0.x0] 0-2,4 used */	\
			__asm	unpcklpd	xmm1,xmm4	/* [c0.y0,c1.y1] */	\
			__asm	unpckhpd	xmm2,xmm4	/* [s0.x0,s1.x1] */	\
			__asm	addpd	xmm1,xmm2	/* XMM1 = [wt_i0,wt_i1] 0-1 used */	\
			/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */	\
			__asm	mov		ecx, half_arr	/* No longer need __NRT_BITS, so reuse ecx */	\
			/* Multiply the complex transform output [x,y] = [re,im] by any scale factor: [x,y] *= scale: */	\
			__asm	mov		edx, __data	\
			__asm	movaps	xmm4,[edx     ] /* x = [a.re,b.re] */	\
			__asm	movaps	xmm2,[edx+0x10]	/* y = [a.im,b.im] */	\
			__asm	movaps	xmm5,[ecx+0x20]	/* [scale,scale] */	\
			__asm	mulpd	xmm4,xmm5	\
			__asm	mulpd	xmm2,xmm5	\
			__asm	movaps	xmm5,xmm4	/* x copy */	\
			__asm	movaps	xmm3,xmm2	/* y copy */	\
			/* Inverse weight is (wt_re, -wt_im): */	\
			__asm	mulpd	xmm4,xmm0	/* [x     ]*wt_re */	\
			__asm	mulpd	xmm3,xmm1	/* [y copy]*wt_im */	\
			__asm	mulpd	xmm2,xmm0	/* [y     ]*wt_re */	\
			__asm	mulpd	xmm5,xmm1	/* [x copy]*wt_im */	\
			__asm	addpd	xmm4,xmm3	/* [a.re,b.re] = x*wt_re + y*wt_im */	\
			__asm	subpd	xmm2,xmm5	/* [a.im,b.im] = y*wt_re - x*wt_im */	\
			__asm	mov		ebx, __cy	\
			__asm	movaps	xmm5,xmm4	/* [a.re,b.re] copy */	\
			__asm	shufpd	xmm4,xmm2,0	/* XMM4 = x = [a.re,a.im] */	\
			__asm	shufpd	xmm5,xmm2,3	/* XMM5 = y = [b.re,b.im] 0,1,4,5 uaed */	\
			/* normalize a-pair, compute carryout, compute ROE: */	\
			__asm	addpd	xmm4,[ebx]	/* [a.re,a.im] + [cx,cy] */	\
			__asm	movaps	xmm6,[ecx-0x20]		/* XMM6 = maxerr */	\
			__asm	movaps	xmm7,[ecx-0x10]		/* XMM7 = rnd_const */	\
			__asm	movaps	xmm2,xmm4			/* copy x */	\
			__asm	addpd	xmm4,xmm7	\
			__asm	subpd	xmm4,xmm7			/* temp = DNINT(x) */	\
			__asm	mov		eax, sign_mask	\
			__asm	subpd	xmm2,xmm4			/* frac = [x - temp] */	\
			__asm	andpd	xmm2,[eax]			/* frac = fabs(frac) */	\
			__asm	maxpd	xmm2,xmm6		/* if(frac > maxerr) maxerr=frac */	\
			__asm	movaps	xmm6,xmm2		/* Note serialization here! */	\
			__asm	movaps	xmm2,xmm4			/* cpy temp */	\
			__asm	mulpd	xmm2,[ecx+0x10]		/* temp*baseinv[0] */	\
			__asm	addpd	xmm2,xmm7	\
			__asm	subpd	xmm2,xmm7			/* [cx,cy] = DNINT(temp*baseinv[0]) */	\
			__asm	movaps	xmm3,xmm2			/* cpy [cx,cy] */	\
			__asm	mulpd	xmm3,[ecx     ]		/* [cx,cy]*base[0] */	\
			__asm	subpd	xmm4,xmm3			/* XMM4 = [a.re,a.im] = temp-[cx,cy]*base[0] */	\
			/* Now do b-pair: [b.re,b.im] in xmm5, carry in xmm2, xmm3 free, wt_[re,im] in xmmA,B, xmm6 free, rnd_const in xmm7: */	\
			__asm	addpd	xmm5,xmm2	/* [b.re,b.im] + [cx,cy] */	\
			__asm	movaps	xmm2,xmm5			/* copy y */	\
			__asm	addpd	xmm5,xmm7	\
			__asm	subpd	xmm5,xmm7			/* temp = DNINT(y) */	\
			__asm	subpd	xmm2,xmm5			/* frac = [y - temp] */	\
			__asm	andpd	xmm2,[eax]			/* frac = fabs(frac) */	\
			__asm	maxpd	xmm2,xmm6		/* if(frac > maxerr) maxerr=frac */	\
			__asm	movaps	xmm6,xmm2		/* Note serialization here! */	\
			__asm	movaps	xmm2,xmm5			/* cpy temp */	\
			__asm	mulpd	xmm2,[ecx+0x10]		/* temp*baseinv[0] */	\
			__asm	addpd	xmm2,xmm7	\
			__asm	subpd	xmm2,xmm7			/* [cx,cy] = DNINT(temp*baseinv[0]) */	\
			__asm	movaps	xmm3,xmm2			/* cpy [cx,cy] */	\
			__asm	mulpd	xmm3,[ecx     ]		/* [cx,cy]*base[0] */	\
			__asm	subpd	xmm5,xmm3			/* XMM5 = [b.re,b.im] = temp-[cx,cy]*base[0] */	\
			__asm	movaps	[ebx],xmm2			/* store cy_out */	\
			__asm	movaps	xmm2,xmm4	/* [a.re,a.im] copy */	\
			__asm	shufpd	xmm4,xmm5,0	/* x = [a.re,b.re] */	\
			__asm	shufpd	xmm2,xmm5,3	/* y = [a.im,b.im] */	\
			__asm	movaps	xmm5,xmm4	/* x copy */	\
			__asm	movaps	xmm3,xmm2	/* y copy */	\
			__asm	movaps	[ecx-0x20],xmm6		/* Store maxerr */	\
			/* Forward weight is (wt_re, +wt_im): */	\
			__asm	mulpd	xmm4,xmm0	/* [x     ]*wt_re */	\
			__asm	mulpd	xmm3,xmm1	/* [y copy]*wt_im */	\
			__asm	mulpd	xmm2,xmm0	/* [y     ]*wt_re */	\
			__asm	mulpd	xmm5,xmm1	/* [x copy]*wt_im */	\
			__asm	subpd	xmm4,xmm3	/* rt = x*wt_re - y*wt_im */	\
			__asm	addpd	xmm5,xmm2	/* it = x*wt_im + y*wt_re */	\
			__asm	movaps	[edx     ],xmm4 /* store rt = ~[a.re,b.re] */	\
			__asm	movaps	[edx+0x10],xmm5	/* store it = ~[a.im,b.im] */	\
			/* Prepare for next pair of complex data: */	\
			__asm	add		esi, __idx_incr	/* idx_offset += idx_incr */	\
			__asm	mov		__idx_offset, esi
		}

		#define SSE2_fermat_carry_norm_errcheck(__data,__cy,__idx_offset,__idx_incr,__odd_radix,__offset0,__offset1,__NRTM1,__NRT_BITS)\
		{\
			__asm	mov		esi, __idx_offset	/* esi stores [j + idx_offset], idx_offset starts = 0, gets incremented by idx_incr each macro invocation */\
			__asm	mov		edi, __odd_radix	/* [1,2,3]*odd_radix are the index offsets to the wtinv, base, and base_inv values, respectively. */\
			__asm	mov		ecx, __NRT_BITS	\
			__asm	mov		eax, esi	/* j + idx_offset */	\
			__asm	shr		eax, 1	/* l = ((j + idx_offset) >> 1) */	\
			__asm	mov		ebx, eax	\
			__asm	and		eax, __NRTM1		/* k1 = (l & __NRTM1) */	\
			__asm	shr		ebx, cl	/* k2=(l >> __NRT_BITS) */	\
			__asm	shl		eax, 4	/* 16 bytes for array-of-complex */	\
			__asm	shl		ebx, 4	/* 16 bytes for array-of-complex */	\
			__asm	shl		edi, 4	/* 16 bytes for array-of-complex */	\
			__asm	add		eax, add1	/* rn0[k1] */	\
			__asm	add		ebx, add2	/* rn1[k2] */	\
			__asm	movaps	xmm0,[eax]	/* [c0,s0] */	\
			__asm	movaps	xmm1,[ebx]	/* [x0,y0] */	\
			__asm	mov		eax, esi	\
			__asm	movaps	xmm2,xmm1	/* [x0,y0] copy */	\
			__asm	shufpd	xmm2,xmm2,1	/* [y0,x0] (swap re <--> im) */	\
			__asm	mulpd	xmm1,xmm0	/* [c0.x0,s0.y0] */	\
			__asm	mulpd	xmm2,xmm0	/* [c0.y0,s0.x0] 1,2 used */	\
			/* Get next root for interleaving with the first: */	\
			__asm	add		eax, 2	\
			__asm	shr		eax, 1	/* l = ((j + idx_offset) >> 1) */	\
			__asm	mov		ebx, eax	\
			__asm	and		eax, __NRTM1		/* k1 = (l & __NRTM1) */	\
			__asm	shr		ebx, cl	/* k2=(l >> __NRT_BITS) */	\
			__asm	shl		eax, 4	/* 16 bytes for array-of-complex */	\
			__asm	shl		ebx, 4	/* 16 bytes for array-of-complex */	\
			__asm	add		eax, add1	/* rn0[k1] */	\
			__asm	add		ebx, add2	/* rn1[k2] */	\
			__asm	movaps	xmm0,[eax]	/* [c1,s1] */	\
			__asm	movaps	xmm3,[ebx]	/* [x1,y1] 0-3 used*/	\
			__asm	mov		eax, esi	\
			__asm	movaps	xmm4,xmm3	/* [x1,y1] copy */	\
			__asm	shufpd	xmm4,xmm4,1	/* [y1,x1] (swap re <--> im) */	\
			__asm	mulpd	xmm3,xmm0	/* [c1.x1,s1.y1] */	\
			__asm	mulpd	xmm4,xmm0	/* [c1.y1,s1.x1] 1-4 used */	\
			__asm	movaps	xmm0,xmm1	/* xmm0 <- copy [c0.x0,s0.y0] */	\
			__asm	unpcklpd	xmm0,xmm3	/* [c0.x0,c1.x1] */	\
			__asm	unpckhpd	xmm1,xmm3	/* [s0.y0,s1.y1], 0-2,4 used */	\
			__asm	subpd	xmm0,xmm1	/* XMM0 = [wt_r0,wt_r1] 0,2,4 used */	\
			__asm	movaps	xmm1,xmm2	/* xmm1 <- copy [c0.y0,s0.x0] 0-2,4 used */	\
			__asm	unpcklpd	xmm1,xmm4	/* [c0.y0,c1.y1] */	\
			__asm	unpckhpd	xmm2,xmm4	/* [s0.x0,s1.x1] */	\
			__asm	addpd	xmm1,xmm2	/* XMM1 = [wt_i0,wt_i1] 0-1 used */	\
			/* half_arr[0,1,2,3] = [base*2, baseinv*2,wt_re*2,wt_im*2] */	\
			__asm	mov		ecx, half_arr	/* No longer need __NRT_BITS, so reuse ecx */	\
			/* Multiply the complex transform output [x,y] = [re,im] by the inverse IBDWT weight, which includes the scale factor: [x,y] *= wtinv: */	\
			__asm	mov		edx, __data	\
			__asm	movaps	xmm4,[edx     ] /* x = [a.re,b.re] */	\
			__asm	movaps	xmm2,[edx+0x10]	/* y = [a.im,b.im] */	\
			__asm	add		ecx, __offset0	\
			__asm	movaps	xmm5,[ecx+edi]	/* [wtinv0,wtinv1] */	\
			__asm	sub		ecx, __offset0	\
			__asm	mulpd	xmm4,xmm5	\
			__asm	mulpd	xmm2,xmm5	\
			__asm	movaps	xmm5,xmm4	/* x copy */	\
			__asm	movaps	xmm3,xmm2	/* y copy */	\
			/* Inverse weight is (wt_re, -wt_im): */	\
			__asm	mulpd	xmm4,xmm0	/* [x     ]*wt_re */	\
			__asm	mulpd	xmm3,xmm1	/* [y copy]*wt_im */	\
			__asm	mulpd	xmm2,xmm0	/* [y     ]*wt_re */	\
			__asm	mulpd	xmm5,xmm1	/* [x copy]*wt_im */	\
			__asm	addpd	xmm4,xmm3	/* [a.re,b.re] = x*wt_re + y*wt_im */	\
			__asm	subpd	xmm2,xmm5	/* [a.im,b.im] = y*wt_re - x*wt_im */	\
			__asm	mov		ebx, __cy	\
			__asm	movaps	xmm5,xmm4	/* [a.re,b.re] copy */	\
			__asm	shufpd	xmm4,xmm2,0	/* XMM4 = x = [a.re,a.im] */	\
			__asm	shufpd	xmm5,xmm2,3	/* XMM5 = y = [b.re,b.im] 0,1,4,5 uaed */	\
			/* normalize a-pair, compute carryout, compute ROE: */	\
			__asm	addpd	xmm4,[ebx]	/* [a.re,a.im] + [cx,cy] */	\
			__asm	movaps	xmm6,[ecx-0x20]		/* XMM6 = maxerr */	\
			__asm	movaps	xmm7,[ecx-0x10]		/* XMM7 = rnd_const */	\
			__asm	add		ecx, __offset0	\
			__asm	movaps	xmm2,xmm4			/* copy x */	\
			__asm	shl		edi, 1	\
			__asm	addpd	xmm4,xmm7	\
			__asm	subpd	xmm4,xmm7			/* temp = DNINT(x) */	\
			__asm	mov		eax, sign_mask	\
			__asm	subpd	xmm2,xmm4			/* frac = [x - temp] */	\
			__asm	andpd	xmm2,[eax]			/* frac = fabs(frac) */	\
			__asm	maxpd	xmm2,xmm6		/* if(frac > maxerr) maxerr=frac */	\
			__asm	movaps	xmm6,xmm2		/* Note serialization here! */	\
			__asm	add		ecx, edi	\
			__asm	shr		edi, 1	\
			__asm	movaps	xmm2,xmm4			/* cpy temp */	\
			__asm	mulpd	xmm2,[ecx+edi]	/* temp*baseinv[0] */	\
			__asm	addpd	xmm2,xmm7	\
			__asm	subpd	xmm2,xmm7			/* [cx,cy] = DNINT(temp*baseinv[0]) */	\
			__asm	movaps	xmm3,xmm2			/* cpy [cx,cy] */	\
			__asm	mulpd	xmm3,[ecx    ]	/* [cx,cy]*base[0] */	\
			__asm	sub		ecx, __offset0	\
			__asm	subpd	xmm4,xmm3			/* XMM4 = [a.re,a.im] = temp-[cx,cy]*base[0] */	\
			/* Now do b-pair: [b.re,b.im] in xmm5, carry in xmm2, xmm3 free, wt_[re,im] in xmmA,B, xmm6 free, rnd_const in xmm7: */	\
			__asm	addpd	xmm5,xmm2	/* [b.re,b.im] + [cx,cy] */	\
			__asm	movaps	xmm2,xmm5			/* copy y */	\
			__asm	addpd	xmm5,xmm7	\
			__asm	subpd	xmm5,xmm7			/* temp = DNINT(y) */	\
			__asm	subpd	xmm2,xmm5			/* frac = [y - temp] */	\
			__asm	andpd	xmm2,[eax]			/* frac = fabs(frac) */	\
			__asm	maxpd	xmm2,xmm6		/* if(frac > maxerr) maxerr=frac */	\
			__asm	movaps	xmm6,xmm2		/* Note serialization here! */	\
			__asm	movaps	xmm2,xmm5			/* cpy temp */	\
			__asm	add		ecx, __offset1	\
			__asm	mulpd	xmm2,[ecx+edi]	/* temp*baseinv[1] */	\
			__asm	addpd	xmm2,xmm7	\
			__asm	subpd	xmm2,xmm7			/* [cx,cy] = DNINT(temp*baseinv[1]) */	\
			__asm	shl		edi, 1			/* prepare to re-subtract 2*odd_radix from local-store pointer */\
			__asm	movaps	xmm3,xmm2			/* cpy [cx,cy] */	\
			__asm	mulpd	xmm3,[ecx    ]	/* [cx,cy]*base[1] */	\
			__asm	sub		ecx, __offset1	\
			__asm	subpd	xmm5,xmm3			/* XMM5 = [b.re,b.im] = temp-[cx,cy]*base[1] */	\
			__asm	movaps	[ebx],xmm2			/* store cy_out */	\
			__asm	movaps	xmm2,xmm4	/* [a.re,a.im] copy */	\
			__asm	shufpd	xmm4,xmm5,0	/* x = [a.re,b.re] */	\
			__asm	shufpd	xmm2,xmm5,3	/* y = [a.im,b.im] */	\
			__asm	movaps	xmm5,xmm4	/* x copy */	\
			__asm	movaps	xmm3,xmm2	/* y copy */	\
			/* Forward acyclic-convo weight is (wt_re, +wt_im): */	\
			__asm	sub		ecx, edi	\
			__asm	mulpd	xmm4,xmm0	/* [x     ]*wt_re */	\
			__asm	mulpd	xmm3,xmm1	/* [y copy]*wt_im */	\
			__asm	movaps	[ecx-0x20],xmm6		/* Store maxerr */	\
			__asm	add		ecx, __offset0	\
			__asm	mulpd	xmm2,xmm0	/* [y     ]*wt_re */	\
			__asm	mulpd	xmm5,xmm1	/* [x copy]*wt_im */	\
			__asm	movaps	xmm0,[ecx]	/* [wt0,wt1] */	\
			__asm	subpd	xmm4,xmm3	/* rt = x*wt_re - y*wt_im */	\
			__asm	addpd	xmm5,xmm2	/* it = x*wt_im + y*wt_re */	\
			/* Forward IBDWT weight: */\
			__asm	mulpd	xmm4,xmm0	\
			__asm	mulpd	xmm5,xmm0	\
			__asm	movaps	[edx     ],xmm4 /* store rt = ~[a.re,b.re] */	\
			__asm	movaps	[edx+0x10],xmm5	/* store it = ~[a.im,b.im] */	\
			/* Prepare for next pair of complex data: */	\
			__asm	add		esi, __idx_incr	/* idx_offset += idx_incr */	\
			__asm	mov		__idx_offset, esi
		}

	/*************************************************************/
	/**************** MERSENNE-MOD CARRY MACROS ******************/
	/*************************************************************/

	/* These are versions specialized for power-of-2 runlengths: */

	/* SSE2 version assumes the following:

		- x and y are in xmm0 and xmm1 on entry, xmm2 and xmm3 hold the next-higher words in the chain, i.e. xmm0-3 are reserved;
		- wtA address [assumed 16-byte aligned] points to wtA[col], i.e. the low word of the resulting xmm load = set=0 array location
		- wtA,B addresses [unaligned] point to wtB,C[c02,co3-1]   , i.e. the hi  word of the resulting xmm load = set=0 array location
			[that means we need an unaligned load and a shufpd to swap lo,hi prior to using]
		- All four of the 32-bit address registers eax,ebx,ecx,edx are available;
		- Doubled rnd_const is in memloc half_arr-1;
		- Doubled wtl,wtn,wtlp1,wtnm1 pairs are in memlocs half_arr+16,17,18,19, respectively.

	The SSE2 version processes 4 complex data per macro invocation, e.g. carries among

		j:			re0->im0->re1->im1
		j+stride:	re0->im0->re1->im1

	In other words we replace 2 passes through the non-SSE2 sequence [e.g for radix = 16]:

	 cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0    );
		cmplx_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,bjmodn1,0x1);
		cmplx_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,bjmodn2,0x2);
		cmplx_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,bjmodn3,0x3);
		cmplx_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,bjmodn4,0x4);
		cmplx_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,bjmodn5,0x5);
		cmplx_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,bjmodn6,0x6);
		cmplx_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,bjmodn7,0x7);
		cmplx_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,bjmodn8,0x8);
		cmplx_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,bjmodn9,0x9);
		cmplx_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,bjmodnA,0xA);
		cmplx_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,bjmodnB,0xB);
		cmplx_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,bjmodnC,0xC);
		cmplx_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,bjmodnD,0xD);
		cmplx_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,bjmodnE,0xE);
		cmplx_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,bjmodnF,0xF);

	With one pass through the SSE2ified versions of the above macros:

		SSE2_cmplx_carry_norm_pow2_errcheck0(r1 ,add0,add1,add2,cy_r01,bjmodn0,bjmodn1);
		SSE2_cmplx_carry_norm_pow2_errcheck (r5 ,add0,add1,add2,cy_r23,bjmodn2,bjmodn3);
		SSE2_cmplx_carry_norm_pow2_errcheck (r9 ,add0,add1,add2,cy_r45,bjmodn4,bjmodn5);
		SSE2_cmplx_carry_norm_pow2_errcheck (r13,add0,add1,add2,cy_r67,bjmodn6,bjmodn7);
		SSE2_cmplx_carry_norm_pow2_errcheck (r17,add0,add1,add2,cy_r89,bjmodn8,bjmodn9);
		SSE2_cmplx_carry_norm_pow2_errcheck (r21,add0,add1,add2,cy_rAB,bjmodnA,bjmodnB);
		SSE2_cmplx_carry_norm_pow2_errcheck (r25,add0,add1,add2,cy_rCD,bjmodnC,bjmodnD);
		SSE2_cmplx_carry_norm_pow2_errcheck (r29,add0,add1,add2,cy_rEF,bjmodnE,bjmodnF);
		[address-calc stuff]
		SSE2_cmplx_carry_norm_pow2_errcheck2(r1 ,add0,add1     ,cy_r01,bjmodn0,bjmodn1);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r5 ,add0,add1     ,cy_r23,bjmodn2,bjmodn3);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r9 ,add0,add1     ,cy_r45,bjmodn4,bjmodn5);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r13,add0,add1     ,cy_r67,bjmodn6,bjmodn7);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r17,add0,add1     ,cy_r89,bjmodn8,bjmodn9);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r21,add0,add1     ,cy_rAB,bjmodnA,bjmodnB);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r25,add0,add1     ,cy_rCD,bjmodnC,bjmodnD);
		SSE2_cmplx_carry_norm_pow2_errcheck2(r29,add0,add1     ,cy_rEF,bjmodnE,bjmodnF);
	*/
		#define SSE2_cmplx_carry_norm_pow2_errcheck0(__data,__wtA,__wtB,__wtC,__cy,__bjmod_0,__bjmod_1)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */\
				__asm	movaps		xmm0,[eax     ]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm0,[eax+0x20]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */\
				__asm	movaps		xmm1,[eax+0x10]\
				__asm	movaps		xmm3,[eax+0x10]\
				__asm	unpcklpd	xmm1,[eax+0x30]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm1	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		edx, i	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_sil\
				__asm	mov		edi, n_minus_sil\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_sil-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_sil-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwt) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwt) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax     ]	// R0~\
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtB[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	sub		ecx, 10H	/* add1 -= 2 */\
				__asm	mov		__wtB, ecx\
				__asm	mulpd		xmm1,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* x */\
				__asm	movaps		[eax     ],xmm0	/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sw\
				__asm	sub		edx, ebx\
				__asm	shr		edx, 31	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_silp1\
				__asm	mov		edi, n_minus_silp1\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_silp1-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_silp1-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwtm1) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x10]	/* I0~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtC[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	add		ebx, 10H	/* add0 += 2 */\
				__asm	sub		ecx, 10H	/* add2 -= 2 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	mulpd		xmm1,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* y */\
				__asm	movaps		[eax+0x10],xmm0	/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
		}

		/************* Use this for all first-pass carry pairs but the first: **************/

		#define SSE2_cmplx_carry_norm_pow2_errcheck1(__data,__wtA,__wtB,__wtC,__cy,__bjmod_0,__bjmod_1)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */\
				__asm	movaps		xmm0,[eax     ]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm0,[eax+0x20]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */\
				__asm	movaps		xmm1,[eax+0x10]\
				__asm	movaps		xmm3,[eax+0x10]\
				__asm	unpcklpd	xmm1,[eax+0x30]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm1	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		edx, sw\
				__asm	sub		edx, ebx\
				__asm	shr		edx, 31	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_sil\
				__asm	mov		edi, n_minus_sil\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_sil-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_sil-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwt) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwt) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax     ]	// R0~\
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtB[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	sub		ecx, 10H	/* add1 -= 2 */\
				__asm	mov		__wtB, ecx\
				__asm	mulpd		xmm1,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* x */\
				__asm	movaps		[eax     ],xmm0	/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sw\
				__asm	sub		edx, ebx\
				__asm	shr		edx, 31	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_silp1\
				__asm	mov		edi, n_minus_silp1\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_silp1-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_silp1-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwtm1) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x10]	/* I0~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtC[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	add		ebx, 10H	/* add0 += 2 */\
				__asm	sub		ecx, 10H	/* add2 -= 2 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	mulpd		xmm1,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* y */\
				__asm	movaps		[eax+0x10],xmm0	/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
		}

		/********* 2nd-pass version of the above; no special 0-index case needed here:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck2(__data,__wtA,__wtB,__cy,__bjmod_0,__bjmod_1)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		edx, sw\
				__asm	sub		edx, ebx\
				__asm	shr		edx, 31	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_sil\
				__asm	mov		edi, n_minus_sil\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_sil-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_sil-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwt) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwt) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x20]	/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtB[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	mulpd		xmm1,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* x */\
				__asm	movaps		[eax+0x20],xmm0	/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sw\
				__asm	sub		edx, ebx\
				__asm	shr		edx, 31	/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	mov		esi, sw\
				__asm	sub		esi, ecx\
				__asm	shr		esi, 31	/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shl		esi, 1	/* (i1 << 1) */\
				__asm	add		esi, edx\
				__asm	shl		esi, 28	/* i0 = (i0 + (i1 << 1)) << 4; Address offset into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		eax, n_minus_silp1\
				__asm	mov		edi, n_minus_silp1\
				__asm	sub		eax, ebx\
				__asm	sub		edi, ecx\
				__asm	shr		eax, 31		/* m0=((uint32)(n_minus_silp1-__bjmod_0) >> 31); */\
				__asm	shr		edi, 31		/* m1=((uint32)(n_minus_silp1-__bjmod_1) >> 31); */\
				__asm	shl		edi, 1	/* (m1 << 1) */\
				__asm	add		eax, edi\
				__asm	shl		eax, 20	/*	m0 = (m0 + (m1 << 1)) << 4;	Address offset into one_half[m01] table; move into byte[2] of eax... */\
				__asm	add		esi, eax	/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				__asm	sub		ebx, edi	/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ecx, edi	/* (__bjmod_1 - sinwtm1) */\
				__asm	shr		ebx, 31		/* m2=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31		/* m3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ecx, 1	/* (m3 << 1) */\
				__asm	add		ebx, ecx\
				__asm	shl		ebx, 12	/*	m2 = (m2 + (m3 << 1)) << 4;	Address offset into one_half[m23] table; move into byte[1] of ebx... */\
				__asm	add		esi, ebx	/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x30]	/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr	/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	shld	edi,esi,16\
				__asm	and		edi, 000000ffH	/* m0 */\
				__asm	shld	edx,esi,24\
				__asm	and		edx, 000000ffH	/* m2 */\
				__asm	add		edi, eax\
				__asm	add		edx, eax\
				__asm	movaps		xmm1,[ebx]	/* wtA[j  ]; ebx FREE */\
				__asm	movupd		xmm2,[ecx]	/* wtC[j-1]; ecx FREE */\
				__asm	shufpd		xmm2,xmm2,1\
				__asm	add		ebx, 10H	/* add0 += 2 */\
				__asm	sub		ecx, 10H	/* add1 -= 2 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	mulpd		xmm1,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]	/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]	/* wtinv=wtinv*one_half[4+m23] */\
				__asm	mov		ecx, __cy		/* cy_in */\
				__asm	mulpd		xmm0,xmm2	/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]	/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0	/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask\
				__asm	subpd		xmm0,xmm3	/* x - temp */\
				__asm	andpd		xmm0,[ebx]	/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,[eax-0x20]	/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edx,esi, 8\
				__asm	and		edx, 000000ffH	/* i0 */\
				__asm	add		edx, eax\
				__asm	movaps		xmm0,xmm3	/* cpy temp */\
				__asm	mulpd		xmm3,[edx+0xc0]	/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]\
				__asm	subpd		xmm3,[eax-0x10]	/* cy_out */\
				__asm	movaps		[ecx],xmm3	/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edx+0x80]	/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3	/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1	/* y */\
				__asm	movaps		[eax+0x30],xmm0	/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, __bjmod_0\
				__asm	mov		ecx, __bjmod_1\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	and		ebx, edi	/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ecx, edi	/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	mov		__bjmod_0, ebx\
				__asm	mov		__bjmod_1, ecx\
		/***************Repack the data:*************************/\
				__asm	movaps	xmm1,[eax+0x10]	/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]	/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1	/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0	/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]\
				__asm	unpcklpd	xmm1,[eax+0x30]\
				__asm	movaps	[eax+0x30],xmm3	/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]\
				__asm	unpcklpd	xmm0,[eax+0x20]\
				__asm	movaps	[eax+0x20],xmm2	/* Store hi real in aj2 */\
				__asm	movaps	[eax+0x10],xmm1	/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0	/* a[jt+p0 ] */\
		}

		/********* Double-wide version of SSE2_cmplx_carry_norm_pow2_errcheck0:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck0_2x(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0,__bjmod_1,__bjmod_2,__bjmod_3)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm0,[eax     ]		__asm	movaps		xmm4,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		edi, sw\
				/* Load bjmodn pointers: */\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	mov		eax, i					/*	i0=i for first block */\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_sil\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_sil */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_sil */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_sil */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_sil */\
				__asm	neg		eax						/*	n_minus_sil - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_sil - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_sil - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_sil - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_sil - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_sil - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_sil - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_sil - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwt) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwt) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwt) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwt) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax     ]		__asm	movaps		xmm4,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x100]	__asm	mulpd		xmm5,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	__asm	mulpd		xmm6,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* x - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* x */\
				__asm	movaps		[eax     ],xmm0		__asm	movaps		[eax+0x40],xmm4		/* store x */\
				/* Get ready for next set [IM0~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edi, sw\
				__asm	sub		eax, edi				/*	__bjmod_0 - sw */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		eax						/*	sw - __bjmod_0 */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	shr		eax, 31					/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_silp1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_silp1 */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_silp1 */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_silp1 */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_silp1 */\
				__asm	neg		eax						/*	n_minus_silp1 - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_silp1 - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_silp1 - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_silp1 - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_silp1 - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_silp1 - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_silp1 - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_silp1 - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwtm1) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwtm1) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwtm1) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x10]		__asm	movaps		xmm4,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x120]	__asm	mulpd		xmm5,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	__asm	mulpd		xmm6,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(y) */\
										/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* y - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* y */\
				__asm	movaps		[eax+0x10],xmm0		__asm	movaps		[eax+0x50],xmm4		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
		}

		/********* Double-wide version of SSE2_cmplx_carry_norm_pow2_errcheck1:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck1_2x(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0,__bjmod_1,__bjmod_2,__bjmod_3)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm0,[eax     ]		__asm	movaps		xmm4,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		edi, sw\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - sw */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		eax						/*	sw - __bjmod_0 */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	shr		eax, 31					/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_sil\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_sil */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_sil */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_sil */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_sil */\
				__asm	neg		eax						/*	n_minus_sil - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_sil - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_sil - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_sil - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_sil - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_sil - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_sil - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_sil - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwt) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwt) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwt) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwt) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax     ]		__asm	movaps		xmm4,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x100]	__asm	mulpd		xmm5,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	__asm	mulpd		xmm6,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* x - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* x */\
				__asm	movaps		[eax     ],xmm0		__asm	movaps		[eax+0x40],xmm4		/* store x */\
				/* Get ready for next set [IM0~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edi, sw\
				__asm	sub		eax, edi				/*	__bjmod_0 - sw */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		eax						/*	sw - __bjmod_0 */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	shr		eax, 31					/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_silp1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_silp1 */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_silp1 */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_silp1 */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_silp1 */\
				__asm	neg		eax						/*	n_minus_silp1 - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_silp1 - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_silp1 - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_silp1 - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_silp1 - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_silp1 - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_silp1 - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_silp1 - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwtm1) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwtm1) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwtm1) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x10]		__asm	movaps		xmm4,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x120]	__asm	mulpd		xmm5,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	__asm	mulpd		xmm6,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(y) */\
										/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* y - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* y */\
				__asm	movaps		[eax+0x10],xmm0		__asm	movaps		[eax+0x50],xmm4		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
		}

		/********* Double-wide version of SSE2_cmplx_carry_norm_pow2_errcheck2:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck2_2x(__data,__wtA,__wtB,__cyA,__cyB,__bjmod_0,__bjmod_1,__bjmod_2,__bjmod_3)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, sw	/* Need this starting e*x register ref to work around MSVC state-save bug */\
				__asm	mov		edi, sw\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - sw */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		eax						/*	sw - __bjmod_0 */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	shr		eax, 31					/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_sil\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_sil */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_sil */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_sil */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_sil */\
				__asm	neg		eax						/*	n_minus_sil - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_sil - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_sil - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_sil - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_sil - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_sil - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_sil - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_sil - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwt\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwt) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwt) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwt) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwt) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwt) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwt) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x20]		__asm	movaps		xmm4,[eax+0x60]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x100]	__asm	mulpd		xmm5,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm2,[eax+0x110]	__asm	mulpd		xmm6,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* x = x*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = x */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* x - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* x */\
				__asm	movaps		[eax+0x20],xmm0		__asm	movaps		[eax+0x60],xmm4		/* store x */\
				/* Get ready for next set [IM0~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edi, sw\
				__asm	sub		eax, edi				/*	__bjmod_0 - sw */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - sw */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - sw */\
				__asm	sub		edx, edi				/*	__bjmod_3 - sw */\
				__asm	neg		eax						/*	sw - __bjmod_0 */\
				__asm	neg		ebx						/*	sw - __bjmod_1 */\
				__asm	neg		ecx						/*	sw - __bjmod_2 */\
				__asm	neg		edx						/*	sw - __bjmod_3 */\
				__asm	shr		eax, 31					/*	i0=((uint32)(sw - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	i1=((uint32)(sw - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	i2=((uint32)(sw - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	i3=((uint32)(sw - __bjmod_3) >> 31);	*/\
				__asm	mov		esi, eax				/*	for result */\
				__asm	shl		ebx, 1					/* (i1 << 1) */\
				__asm	shl		ecx, 2					/* (i2 << 2) */\
				__asm	shl		edx, 3					/* (i3 << 3) */\
				__asm	add		esi, ebx\
				__asm	add		esi, ecx\
				__asm	add		esi, edx\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	mov		edi, n_minus_silp1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/*	__bjmod_0 - n_minus_silp1 */\
				__asm	sub		ebx, edi				/*	__bjmod_1 - n_minus_silp1 */\
				__asm	sub		ecx, edi				/*	__bjmod_2 - n_minus_silp1 */\
				__asm	sub		edx, edi				/*	__bjmod_3 - n_minus_silp1 */\
				__asm	neg		eax						/*	n_minus_silp1 - __bjmod_0 */\
				__asm	neg		ebx						/*	n_minus_silp1 - __bjmod_1 */\
				__asm	neg		ecx						/*	n_minus_silp1 - __bjmod_2 */\
				__asm	neg		edx						/*	n_minus_silp1 - __bjmod_3 */\
				__asm	shr		eax, 31					/*	m0=((uint32)(n_minus_silp1 - __bjmod_0) >> 31);	*/\
				__asm	shr		ebx, 31					/*	m1=((uint32)(n_minus_silp1 - __bjmod_1) >> 31);	*/\
				__asm	shr		ecx, 31					/*	m2=((uint32)(n_minus_silp1 - __bjmod_2) >> 31);	*/\
				__asm	shr		edx, 31					/*	m3=((uint32)(n_minus_silp1 - __bjmod_3) >> 31);	*/\
				__asm	shl		ebx, 1					/* (m1 << 1) */\
				__asm	shl		ecx, 2					/* (m2 << 2) */\
				__asm	shl		edx, 3					/* (m3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		edi, sinwtm1\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	sub		eax, edi				/* (__bjmod_0 - sinwtm1) */\
				__asm	sub		ebx, edi				/* (__bjmod_1 - sinwtm1) */\
				__asm	sub		ecx, edi				/* (__bjmod_2 - sinwtm1) */\
				__asm	sub		edx, edi				/* (__bjmod_3 - sinwtm1) */\
				__asm	shr		eax, 31					/* n0=1 + ((uint32)(__bjmod_0 - sinwtm1) >> 31); */\
				__asm	shr		ebx, 31					/* n1=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		ecx, 31					/* n2=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shr		edx, 31					/* n3=1 + ((uint32)(__bjmod_1 - sinwtm1) >> 31); */\
				__asm	shl		ebx, 1					/* (n1 << 1) */\
				__asm	shl		ecx, 2					/* (n2 << 2) */\
				__asm	shl		edx, 3					/* (n3 << 3) */\
				__asm	add		eax, ebx\
				__asm	add		eax, ecx\
				__asm	add		eax, edx\
				__asm	shl		eax, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of eax... */\
				__asm	add		esi, eax				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm0,[eax+0x30]		__asm	movaps		xmm4,[eax+0x70]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	movaps		xmm1,[ebx]			__asm	movaps		xmm5,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm2,[ecx]			__asm	movhpd		xmm6,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm2,[ecx+0x08]		__asm	movlpd		xmm6,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm1,[eax+0x120]	__asm	mulpd		xmm5,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm2,[eax+0x130]	__asm	mulpd		xmm6,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm1,[edi     ]		__asm	mulpd		xmm5,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm2,[edx+0x40]		__asm	mulpd		xmm6,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm0,xmm2			__asm	mulpd		xmm4,xmm6			/* y = y*wtinv */\
				__asm	addpd		xmm0,[ecx]			__asm	addpd		xmm4,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm0			__asm	movaps		xmm7,xmm4			/* temp = y */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* temp = DNINT(y) */\
										/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* y - temp */\
				__asm	andpd		xmm0,[ebx]			__asm	andpd		xmm4,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm0,xmm4					__asm	maxpd		xmm0,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm0		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm0,xmm3			__asm	movaps		xmm4,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,[eax-0x10]		__asm	addpd		xmm7,[eax-0x10]				__asm	subpd		xmm3,[eax-0x10]		__asm	subpd		xmm7,[eax-0x10]		/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* y = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm0,xmm3			__asm	subpd		xmm4,xmm7			/* (temp-cy*base[i1]) */\
				__asm	mulpd		xmm0,xmm1			__asm	mulpd		xmm4,xmm5			/* y */\
				__asm	movaps		[eax+0x30],xmm0		__asm	movaps		[eax+0x70],xmm4		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				/* Load bjmodn pointers: */\
				__asm	mov		eax, __bjmod_0\
				__asm	mov		ebx, __bjmod_1\
				__asm	mov		ecx, __bjmod_2\
				__asm	mov		edx, __bjmod_3\
				/* Dereference 'em: */\
				__asm	mov		eax, [eax]\
				__asm	mov		ebx, [ebx]\
				__asm	mov		ecx, [ecx]\
				__asm	mov		edx, [edx]\
				__asm	mov		esi, bw\
				__asm	mov		edi, nm1\
				__asm	add		eax, esi\
				__asm	add		ebx, esi\
				__asm	add		ecx, esi\
				__asm	add		edx, esi\
				__asm	and		eax, edi				/* __bjmod_0 = (__bjmod_0 + bw) & nm1; */\
				__asm	and		ebx, edi				/* __bjmod_1 = (__bjmod_1 + bw) & nm1; */\
				__asm	and		ecx, edi				/* __bjmod_2 = (__bjmod_2 + bw) & nm1; */\
				__asm	and		edx, edi				/* __bjmod_3 = (__bjmod_3 + bw) & nm1; */\
				/* Rereference 'em: */\
				__asm	mov		edi, __bjmod_0\
				__asm	mov		[edi], eax\
				__asm	mov		edi, __bjmod_1\
				__asm	mov		[edi], ebx\
				__asm	mov		edi, __bjmod_2\
				__asm	mov		[edi], ecx\
				__asm	mov		edi, __bjmod_3\
				__asm	mov		[edi], edx\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __data\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[eax+0x50]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[eax+0x40]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]				__asm	movaps	[eax+0x30],xmm3			__asm	movaps	[eax+0x70],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]				__asm	movaps	[eax+0x20],xmm2			__asm	movaps	[eax+0x60],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[eax+0x50],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[eax+0x40],xmm4			/* a[jt+p0 ] */\
		}

		/******************************************************************************************/
		/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck0_2x:***********/
		/******************************************************************************************/

		#define SSE2_cmplx_carry_norm_pow2_errcheck0_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				/* i0=i for first block: */\
				__asm	mov		ecx, i							__asm	and		esi, 0xfffffffe			/* Mask off lowest bit */\
				__asm	add		esi, ecx						__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw  */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1 */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck1_2x:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck1_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw  */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1 */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/

		#define SSE2_cmplx_carry_norm_pow2_errcheck2_2B(__data,__wtA,__wtB,__cyA,__cyB,__bjmod_0)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm1,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x20]		__asm	movaps		xmm5,[eax+0x60]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax+0x20],xmm1		__asm	movaps		[eax+0x60],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x30]		__asm	movaps		xmm5,[eax+0x70]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x30],xmm1		__asm	movaps		[eax+0x70],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw  */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1 */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __data\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[eax+0x50]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[eax+0x40]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]				__asm	movaps	[eax+0x30],xmm3			__asm	movaps	[eax+0x70],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]				__asm	movaps	[eax+0x20],xmm2			__asm	movaps	[eax+0x60],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[eax+0x50],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[eax+0x40],xmm4			/* a[jt+p0 ] */\
		}

		/******************************************************************************************************************************************************************/
		/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
		/******************************************************************************************************************************************************************/

		/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck1_2x:***********/

		#define SSE2_cmplx_carry_norm_pow2_nocheck1_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw  */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1 */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/

		#define SSE2_cmplx_carry_norm_pow2_nocheck2_2B(__data,__wtA,__wtB,__cyA,__cyB,__bjmod_0)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm1,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x20]		__asm	movaps		xmm5,[eax+0x60]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax+0x20],xmm1		__asm	movaps		[eax+0x60],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1; & doesn't care whether integer [pand] or floating [andpd], but data are int, so use pand for form's sake */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x30]		__asm	movaps		xmm5,[eax+0x70]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x30],xmm1		__asm	movaps		[eax+0x70],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		eax, sse_bw\
				__asm	mov		ebx, sse_nm1\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw  */\
				__asm	pand	xmm0,[ebx]				/* bjmod[0:3] &= nm1 */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __data\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[eax+0x50]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[eax+0x40]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]				__asm	movaps	[eax+0x30],xmm3			__asm	movaps	[eax+0x70],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]				__asm	movaps	[eax+0x20],xmm2			__asm	movaps	[eax+0x60],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[eax+0x50],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[eax+0x40],xmm4			/* a[jt+p0 ] */\
		}



		/******************************************************************************************/
		/********* Non-power-of-2-FFT version of SSE2_cmplx_carry_norm_pow2_errcheck0_2B:**********/
		/******************************************************************************************/

		#define SSE2_cmplx_carry_norm_errcheck0_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw						__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
														/* i0=i for first block: */\
				__asm	mov		ecx, i							__asm	and		esi, 0xfffffffe			/* Mask off lowest bit */\
				__asm	add		esi, ecx						__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr			/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Non-power-of-2-FFT version of SSE2_cmplx_carry_norm_pow2_errcheck1_2B:**********/

		#define SSE2_cmplx_carry_norm_errcheck1_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Non-power-of-2-FFT version of SSE2_cmplx_carry_norm_pow2_errcheck2_2B:**********/

		#define SSE2_cmplx_carry_norm_errcheck2_2B(__data,__wtA,__wtB,__cyA,__cyB,__bjmod_0)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm1,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x20]		__asm	movaps		xmm5,[eax+0x60]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax+0x20],xmm1		__asm	movaps		[eax+0x60],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x30]		__asm	movaps		xmm5,[eax+0x70]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x30],xmm1		__asm	movaps		[eax+0x70],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __data\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[eax+0x50]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[eax+0x40]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]				__asm	movaps	[eax+0x30],xmm3			__asm	movaps	[eax+0x70],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]				__asm	movaps	[eax+0x20],xmm2			__asm	movaps	[eax+0x60],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[eax+0x50],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[eax+0x40],xmm4			/* a[jt+p0 ] */\
		}

		/******************************************************************************************************************************************************************/
		/********** No-ROE-Check versions of the latter 2 of the above 3 macros - the first is called too infrequently to bother with a special non-ROE version: **********/
		/******************************************************************************************************************************************************************/

		#define SSE2_cmplx_carry_norm_nocheck1_2B(__data,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __data\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[eax+0x40]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[eax+0x20]		__asm	unpcklpd	xmm5,[eax+0x60]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[eax+0x20],xmm2		__asm	movaps		[eax+0x60],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[eax+0x50]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[eax+0x50]	\
				__asm	unpcklpd	xmm2,[eax+0x30]		__asm	unpcklpd	xmm6,[eax+0x70]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[eax+0x50],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[eax+0x30],xmm3		__asm	movaps		[eax+0x70],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[eax+0x40],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[eax+0x50]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[eax+0x50],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		/********* Non-power-of-2-FFT version of SSE2_cmplx_carry_norm_pow2_nocheck2_2B:**********/

		#define SSE2_cmplx_carry_norm_nocheck2_2B(__data,__wtA,__wtB,__cyA,__cyB,__bjmod_0)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm1,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x20]		__asm	movaps		xmm5,[eax+0x60]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax+0x20],xmm1		__asm	movaps		[eax+0x60],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax+0x30]		__asm	movaps		xmm5,[eax+0x70]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __data\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x30],xmm1		__asm	movaps		[eax+0x70],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __data\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[eax+0x50]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[eax+0x40]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[eax+0x30]		__asm	unpckhpd	xmm7,[eax+0x70]				__asm	unpcklpd	xmm1,[eax+0x30]		__asm	unpcklpd	xmm5,[eax+0x70]				__asm	movaps	[eax+0x30],xmm3			__asm	movaps	[eax+0x70],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[eax+0x20]		__asm	unpckhpd	xmm6,[eax+0x60]				__asm	unpcklpd	xmm0,[eax+0x20]		__asm	unpcklpd	xmm4,[eax+0x60]				__asm	movaps	[eax+0x20],xmm2			__asm	movaps	[eax+0x60],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[eax+0x50],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[eax+0x40],xmm4			/* a[jt+p0 ] */\
		}

		/******************************************************************************************************************************************************************/
		/********** In-place-suitable versions of the 3 key carry routines above. *****************************************************************************************/
		/******************************************************************************************************************************************************************/

		#define SSE2_cmplx_carry_norm_errcheck0_2C(__in0,__in1,__in2,__in3,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __in0\
				__asm	mov		ebx, __in1\
				__asm	mov		ecx, __in2\
				__asm	mov		edx, __in3\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[ecx     ]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[ecx     ]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[ebx     ]		__asm	unpcklpd	xmm5,[edx     ]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[ebx     ]		__asm	unpckhpd	xmm6,[edx     ]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[ebx     ],xmm2		__asm	movaps		[edx     ],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[ecx+0x10]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[ecx+0x10]	\
				__asm	unpcklpd	xmm2,[ebx+0x10]		__asm	unpcklpd	xmm6,[edx+0x10]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[ebx+0x10]		__asm	unpckhpd	xmm7,[edx+0x10]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[ecx+0x10],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[ebx+0x10],xmm3		__asm	movaps		[edx+0x10],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				/* i0=i for first block: */\
				__asm	mov		ecx, i							__asm	and		esi, 0xfffffffe			/* Mask off lowest bit */\
				__asm	add		esi, ecx						__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(x-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[edx     ],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[edx+0x10]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				__asm	mov		ebx, sign_mask					__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y - temp */\
				__asm	andpd		xmm1,[ebx]			__asm	andpd		xmm5,[ebx]			/* frac = fabs(y-temp) */\
				__asm	maxpd		xmm1,xmm5					__asm	maxpd		xmm1,[eax-0x20]		/* if(frac > maxerr) maxerr=frac */\
				__asm	movaps		[eax-0x20],xmm1		/* Note serialization here! */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[edx+0x10],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		#define SSE2_cmplx_carry_norm_nocheck1_2C(__in0,__in1,__in2,__in3,__wtA,__wtB,__wtC,__cyA,__cyB,__bjmod_0)\
		{\
		/***************Unpack the data:*************************/\
				__asm	mov		eax, __in0\
				__asm	mov		ebx, __in1\
				__asm	mov		ecx, __in2\
				__asm	mov		edx, __in3\
				/* Real parts: */							__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[ecx     ]	/* r1, this is the active  xmm register */\
				__asm	movaps		xmm2,[eax     ]		__asm	movaps		xmm6,[ecx     ]	/* r1, this is the scratch xmm register */\
				__asm	unpcklpd	xmm1,[ebx     ]		__asm	unpcklpd	xmm5,[edx     ]	/* r1 -x- r3 (lo halves) ==> R0~ */\
				__asm	unpckhpd	xmm2,[ebx     ]		__asm	unpckhpd	xmm6,[edx     ]	/* r1 -x- r3 (hi halves) ==> R1~ */\
				__asm	movaps		[ebx     ],xmm2		__asm	movaps		[edx     ],xmm6	/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Imag parts: */							__asm	movaps		xmm2,[eax+0x10]		__asm	movaps		xmm6,[ecx+0x10]	\
				__asm	movaps		xmm3,[eax+0x10]		__asm	movaps		xmm7,[ecx+0x10]	\
				__asm	unpcklpd	xmm2,[ebx+0x10]		__asm	unpcklpd	xmm6,[edx+0x10]	/* r2 -x- r4 (lo halves) ==> I0~ */\
				__asm	unpckhpd	xmm3,[ebx+0x10]		__asm	unpckhpd	xmm7,[edx+0x10]	/* r2 -x- r4 (hi halves) ==> I1~ */\
				__asm	movaps		[eax+0x10],xmm2		__asm	movaps		[ecx+0x10],xmm6	/* Tmp store I0~ until needed by imaginary-part-processing section */\
				__asm	movaps		[ebx+0x10],xmm3		__asm	movaps		[edx+0x10],xmm7	/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
				/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm7,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm7,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm7,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm7,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm7				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
		/*\
				__asm	mov		eax, __data\
				__asm	movaps		xmm1,[eax     ]		__asm	movaps		xmm5,[eax+0x40]		// R1~ \
		*/\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]		/* [NOTE: movhpd/movlpd preferable to movupd/shufpd] */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[edx     ],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[edx+0x10]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtC				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add2 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtC, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in0\
				__asm	mov		edx, __in2\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[edx+0x10],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		}

		#define SSE2_cmplx_carry_norm_nocheck2_2C(__in0,__in1,__in2,__in3,__wtA,__wtB,__cyA,__cyB,__bjmod_0)\
		{\
			/**********************************************/\
			/*          Real      parts                   */\
			/**********************************************/\
				__asm	mov		eax, __bjmod_0			/* Pointer to bjmodn data */\
				__asm	movaps	xmm0,[eax]				/* bjmod[0:3] */\
				__asm	mov		ebx, sse_sw\
				__asm	movaps	xmm1,[ebx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_sil\
				__asm	movd	xmm2,ecx				/* n_minus_sil in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_sil - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwt\
				__asm	movd	xmm3,edx				/* sinwt in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwt */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		ebx, __in1\
				__asm	mov		edx, __in3\
				__asm	movaps		xmm1,[ebx     ]		__asm	movaps		xmm5,[edx     ]		/* R1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB\
				__asm	movaps		xmm4,[eax-0x10]		/* RND_CONST */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x100]	__asm	mulpd		xmm6,[eax+0x100]	/* wt   =wtA*wtl */\
				__asm	mulpd		xmm3,[eax+0x110]	__asm	mulpd		xmm7,[eax+0x110]	/* wtinv=wtB*wtn */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* x = x*wtinv; xmm3,xmm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* x = x*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = x */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(x) */\
				/*\
				frac = fabs(x-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in1\
				__asm	mov		edx, __in3\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* x = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* x*= wt */\
				__asm	movaps		[eax     ],xmm1		__asm	movaps		[edx     ],xmm5		/* store x */\
				/* Get ready for next set [IM0~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
			/**********************************************/\
			/*          Imaginary parts                   */\
			/**********************************************/\
				__asm	mov		edx, sse_sw\
				__asm	movaps	xmm1,[edx]				/* sw[0:3] */\
				__asm	psubd	xmm1,xmm0				/* sw[0:3] - bjmod[0:3] */\
				__asm	movmskps esi,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		esi, 24					/* <i3|i2|i1|i0>; Packed indices into base,base_inv tables; move into leftmost byte[3] */\
				__asm	movaps	xmm1,xmm0				/* bjmod[0:3] COPY */\
				__asm	mov		ecx, n_minus_silp1\
				__asm	movd	xmm2,ecx				/* n_minus_silp1 in low 32 bits of xmm2 */\
				__asm	pshufd	xmm2,xmm2,0				/* Broadcast low 32 bits of xmm2 to all 4 slots of xmm2 */\
				__asm	psubd	xmm2,xmm0				/* n_minus_silp1 - bjmod[0:3] */\
				__asm	movmskps ecx,xmm2				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		ecx, 16					/* <m3|m2|m1|m0>; Packed indices into base,base_inv tables; move into leftmost byte[2] of ecx... */\
				__asm	add		esi, ecx				/* ....and fold into esi. PERSISTENT COPY OF BJMOD[0:3] REMAINS IN XMM0. */\
				__asm	mov		edx, sinwtm1\
				__asm	movd	xmm3,edx				/* sinwtm1 in low 32 bits of xmm3 */\
				__asm	pshufd	xmm3,xmm3,0				/* Broadcast low 32 bits of xmm3 to all 4 slots of xmm3 */\
				__asm	psubd	xmm1,xmm3				/* bjmod[0:3] - sinwtm1 */\
				__asm	movmskps edx,xmm1				/* Extract sign bits into 4-bit signmask */\
				__asm	shl		edx, 8					/* <n3|n2|n1|n0>; Packed indices into base,base_inv tables; move into leftmost byte[1] of edx... */\
				__asm	add		esi, edx				/* ....and fold into esi. */\
				__asm	mov		eax, __in1\
				__asm	mov		edx, __in3\
				__asm	movaps		xmm1,[eax+0x10]		__asm	movaps		xmm5,[edx+0x10]		/* I1~ */\
				/* Don't explicitly load address of sse2_rnd, since we know it's in [half_arr - 0x10]. */\
				__asm	mov		eax, half_arr							/* This is a real array address from the calling routine, hence no prepended __ . */\
				__asm	mov		ebx, __wtA\
				__asm	mov		ecx, __wtB				/* wtB == wtC for this latter set of carries */\
				__asm	movaps		xmm2,[ebx]			__asm	movaps		xmm6,[ebx+0x10]		/* wtA[j  ]; ebx FREE */\
				__asm	movhpd		xmm3,[ecx]			__asm	movhpd		xmm7,[ecx-0x10]		/* wtC[j-1]; ecx FREE */\
				__asm	movlpd		xmm3,[ecx+0x08]		__asm	movlpd		xmm7,[ecx-0x08]				__asm	add		ebx, 20H	/* add0 += 4 */\
				__asm	sub		ecx, 20H	/* add1 -= 4 */\
				__asm	mov		__wtA, ebx\
				__asm	mov		__wtB, ecx\
				__asm	shld	edi,esi,20				__asm	shld	ebx,esi,18\
				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* m0 */\
				__asm	shld	edx,esi,28				__asm	shld	ecx,esi,26\
				__asm	and		edx, 00000030H			__asm	and		ecx, 00000030H			/* m2 */\
				__asm	add		edi, eax				__asm	add		ebx, eax\
				__asm	add		edx, eax				__asm	add		ecx, eax\
				__asm	mulpd		xmm2,[eax+0x120]	__asm	mulpd		xmm6,[eax+0x120]	/* wt   =wtA*wtlp1 */\
				__asm	mulpd		xmm3,[eax+0x130]	__asm	mulpd		xmm7,[eax+0x130]	/* wtinv=wtC*wtnm1 */\
				__asm	mulpd		xmm2,[edi     ]		__asm	mulpd		xmm6,[ebx     ]		/* wt   =wt   *one_half[m01] */\
				__asm	mulpd		xmm3,[edx+0x40]		__asm	mulpd		xmm7,[ecx+0x40]		/* wtinv=wtinv*one_half[4+m23] */\
																										__asm	mov		ecx, __cyA				__asm	mov		edx, __cyB				/* cy_in */\
				__asm	mulpd		xmm1,xmm3			__asm	mulpd		xmm5,xmm7			/* y = y*wtinv; ymm3,ymm7 FREE */\
				__asm	addpd		xmm1,[ecx]			__asm	addpd		xmm5,[edx]			/* y = y*wtinv + cy */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* temp = y */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* temp = DNINT(y) */\
				/*\
				frac = fabs(y-temp);\
				if(frac > maxerr) maxerr=frac;\
				*/\
				/* NO ROE HERE */\
				/* cy   = DNINT(temp*baseinv[i1]) */\
				__asm	shld	edi,esi,12				__asm	shld	ebx,esi,10				__asm	and		edi, 00000030H			__asm	and		ebx, 00000030H			/* i0 */\
				__asm	add		edi, eax				__asm	add		ebx, eax				__asm	movaps		xmm1,xmm3			__asm	movaps		xmm5,xmm7			/* cpy temp */\
				__asm	mulpd		xmm3,[edi+0xc0]		__asm	mulpd		xmm7,[ebx+0xc0]		/* temp*baseinv[i1] */\
				__asm	addpd		xmm3,xmm4			__asm	addpd		xmm7,xmm4					__asm	subpd		xmm3,xmm4			__asm	subpd		xmm7,xmm4			/* cy_out */\
				__asm	movaps		[ecx],xmm3			__asm	movaps		[edx],xmm7			/* store cy_out */\
				/* x = (temp-cy*base[i1])*wt */\
				__asm	mov		eax, __in1\
				__asm	mov		edx, __in3\
				__asm	mulpd		xmm3,[edi+0x80]		__asm	mulpd		xmm7,[ebx+0x80]		/* cy*base[i1] */\
				__asm	subpd		xmm1,xmm3			__asm	subpd		xmm5,xmm7			/* y = (temp-cy*base[i1]) */\
				__asm	mulpd		xmm1,xmm2			__asm	mulpd		xmm5,xmm6			/* y*= wt */\
				__asm	movaps		[eax+0x10],xmm1		__asm	movaps		[edx+0x10],xmm5		/* store y */\
				/* Get ready for next set [RE1~, IM1~] : */\
				__asm	mov		ebx, sse_n\
				__asm	movaps	xmm2,[ebx]\
				__asm	mov		eax, sse_bw\
				__asm	paddd	xmm0,[eax]				/* bjmod[0:3] += bw ; must use packed-INTEGER add [not addpd!] here, severe performance penalty from using addpd. */\
				__asm	movaps	xmm1,xmm0				/* bjmodn COPY */\
				__asm	pcmpgtd	xmm1,xmm2				/* (bj > n)? Gives all 1s in the 32 bit slot of xmm1 if (n < bj) for the bj in that slot of xmm0 */\
				__asm	pand	xmm1,[ebx]				/* n in each slot of xmm1 for which bj += bw overflowed n, 0 otherwise */\
				__asm	psubd	xmm0,xmm1				/* bjmod[0:3] % n */\
				__asm	mov		ecx, __bjmod_0\
				__asm	movaps	[ecx],xmm0				/* Write bjmod[0:3] */\
		/***************Repack the data:*************************/\
				__asm	mov		eax, __in0\
				__asm	mov		ebx, __in1\
				__asm	mov		ecx, __in2\
				__asm	mov		edx, __in3\
				__asm	movaps	xmm1,[eax+0x10]			__asm	movaps	xmm5,[ecx+0x10]			/* reload a[jp+p0 ] */\
				__asm	movaps	xmm0,[eax     ]			__asm	movaps	xmm4,[ecx     ]			/* reload a[jt+p0 ] */\
				__asm	movaps		xmm3,xmm1			__asm	movaps		xmm7,xmm5			/* cpy a[jp    ] */\
				__asm	movaps		xmm2,xmm0			__asm	movaps		xmm6,xmm4			/* cpy a[jt    ] */\
				__asm	unpckhpd	xmm3,[ebx+0x10]		__asm	unpckhpd	xmm7,[edx+0x10]				__asm	unpcklpd	xmm1,[ebx+0x10]		__asm	unpcklpd	xmm5,[edx+0x10]				__asm	movaps	[ebx+0x10],xmm3			__asm	movaps	[edx+0x10],xmm7			/* Store hi imag in aj2 */\
				__asm	unpckhpd	xmm2,[ebx     ]		__asm	unpckhpd	xmm6,[edx     ]				__asm	unpcklpd	xmm0,[ebx     ]		__asm	unpcklpd	xmm4,[edx     ]				__asm	movaps	[ebx     ],xmm2			__asm	movaps	[edx     ],xmm6			/* Store hi real in aj2 */\
																										__asm	movaps	[eax+0x10],xmm1			__asm	movaps	[ecx+0x10],xmm5			/* a[jp+p0 ] */\
				__asm	movaps	[eax     ],xmm0			__asm	movaps	[ecx     ],xmm4			/* a[jt+p0 ] */\
		}

	/******************************************************************************************/
	/******************************************************************************************/
	/******************************************************************************************/

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "carry_gcc32.h"

		#else

			#include "carry_gcc64.h"

		#endif	/* #if(OS_BITS == 32) */

	#endif	/* MSVC or GCC */

#endif	/* USE_SSE2 */

#endif	/* carry_h_included */
