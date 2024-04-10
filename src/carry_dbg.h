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
#ifndef carry_dbg_h_included
#define carry_dbg_h_included

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
#define DBG_fermat_carry_norm_pow2_errcheck(x,y,cx,cy,idx_offset,NRTM1,NRT_BITS)\
{\
	/* Multiply the current transform output by any scale factor: */\
		x *= scale;\
		y *= scale;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		t1=rn0[k1].re;		t2=rn0[k1].im;\
		rt=rn1[k2].re;		it=rn1[k2].im;\
		wt_re =t1*rt-t2*it;	wt_im =t1*it+t2*rt;\
		\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp  = (rt + RND_A) - RND_B;\
		frac = fabs(rt - temp);\
		if(frac > maxerr) maxerr=frac;\
		cx = (temp*baseinv[0] + RND_A) - RND_B;\
		rt = temp-cx*base[0];\
		temp  = (it + RND_A) - RND_B;\
		frac = fabs(it - temp);\
		if(frac > maxerr) maxerr=frac;\
		cy = (temp*baseinv[0] + RND_A) - RND_B;\
		it = temp-cy*base[0];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
}

#define DBG_fermat_carry_norm_pow2_nocheck(x,y,cx,cy,idx_offset,NRTM1,NRT_BITS)\
{\
	/* Multiply the current transform output by any scale factor: */\
		x *= scale;\
		y *= scale;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		t1=rn0[k1].re;		t2=rn0[k1].im;\
		rt=rn1[k2].re;		it=rn1[k2].im;\
		wt_re =t1*rt-t2*it;	wt_im =t1*it+t2*rt;\
		\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp  = (rt + RND_A) - RND_B;\
		cx = (temp*baseinv[0] + RND_A) - RND_B;\
		rt = temp-cx*base[0];\
		temp  = (it + RND_A) - RND_B;\
		cy = (temp*baseinv[0] + RND_A) - RND_B;\
		it = temp-cy*base[0];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
}

/*
Non-power-of-2 runlengths combine the acyclic sincos weights
with the Mersenne-mod-style IBDWT roots-of-2 weights:
*/
#define DBG_fermat_carry_norm_errcheck(x,y,cx,cy,ii,bjmodn,idx_offset,NRTM1,NRT_BITS)\
{\
	/* For Fermat-mod case, combine inverse weight (same for real and imaginary */\
	/* parts of the output) with inverse-FFT scale factor: */\
		wt    =       wt0[ii];\
		wtinv = scale*wt1[ii];\
		ii += SW_DIV_N - nwt;\
		ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
		\
		x *= wtinv;\
		y *= wtinv;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		t1=rn0[k1].re;		t2=rn0[k1].im;\
		rt=rn1[k2].re;		it=rn1[k2].im;\
		wt_re =t1*rt-t2*it;	wt_im =t1*it+t2*rt;\
		\
/*\
if(iter==11 && outer==0 && (j <= 256) && idx_offset==0)\
{\
	printf("j=%3d x,y = %20.3f %20.3f s,c = %18.16f %18.16f\n",j,x,y,wt_re,wt_im);\
}\
*/\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp  = (rt + RND_A) - RND_B;\
		frac = fabs(rt - temp);\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Re(a) at j = %10d\n",frac,j);\
*/\
		if(frac > maxerr) maxerr=frac;\
i = (bjmodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */\
/*if(!idx_offset) printf("k = %8d, bjmodn = %10d, sw = %10d, base[ %1d ]\n",l,bjmodn,sw,i);*/\
bjmodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */\
bjmodn += ( ((int)i-1) & n);		/*       add 0 if a bigword,   N if a smallword */\
		cx = (temp*baseinv[i] + RND_A) - RND_B;\
		rt = temp-cx*base[i];\
		temp  = (it + RND_A) - RND_B;\
		frac = fabs(it - temp);\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Im(a) at j = %10d\n",frac,j);\
*/\
		if(frac > maxerr) maxerr=frac;\
		cy = (temp*baseinv[i] + RND_A) - RND_B;\
		it = temp-cy*base[i];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
	/* Forward IBDWT weight: */\
		x *= wt;\
		y *= wt;\
}

#define DBG_fermat_carry_norm_nocheck(x,y,cx,cy,ii,bjmodn,idx_offset,NRTM1,NRT_BITS)\
{\
	/* For Fermat-mod case, combine inverse weight (same for real and imaginary */\
	/* parts of the output) with inverse-FFT scale factor: */\
		wt    =       wt0[ii];\
		wtinv = scale*wt1[ii];\
		ii += SW_DIV_N - nwt;\
		ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
		\
		x *= wtinv;\
		y *= wtinv;\
	/* Get the needed Nth root of -1: */\
		l = ((j + idx_offset) >> 1);\
		k1=(l & NRTM1);\
		k2=(l >> NRT_BITS);\
		t1=rn0[k1].re;		t2=rn0[k1].im;\
		rt=rn1[k2].re;		it=rn1[k2].im;\
		wt_re =t1*rt-t2*it;	wt_im =t1*it+t2*rt;\
		\
	/* Inverse weight is (wt_re, -wt_im): */\
		rt = x*wt_re + y*wt_im + cx;\
		it = y*wt_re - x*wt_im + cy;\
		temp  = (rt + RND_A) - RND_B;\
i = (bjmodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */\
/*if(!idx_offset) printf("k = %8d, bjmodn = %10d, sw = %10d, base[ %1d ]\n",l,bjmodn,sw,i);*/\
bjmodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */\
bjmodn += ( ((int)i-1) & n);		/*       add 0 if a bigword,   N if a smallword */\
		cx = (temp*baseinv[i] + RND_A) - RND_B;\
		rt = temp-cx*base[i];\
		temp  = (it + RND_A) - RND_B;\
		cy = (temp*baseinv[i] + RND_A) - RND_B;\
		it = temp-cy*base[i];\
	/* Forward weight is (wt_re, +wt_im): */\
		x = rt*wt_re - it*wt_im;\
		y = rt*wt_im + it*wt_re;\
	/* Forward IBDWT weight: */\
		x *= wt;\
		y *= wt;\
}

/*************************************************************/
/**************** MERSENNE-MOD CARRY MACROS ******************/
/*************************************************************/

/* These are versions specialized for power-of-2 runlengths: */

#define DBG_cmplx_carry_norm_pow2_errcheck0(x,y,cy,bjmodn)\
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
/*	\
fprintf(dbg_file, "A: j,col,co2,co3 = %d, %d %d %d\n",j,col,co2,co3);\
fprintf(dbg_file, "A: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "A: x, y, cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
fprintf(dbg_file, "A: wt,inv = %20.15f %20.15f\n",wt,wtinv);\
*/	\
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
/*if(x != 0 && j > len_a) len_a = j;*/\
	\
if(frac > 0.1)\
fprintf(dbg_file, "WARN: frac = %10.8f occurred in Re(a[%2u]) at j = %10d\n",frac,j,0);\
	\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
/*if(y != 0 && j > len_a) len_a = j;*/\
	\
if(frac > 0.1)\
fprintf(dbg_file, "WARN: frac = %10.8f occurred in Im(a[%2u]) at j = %10d\n",frac,j,0);\
	\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
/*	\
fprintf(dbg_file, "A: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "A:~x,~y,~cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
*/	\
}

#define DBG_cmplx_carry_norm_pow2_errcheck(x,y,cy,bjmodn,set)\
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
/*	\
fprintf(dbg_file, "B: j, set = %d, %d\n",j,set);\
fprintf(dbg_file, "B: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "B: x, y, cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
fprintf(dbg_file, "B: wt,inv = %20.15f %20.15f\n",wt,wtinv);\
*/	\
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
/*if(x != 0 && (j+set*n16) > len_a) len_a = (j+set*n16);*/\
	\
if(frac > 0.1)\
fprintf(dbg_file, "WARN: frac = %10.8f occurred in Re(a[%2u]) at j = %10d\n",frac,j,set);\
	\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
/*if(y != 0 && (j+set*n16) > len_a) len_a = (j+set*n16);*/\
	\
if(frac > 0.1)\
fprintf(dbg_file, "WARN: frac = %10.8f occurred in Im(a[%2u]) at j = %10d\n",frac,j,set);\
	\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
/*	\
fprintf(dbg_file, "B: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "B:~x,~y,~cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
*/	\
}

#define DBG_cmplx_carry_norm_pow2_nocheck0(x,y,cy,bjmodn)\
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
/*	\
fprintf(dbg_file, "C: col,co2,co3 = %d %d %d\n",col,co2,co3);\
fprintf(dbg_file, "C: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "C: x, y, cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
fprintf(dbg_file, "C: wt,inv = %20.15f %20.15f\n",wt,wtinv);\
*/	\
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
/*	\
fprintf(dbg_file, "C: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "C:~x,~y,~cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
*/	\
}

#define DBG_cmplx_carry_norm_pow2_nocheck(x,y,cy,bjmodn,set)\
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
/*	\
fprintf(dbg_file, "D: set = %d\n",set);\
fprintf(dbg_file, "D: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "D: x, y, cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
fprintf(dbg_file, "D: wt,inv = %20.15f %20.15f\n",wt,wtinv);\
*/	\
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn = (bjmodn + bw) & nm1;\
/*	\
fprintf(dbg_file, "D: i  ,m  ,m2  = %d %d %d\n",i  ,m  ,m2 );\
fprintf(dbg_file, "D:~x,~y,~cy = %20.15f %20.15f %20.15f\n",x,y,cy);\
*/	\
}


/*
Non-power-of-2 runlengths:
*/
#define DBG_cmplx_carry_norm_errcheck0(x,y,cy,bjmodn)\
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
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Re(a[%2u]) at j = %10d\n",frac,j,0);\
*/\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i]);\
ASSERT(fabs(x+x) <= base[i], "X-output out of range!");\
		x *= wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Im(a[%2u]) at j = %10d\n",frac,j,0);\
*/\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i]);\
ASSERT(fabs(y+y) <= base[i], "Y-output out of range!");\
		y *= wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define DBG_cmplx_carry_norm_errcheck(x,y,cy,bjmodn,set)\
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
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Re(a[%2u]) at j = %10d\n",frac,j,set);\
*/\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i]);\
ASSERT(fabs(x+x) <= base[i], "X-output out of range!");\
		x *= wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
/*\
if(frac > 0.1)\
printf("WARN: frac = %10.8f occurred in Im(a[%2u]) at j = %10d\n",frac,j,set);\
*/\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i]);\
ASSERT(fabs(y+y) <= base[i], "Y-output out of range!");\
		y *= wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

/* A version with a slightly different sequencing, attempting to achieve better pipelining: */
#define DBG_cmplx_carry_norm_errcheckB(x,y,cy,bjmodn,set)\
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
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt =     wt*one_half[m];\
		y  = cy + y*one_half[m2];\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define DBG_cmplx_carry_norm_nocheck0(x,y,cy,bjmodn)\
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
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
	  bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
	  bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define DBG_cmplx_carry_norm_nocheck(x,y,cy,bjmodn,set)\
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
		\
		x = cy+ x*wtinv;\
		temp  = (x + RND_A) - RND_B;\
check_nint(temp, x);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		x = (temp-cy*base[i])*wt;\
		\
		bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
		y = cy+ y*wtinv;\
		temp  = (y + RND_A) - RND_B;\
check_nint(temp, y);\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
check_nint(cy, temp*baseinv[i]);\
		y = (temp-cy*base[i])*wt;\
		\
		bjmodn -= sw;					/* result >= 0 if a bigword, < 0 if a smallword */\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}


/*** Integer-based carry macros (not currently used, as they're slow): ***/

#define DBG_Icmplx_carry_norm_errcheck0(x,y,cy,bjmodn)\
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
		\
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
		temp  = (x + RND_A) - RND_B;\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		x = (temp-cy*base[i])*wt;\
		*/\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
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
		temp  = (y + RND_A) - RND_B;\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		y = (temp-cy*base[i])*wt;\
		*/\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}

#define DBG_Icmplx_carry_norm_errcheck(x,y,cy,bjmodn,set)\
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
		\
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
		temp  = (x + RND_A) - RND_B;\
		frac = fabs(x-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		x = (temp-cy*base[i])*wt;\
		*/\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
		\
		wt   =wtlp1*wtA;\
		wtinv=wtnm1*wtC;\
		i =((uint32)(sw - bjmodn) >> 31);\
		m =((uint32)(n_minus_silp1-bjmodn) >> 31);\
		m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);\
		wt   =wt   *one_half[m];\
		wtinv=wtinv*one_half[m2];\
		\
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
		temp  = (y + RND_A) - RND_B;\
		frac = fabs(y-temp);\
		if(frac > maxerr) maxerr=frac;\
		cy   = (temp*baseinv[i] + RND_A) - RND_B;\
		y = (temp-cy*base[i])*wt;\
		*/\
	  bjmodn -= sw;\
	  bjmodn += ( (-(int)((uint32)bjmodn >> 31)) & n);\
}


#endif	/* carry_dbg_h_included */
