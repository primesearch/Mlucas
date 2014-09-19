/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...The radix-224 DIT pass is here:	*/
	#ifdef USE_SSE2

	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 7 radix-32 transforms...*/
	  #if USE_COMPACT_OBJ_CODE
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			add0 = &a[j1+dit_phi[l]];	itmp = (int *)dit_offsets+(l<<5);
			SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	tmp += 64;
		}
	  #else
		tmp = r00; l = 0;
		add0 = &a[j1    ];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p60];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+pc0];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p40];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+pa0];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p20];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p80];	itmp = (int *)dit_offsets+l;	SSE2_RADIX32_DIT_NOTWIDDLE(add0,itmp, tmp, isrt2);
	  #endif
	/*...and now do 32 radix-7 transforms: */
		tmp = r00;
		for(l = 0; l < 32; l++) {
			// Input-ptrs are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x40;
			va2 = tmp + 0x80;
			va3 = tmp + 0xc0;
			va4 = tmp + 0x100;
			va5 = tmp + 0x140;
			va6 = tmp + 0x180;
			// Output pointers are into s1p** memblock:
			int kk = dit_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of kk; the "which length-13 half of the dit_p20_cperms array?" selector is via (kk < 0):
			int jj = ((-(kk < 0)) & 13)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/13)
						+ (kk & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[jj], k1 = dit_p20_cperms[jj+1], k2 = dit_p20_cperms[jj+2], k3 = dit_p20_cperms[jj+3], k4 = dit_p20_cperms[jj+4], k5 = dit_p20_cperms[jj+5], k6 = dit_p20_cperms[jj+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			kk = (kk & 0x7fffffff) >> 3;
			vb0 = s1p00 + (k0+kk);
			vb1 = s1p00 + (k1+kk);
			vb2 = s1p00 + (k2+kk);
			vb3 = s1p00 + (k3+kk);
			vb4 = s1p00 + (k4+kk);
			vb5 = s1p00 + (k5+kk);
			vb6 = s1p00 + (k6+kk);
			SSE2_RADIX_07_DFT(
				va0,va1,va2,va3,va4,va5,va6,
				dc0,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6
			);	tmp += 2;
		}

	#else	// USE_SSE2 = False:

	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 7 radix-32 transforms...*/
		tptr = t; l = 0;
		jt = j1    ; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+p60; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+pc0; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+p40; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+pa0; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+p20; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);	l += 32;
		jt = j1+p80; RADIX_32_DIT((a+jt),dit_offsets+l,RE_IM_STRIDE, (double *)(tptr+l),t_offsets,1);

	/*...and now do 32 radix-7 transforms: */
		tptr = t;
	  #if USE_COMPACT_OBJ_CODE
		for(l = 0; l < 32; l++) {
			int kk = dit_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of kk; the "which length-13 half of the dit_p20_cperms array?" selector is via (kk < 0):
			int jj = ((-(kk < 0)) & 13)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/13)
						+ (kk & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[jj], k1 = dit_p20_cperms[jj+1], k2 = dit_p20_cperms[jj+2], k3 = dit_p20_cperms[jj+3], k4 = dit_p20_cperms[jj+4], k5 = dit_p20_cperms[jj+5], k6 = dit_p20_cperms[jj+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			kk = (kk & 0x7fffffff) >> 3;
			jt = j1+kk; jp = j2+kk;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		}
	  #else
		jt = j1   ; jp = j2   ;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pf; jp = j2+pf;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pe; jp = j2+pe;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pd; jp = j2+pd;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pc; jp = j2+pc;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pb; jp = j2+pb;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pa; jp = j2+pa;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p9; jp = j2+p9;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p8; jp = j2+p8;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p7; jp = j2+p7;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p6; jp = j2+p6;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p5; jp = j2+p5;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p4; jp = j2+p4;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p3; jp = j2+p3;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p2; jp = j2+p2;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p1; jp = j2+p1;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1   ; jp = j2   ;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p70],a[jp+p70],a[jt+pb0],a[jp+pb0],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p90],a[jp+p90],a[jt+pd0],a[jp+pd0],a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pf; jp = j2+pf;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pe; jp = j2+pe;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pd; jp = j2+pd;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pc; jp = j2+pc;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pb; jp = j2+pb;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+pa; jp = j2+pa;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p9; jp = j2+p9;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p8; jp = j2+p8;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p7; jp = j2+p7;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p6; jp = j2+p6;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p5; jp = j2+p5;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p4; jp = j2+p4;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p3; jp = j2+p3;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p2; jp = j2+p2;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p60],a[jp+p60],a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
		jt = j1+p1; jp = j2+p1;	RADIX_07_DFT(tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+pa0],a[jp+pa0],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p80],a[jp+p80],a[jt+pc0],a[jp+pc0],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	tptr++;
	  #endif

	#endif	// USE_SSE2?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

	/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p00; tmp = cy_r; itmp = bjmodn;
		AVX_cmplx_carry_norm_errcheck0_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		tm1 += 8; tmp += 1; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			AVX_cmplx_carry_norm_errcheck1_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 1; itmp += 4;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

		l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

	/*	i =((uint32)(sw - bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);	tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#else	// Scalar-double mode:

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */		
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p4,p8,...
			cmplx_carry_norm_errcheck(a[jt   ],a[jp   ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p1],a[jp+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p2],a[jp+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
			cmplx_carry_norm_errcheck(a[jt+p3],a[jp+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		}

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}
	else	/* MODULUS_TYPE_FERMAT */
	{

	#ifdef USE_AVX
		int k3,k4,k5,k6,k7;
		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		tmp = base_negacyclic_root;	tm2 = tmp+1;

	  #if HIACC
		// Hi-accuracy version needs RADIX/4 copies of each base root:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		for(i = 0; i < (RADIX << 1); i += 8) {
			VEC_DBL_INIT(tmp+ i,wt_re);	VEC_DBL_INIT(tm2+ i,wt_im);
		}

	  #else	// HIACC = false:

		// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
		// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

		// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
		SSE2_fermat_carry_init_loacc(base_negacyclic_root);

	  #endif

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indices j,j+2,j+4,j+6:
	  #if HIACC

		/* The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls: */

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0x3800;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < isrt2) {	// Can't use l for loop index here since need it for byte offset in carry macro call
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x700, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #else	/* HIACC = false: */

		tm0 = s1p00; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r;
		ic = 0; jc = 1; kc = 2; lc = 3;
		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x700, 0xe0,0x1c0,0x2a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #if (OS_BITS == 64)	// Run out of registers here in serial-build mode, so use (threaded or not?) to toggle carry-macro version selection here:

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < isrt2) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];
			k2 = jcycle[ic];
			int k3 = icycle[jc];
			int k4 = jcycle[jc];
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		l = ODD_RADIX << 4;	// 32-bit version needs preshifted << 4 input value
		while(tm1 < isrt2) {
			//Sep 2014: Even with reduced-register version of the 32-bit Fermat-mod carry macro,
			// GCC runs out of registers on this one, without some playing-around-with-alternate code-sequences ...
			// Pulling the array-refs out of the carry-macro call like so solves the problem:
			k1 = icycle[ic];
			k2 = jcycle[ic];
			SSE2_fermat_carry_norm_errcheck(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2);
			tm1 += 2; tmp++;
			MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

	  #endif

	#else	// Scalar-double mode:

		// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
		ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into icycle mini-array, gets incremented (mod ODD_RADIX) between macro calls
		for(m = 0; m < RADIX>>2; m++) {
			jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p4,p8,...,p56
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}

		icycle[ 0] += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
		icycle[ 1] += wts_idx_incr;
		icycle[ 2] += wts_idx_incr;
		icycle[ 3] += wts_idx_incr;
		icycle[ 4] += wts_idx_incr;
		icycle[ 5] += wts_idx_incr;
		icycle[ 6] += wts_idx_incr;
		icycle[ 0] += ( (-(int)((uint32)icycle[ 0] >> 31)) & nwt);
		icycle[ 1] += ( (-(int)((uint32)icycle[ 1] >> 31)) & nwt);
		icycle[ 2] += ( (-(int)((uint32)icycle[ 2] >> 31)) & nwt);
		icycle[ 3] += ( (-(int)((uint32)icycle[ 3] >> 31)) & nwt);
		icycle[ 4] += ( (-(int)((uint32)icycle[ 4] >> 31)) & nwt);
		icycle[ 5] += ( (-(int)((uint32)icycle[ 5] >> 31)) & nwt);
		icycle[ 6] += ( (-(int)((uint32)icycle[ 6] >> 31)) & nwt);

	#endif	/* #ifdef USE_SSE2 */

	// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
	#ifdef USE_SSE2

		icycle[ 0] += wts_idx_inc2;		icycle[ 0] += ( (-(icycle[ 0] < 0)) & nwt16);
		icycle[ 1] += wts_idx_inc2;		icycle[ 1] += ( (-(icycle[ 1] < 0)) & nwt16);
		icycle[ 2] += wts_idx_inc2;		icycle[ 2] += ( (-(icycle[ 2] < 0)) & nwt16);
		icycle[ 3] += wts_idx_inc2;		icycle[ 3] += ( (-(icycle[ 3] < 0)) & nwt16);
		icycle[ 4] += wts_idx_inc2;		icycle[ 4] += ( (-(icycle[ 4] < 0)) & nwt16);
		icycle[ 5] += wts_idx_inc2;		icycle[ 5] += ( (-(icycle[ 5] < 0)) & nwt16);
		icycle[ 6] += wts_idx_inc2;		icycle[ 6] += ( (-(icycle[ 6] < 0)) & nwt16);

		jcycle[ 0] += wts_idx_inc2;		jcycle[ 0] += ( (-(jcycle[ 0] < 0)) & nwt16);
		jcycle[ 1] += wts_idx_inc2;		jcycle[ 1] += ( (-(jcycle[ 1] < 0)) & nwt16);
		jcycle[ 2] += wts_idx_inc2;		jcycle[ 2] += ( (-(jcycle[ 2] < 0)) & nwt16);
		jcycle[ 3] += wts_idx_inc2;		jcycle[ 3] += ( (-(jcycle[ 3] < 0)) & nwt16);
		jcycle[ 4] += wts_idx_inc2;		jcycle[ 4] += ( (-(jcycle[ 4] < 0)) & nwt16);
		jcycle[ 5] += wts_idx_inc2;		jcycle[ 5] += ( (-(jcycle[ 5] < 0)) & nwt16);
		jcycle[ 6] += wts_idx_inc2;		jcycle[ 6] += ( (-(jcycle[ 6] < 0)) & nwt16);

	  #ifdef USE_AVX
		kcycle[ 0] += wts_idx_inc2;		kcycle[ 0] += ( (-(kcycle[ 0] < 0)) & nwt16);
		kcycle[ 1] += wts_idx_inc2;		kcycle[ 1] += ( (-(kcycle[ 1] < 0)) & nwt16);
		kcycle[ 2] += wts_idx_inc2;		kcycle[ 2] += ( (-(kcycle[ 2] < 0)) & nwt16);
		kcycle[ 3] += wts_idx_inc2;		kcycle[ 3] += ( (-(kcycle[ 3] < 0)) & nwt16);
		kcycle[ 4] += wts_idx_inc2;		kcycle[ 4] += ( (-(kcycle[ 4] < 0)) & nwt16);
		kcycle[ 5] += wts_idx_inc2;		kcycle[ 5] += ( (-(kcycle[ 5] < 0)) & nwt16);
		kcycle[ 6] += wts_idx_inc2;		kcycle[ 6] += ( (-(kcycle[ 6] < 0)) & nwt16);

		lcycle[ 0] += wts_idx_inc2;		lcycle[ 0] += ( (-(lcycle[ 0] < 0)) & nwt16);
		lcycle[ 1] += wts_idx_inc2;		lcycle[ 1] += ( (-(lcycle[ 1] < 0)) & nwt16);
		lcycle[ 2] += wts_idx_inc2;		lcycle[ 2] += ( (-(lcycle[ 2] < 0)) & nwt16);
		lcycle[ 3] += wts_idx_inc2;		lcycle[ 3] += ( (-(lcycle[ 3] < 0)) & nwt16);
		lcycle[ 4] += wts_idx_inc2;		lcycle[ 4] += ( (-(lcycle[ 4] < 0)) & nwt16);
		lcycle[ 5] += wts_idx_inc2;		lcycle[ 5] += ( (-(lcycle[ 5] < 0)) & nwt16);
		lcycle[ 6] += wts_idx_inc2;		lcycle[ 6] += ( (-(lcycle[ 6] < 0)) & nwt16);
	  #endif
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-56 DIF pass is here:	*/

	#ifdef USE_SSE2

	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 32 radix-7 transforms...*/
		tmp = r00;
		for(l = 0; l < 32; l++) {
			// Input pointers are into s1p** memblock:
			int kk = dif_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of kk; the "which length-13 half of the dif_p20_cperms array?" selector is via (kk < 0):
			int jj = ((-(kk < 0)) & 13)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/13)
						+ (kk & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[jj], k1 = dif_p20_cperms[jj+1], k2 = dif_p20_cperms[jj+2], k3 = dif_p20_cperms[jj+3], k4 = dif_p20_cperms[jj+4], k5 = dif_p20_cperms[jj+5], k6 = dif_p20_cperms[jj+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			kk = (kk & 0x7fffffff) >> 3;
			vb0 = s1p00 + (k0+kk);
			vb1 = s1p00 + (k1+kk);
			vb2 = s1p00 + (k2+kk);
			vb3 = s1p00 + (k3+kk);
			vb4 = s1p00 + (k4+kk);
			vb5 = s1p00 + (k5+kk);
			vb6 = s1p00 + (k6+kk);
			// Output-ptrs [va/vb swap roles here vs DIT] are regular-stride offsets of r00:
			va0 = tmp;
			va1 = tmp + 0x40;
			va2 = tmp + 0x80;
			va3 = tmp + 0xc0;
			va4 = tmp + 0x100;
			va5 = tmp + 0x140;
			va6 = tmp + 0x180;
			SSE2_RADIX_07_DFT(
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,
				dc0,
				va0,va1,va2,va3,va4,va5,va6
			);	tmp += 2;
		}
	/*...and now do 7 radix-32 transforms: */
	  #if USE_COMPACT_OBJ_CODE
		tmp = r00;
		for(l = 0; l < ODD_RADIX; l++) {
			add0 = &a[j1+dif_phi[l]];	itmp = (int *)dif_offsets+(l<<5);
			SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	tmp += 64;
		}
	  #else
		tmp = r00; l = 0;
		add0 = &a[j1    ];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+pc0];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+pa0];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p80];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p60];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p40];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);	l += 32; tmp += 64;
		add0 = &a[j1+p20];	itmp = (int *)dif_offsets+l;	SSE2_RADIX32_DIF_NOTWIDDLE(add0,itmp, tmp, isrt2);
	  #endif

	#else	// USE_SSE2 = False:

	/*...gather the needed data (224 64-bit complex, i.e. 448 64-bit reals) and do 32 radix-7 transforms: */
		tptr = t;
	  #if USE_COMPACT_OBJ_CODE
		// Cf. USE_COMPACT_OBJ_CODE bulld path in radix224_dif_pass1():
		for(l = 0; l < 32; l++) {
			int kk = dif_p20_lo_offset[l];
			// Extract index (in [0-6]) into circ-shift array used for high parts of p-mults. The [0-6] value is
			// in low 3 bits of kk; the "which length-13 half of the dif_p20_cperms array?" selector is via (kk < 0):
			int jj = ((-(kk < 0)) & 13)	// +/- sign on kk puts us into lower/upper half of the cshift array (base index 0/13)
						+ (kk & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[jj], k1 = dif_p20_cperms[jj+1], k2 = dif_p20_cperms[jj+2], k3 = dif_p20_cperms[jj+3], k4 = dif_p20_cperms[jj+4], k5 = dif_p20_cperms[jj+5], k6 = dif_p20_cperms[jj+6];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			kk = (kk & 0x7fffffff) >> 3;
			jt = j1+kk; jp = j2+kk;
			RADIX_07_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im,
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++;
		}
	  #else
		jt = j1   ; jp = j2   ;	RADIX_07_DFT(a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p9; jp = j2+p9;	RADIX_07_DFT(a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p2; jp = j2+p2;	RADIX_07_DFT(a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pb; jp = j2+pb;	RADIX_07_DFT(a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p4; jp = j2+p4;	RADIX_07_DFT(a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pd; jp = j2+pd;	RADIX_07_DFT(a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p6; jp = j2+p6;	RADIX_07_DFT(a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pf; jp = j2+pf;	RADIX_07_DFT(a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p8; jp = j2+p8;	RADIX_07_DFT(a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p1; jp = j2+p1;	RADIX_07_DFT(a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pa; jp = j2+pa;	RADIX_07_DFT(a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p3; jp = j2+p3;	RADIX_07_DFT(a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pc; jp = j2+pc;	RADIX_07_DFT(a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p5; jp = j2+p5;	RADIX_07_DFT(a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pe; jp = j2+pe;	RADIX_07_DFT(a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p7; jp = j2+p7;	RADIX_07_DFT(a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1   ; jp = j2   ;	RADIX_07_DFT(a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p9; jp = j2+p9;	RADIX_07_DFT(a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p2; jp = j2+p2;	RADIX_07_DFT(a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pb; jp = j2+pb;	RADIX_07_DFT(a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p4; jp = j2+p4;	RADIX_07_DFT(a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pd; jp = j2+pd;	RADIX_07_DFT(a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p6; jp = j2+p6;	RADIX_07_DFT(a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pf; jp = j2+pf;	RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p8; jp = j2+p8;	RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p1; jp = j2+p1;	RADIX_07_DFT(a[jt+p30],a[jp+p30],a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pa; jp = j2+pa;	RADIX_07_DFT(a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p3; jp = j2+p3;	RADIX_07_DFT(a[jt+p20],a[jp+p20],a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pc; jp = j2+pc;	RADIX_07_DFT(a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p5; jp = j2+p5;	RADIX_07_DFT(a[jt+p10],a[jp+p10],a[jt+pd0],a[jp+pd0],a[jt+pb0],a[jp+pb0],a[jt+p90],a[jp+p90],a[jt+p70],a[jp+p70],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+pe; jp = j2+pe;	RADIX_07_DFT(a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		jt = j1+p7; jp = j2+p7;	RADIX_07_DFT(a[jt    ],a[jp    ],a[jt+pc0],a[jp+pc0],a[jt+pa0],a[jp+pa0],a[jt+p80],a[jp+p80],a[jt+p60],a[jp+p60],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+32)->re,(tptr+32)->im,(tptr+64)->re,(tptr+64)->im,(tptr+96)->re,(tptr+96)->im,(tptr+128)->re,(tptr+128)->im,(tptr+160)->re,(tptr+160)->im,(tptr+192)->re,(tptr+192)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
	  #endif
	/*...and now do 7 radix-32 transforms: */
		tptr = t; l = 0;
		jt = j1    ; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+pc0; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+pa0; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+p80; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+p60; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+p40; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);	l += 32;
		jt = j1+p20; RADIX_32_DIF((double *)(tptr+l),t_offsets,1, (a+jt),dif_offsets+l,RE_IM_STRIDE);

	#endif	/* !USE_SSE2 */
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */
