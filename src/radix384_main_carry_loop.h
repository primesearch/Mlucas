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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(int k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
	#endif

	/*...The radix-384 DIT pass is here:	*/

	#ifdef USE_SSE2

		// Use s1p000-07f for scratch for these 3 DFTs:
		SSE2_RADIX128_DIT(
			(a+j1     ),(int *)(dit_i_offsets+0x000),
			s1p000, isrt2,two,twid0,
			r000
		);	// Outputs in r0[ 00- 7f]
		SSE2_RADIX128_DIT(
			(a+j1+p080),(int *)(dit_i_offsets+0x080),
			s1p000, isrt2,two,twid0,
			r080
		);	// Outputs in r1[ 80- ff]
		SSE2_RADIX128_DIT(
			(a+j1+p100),(int *)(dit_i_offsets+0x100),
			s1p000, isrt2,two,twid0,
			r100
		);	// Outputs in r2[100-17f]
/******************* AVX debug stuff: *******************/
#if 0
if(tid == 1) {
	printf("DIT-intermediates:\n");
	tmp = r000;
	for(l = 0; l < 2*RADIX; l += 2) {
		printf("%3d: %20.10e,%20.10e\n",l>>1,(tmp+l)->d0,(tmp+l+1)->d0);
	}
	exit(0);	// If printing the above inputs, exit immediately.
}
#endif
/********************************************************/

	/*...and now do 256 radix-3 transforms: */
	  #if 0//OS_BITS == 64	// Dec 2014: Data-doubled-radix-3 versions slightly faster if enough register available

		vec_dbl *vb0,*vb1,*vb2, *vc0,*vc1,*vc2, *vd0,*vd1,*vd2;
		// Data-doubled-macro-call version of 32-bit loop below:
		l = 2; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	// Skip 0-term, which gets saved for wraparound
		l = 4; vc0 = r000+l; vc1 = vc0+512; vc2 = vc0+1024;	// 2nd set of input-ptrs
		l -= 2;	// Need l = 4 for init of vc-ptrs above, but need l = 2 to match loop control of template 32-bit loop
		for(l1 = 0; l1 < 48; l1 += 3) {
			k0 = dit_triplets[l1]; k1 = dit_triplets[l1+1]; k2 = dit_triplets[l1+2];
			l2 = 0x1e;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x1a;	vb0 = s1p000+l2; vb1 = vb0+k0; vb2 = vb0+k1; vb0 += k2;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x16;	vb0 = s1p000+l2; vb1 = vb0+k2; vb2 = vb0+k0; vb0 += k1;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x12;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x0e;	vb0 = s1p000+l2; vb1 = vb0+k0; vb2 = vb0+k1; vb0 += k2;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x0a;	vb0 = s1p000+l2; vb1 = vb0+k2; vb2 = vb0+k0; vb0 += k1;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k0; vd2 = vd0+k1; vd0 += k2;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x06;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k2; vd2 = vd0+k0; vd0 += k1;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4;             tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	/* needed for final-loop-pass wraparound of 0-term [vc-ptrs in next macro call] */
			l+=2; l &= 0x1ff; vc0 = r000+l; vc1 = vc0+512; vc2 = vc0+1024;
			l-=2;
			l2 = 0x02;	vb0 = s1p000+l2; vb1 = vb0+k0; vb2 = vb0+k1; vb0 += k2;
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;

			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */tmp,tm1,tm2, vb0,vb1,vb2, /* DFT #2: */vc0,vc1,vc2, vd0,vd1,vd2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;
		}

	  #else	// 32-bit:

		vec_dbl *vd0,*vd1,*vd2;
		// Loop-based compact-obj-code impl exploits above index pattern to group interior [sandwiched
		// between single leading and 15 trailing] macro calls into sets of 16 with neatly cutoff index groupings:
		l = 2; tmp = r000+l; tm1 = tmp+256; tm2 = tmp+512;	// Skip 0-term, which gets saved for wraparound
		for(l1 = 0; l1 < 24; l1 += 3) {
			k0 = dit_triplets[l1]; k1 = dit_triplets[l1+1]; k2 = dit_triplets[l1+2];	// vd*-ordering here mimics the a[]-index-offsets in the scalar version:
																						//							vvvvvvvvvvv
			l2 = 0x1e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd2,vd1); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x1a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd2,vd1,vd0); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x16;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd2,vd1); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd2,vd1,vd0); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x12;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd2,vd1); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd2,vd1,vd0); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd2,vd1); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd2,vd1,vd0); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x06;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd0,vd2,vd1); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x02;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd2,vd1,vd0); l+=2; l &= 0xff; tmp = r000+l; tm1 = tmp+256; tm2 = tmp+512;	/* needed for final-loop-pass wraparound of 0-term */
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;      SSE2_RADIX_03_DFT(tmp,tm1,tm2, cc1, vd1,vd0,vd2); l+=2; tmp+=2; tm1+=2; tm2+=2;
		}

	  #endif	// 32 / 64-bit?

	#else	/* !USE_SSE2 */

	/*...gather the needed data (384 64-bit complex) and do 3 radix-128 transforms,	*/
		RADIX_128_DIT((a+j1     ),dit_i_offsets      ,RE_IM_STRIDE, (double *)(t+0x000),dit_o_offsets,1);	// Outputs in t[ 00- 7f]
		RADIX_128_DIT((a+j1+p080),dit_i_offsets+0x080,RE_IM_STRIDE, (double *)(t+0x080),dit_o_offsets,1);	// Outputs in t[ 80- ff]
		RADIX_128_DIT((a+j1+p100),dit_i_offsets+0x100,RE_IM_STRIDE, (double *)(t+0x100),dit_o_offsets,1);	// Outputs in t[100-17f]
/******************* AVX debug stuff: *******************/
#if 0
	printf("DIT-intermediates:\n");
	tptr = t;
	for(l = 0; l < RADIX; l++, tptr++) {
		printf("%3d: %20.10e,%20.10e\n",l  ,tptr->re,tptr->im);
	}
	exit(0);	// If printing the above inputs, exit immediately.
#endif
/********************************************************/

	/*...and now do 128 radix-3 transforms: */
		l = 1; l1 = l+128; l2 = l+256;	// Skip 0-term, which gets saved for wraparound
		for(m = 0; m < 24; m += 3) {
			k0 = dit_triplets[m]; k1 = dit_triplets[m+1]; k2 = dit_triplets[m+2];
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1]); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k2],a[jp+k2],a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0]); ++l; l &= 0x07f; l1 = l+128; l2 = l+256;	// <*** needed for final-loop-pass wraparound of 0-term
			jt = j1     ; jp = j2     ;	RADIX_03_DFT(s,c3m1, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im, t00,t01,t02,t03,t04,t05, a[jt+k1],a[jp+k1],a[jt+k0],a[jp+k0],a[jt+k2],a[jp+k2]); ++l; ++l1; ++l2;
		}

	#endif

/******************* AVX debug stuff: *******************/
#if 0
if(tid == 1) {
	printf("DIT-outputs/CY-inputs:\n");
  #ifdef USE_SSE2
	tmp = s1p000;
  #else
	int ipad;
  #endif
	for(l = 0; l < RADIX; l+=4) {
	#ifdef USE_SSE2
		printf("A[%3d] = %20.10e,%20.10e\n",l  ,(tmp+2*l  )->d0,(tmp+2*l+1)->d0);
		printf("A[%3d] = %20.10e,%20.10e\n",l+1,(tmp+2*l+2)->d0,(tmp+2*l+3)->d0);
		printf("A[%3d] = %20.10e,%20.10e\n",l+2,(tmp+2*l+4)->d0,(tmp+2*l+5)->d0);
		printf("A[%3d] = %20.10e,%20.10e\n",l+3,(tmp+2*l+6)->d0,(tmp+2*l+7)->d0);
	#else
		ipad = poff[l>>2];	// Get padded p04-multiple for current print-quartet
		printf("A[%3d] = %20.10e,%20.10e\n",l  ,a[ipad   ],a[ipad   +1]);
		printf("A[%3d] = %20.10e,%20.10e\n",l+1,a[ipad+p1],a[ipad+p1+1]);
		printf("A[%3d] = %20.10e,%20.10e\n",l+2,a[ipad+p2],a[ipad+p2+1]);
		printf("A[%3d] = %20.10e,%20.10e\n",l+3,a[ipad+p3],a[ipad+p3+1]);
	#endif
	}
	exit(0);	// If printing the above inputs, exit immediately.
}
#endif
/********************************************************/

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 48 separate blocks of the A-array, we need 48 separate carries.	*/

		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			double *addr_ = (double *)s1p000 + target_set;
			*addr_ += target_cy*(n>>1);	// target_cy = [-2 << within-word-shift]*[DWT weight]*n/2, i.e. includes fwd DWT weight and n/2 factor
		#else
			// target_set in [0,2*RADIX); tidx_mod_stride [even|odd] means shifted-carry goes into [Re|Im] part of the complex FFT datum:
			l = target_set&1;	target_set >>= 1;
			a[j1+poff[target_set>>2]+p0123[target_set&3]+l] += target_cy*(n>>1);
		#endif
			target_idx = -1;
		}

	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];
		/* ptr to local storage for the doubled wtl,wtn terms: */
	  #ifdef USE_AVX512
		tmp = half_arr +  64;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling;
								// 1st 64 slots hold outputs of wtsinit call. Only half of said slots used in 8-way-init mode.
	  #else
		tmp = half_arr + 128;	// 1st 64 slots are basic-4 LUTs, next 32 are the additional 2 LOACC LUTs, next 32 hold outputs of wtsinit call
	  #endif
		l= j & (nwt-1);						// These rcol wts-terms are for individual-double-broadcast-to-full-vector-width,
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];	// hence the mixing of fwd/inv wts, which is normally taboo.
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;
		n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;
	  #ifdef USE_AVX512
		l= (j+8) & (nwt-1);					tmp -= 3;	// Reset to same tmp-startval as above, now copy data into d4-7 slots of vec_dbl
		n_minus_sil  ->d4 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d4 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d4 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d4 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+10) & (nwt-1);				++tmp;
		n_minus_sil  ->d5 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d5 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d5 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d5 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+12) & (nwt-1);				++tmp;
		n_minus_sil  ->d6 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d6 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d6 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d6 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+14) & (nwt-1);				++tmp;
		n_minus_sil  ->d7 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d7 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d7 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d7 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;
	  #endif

	  if(incr) {	// Have no specialized HIACC carry macro in AVX-512, so use 0-or-not-ness of the low
						// inc_arr element to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 ii,loop, co2save = co2;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every 8th pass or better:
		// incr must divide nloop [RADIX/8 = 96 or RADIX/16 = 48, depending on whether we use 8-or-16-way carry macros]!
	  #ifdef CARRY_16_WAY
		const uint32 nloop = RADIX>>4;
	  #else
		const uint32 nloop = RADIX>>3;
	  #endif
		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tmp = s1p000; tm1 = cy; itmp = bjmodn;
	  #ifndef USE_AVX512
		tm2 = cy+1; itm2 = bjmodn+4;	// tm2,itm2 not used in AVX-512 mode
	  #endif
		for(loop = 0; loop < nloop; loop += incr)
		{
			co2 = co2save;	// Need this for all wts-inits beynd the initial set, due to the co2 = co3 preceding the (j+2) data
		  #ifdef CARRY_16_WAY
			ii = loop << 4;	// Reflects 16 independent carry chains being done in each AVX_cmplx_carry_fast_errcheck_X8 call
		  #else
			ii = loop << 3;	// Reflects  8 independent carry chains being done in each AVX_cmplx_carry_fast_errcheck_X8 call
		  #endif
			add1 = &wt1[col  +ii];	/* Don't use add0 here, to avoid need to reload main-array address */
			add2 = &wt1[co2-1-ii];
			add3 = &wt1[co3-1-ii];

			// Since use wt1-array in the wtsinit macro, need to fiddle this here:
			co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
			// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
		  #ifdef CARRY_16_WAY
			AVX_cmplx_carry_fast_wtsinit_X16(add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		  #else
			AVX_cmplx_carry_fast_wtsinit_X8 (add1,add2,add3, itmp, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_n)
		  #endif
			for(l = loop; l < loop+incr; l++) {
				// Each AVX carry macro call also processes 8 prefetches of main-array data
				add0 = a + j1 + pfetch_dist + poff[l+l];
			  // In AVX-512 mode, the 4 doubles base[0],baseinv[1],wts_mult[1],inv_mult[0] are in the d0-3 slots of the otherwise-unused sse2_rnd vec_dbl:
			  #ifdef USE_AVX512
			   #ifdef CARRY_16_WAY
				AVX_cmplx_carry_fast_errcheck_X16(tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 32; tm1 += 2;           itmp += 16;           i = 0;
			   #else
				AVX_cmplx_carry_fast_errcheck_X8 (tmp, tm1    , itmp     , half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 1;           itmp +=  8;           i = 0;	// CY-ptr only advances 1 in AVX-512/CARRY_8_WAY mode, since all 8 dbl-carries fit in a single vec_dbl
			   #endif
			  #else	// USE_AVX:
				AVX_cmplx_carry_fast_errcheck_X8(tmp, tm1,tm2, itmp,itm2, half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3,p4, addr);
				tmp += 16; tm1 += 2; tm2 += 2; itmp += 8; itm2 += 8; i = 0;
			  #endif
			}
		}

	  } else {	// HiACC:

		/* In AVX mode advance carry-ptrs just 1 for each vector-carry-macro call: */
		tm1 = s1p000; tmp = cy; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each AVX carry macro call also processes 4 prefetches of main-array data
			tm2 = (vec_dbl *)(a + j1 + pfetch_dist + poff[l]);
			AVX_cmplx_carry_norm_errcheck_X4(tm1,add1,add2,add3,tmp,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, tm2,p1,p2,p3, addr);
			tm1 += 8; tmp += 1; itmp += 4; i = 0;
		}

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

	  if(incr) {	// Have no specialized HIACC carry macro in ARM_V8_SIMD, so use 0-or-not-ness of incr
						// in lieu of (USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) is for ARMv8
						// to divert non-AVX512 builds to 'else' clause of if() in HIACC mode.
		uint32 i0,i1,i2,i3, ii,nwtml, loop,nloop = RADIX>>2, co2save = co2;

		i = (!j);	// Need this to force 0-wod to be bigword
		addr = &prp_mult;
		tm1 = s1p000; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		// Beyond chain length 8, the chained-weights scheme becomes too inaccurate, so re-init seed-wts every few passes:
		for(loop = 0; loop < nloop; loop += incr)
		{
			ii = loop << 2;	// Reflects 4 independent carry chains being done in each SSE2_cmplx_carry_fast_pow2_errcheck call
			/*** wt_re,wi_re,wt_im,wi_im inits. Cf. radix16_main_carry_loop.h for scalar-macro prototyping of this: ***/
			l = j & (nwt-1);	nwtml = nwt-l;
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwtml  ];
			sinwtm1 = si[nwtml-1];
			wtl     = wt0[    l  ];
			wtn     = wt0[nwtml  ]*scale;
			wtlp1   = wt0[    l+1];
			wtnm1   = wt0[nwtml-1]*scale;

			co2 = co2save;	// Need this for all wts-inits beynd the initial set, due to the co2 = co3 preceding the (j+2) data
			ctmp = (struct complex *)half_arr + 24;	// ptr to local storage for the doubled wtl,wtn terms:
			// (j)-data occupy the 8 xmm-sized slots above the 16 used by fixed auxiliary-data, and overwrite these inits:
			ctmp->re = ctmp->im = wtl;		ctmp += 2;
			ctmp->re = ctmp->im = wtn;		ctmp += 2;
			ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
			ctmp->re = ctmp->im = wtnm1;

			l = (j+2) & (nwt-1);	nwtml = nwt-l;;
			i0 = n-si[l  ];
			i1 = n-si[l+1];
			i2 = si[nwtml  ];
			i3 = si[nwtml-1];
			wtl     = wt0[    l  ];
			wtn     = wt0[nwtml  ]*scale;
			wtlp1   = wt0[    l+1];
			wtnm1   = wt0[nwtml-1]*scale;

			ctmp = (struct complex *)half_arr + 32;	// (j+2) data start at ctmp + 8
			ctmp->re = ctmp->im = wtl;		ctmp += 2;
			ctmp->re = ctmp->im = wtn;		ctmp += 2;
			ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
			ctmp->re = ctmp->im = wtnm1;

			add1 = &wt1[col  +ii];	/* Don't use add0 here, to avoid need to reload main-array address */
			add2 = &wt1[co2-1-ii];
			add3 = &wt1[co3-1-ii];

			// Since use wt1-array in the wtsinit macro, need to fiddle this here:
			co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
			// *But*: since the init macro does an on-the-fly version of this between j,j+2 portions, external code co2=co3 must come *after* both ctmp-data octets are inited.
			add0 = (double*)(bjmodn+ii);
			SSE2_cmplx_carry_fast_wtsinit(add1,add2,add3, add0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, i0,i1,i2,i3, sse_bw,sse_n)

			for(l = loop; l < loop+incr; l++) {
				// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
				add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
				SSE2_cmplx_carry_fast_errcheck(tm1,tmp,tm2,itmp,half_arr,i,sign_mask,sse_bw,sse_n,sse_sw, add0,p1,p2,p3, addr);
				tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
			}
		}

	  } else {	// HiACC:

		l= j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

	/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		tm1 = s1p000; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		addr = &prp_mult;
		i = (!j);
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p1, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4; i = 0;
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

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p000; tmp = cy; tm2 = cy+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			// Each SSE2 carry macro call also processes 2 prefetches of main-array data
			add0 = a + j1 + pfetch_dist + poff[l];	// poff[] = p0,4,8,...
			add0 += (-(l&0x1)) & p2;	// Base-addr incr by extra p2 on odd-index passes
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw, add0,p2,p3, addr);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		}

	  }	// LOACC or HIACC?

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

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; itmp = bjmodn;
		double *addr_ = cy;
		for(l1 = 0; l1 < RADIX>>2; l1++) {
			jt = j1 + poff[l1]; jp = j2 + poff[l1];	// poff[] = p04,08,...
			// Re-init weights every 4th macro invocation to keep errors under control:
			cmplx_carry_norm_errcheck0(a[jt   ],a[jp   ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p1],a[jp+p1],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p2],a[jp+p2],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_fast_errcheck (a[jt+p3],a[jp+p3],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
		}

	  } else {	// HiACC:

		/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
		l = 0; itmp = bjmodn;
		double *addr_ = cy;
		for(l1 = 0; l1 < RADIX>>2; l1++) {
			jt = j1 + poff[l1]; jp = j2 + poff[l1];	// poff[] = p04,08,...
			cmplx_carry_norm_errcheck0(a[jt   ],a[jp   ],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p1],a[jp+p1],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p2],a[jp+p2],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
			cmplx_carry_norm_errcheck (a[jt+p3],a[jp+p3],*addr_,*itmp,l,prp_mult); ++l; ++addr_; ++itmp;
		}

	  }	// LOACC or HIACC?

		i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	/*...The radix-384 DIF pass is here:	*/

	#ifdef USE_SSE2

	/*...gather the needed data (384 64-bit complex) and do 256 radix-3 transforms - We want unit-strides in the radix384-DFT macro, so use large output strides here: */
	  #if 0//OS_BITS == 64

		// Data-doubled-macro-call version of 32-bit loop below:
		l = 2; tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	// Skip 0-term, which gets saved for wraparound
		l = 4; vc0 = r000+l; vc1 = vc0+512; vc2 = vc0+1024;	// 2nd set of input-ptrs
		l -= 2;	// Need l = 4 for init of vc-ptrs above, but need l = 2 to match loop control of template 32-bit loop
		for(l1 = 0; l1 < 144; l1 += 9) {
			k0 = dif_triplets[l1  ]; k1 = dif_triplets[l1+1]; k2 = dif_triplets[l1+2];
			l2 = 0x1a;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x0e;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x02;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			k0 = dif_triplets[l1+3]; k1 = dif_triplets[l1+4]; k2 = dif_triplets[l1+5];
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x16;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x0a;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			k0 = dif_triplets[l1+6]; k1 = dif_triplets[l1+7]; k2 = dif_triplets[l1+8];
			l2 = 0x1e;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;

			l2 = 0x12;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4;             tmp = r000+l; tm1 = tmp+512; tm2 = tmp+1024;	/* needed for final-loop-pass wraparound of 0-term [vc-ptrs in next macro call] */
			l+=2; l &= 0x1ff; vc0 = r000+l; vc1 = vc0+512; vc2 = vc0+1024;
			l-=2;
			l2 = 0x06;	vb0 = s1p000+l2; vb1 = vb0+k1; vb2 = vb0+k2; vb0 += k0;
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;
			SSE2_RADIX_03_DFT_X2(cc1, /* DFT #1: */vb0,vb1,vb2, tmp,tm1,tm2, /* DFT #2: */vd0,vd1,vd2, vc0,vc1,vc2)
			l+=4; tmp+=4; tm1+=4; tm2+=4;	vc0+=4; vc1+=4; vc2+=4;
		}

	  #else	// 32-bit:

		l = 2; tmp = r000+l; tm1 = tmp+256; tm2 = tmp+512;	// Skip 0-term, which gets saved for wraparound
		for(l1 = 0; l1 < 72; l1 += 9) {
			k0 = dif_triplets[l1  ]; k1 = dif_triplets[l1+1]; k2 = dif_triplets[l1+2];
			l2 = 0x1a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x14;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x08;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x02;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			k0 = dif_triplets[l1+3]; k1 = dif_triplets[l1+4]; k2 = dif_triplets[l1+5];
			l2 = 0x1c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x16;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x10;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0a;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x04;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			k0 = dif_triplets[l1+6]; k1 = dif_triplets[l1+7]; k2 = dif_triplets[l1+8];
			l2 = 0x1e;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x18;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x12;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x0c;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
			l2 = 0x06;	vd0 = s1p000+l2; vd1 = vd0+k1; vd2 = vd0+k2; vd0 += k0; SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; l &= 0xff; tmp = r000+l; tm1 = tmp+256; tm2 = tmp+512;	// <*** needed for final-loop-pass wraparound of 0-term
						vd0 = s1p000+k0; vd1 = s1p000+k1; vd2 = s1p000+k2;      SSE2_RADIX_03_DFT(vd0,vd1,vd2, cc1, tmp,tm1,tm2); l+=2; tmp+=2; tm1+=2; tm2+=2;
		}

	  #endif	// 32 / 64-bit?

	/*...and now do 3 radix-128 transforms, using s1p-storage for intermediates: */

/******************* AVX debug stuff: *******************/
#if 0
	printf("DIF-intermediates:\n");
	tmp = r000;
	for(l = 0; l < 2*RADIX; l += 2) {
		printf("%3d: %20.10e,%20.10e\n",l>>1,(tmp+l)->d0,(tmp+l+1)->d0);
	}
	exit(0);	// If printing the above inputs, exit immediately.
#endif
/********************************************************/
		SSE2_RADIX128_DIF(	// Inputs in r0[ 00- 7f]
			r000,
			s1p000, isrt2,two,twid0,
			(a+j1     ),dif_o_offsets
		);
		SSE2_RADIX128_DIF(	// Inputs in r1[ 80- ff]
			r080,
			s1p000, isrt2,two,twid0,
			(a+j1+p100),dif_o_offsets+0x080
		);
		SSE2_RADIX128_DIF(	// Inputs in r2[100-17f]
			r100,
			s1p000, isrt2,two,twid0,
			(a+j1+p080),dif_o_offsets+0x100
		);

	#else	/* !USE_SSE2 */

	/*...gather the needed data (384 64-bit complex) and do 128 radix-3 transforms - We want unit-strides in the radix384-DFT macro, so use large output strides here: */
		l = 1; l1 = l+128; l2 = l+256;	// Skip 0-term, which gets saved for wraparound
		for(m = 0; m < 72; m += 3) {
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2]; m += 3;
			jt = j1 + pd; jp = j2 + pd;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pa; jp = j2 + pa;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p7; jp = j2 + p7;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p4; jp = j2 + p4;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p1; jp = j2 + p1;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2]; m += 3;
			jt = j1 + pe; jp = j2 + pe;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pb; jp = j2 + pb;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p8; jp = j2 + p8;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p5; jp = j2 + p5;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p2; jp = j2 + p2;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			k0 = dif_triplets[m]; k1 = dif_triplets[m+1]; k2 = dif_triplets[m+2];	// Skip k-incr here since loop control handles this one
			jt = j1 + pf; jp = j2 + pf;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + pc; jp = j2 + pc;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p9; jp = j2 + p9;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p6; jp = j2 + p6;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
			jt = j1 + p3; jp = j2 + p3;	RADIX_03_DFT(s,c3m1, a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; l &= 0x7f; l1 = l+128; l2 = l+256;	// <*** needed for final-loop-pass wraparound of 0-term
										RADIX_03_DFT(s,c3m1, a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2], t00,t01,t02,t03,t04,t05, t[l].re,t[l].im,t[l1].re,t[l1].im,t[l2].re,t[l2].im); ++l; ++l1; ++l2;
		}
/******************* AVX debug stuff: *******************/
#if 0
	printf("DIF-intermediates:\n");
	tptr = t;
	for(l = 0; l < RADIX; l++, tptr++) {
		printf("%3d: %20.10e,%20.10e\n",l  ,tptr->re,tptr->im);
	}
	exit(0);	// If printing the above inputs, exit immediately.
#endif
/********************************************************/
	/*...and now do 3 radix-128 transforms: */
		RADIX_128_DIF((double *)(t+0x000),dif_i_offsets,1, (a+j1     ),dif_o_offsets      ,RE_IM_STRIDE);	// Inputs in t[ 00- 7f]
		RADIX_128_DIF((double *)(t+0x080),dif_i_offsets,1, (a+j1+p100),dif_o_offsets+0x080,RE_IM_STRIDE);	// Inputs in t[ 80- ff]
		RADIX_128_DIF((double *)(t+0x100),dif_i_offsets,1, (a+j1+p080),dif_o_offsets+0x100,RE_IM_STRIDE);	// Inputs in t[100-17f]

	#endif
	}

	jstart += nwt;
	jhi    += nwt;

	col += RADIX;
	co3 -= RADIX;
}	/* end for(int k=1; k <= khi; k++) */
