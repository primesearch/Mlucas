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

	/*...The radix-960 DIT pass is here:	*/

	#ifdef USE_SSE2

	// Gather the needed data (960 64-bit complex) and do 15 radix-64 transforms.
		tmp = r00;
		for(l = 0; l < 15; l++) {
			// Compute dit_i_offsets for current 64-block:
			for(kk = 0; kk < 4; kk++) {
				jp = (l<<2) + kk;
				i64 = dit_perm16[dit_perm16_idx[jp]];
 				jp = j1 + phi[dit_perm16_phidx[jp]];	// Main-array base index plus hi-part p-offset
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (kk<<4);
				k0 = plo[(i64 >> 60)&0xf];		dit_i_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dit_i_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dit_i_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dit_i_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dit_i_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dit_i_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dit_i_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dit_i_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dit_i_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dit_i_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dit_i_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dit_i_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dit_i_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dit_i_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dit_i_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dit_i_offsets[jt+0xf] = jp + kf;
			}
			// Use s1p00-3f for scratch for these 15 DFTs:
			SSE2_RADIX_64_DIT( FALSE, thr_id,
				a,dit_i_offsets,
				s1p00,	/* local scratch storage */
				tmp,t_offsets
			); tmp += 2;
		}

		/*...and now do 64 radix-15 transforms. See radix960_dit_pass1() for details on the oindex patterning.
		In 64-bit mode we use the 16-register doubled-radix-15-DFT macros as detailed in the radix60 carry routine:
		*/
	  #if (OS_BITS == 32)
		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce;

		tmp = r00;	tm2 = s1p00;	kk = 0x3c;	// Use kk as idx into phi[], so drop the trailing hex 0
		for(l = 0; l < 64; l++) {
		// NB: All offsets end up being shifted "one too many" bits leftward since
		// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:
		// [0] here:
			jp = ((64 - l)&0xf)<<1;	// Replace plo[] of scalar-double code with []<<1
		// [1] here:
			mask3 = (-(l > 0));
			k0 = kk & mask3;
		// [2] here:
			// Remember: All indices /= 16 here:		Now get the resulting p* offsets: Replace phi[] of scalar-double code with []<<5
			k1 = k0-0x10; k1 += (-(k1 < 0))&0x3c;		k0 = jp + (k0 << 5);
			k2 = k1-0x10; k2 += (-(k2 < 0))&0x3c;		k1 = jp + (k1 << 5);
			k3 = k2-0x10; k3 += (-(k3 < 0))&0x3c;		k2 = jp + (k2 << 5);
			k4 = k3-0x10; k4 += (-(k4 < 0))&0x3c;		k3 = jp + (k3 << 5);
			k5 = k4-0x10; k5 += (-(k5 < 0))&0x3c;		k4 = jp + (k4 << 5);
			k6 = k5-0x10; k6 += (-(k6 < 0))&0x3c;		k5 = jp + (k5 << 5);
			k7 = k6-0x10; k7 += (-(k7 < 0))&0x3c;		k6 = jp + (k6 << 5);
			k8 = k7-0x10; k8 += (-(k8 < 0))&0x3c;		k7 = jp + (k7 << 5);
			k9 = k8-0x10; k9 += (-(k9 < 0))&0x3c;		k8 = jp + (k8 << 5);
			ka = k9-0x10; ka += (-(ka < 0))&0x3c;		k9 = jp + (k9 << 5);
			kb = ka-0x10; kb += (-(kb < 0))&0x3c;		ka = jp + (ka << 5);
			kc = kb-0x10; kc += (-(kc < 0))&0x3c;		kb = jp + (kb << 5);
			kd = kc-0x10; kd += (-(kd < 0))&0x3c;		kc = jp + (kc << 5);
			ke = kd-0x10; ke += (-(ke < 0))&0x3c;		kd = jp + (kd << 5);
														ke = jp + (ke << 5);
			// Set up for next loop execution:
			kk -= 0x2c + ((l&0xf) == 0); //<*** further decr leading term whenever lo part passes thru 0
			kk += (-(kk < 0))&0x3c;
			// 1st set of Output ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;

			// Input ptrs:		
			va0 = tmp     ;
			va1 = tmp+0x02;
			va2 = tmp+0x04;
			va3 = tmp+0x06;
			va4 = tmp+0x08;
			va5 = tmp+0x0a;
			va6 = tmp+0x0c;
			va7 = tmp+0x0e;
			va8 = tmp+0x10;
			va9 = tmp+0x12;
			vaa = tmp+0x14;
			vab = tmp+0x16;
			vac = tmp+0x18;
			vad = tmp+0x1a;
			vae = tmp+0x1c;

			SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1,
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,	/* Ins: r00 +... */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* Scratch */
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce		/* Outs: s1p00 +... */
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}

	  #else

		vec_dbl
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,*va9,*vaa,*vab,*vac,*vad,*vae,
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,*vb9,*vba,*vbb,*vbc,*vbd,*vbe,
		*vc0,*vc1,*vc2,*vc3,*vc4,*vc5,*vc6,*vc7,*vc8,*vc9,*vca,*vcb,*vcc,*vcd,*vce,
		*vd0,*vd1,*vd2,*vd3,*vd4,*vd5,*vd6,*vd7,*vd8,*vd9,*vda,*vdb,*vdc,*vdd,*vde;

		tmp = r00;	tm2 = s1p00;	kk = 0x3c;	// Use kk as idx into phi[], so drop the trailing hex 0
		for(l = 0; l < 64; l += 2) {
		// NB: All offsets end up being shifted "one too many" bits leftward since
		// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:
		// [0] here:
			jp = ((64 - l)&0xf)<<1;	// Replace plo[] of scalar-double code with []<<1
		// [1] here:
			mask3 = (-(l > 0));
			k0 = kk & mask3;
		// [2] here:
			// Remember: All indices /= 16 here:		Now get the resulting p* offsets: Replace phi[] of scalar-double code with []<<5
			k1 = k0-0x10; k1 += (-(k1 < 0))&0x3c;		k0 = jp + (k0 << 5);
			k2 = k1-0x10; k2 += (-(k2 < 0))&0x3c;		k1 = jp + (k1 << 5);
			k3 = k2-0x10; k3 += (-(k3 < 0))&0x3c;		k2 = jp + (k2 << 5);
			k4 = k3-0x10; k4 += (-(k4 < 0))&0x3c;		k3 = jp + (k3 << 5);
			k5 = k4-0x10; k5 += (-(k5 < 0))&0x3c;		k4 = jp + (k4 << 5);
			k6 = k5-0x10; k6 += (-(k6 < 0))&0x3c;		k5 = jp + (k5 << 5);
			k7 = k6-0x10; k7 += (-(k7 < 0))&0x3c;		k6 = jp + (k6 << 5);
			k8 = k7-0x10; k8 += (-(k8 < 0))&0x3c;		k7 = jp + (k7 << 5);
			k9 = k8-0x10; k9 += (-(k9 < 0))&0x3c;		k8 = jp + (k8 << 5);
			ka = k9-0x10; ka += (-(ka < 0))&0x3c;		k9 = jp + (k9 << 5);
			kb = ka-0x10; kb += (-(kb < 0))&0x3c;		ka = jp + (ka << 5);
			kc = kb-0x10; kc += (-(kc < 0))&0x3c;		kb = jp + (kb << 5);
			kd = kc-0x10; kd += (-(kd < 0))&0x3c;		kc = jp + (kc << 5);
			ke = kd-0x10; ke += (-(ke < 0))&0x3c;		kd = jp + (kd << 5);
														ke = jp + (ke << 5);
			// Set up for next loop execution:
			kk -= 0x2c + ((l&0xf) == 0); //<*** further decr leading term whenever lo part passes thru 0
			kk += (-(kk < 0))&0x3c;
			// 1st set of Output ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
		// Now emulate the elided odd-index loop pass (replace l with l+1 in the above sequence) to get next set of ks:
			jp = ((63 - l)&0xf)<<1;	// Replace plo[] of scalar-double code with []<<1
		// [1] here: Know l > 0 so no need to compute mask, just assume = all-ones:
			k0 = kk;
		// [2] here:
			// Remember: All indices /= 16 here:		Now get the resulting p* offsets: Replace phi[] of scalar-double code with []<<5
			k1 = k0-0x10; k1 += (-(k1 < 0))&0x3c;		k0 = jp + (k0 << 5);
			k2 = k1-0x10; k2 += (-(k2 < 0))&0x3c;		k1 = jp + (k1 << 5);
			k3 = k2-0x10; k3 += (-(k3 < 0))&0x3c;		k2 = jp + (k2 << 5);
			k4 = k3-0x10; k4 += (-(k4 < 0))&0x3c;		k3 = jp + (k3 << 5);
			k5 = k4-0x10; k5 += (-(k5 < 0))&0x3c;		k4 = jp + (k4 << 5);
			k6 = k5-0x10; k6 += (-(k6 < 0))&0x3c;		k5 = jp + (k5 << 5);
			k7 = k6-0x10; k7 += (-(k7 < 0))&0x3c;		k6 = jp + (k6 << 5);
			k8 = k7-0x10; k8 += (-(k8 < 0))&0x3c;		k7 = jp + (k7 << 5);
			k9 = k8-0x10; k9 += (-(k9 < 0))&0x3c;		k8 = jp + (k8 << 5);
			ka = k9-0x10; ka += (-(ka < 0))&0x3c;		k9 = jp + (k9 << 5);
			kb = ka-0x10; kb += (-(kb < 0))&0x3c;		ka = jp + (ka << 5);
			kc = kb-0x10; kc += (-(kc < 0))&0x3c;		kb = jp + (kb << 5);
			kd = kc-0x10; kd += (-(kd < 0))&0x3c;		kc = jp + (kc << 5);
			ke = kd-0x10; ke += (-(ke < 0))&0x3c;		kd = jp + (kd << 5);
														ke = jp + (ke << 5);
			// Set up for next loop execution: Odd l here, so no need for "pass thru 0?" check:
			kk -= 0x2c;
			kk += (-(kk < 0))&0x3c;
			// 2nd set of Output ptrs:
			vd0 = tm2+k0;
			vd1 = tm2+k1;
			vd2 = tm2+k2;
			vd3 = tm2+k3;
			vd4 = tm2+k4;
			vd5 = tm2+k5;
			vd6 = tm2+k6;
			vd7 = tm2+k7;
			vd8 = tm2+k8;
			vd9 = tm2+k9;
			vda = tm2+ka;
			vdb = tm2+kb;
			vdc = tm2+kc;
			vdd = tm2+kd;
			vde = tm2+ke;

			// 2 sets of Input ptrs:		
							tm1 = tmp + 0x1e;
			va0 = tmp     ;	vb0 = tm1     ;
			va1 = tmp+0x02;	vb1 = tm1+0x02;
			va2 = tmp+0x04;	vb2 = tm1+0x04;
			va3 = tmp+0x06;	vb3 = tm1+0x06;
			va4 = tmp+0x08;	vb4 = tm1+0x08;
			va5 = tmp+0x0a;	vb5 = tm1+0x0a;
			va6 = tmp+0x0c;	vb6 = tm1+0x0c;
			va7 = tmp+0x0e;	vb7 = tm1+0x0e;
			va8 = tmp+0x10;	vb8 = tm1+0x10;
			va9 = tmp+0x12;	vb9 = tm1+0x12;
			vaa = tmp+0x14;	vba = tm1+0x14;
			vab = tmp+0x16;	vbb = tm1+0x16;
			vac = tmp+0x18;	vbc = tm1+0x18;
			vad = tmp+0x1a;	vbd = tm1+0x1a;
			vae = tmp+0x1c;	vbe = tm1+0x1c;

			SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,	/* In1: r00 +... */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* Scratch1 */
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,	/* Out1: s1p00 +... */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe,	/* In2: r0f +... */
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,	/* Scratch2 */
				vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde		/* Out2: s1p00+jt +... */
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	  #endif	// (OS_BITS == 32)

	#else	/* !USE_SSE2 */

	// Gather the needed data (960 64-bit complex) and do 15 radix-64 transforms.
		for(l = 0; l < 15; l++) {
			// Compute dit_i_offsets for current 64-block:
			for(kk = 0; kk < 4; kk++) {
				jp = (l<<2) + kk;
				i64 = dit_perm16[dit_perm16_idx[jp]];
 				jp = j1 + phi[dit_perm16_phidx[jp]];	// Main-array base index plus hi-part p-offset
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (kk<<4);
				k0 = plo[(i64 >> 60)&0xf];		dit_i_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dit_i_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dit_i_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dit_i_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dit_i_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dit_i_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dit_i_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dit_i_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dit_i_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dit_i_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dit_i_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dit_i_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dit_i_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dit_i_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dit_i_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dit_i_offsets[jt+0xf] = jp + kf;
			}
			RADIX_64_DIT(a,dit_i_offsets,RE_IM_STRIDE, (double *)(t+l),t_offsets,1);
		}

		//...and now do 64 radix-15 transforms:
		tptr = t;	kk = 0x3c;	// Use kk as idx into phi[], so drop the trailing hex 0
		for(l = 0; l < 64; l++) {
		// [0] here:
			jp = plo[(64 - l)&0xf];
		// [1] here:
			mask3 = (-(l > 0));
			k0 = kk & mask3;
		// [2] here:
			// Remember: All indices /= 16 here:		Now get the resulting p* offsets:
			k1 = k0-0x10; k1 += (-(k1 < 0))&0x3c;		k0 = jp + phi[k0];
			k2 = k1-0x10; k2 += (-(k2 < 0))&0x3c;		k1 = jp + phi[k1];
			k3 = k2-0x10; k3 += (-(k3 < 0))&0x3c;		k2 = jp + phi[k2];
			k4 = k3-0x10; k4 += (-(k4 < 0))&0x3c;		k3 = jp + phi[k3];
			k5 = k4-0x10; k5 += (-(k5 < 0))&0x3c;		k4 = jp + phi[k4];
			k6 = k5-0x10; k6 += (-(k6 < 0))&0x3c;		k5 = jp + phi[k5];
			k7 = k6-0x10; k7 += (-(k7 < 0))&0x3c;		k6 = jp + phi[k6];
			k8 = k7-0x10; k8 += (-(k8 < 0))&0x3c;		k7 = jp + phi[k7];
			k9 = k8-0x10; k9 += (-(k9 < 0))&0x3c;		k8 = jp + phi[k8];
			ka = k9-0x10; ka += (-(ka < 0))&0x3c;		k9 = jp + phi[k9];
			kb = ka-0x10; kb += (-(kb < 0))&0x3c;		ka = jp + phi[ka];
			kc = kb-0x10; kc += (-(kc < 0))&0x3c;		kb = jp + phi[kb];
			kd = kc-0x10; kd += (-(kd < 0))&0x3c;		kc = jp + phi[kc];
			ke = kd-0x10; ke += (-(ke < 0))&0x3c;		kd = jp + phi[kd];
														ke = jp + phi[ke];
			// Set up for next loop execution:
			kk -= 0x2c + ((l&0xf) == 0); //<*** further decr leading term whenever lo part passes thru 0
			kk += (-(kk < 0))&0x3c;
			RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],a[j1+k7],a[j2+k7],a[j1+k8],a[j2+k8],a[j1+k9],a[j2+k9],a[j1+ka],a[j2+ka],a[j1+kb],a[j2+kb],a[j1+kc],a[j2+kc],a[j1+kd],a[j2+kd],a[j1+ke],a[j2+ke]
			);	tptr += 0xf;
		}

	#endif	// USE_SSE2?

	/*...Now do the carries. Since the outputs would
	normally be getting dispatched to RADIX separate blocks of the A-array, we need 28 separate carries.	*/

/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
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
		SSE2_cmplx_carry_norm_errcheck0_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
		tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
		for(l = 1; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck1_2B(tm1,add1,add2,add3,tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

	/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		tm1 = s1p00; tmp = cy_r; tm2 = cy_r+0x01; itmp = bjmodn;
		for(l = 0; l < RADIX>>2; l++) {
			SSE2_cmplx_carry_norm_errcheck2_2B(tm1,add1,add2,     tmp,tm2,itmp,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			tm1 += 8; tmp += 2; tm2 += 2; itmp += 4;
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

		/*...set0 is slightly different from others; divide work into blocks of RADIX/4 macro calls, 1st set of which gets pulled out of loop: */		
		l = 0; addr = cy_r; itmp = bjmodn;
	   cmplx_carry_norm_errcheck0(a[j1   ],a[j2   ],*addr,*itmp  ); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p1],a[j2+p1],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p2],a[j2+p2],*addr,*itmp,l); ++l; ++addr; ++itmp;
		cmplx_carry_norm_errcheck(a[j1+p3],a[j2+p3],*addr,*itmp,l); ++l; ++addr; ++itmp;
		// Remaining quartets of macro calls done in loop:
		for(ntmp = 1; ntmp < RADIX>>2; ntmp++) {
			jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];
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

		tm0 = s1p00; tmp = base_negacyclic_root; l = 0xf000;
		tm1 = cy_r; // *cycle[] indices increment by +4 (mod ODD_RADIX) between macro calls
		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],icycle[ic+1],icycle[ic+2],icycle[ic+3], jcycle[ic],kcycle[ic],lcycle[ic] of the non-looped version with
		// icycle[ic],icycle[jc],icycle[kc],icycle[lc], jcycle[ic],kcycle[ic],lcycle[ic] :
		ic = 0; jc = 1; kc = 2; lc = 3;
		while(tm0 < x00)	// Can't use l for loop index here since need it for byte offset in carry macro call
		{
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_hiacc(tm0,tmp,l,tm1,0x1e00, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++; tmp += 8; l -= 0xc0;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #else	/* HIACC = false: */

		tm0 = s1p00; tmp = base_negacyclic_root;	// tmp *not* incremented between macro calls in loacc version
		tm1 = cy_r; // tm2 = cy_i;	*** replace with literal-byte-offset in macro call to save a reg
		ic = 0; jc = 1; kc = 2; lc = 3;
		for(l = 0; l < RADIX>>2; l++) {	// RADIX/4 loop passes
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];	k5 = jcycle[ic];	k6 = kcycle[ic];	k7 = lcycle[ic];
			k2 = icycle[ic];
			k3 = icycle[kc];
			k4 = icycle[lc];
																		/* vvvvvvvvvvvvvvv [1,2,3]*ODD_RADIX; assumed << l2_sz_vd on input: */
			SSE2_fermat_carry_norm_errcheck_X4_loacc(tm0,tmp,tm1,0x1e00, 0x1e0,0x3c0,0x5a0, half_arr,sign_mask,k1,k2,k3,k4,k5,k6,k7);
			tm0 += 8; tm1++;
			MOD_ADD32(ic, 4, ODD_RADIX, ic);
			MOD_ADD32(jc, 4, ODD_RADIX, jc);
			MOD_ADD32(kc, 4, ODD_RADIX, kc);
			MOD_ADD32(lc, 4, ODD_RADIX, lc);
		}

	  #endif	/* HIACC? */

	#elif defined(USE_SSE2)

		/* For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c. */

		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #if (OS_BITS == 64)

		// [ijkl]c = indices into icycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls; replace the
		// icycle[ic],jcycle[ic],icycle[ic+1],jcycle[ic+1] of the non-looped version with icycle[ic],jcycle[ic],icycle[jc],jcycle[jc]:
		ic = 0; jc = 1;
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX;	// Need to stick this #def into an intvar to work around [error: invalid lvalue in asm input for constraint 'm']
		while(tm1 < x00) {
			//See "Sep 2014" note in 32-bit SSE2 version of this code below
			k1 = icycle[ic];
			k2 = jcycle[ic];
			k3 = icycle[jc];
			k4 = jcycle[jc];
			SSE2_fermat_carry_norm_errcheck_X2(tm1,tmp,NRT_BITS,NRTM1,idx_offset,idx_incr,l,half_arr,sign_mask,add1,add2,k1,k2,k3,k4);
			tm1 += 4; tmp += 2;
			MOD_ADD32(ic, 2, ODD_RADIX, ic);
			MOD_ADD32(jc, 2, ODD_RADIX, jc);
		}

	  #else // Mar 2014: Worked around the out-of-regs compiler issues with the _X2 version of this macro (the
			// code in carry_gcc64.h has details), but keep non-X2 version in case hit out-of-regs again at some point

		ic = 0;	// ic = idx into [i|j]cycle mini-arrays, gets incremented (mod ODD_RADIX) between macro calls
		tm1 = s1p00; tmp = cy_r;	// <*** Again rely on contiguity of cy_r,i here ***
		l = ODD_RADIX << 4;	// 32-bit version needs preshifted << 4 input value
		while(tm1 < x00) {
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
			jt = j1 + poff[m]; jp = j2 + poff[m];
			fermat_carry_norm_errcheckB(a[jt   ],a[jp   ],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p1],a[jp+p1],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p2],a[jp+p2],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
			fermat_carry_norm_errcheckB(a[jt+p3],a[jp+p3],*addr,*addi,icycle[ic],ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR; ++addr; ++addi; MOD_ADD32(ic, 1, ODD_RADIX, ic);
		}
		for(ntmp = 0; ntmp < ODD_RADIX; ntmp++) {
			icycle[ntmp] += wts_idx_incr;	// Inside the loop use this, as it is faster than general-mod '% nwt'
			icycle[ntmp] += ( (-(int)((uint32)icycle[ntmp] >> 31)) & nwt);
		}

	#endif	/* #ifdef USE_SSE2 */

	// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
	#ifdef USE_SSE2
		for(ntmp = 0; ntmp < ODD_RADIX; ntmp++)
		{
			icycle[ntmp] += wts_idx_inc2;		icycle[ntmp] += ( (-(icycle[ntmp] < 0)) & nwt16);
			jcycle[ntmp] += wts_idx_inc2;		jcycle[ntmp] += ( (-(jcycle[ntmp] < 0)) & nwt16);
		#ifdef USE_AVX
			kcycle[ntmp] += wts_idx_inc2;		kcycle[ntmp] += ( (-(kcycle[ntmp] < 0)) & nwt16);
			lcycle[ntmp] += wts_idx_inc2;		lcycle[ntmp] += ( (-(lcycle[ntmp] < 0)) & nwt16);
		#endif
		}
	#endif

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-960 DIF pass is here:	*/

	#ifdef USE_SSE2

	//...gather the needed data (960 64-bit complex) and do 64 radix-15 transforms:
	  #if (OS_BITS == 32)	//Aug 2014: Too much work to 32-bit-ify the register-intensive large-pow2 macros, so punt on those
	
		tmp = r00;	tm2 = s1p00;	kk = 0;
		for(l = 0; l < 64; l++) {
		// NB: All offsets end up being shifted "one too many" bits leftward since
		// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:
		// [0] here:
			jp = (l&0xf)<<1;	// p0,..,f, repeated 4x (Replace plo[] of scalar-double code with []<<1)
		// [1],[2] here:
			// Now get the remaining row terms and the resulting p* offsets:
			k0 = kk;							// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-4; k1 += (-(k1 < 0))&60;		k0 = (k0 << 5) + jp;
			k2 = k1-4; k2 += (-(k2 < 0))&60;		k1 = (k1 << 5) + jp;
			k3 = k2-4; k3 += (-(k3 < 0))&60;		k2 = (k2 << 5) + jp;
			k4 = k3-4; k4 += (-(k4 < 0))&60;		k3 = (k3 << 5) + jp;
			k5 = k4-4; k5 += (-(k5 < 0))&60;		k4 = (k4 << 5) + jp;
			k6 = k5-4; k6 += (-(k6 < 0))&60;		k5 = (k5 << 5) + jp;
			k7 = k6-4; k7 += (-(k7 < 0))&60;		k6 = (k6 << 5) + jp;
			k8 = k7-4; k8 += (-(k8 < 0))&60;		k7 = (k7 << 5) + jp;
			k9 = k8-4; k9 += (-(k9 < 0))&60;		k8 = (k8 << 5) + jp;
			ka = k9-4; ka += (-(ka < 0))&60;		k9 = (k9 << 5) + jp;
			kb = ka-4; kb += (-(kb < 0))&60;		ka = (ka << 5) + jp;
			kc = kb-4; kc += (-(kc < 0))&60;		kb = (kb << 5) + jp;
			kd = kc-4; kd += (-(kd < 0))&60;		kc = (kc << 5) + jp;
			ke = kd-4; ke += (-(ke < 0))&60;		kd = (kd << 5) + jp;
													ke = (ke << 5) + jp;
			// Set up for next loop execution:
			kk -= ((l & 0xf) != 0xf);
			kk += (-(kk < 0))&60;
			// 1st set of input ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;

			// Output ptrs:
			va0 = tmp     ;
			va1 = tmp+0x02;
			va2 = tmp+0x04;
			va3 = tmp+0x06;
			va4 = tmp+0x08;
			va5 = tmp+0x0a;
			va6 = tmp+0x0c;
			va7 = tmp+0x0e;
			va8 = tmp+0x10;
			va9 = tmp+0x12;
			vaa = tmp+0x14;
			vab = tmp+0x16;
			vac = tmp+0x18;
			vad = tmp+0x1a;
			vae = tmp+0x1c;

			SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,	/* Ins: s1p00 +... */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* Scratch */
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae		/* Outs: r00 +... */
			);
			tmp += 0x1e;	// advance ptr by 15 vec_cmplx [= 30 vec_dbl] elements
		}

	  #else

		tmp = r00;	tm2 = s1p00;	kk = 0;
		for(l = 0; l < 64; l += 2) {
		// NB: All offsets end up being shifted "one too many" bits leftward since
		// we need vec_complex pointer offsets, whereas underlying pointers are vec_dbl:
		// [0] here:
			jp = (l&0xf)<<1;	// p0,..,f, repeated 4x (Replace plo[] of scalar-double code with []<<1)
		// [1],[2] here:
			// Now get the remaining row terms and the resulting p* offsets:
			k0 = kk;							// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-4; k1 += (-(k1 < 0))&60;		k0 = (k0 << 5) + jp;
			k2 = k1-4; k2 += (-(k2 < 0))&60;		k1 = (k1 << 5) + jp;
			k3 = k2-4; k3 += (-(k3 < 0))&60;		k2 = (k2 << 5) + jp;
			k4 = k3-4; k4 += (-(k4 < 0))&60;		k3 = (k3 << 5) + jp;
			k5 = k4-4; k5 += (-(k5 < 0))&60;		k4 = (k4 << 5) + jp;
			k6 = k5-4; k6 += (-(k6 < 0))&60;		k5 = (k5 << 5) + jp;
			k7 = k6-4; k7 += (-(k7 < 0))&60;		k6 = (k6 << 5) + jp;
			k8 = k7-4; k8 += (-(k8 < 0))&60;		k7 = (k7 << 5) + jp;
			k9 = k8-4; k9 += (-(k9 < 0))&60;		k8 = (k8 << 5) + jp;
			ka = k9-4; ka += (-(ka < 0))&60;		k9 = (k9 << 5) + jp;
			kb = ka-4; kb += (-(kb < 0))&60;		ka = (ka << 5) + jp;
			kc = kb-4; kc += (-(kc < 0))&60;		kb = (kb << 5) + jp;
			kd = kc-4; kd += (-(kd < 0))&60;		kc = (kc << 5) + jp;
			ke = kd-4; ke += (-(ke < 0))&60;		kd = (kd << 5) + jp;
													ke = (ke << 5) + jp;
			// Set up for next loop execution:
			kk -= ((l & 0xf) != 0xf);
			kk += (-(kk < 0))&60;
			// 1st set of input ptrs:
			vc0 = tm2+k0;
			vc1 = tm2+k1;
			vc2 = tm2+k2;
			vc3 = tm2+k3;
			vc4 = tm2+k4;
			vc5 = tm2+k5;
			vc6 = tm2+k6;
			vc7 = tm2+k7;
			vc8 = tm2+k8;
			vc9 = tm2+k9;
			vca = tm2+ka;
			vcb = tm2+kb;
			vcc = tm2+kc;
			vcd = tm2+kd;
			vce = tm2+ke;
			
		// Now emulate the elided odd-index loop pass (replace l with l+1 in the above sequence) to get next set of ks:
		// [0] here:
			jp += 2;
		// [1],[2] here:
			k0 = kk;							// Now get the resulting shifted offsets, adding the low bits to each:
			k1 = k0-4; k1 += (-(k1 < 0))&60;		k0 = (k0 << 5) + jp;
			k2 = k1-4; k2 += (-(k2 < 0))&60;		k1 = (k1 << 5) + jp;
			k3 = k2-4; k3 += (-(k3 < 0))&60;		k2 = (k2 << 5) + jp;
			k4 = k3-4; k4 += (-(k4 < 0))&60;		k3 = (k3 << 5) + jp;
			k5 = k4-4; k5 += (-(k5 < 0))&60;		k4 = (k4 << 5) + jp;
			k6 = k5-4; k6 += (-(k6 < 0))&60;		k5 = (k5 << 5) + jp;
			k7 = k6-4; k7 += (-(k7 < 0))&60;		k6 = (k6 << 5) + jp;
			k8 = k7-4; k8 += (-(k8 < 0))&60;		k7 = (k7 << 5) + jp;
			k9 = k8-4; k9 += (-(k9 < 0))&60;		k8 = (k8 << 5) + jp;
			ka = k9-4; ka += (-(ka < 0))&60;		k9 = (k9 << 5) + jp;
			kb = ka-4; kb += (-(kb < 0))&60;		ka = (ka << 5) + jp;
			kc = kb-4; kc += (-(kc < 0))&60;		kb = (kb << 5) + jp;
			kd = kc-4; kd += (-(kd < 0))&60;		kc = (kc << 5) + jp;
			ke = kd-4; ke += (-(ke < 0))&60;		kd = (kd << 5) + jp;
													ke = (ke << 5) + jp;
			// Set up for next loop execution:
			kk -= (((l+1) & 0xf) != 0xf);
			kk += (-(kk < 0))&60;
			// 2nd set of input ptrs:
			vd0 = tm2+k0;
			vd1 = tm2+k1;
			vd2 = tm2+k2;
			vd3 = tm2+k3;
			vd4 = tm2+k4;
			vd5 = tm2+k5;
			vd6 = tm2+k6;
			vd7 = tm2+k7;
			vd8 = tm2+k8;
			vd9 = tm2+k9;
			vda = tm2+ka;
			vdb = tm2+kb;
			vdc = tm2+kc;
			vdd = tm2+kd;
			vde = tm2+ke;

			// 2 sets of output ptrs:		
							tm1 = tmp + 0x1e;
			va0 = tmp     ;	vb0 = tm1     ;
			va1 = tmp+0x02;	vb1 = tm1+0x02;
			va2 = tmp+0x04;	vb2 = tm1+0x04;
			va3 = tmp+0x06;	vb3 = tm1+0x06;
			va4 = tmp+0x08;	vb4 = tm1+0x08;
			va5 = tmp+0x0a;	vb5 = tm1+0x0a;
			va6 = tmp+0x0c;	vb6 = tm1+0x0c;
			va7 = tmp+0x0e;	vb7 = tm1+0x0e;
			va8 = tmp+0x10;	vb8 = tm1+0x10;
			va9 = tmp+0x12;	vb9 = tm1+0x12;
			vaa = tmp+0x14;	vba = tm1+0x14;
			vab = tmp+0x16;	vbb = tm1+0x16;
			vac = tmp+0x18;	vbc = tm1+0x18;
			vad = tmp+0x1a;	vbd = tm1+0x1a;
			vae = tmp+0x1c;	vbe = tm1+0x1c;

			SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1,
				vc0,vc1,vc2,vc3,vc4,vc5,vc6,vc7,vc8,vc9,vca,vcb,vcc,vcd,vce,	/* Ins1: s1p00 +... */
				x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,	/* Scratch1 */
				va0,va1,va2,va3,va4,va5,va6,va7,va8,va9,vaa,vab,vac,vad,vae,	/* Out1: r00 +... */
				vd0,vd1,vd2,vd3,vd4,vd5,vd6,vd7,vd8,vd9,vda,vdb,vdc,vdd,vde,	/* Ins2: s1p00 + jt... */
				y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,	/* Scratch2 */
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8,vb9,vba,vbb,vbc,vbd,vbe		/* Out2: r0f +... */
			);
			tmp += 0x3c;	// advance ptr by 30 vec_cmplx [= 60 vec_dbl] elements
		}

	  #endif	// (OS_BITS == 32)

		/*...and now do 15 radix-64 transforms, with index-perms as described in radix960_dif_pass1().
		inputs of SSE2_RADIX16_DIT_0TWIDDLE from 30*vec_dbl - separated memlocs, same offsets as already set for DIT: */

		tmp = r00;
		for(l = 0; l < 15; l++) {
			// Compute dif_o_offsets for current 64-block:
			for(kk = 0; kk < 4; kk++) {
				jp = (l<<2) + kk;
				i64 = dif_perm16[jp];
 				jp = p_out_hi[jp];	// Add the p[0123]0 term stored in the current p_out_hi quartet
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (kk<<4);
				k0 = plo[(i64 >> 60)&0xf];		dif_o_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dif_o_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dif_o_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dif_o_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dif_o_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dif_o_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dif_o_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dif_o_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dif_o_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dif_o_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dif_o_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dif_o_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dif_o_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dif_o_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dif_o_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dif_o_offsets[jt+0xf] = jp + kf;
			}
			// Compute the p40 multiplier which acts as the base offset:
			mask3 = -(l<3);
			jp = ( ( (3-l)&(-(l>0)) ) & mask3 ) + ( (17-l) & ~mask3 );	// j = 0,2,1,14,13,12,11,10,9,8,7,6,5,4,3
			jt = j1 + phi[jp<<2];		//  Add main-array index  to j*p40 - dif_o_offsets are all relative to this
		// Use s1p00-3f for scratch for these 3 DFTs ... since transform length N = odd*64,
		// the leading pow2-shift arg = trailz(N) - trailz(64) = 0:
			SSE2_RADIX_64_DIF( FALSE, thr_id,
				0,
				tmp,t_offsets,
				s1p00,	// tmp-storage
				a+jt,dif_o_offsets
			); tmp += 2;
		}

	#else	/* !USE_SSE2 */

	//...gather the needed data (960 64-bit complex) and do 64 radix-15 transforms:
		tptr = t; kk = 0;
		for(l = 0; l < 64; l++) {
		// [0] here:
			jp = plo[l&0xf];	// p0,..,f, repeated 4x
			jt = j1 + jp; jp += j2;
		// [1],[2] here:
			// Now get the remaining row terms and the resulting p* offsets:
			k0 = kk;	// Need a separate master index kk because k0 will get overwritten with phi[k0]
			k1 = k0-4; k1 += (-(k1 < 0))&60;		k0 = phi[k0];
			k2 = k1-4; k2 += (-(k2 < 0))&60;		k1 = phi[k1];
			k3 = k2-4; k3 += (-(k3 < 0))&60;		k2 = phi[k2];
			k4 = k3-4; k4 += (-(k4 < 0))&60;		k3 = phi[k3];
			k5 = k4-4; k5 += (-(k5 < 0))&60;		k4 = phi[k4];
			k6 = k5-4; k6 += (-(k6 < 0))&60;		k5 = phi[k5];
			k7 = k6-4; k7 += (-(k7 < 0))&60;		k6 = phi[k6];
			k8 = k7-4; k8 += (-(k8 < 0))&60;		k7 = phi[k7];
			k9 = k8-4; k9 += (-(k9 < 0))&60;		k8 = phi[k8];
			ka = k9-4; ka += (-(ka < 0))&60;		k9 = phi[k9];
			kb = ka-4; kb += (-(kb < 0))&60;		ka = phi[ka];
			kc = kb-4; kc += (-(kc < 0))&60;		kb = phi[kb];
			kd = kc-4; kd += (-(kd < 0))&60;		kc = phi[kc];
			ke = kd-4; ke += (-(ke < 0))&60;		kd = phi[kd];
													ke = phi[ke];
			// Set up for next loop execution:
			kk -= ((l & 0xf) != 0xf);
			kk += (-(kk < 0))&60;
			RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
			);	tptr += 0xf;
		}

	//...and now do 15 radix-64 transforms:
		for(l = 0; l < 15; l++) {
			// Compute dif_o_offsets for current 64-block:
			for(kk = 0; kk < 4; kk++) {
				jp = (l<<2) + kk;
				i64 = dif_perm16[jp];
 				jp = p_out_hi[jp];	// Add the p[0123]0 term stored in the current p_out_hi quartet
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (kk<<4);
				k0 = plo[(i64 >> 60)&0xf];		dif_o_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dif_o_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dif_o_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dif_o_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dif_o_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dif_o_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dif_o_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dif_o_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dif_o_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dif_o_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dif_o_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dif_o_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dif_o_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dif_o_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dif_o_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dif_o_offsets[jt+0xf] = jp + kf;
			}
			// Compute the p40 multiplier which acts as the base offset:
			mask3 = -(l<3);
			jp = ( ( (3-l)&(-(l>0)) ) & mask3 ) + ( (17-l) & ~mask3 );	// j = 0,2,1,14,13,12,11,10,9,8,7,6,5,4,3
			jt = j1 + phi[jp<<2];		//  Add main-array index  to j*p40 - dif_o_offsets are all relative to this
			RADIX_64_DIF((double *)(t+l),t_offsets,1, (a+jt),dif_o_offsets,RE_IM_STRIDE);
		}

	#endif	// USE_SSE2 ?
	}

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(k=1; k <= khi; k++) */

