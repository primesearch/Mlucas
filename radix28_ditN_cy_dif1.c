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

#include "Mlucas.h"

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#define FFT_DEBUG	0
#if FFT_DEBUG
	char dbg_fname[] = "a.txt";
#endif

#ifdef CTIME	// define at compile time to enable internal timing diagnostics
	double dt_fwd, dt_inv, dt_cy, dt_tot;
	clock_t clock1, clock2, clock3;
#endif

#ifdef USE_SSE2

	const int radix28_creals_in_local_store = 128;

	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

//	#define	USE_SCALAR_CARRY	// Uncomment if want to use original non-SSE carry macros

		/*
		Here more details about the small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.
		Here are sample data for:

		F24 using a length-7*2^17 transform, complex radices	F24 using a length-7*2^18 transform, complex radices
		= 28,16,32,32 [only the leading 28 matters here]:		= 28,32,32,32:

			J = 0:												J = 0:
			ii00 = 0, bjmodn00 = 917504, i = 1					ii00 = 0, bjmodn00 = 1835008, i = 1
			ii01 = 6, bjmodn01 = 131072, i = 0					ii01 = 5, bjmodn01 =  524288, i = 0
			ii02 = 5, bjmodn02 = 262144, i = 0					ii02 = 3, bjmodn02 = 1048576, i = 0
			ii03 = 4, bjmodn03 = 393216, i = 0					ii03 = 1, bjmodn03 = 1572864, i = 1
			ii04 = 3, bjmodn04 = 524288, i = 0					ii04 = 6, bjmodn04 =  262144, i = 0
			ii05 = 2, bjmodn05 = 655360, i = 0					ii05 = 4, bjmodn05 =  786432, i = 0
			ii06 = 1, bjmodn06 = 786432, i = 1					ii06 = 2, bjmodn06 = 1310720, i = 0
			J = 2:												J = 2:
			ii00 = 5, bjmodn00 = 262144, i = 0					ii00 = 5, bjmodn00 =  524288, i = 0
			ii01 = 4, bjmodn01 = 393216, i = 0					ii01 = 3, bjmodn01 = 1048576, i = 0
			ii02 = 3, bjmodn02 = 524288, i = 0					ii02 = 1, bjmodn02 = 1572864, i = 1
			ii03 = 2, bjmodn03 = 655360, i = 0					ii03 = 6, bjmodn03 =  262144, i = 0
			ii04 = 1, bjmodn04 = 786432, i = 1					ii04 = 4, bjmodn04 =  786432, i = 0
			ii05 = 0, bjmodn05 = 917504, i = 1					ii05 = 2, bjmodn05 = 1310720, i = 0
			ii06 = 6, bjmodn06 = 131072, i = 0					ii06 = 0, bjmodn06 = 1835008, i = 1
			J = 4:												J = 4:
			ii00 = 3, bjmodn00 = 524288, i = 0					ii00 = 3, bjmodn00 = 1048576, i = 0
			ii01 = 2, bjmodn01 = 655360, i = 0					ii01 = 1, bjmodn01 = 1572864, i = 1
			ii02 = 1, bjmodn02 = 786432, i = 1					ii02 = 6, bjmodn02 =  262144, i = 0
			ii03 = 0, bjmodn03 = 917504, i = 1					ii03 = 4, bjmodn03 =  786432, i = 0
			ii04 = 6, bjmodn04 = 131072, i = 0					ii04 = 2, bjmodn04 = 1310720, i = 0
			ii05 = 5, bjmodn05 = 262144, i = 0					ii05 = 0, bjmodn05 = 1835008, i = 1
			ii06 = 4, bjmodn06 = 393216, i = 0					ii06 = 5, bjmodn06 =  524288, i = 0
			J = 6:												J = 6:
			ii00 = 1, bjmodn00 = 786432, i = 1					ii00 = 1, bjmodn00 = 1572864, i = 1
			ii01 = 0, bjmodn01 = 917504, i = 1					ii01 = 6, bjmodn01 =  262144, i = 0
			ii02 = 6, bjmodn02 = 131072, i = 0					ii02 = 4, bjmodn02 =  786432, i = 0
			ii03 = 5, bjmodn03 = 262144, i = 0					ii03 = 2, bjmodn03 = 1310720, i = 0
			ii04 = 4, bjmodn04 = 393216, i = 0					ii04 = 0, bjmodn04 = 1835008, i = 1
			ii05 = 3, bjmodn05 = 524288, i = 0					ii05 = 5, bjmodn05 =  524288, i = 0
			ii06 = 2, bjmodn06 = 655360, i = 0					ii06 = 3, bjmodn06 = 1048576, i = 0
			J = 8:												J = 8:
			ii00 = 6, bjmodn00 = 131072, i = 0					ii00 = 6, bjmodn00 =  262144, i = 0
			ii01 = 5, bjmodn01 = 262144, i = 0					ii01 = 4, bjmodn01 =  786432, i = 0
			ii02 = 4, bjmodn02 = 393216, i = 0					ii02 = 2, bjmodn02 = 1310720, i = 0
			ii03 = 3, bjmodn03 = 524288, i = 0					ii03 = 0, bjmodn03 = 1835008, i = 1
			ii04 = 2, bjmodn04 = 655360, i = 0					ii04 = 5, bjmodn04 =  524288, i = 0
			ii05 = 1, bjmodn05 = 786432, i = 1					ii05 = 3, bjmodn05 = 1048576, i = 0
			ii06 = 0, bjmodn06 = 917504, i = 1					ii06 = 1, bjmodn06 = 1572864, i = 1
			J = 10:												J = 10:
			ii00 = 4, bjmodn00 = 393216, i = 0					ii00 = 4, bjmodn00 =  786432, i = 0
			ii01 = 3, bjmodn01 = 524288, i = 0					ii01 = 2, bjmodn01 = 1310720, i = 0
			ii02 = 2, bjmodn02 = 655360, i = 0					ii02 = 0, bjmodn02 = 1835008, i = 1
			ii03 = 1, bjmodn03 = 786432, i = 1					ii03 = 5, bjmodn03 =  524288, i = 0
			ii04 = 0, bjmodn04 = 917504, i = 1					ii04 = 3, bjmodn04 = 1048576, i = 0
			ii05 = 6, bjmodn05 = 131072, i = 0					ii05 = 1, bjmodn05 = 1572864, i = 1
			ii06 = 5, bjmodn06 = 262144, i = 0					ii06 = 6, bjmodn06 =  262144, i = 0
			J = 12:												J = 12:
			ii00 = 2, bjmodn00 = 655360, i = 0					ii00 = 2, bjmodn00 = 1310720, i = 0
			ii01 = 1, bjmodn01 = 786432, i = 1					ii01 = 0, bjmodn01 = 1835008, i = 1
			ii02 = 0, bjmodn02 = 917504, i = 1					ii02 = 5, bjmodn02 =  524288, i = 0
			ii03 = 6, bjmodn03 = 131072, i = 0					ii03 = 3, bjmodn03 = 1048576, i = 0
			ii04 = 5, bjmodn04 = 262144, i = 0					ii04 = 1, bjmodn04 = 1572864, i = 1
			ii05 = 4, bjmodn05 = 393216, i = 0					ii05 = 6, bjmodn05 =  262144, i = 0
			ii06 = 3, bjmodn06 = 524288, i = 0					ii06 = 4, bjmodn06 =  786432, i = 0
	`	...And now (i.e. every [nwt]th pass) repeat the j=0 pattern: ...
			J = 14:												J = 14:
			ii00 = 0, bjmodn00 = 917504, i = 1					ii00 = 0, bjmodn00 = 1835008, i = 1
			ii01 = 6, bjmodn01 = 131072, i = 0					ii01 = 5, bjmodn01 =  524288, i = 0
			ii02 = 5, bjmodn02 = 262144, i = 0					ii02 = 3, bjmodn02 = 1048576, i = 0
			ii03 = 4, bjmodn03 = 393216, i = 0					ii03 = 1, bjmodn03 = 1572864, i = 1
			ii04 = 3, bjmodn04 = 524288, i = 0					ii04 = 6, bjmodn04 =  262144, i = 0
			ii05 = 2, bjmodn05 = 655360, i = 0					ii05 = 4, bjmodn05 =  786432, i = 0
			ii06 = 1, bjmodn06 = 786432, i = 1					ii06 = 2, bjmodn06 = 1310720, i = 0

		For the F24 case, the cyclical per-loop-pass index-pattern shift = 2; for F25 it = 1. How to compute it in advance?

		The per-loop-pass update of the bjmodn terms is this:

			i = (bjmodn > sw);					//       i = 1 if a bigword,   0 if a smallword
			bjmodn -= sw;						// result >= 0 if a bigword, < 0 if a smallword
			bjmodn += ( ((int)i-1) & n);		//       add 0 if a bigword,   N if a smallword

		So need to take bjmodn00 start value,which = n, subtract sw, see how many bjmodn's further on we need to go to get the resulting bjmodn00-sw value.



		Take advantage of this as follows:
		Init length-[nwt] integer arrays ii_arr = (initial values of ii00-ii06), i_arr = (initial values of (bjmodn00-06 > sw)),
		treat these as circular arrays, for j = 0 starting index into these arrays = 0,
		on each loop execution we advance the starting index in these arrays by wts_idx_incr := (bw*radix0/n) places.
		These ii and i-array values are used as lookups into suitably initialized weights and base arrays and arrays of their inverses.
		In fact if we use the initial pattersn of the ii adn i-indices (i.e. those corresponding the j=0 initial loop interation)
		to init correspondingly permuted weights and base arrays, we can replace actual in-loop use of separate ii and i-arrays
		with a simple set of [nwt] indices which start with values [0,1,2,...,nwt-1] and get incremented by wts_idx_incr (modulo nwt) each pass.

		Example 1: For p=2^24, n = 917504 = 7*2^17, bjmodn00 - sw = bw = p%n = 262144 which corr. to bjmodn02 in the initial-bjmodn-values list,
		==> per-pass circular shift amount of the wts-array [and related] index = 2 := wts_idx_incr.

		Example 2: For p=2^25, n = 1835008 = 7*2^18, bjmodn00 - sw = bw = p%n = 524288 which corr. to bjmodn01 in the initial-bjmodn-values list,
		==> per-pass circular shift amount of the wts-array [and related] index = 1 := wts_idx_incr.

		The 2 indexing subarrays that come into play are:

		[1] ii: Define SW_DIV_N = sw*nwt/n = (n-bw)*nwt/n = 655360*7/917504 = [5*2^17]*7/[7*2^17] = nwt-wts_idx_incr = 5
			We need (nwt) starting values of ii, with the (j)th of these defined by

			For the above example, with radix0 = 28, we get

			ii[j] = j*(SW_DIV_N*(n/[FFT radix used for the carry step]/2) % nwt
				  = j*[nwt-wts_idx_incr]*[n/28/2] % nwt
				  = j*5*2^14 % 7
				  = j*20 % 7
				  = -j%7
				  = 0,6,5,4,3,2,1 . This determines the order in which we access the wt[] and wtinv[] arrays, so should simply init length-nwt local versions of these arrays which contain the weights ordered according to that index pattern.

		[2] i : See any of the fused final-iFFT-radix/carry/initial-FFT-radix pass routines to see how the parameter bjmodnini is inited...define it as

				bjmodnini := bw*(n/2)/[FFT radix used for the carry step] % n .

			For the above example, with radix0 = 28, we get

				bjmodnini = [2*2^17]*[7*2^16]/28 % n
				= 2^32 % n
				= 131072 ,

			with the (j)th value of this

				bjmodn[j] := j*bjmodini % n,

			and we start this value at n (rather that 0) for j = 0, so the resulting index i into the base and baseinv arrays comes out right:

			i [j] = (bjmodn[j] > sw) .

			For the above example, with radix0 = 28, we get

			i [j] = [1, (j*131072) for j=1...6] > sw = 1,0,0,0,0,0,1, so should simply init length-nwt local versions of the base[] and baseinv[] arrays which contain the entries ordered according to this index pattern.
		*/

	#if defined(COMPILER_TYPE_MSVC)

		#include "sse2_macro.h"

		/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(__i0,__i1,__i2,__i3, __o0,__o1,__o2,__o3)\
		{\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __i3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */	__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p2] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */	__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */					__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */					__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */					__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */					__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into output-array slots: */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	mov	edx, __o3\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[ecx      ],xmm0	/* <- ~t5 */			__asm	movaps	[edx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[ebx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */			__asm	movaps	[ebx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */\
		}

		/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(__i0,__i1,__i2,__i3, __o0,__o1,__o2,__o3)\
		{\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __i3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ebx     ]	/* a[jt+p1] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p1] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */	__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p1] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */	__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p1] */\
			\
			__asm	addpd	xmm0,[ecx     ]	/* t1 */					__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm1,[ecx+0x10]	/* t2 */					__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm2,[ecx     ]	/* t3 */					__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm3,[ecx+0x10]	/* t4 */					__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into output-array slots: */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	mov	edx, __o3\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */			__asm	movaps	[ecx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[edx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */			__asm	movaps	[edx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[ecx+0x010],xmm6	/* <- ~t8 */\
		}

		/*...Radix-7 DFT: Inputs in memory locations __i0-6, outputs go into memory locations __o0-6, possibly coincident with inputs:\
		*/\
		#define SSE2_RADIX_07_DFT(__i0,__i1,__i2,__i3,__i4,__i5,__i6, __cc, __o0,__o1,__o2,__o3,__o4,__o5,__o6)\
		{\
		/*\
			t1r=A1r+A6r;	\
			t6r=A1r-A6r;	\
							\
			t2r=A2r+A5r;	\
			t5r=A2r-A5r;	\
							\
			t3r=A3r+A4r;	\
			t4r=A3r-A4r;	\
		*/\
			__asm	mov	eax, __i1	\
			__asm	mov	ebx, __i2	\
			__asm	mov	ecx, __i3	\
			__asm	mov	edx, __i4	\
			__asm	mov	esi, __i5	\
			__asm	mov	edi, __i6	\
			__asm	movaps	xmm6,[eax     ]	/* A1r */\
			__asm	movaps	xmm1,[edi     ]	/* A6r */\
			__asm	movaps	xmm5,[ebx     ]	/* A2r */\
			__asm	movaps	xmm2,[esi     ]	/* A5r */\
			__asm	movaps	xmm4,[ecx     ]	/* A3r */\
			__asm	movaps	xmm3,[edx     ]	/* A4r */\
			\
			__asm	mov	ebx, __i0	\
			__asm	subpd	xmm6,xmm1	/* t6r = A1r-A6r */\
			__asm	addpd	xmm1,xmm1	/*         2*A6r */\
			__asm	addpd	xmm1,xmm6	/* t1r = A1r+A6r */\
			\
			__asm	subpd	xmm5,xmm2	/* t5r = A2r-A5r */\
			__asm	addpd	xmm2,xmm2	/*         2*A5r */\
			__asm	addpd	xmm2,xmm5	/* t2r = A2r+A5r */\
			\
			__asm	movaps	xmm0,[ebx     ]	/* Ar0 */\
			__asm	subpd	xmm4,xmm3	/* t4r = A3r-A4r */\
			__asm	addpd	xmm3,xmm3	/*         2*A4r */\
			__asm	addpd	xmm3,xmm4	/* t3r = A3r+A4r */\
		/*\
			rt  = t1r+t2r+t3r;	\
			B0r = rt + A0r;		\
			t0r = rt*cx0 + A0r;			t3r=(t6r-t4r+t5r)*sx0;	\
			t1r = t1r-t2r;				t6r= t6r-t5r;			\
			t2r = t3r-t2r;				t5r= t4r+t5r;			\
			t3r =(t1r+t2r)*cx3;			t4r=(t5r-t6r)*sx3;		\
			t1r = t1r*cx1;				t6r= t6r*sx1;			\
			t2r = t2r*cx2;				t5r= t5r*sx2;			\
			tt  = t1r-t3r;				t6r= t4r+t6r;			\
			t2r = t2r-t3r;				t5r= t4r-t5r;			\
																\
			t1r= t0r- tt-t2r;			t4r= t3r-t6r-t5r;		\
			t2r= t0r+t2r;				t5r= t3r+t5r;			\
			t0r= t0r+ tt;				t3r= t3r+t6r;			\
		*/\
			__asm	mov	ecx, __o0	/* Assume that this might be the same address as any of i0-i6 */\
			__asm	mov	esi, __cc	\
			__asm	movaps	[esi+0x80],xmm0	/* cpy t0 into scratch sincos slot */	__asm	movaps	[esi+0x90],xmm6	/* cpy t6 into scratch sincos slot */	\
			__asm	addpd	xmm0,xmm1	/*~A0 = A0+t1 */							__asm	movaps	xmm7,xmm5	/* cpy t5 */			\
			__asm	addpd	xmm3,xmm2	/*~t3 = t3+t2 */							__asm	subpd	xmm5,xmm4	/*~t5 = t5-t4 */		\
			__asm	subpd	xmm1,xmm2	/*~t1 = t1-t2 */							__asm	subpd	xmm6,xmm7	/*~t6 = t6-t5 */		\
			__asm	addpd	xmm2,xmm2	/* 2*t2 */									__asm	addpd	xmm4,xmm7	/*~t5 = t4+t5 */		\
			__asm	addpd	xmm0,xmm3	/* B0 */									__asm	addpd	xmm5,[esi+0x90]	/* t3 = [t5-t4]+t6 */	\
			__asm	subpd	xmm3,xmm2	/*~t2 =  [t2+t3] - 2*t2 = t3-t2 */			__asm	movaps	xmm7,xmm4	/* cpy t5 */			\
			__asm	movaps	[ecx     ],xmm0	/* <-B0, xmm0 FREE */					__asm	subpd	xmm4,xmm6	/* t4 = ~t5-~t6 */		\
			__asm	movaps	xmm2,xmm1	/* cpy ~t1 */																	\
			__asm	subpd	xmm0,[esi+0x80]	/* r = B0 - t0 */						__asm	mulpd	xmm5,[esi+0x10]	/*~t3 = t3*sx0 */	\
			__asm	addpd	xmm2,xmm3	/* ~t1+~t2 */																	\
			__asm	mulpd	xmm3,[esi+0x40]	/* t2 = t2*cx2 */						__asm	mulpd	xmm4,[esi+0x70]	/*~t4 = t4*sx3 */	\
			__asm	mulpd	xmm1,[esi+0x20]	/* t1 = t1*cx1 */						__asm	mulpd	xmm6,[esi+0x30]	/*~t6 = t6*sx1 */	\
			__asm	mulpd	xmm0,[esi]     	/* ~r = r*(cx0-1) */					__asm	mulpd	xmm7,[esi+0x50]	/*~t5 = t5*sx2 */	\
			__asm	mulpd	xmm2,[esi+0x60]	/* t3 =(t1+t2)*cx3 */																		\
			__asm	addpd	xmm0,[ecx     ]	/* t0 =~r + B0 */						__asm	addpd	xmm6,xmm4	/*~t6 = t4+t6 */		\
			__asm	subpd	xmm1,xmm2	/* tt = t1-t3 */							__asm	subpd	xmm4,xmm7	/*~t5 = t4-t5, xmm7 FREE */\
			__asm	subpd	xmm3,xmm2	/* t2 = t2-t3, xmm2 FREE */					\
			__asm	mov	eax, __o1													\
			__asm	mov	ebx, __o2													\
			__asm	mov	ecx, __o3													\
			__asm	mov	edx, __o4													\
			__asm	mov	esi, __o5													\
			__asm	mov	edi, __o6													\
			__asm	movaps	xmm2,xmm0	/* cpy t0 */								__asm	movaps	xmm7,xmm5	/* cpy t3 */		\
			__asm	addpd	xmm0,xmm1	/*~t0 = t0+tt */							__asm	addpd	xmm5,xmm6	/*~t3 = t3+t6 */	\
			__asm	addpd	xmm1,xmm3	/*~tt = tt+t2 */							__asm	addpd	xmm6,xmm4	/*      t6+t5 */	\
			__asm	addpd	xmm3,xmm2	/*~t2 = t2+t0 */							__asm	addpd	xmm4,xmm7	/*~t5 = t5+t3 */	\
			__asm	subpd	xmm2,xmm1	/*~t1 = t0-tt-t2 */							__asm	subpd	xmm7,xmm6	/*~t4 = t3-t6-t5 */	\
			__asm	movaps	[eax     ],xmm0	/* B1 <- t0 */							__asm	movaps	[edi     ],xmm5	/* B6 <- t3 */	\
			__asm	movaps	[ebx     ],xmm2	/* B2 <- t1 */							__asm	movaps	[esi     ],xmm7	/* B5 <- t4 */	\
			__asm	movaps	[ecx     ],xmm3	/* B3 <- t2 */							__asm	movaps	[edx     ],xmm4	/* B4 <- t5 */	\
			\
		/************************** Imaginary Parts: ******************************************/\
			\
			__asm	mov	eax, __i1	\
			__asm	mov	ebx, __i2	\
			__asm	mov	ecx, __i3	\
			__asm	mov	edx, __i4	\
			__asm	mov	esi, __i5	\
			__asm	mov	edi, __i6	\
			__asm	movaps	xmm6,[eax+0x10]	/* A1i */\
			__asm	movaps	xmm1,[edi+0x10]	/* A6i */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A2i */\
			__asm	movaps	xmm2,[esi+0x10]	/* A5i */\
			__asm	movaps	xmm4,[ecx+0x10]	/* A3i */\
			__asm	movaps	xmm3,[edx+0x10]	/* A4i */\
			\
			__asm	mov	ebx, __i0	\
			__asm	subpd	xmm6,xmm1	/* t6i = A1i-A6i */\
			__asm	addpd	xmm1,xmm1	/*         2*A6i */\
			__asm	addpd	xmm1,xmm6	/* t1i = A1i+A6i */\
			\
			__asm	subpd	xmm5,xmm2	/* t5i = A2i-A5i */\
			__asm	addpd	xmm2,xmm2	/*         2*A5i */\
			__asm	addpd	xmm2,xmm5	/* t2i = A2i+A5i */\
			\
			__asm	movaps	xmm0,[ebx+0x10]	/* Ai0 */\
			__asm	subpd	xmm4,xmm3	/* t4i = A3i-A4i */\
			__asm	addpd	xmm3,xmm3	/*         2*A4i */\
			__asm	addpd	xmm3,xmm4	/* t3i = A3i+A4i */\
		/*\
			it  = t1i+t2i+t3i;	\
			B0i = it + A0i;		\
			t0i = it*cx0 + A0i;			t3i=(t6i-t4i+t5i)*sx0;	\
			t1i = t1i-t2i;				t6i= t6i-t5i;			\
			t2i = t2i-t3i;				t5i= t4i+t5i;			\
			t3i =(t1i-t2i)*cx3;			t4i=(t5i-t6i)*sx3;		\
			t1i = t1i*cx1;				t6i= t6i*sx1;			\
			t2i = t2i*cx2;				t5i= t5i*sx2;			\
			it  = t1i-t3i;				t6i= t4i+t6i;			\
			t2i = t2i-t3i;				t5i= t4i-t5i;			\
																\
			t1i= t0i- it-t2i;			t4i= t3i-t6i-t5i;		\
			t2i= t0i+t2i;				t5i= t3i+t5i;			\
			t0i= t0i+ it;				t3i= t3i+t6i;			\
		*/\
			__asm	mov	ecx, __o0	\
			__asm	mov	esi, __cc	\
			__asm	movaps	[esi+0x80],xmm0	/* cpy t0 into scratch sincos slot */	__asm	movaps	[esi+0x90],xmm6	/* cpy t6 into scratch sincos slot */	\
			__asm	addpd	xmm0,xmm1	/*~A0 = A0+t1 */							__asm	movaps	xmm7,xmm5	/* cpy t5 */			\
			__asm	addpd	xmm3,xmm2	/*~t3 = t3+t2 */							__asm	subpd	xmm5,xmm4	/*~t5 = t5-t4 */		\
			__asm	subpd	xmm1,xmm2	/*~t1 = t1-t2 */							__asm	subpd	xmm6,xmm7	/*~t6 = t6-t5 */		\
			__asm	addpd	xmm2,xmm2	/* 2*t2 */									__asm	addpd	xmm4,xmm7	/*~t5 = t4+t5 */		\
			__asm	addpd	xmm0,xmm3	/* B0 */									__asm	addpd	xmm5,[esi+0x90]	/* t3 = [t5-t4]+t6 */	\
			__asm	subpd	xmm3,xmm2	/*~t2 =  [t2+t3] - 2*t2 = t3-t2 */			__asm	movaps	xmm7,xmm4	/* cpy t5 */			\
			__asm	movaps	[ecx+0x10],xmm0	/* <-B0, xmm0 FREE */					__asm	subpd	xmm4,xmm6	/* t4 = ~t5-~t6 */		\
			__asm	movaps	xmm2,xmm1	/* cpy ~t1 */																	\
			__asm	subpd	xmm0,[esi+0x80]	/* r = B0 - t0 */						__asm	mulpd	xmm5,[esi+0x10]	/*~t3 = t3*sx0 */	\
			__asm	addpd	xmm2,xmm3	/* ~t1+~t2 */																	\
			__asm	mulpd	xmm3,[esi+0x40]	/* t2 = t2*cx2 */						__asm	mulpd	xmm4,[esi+0x70]	/*~t4 = t4*sx3 */	\
			__asm	mulpd	xmm1,[esi+0x20]	/* t1 = t1*cx1 */						__asm	mulpd	xmm6,[esi+0x30]	/*~t6 = t6*sx1 */	\
			__asm	mulpd	xmm0,[esi]     	/* ~r = r*(cx0-1) */					__asm	mulpd	xmm7,[esi+0x50]	/*~t5 = t5*sx2 */	\
			__asm	mulpd	xmm2,[esi+0x60]	/* t3 =(t1+t2)*cx3 */																		\
			__asm	addpd	xmm0,[ecx+0x10]	/* t0 =~r + B0 */						__asm	addpd	xmm6,xmm4	/*~t6 = t4+t6 */		\
			__asm	subpd	xmm1,xmm2	/* tt = t1-t3 */							__asm	subpd	xmm4,xmm7	/*~t5 = t4-t5, xmm7 FREE */\
			__asm	subpd	xmm3,xmm2	/* t2 = t2-t3, xmm2 FREE */					\
			__asm	mov	eax, __o1													\
			__asm	mov	ebx, __o2													\
			__asm	mov	ecx, __o3													\
			__asm	movaps	xmm2,xmm0	/* cpy t0 */								__asm	movaps	xmm7,xmm5	/* cpy t3 */		\
			__asm	addpd	xmm0,xmm1	/*~t0 = t0+tt */							__asm	addpd	xmm5,xmm6	/*~t3 = t3+t6 */	\
			__asm	addpd	xmm1,xmm3	/*~tt = tt+t2 */							__asm	addpd	xmm6,xmm4	/*      t6+t5 */	\
			__asm	addpd	xmm3,xmm2	/*~t2 = t2+t0 */							__asm	addpd	xmm4,xmm7	/*~t5 = t5+t3 */	\
			__asm	subpd	xmm2,xmm1	/*~t1 = t0-tt-t2, xmm1 FREE */				__asm	subpd	xmm7,xmm6	/*~t4 = t3-t6-t5, xmm6 FREE */	\
		/*\
			B1r =t0r-t3i;					B1i*=t0i+t3r;\
			B2r =t1r-t4i;					B2i*=t1i+t4r;\
			B3r*=t2r+t5i;					B3i =t2i-t5r;\
			B4r*=t2r-t5i;					B4i =t2i+t5r;\
			B5r =t1r+t4i;					B5i*=t1i-t4r;\
			B6r =t0r+t3i;					B6i*=t0i-t3r;\
		*/\
			__asm	mov	edx, __o4\
			__asm	mov	esi, __o5\
			__asm	mov	edi, __o6\
			/* xmm1,6 FREE */\
			__asm	movaps	xmm1,[eax     ]	/* t0r */					__asm	movaps	xmm6,[edi     ]	/* t3r */					\
			__asm	subpd	xmm1,xmm5	/* B1r =t0r-t3i */				__asm	subpd	xmm0,xmm6	/* B6i =t0i-t3r */				\
			__asm	addpd	xmm5,xmm5	/*        2*t3i */				__asm	addpd	xmm6,xmm6	/*        2*t3r */				\
			__asm	addpd	xmm5,xmm1	/* B6r =t0r+t3i */				__asm	addpd	xmm6,xmm0	/* B1i =t0i+t3r */				\
			__asm	movaps	[eax     ],xmm1	/* <-B1r */					__asm	movaps	[edi+0x10],xmm0	/* <-B6i */		\
			__asm	movaps	[edi     ],xmm5	/* <-B6r */					__asm	movaps	[eax+0x10],xmm6	/* <-B1i */		\
			\
			__asm	movaps	xmm1,[ebx     ]	/* t1r */					__asm	movaps	xmm6,[esi     ]	/* t4r */					\
			__asm	subpd	xmm1,xmm7	/* B2r =t1r-t4i */				__asm	subpd	xmm2,xmm6	/* B5i =t1i-t4r */				\
			__asm	addpd	xmm7,xmm7	/*        2*t4i */				__asm	addpd	xmm6,xmm6	/*        2*t4r */				\
			__asm	addpd	xmm7,xmm1	/* B5r =t1r+t4i */				__asm	addpd	xmm6,xmm2	/* B2i =t1i+t4r */				\
			__asm	movaps	[ebx     ],xmm1	/* <-B2r */					__asm	movaps	[esi+0x10],xmm2	/* <-B5i */		\
			__asm	movaps	[esi     ],xmm7	/* <-B5r*/					__asm	movaps	[ebx+0x10],xmm6	/* <-B2i */		\
			\
			/* Note the order reversal on this pair of outputs: */\
			__asm	movaps	xmm0,[ecx     ]	/* t2r */					__asm	movaps	xmm5,[edx     ]	/* t5r */					\
			__asm	subpd	xmm0,xmm4	/* B4r =t2r-t5i */				__asm	subpd	xmm3,xmm5	/* B3i =t2i-t5r */				\
			__asm	addpd	xmm4,xmm4	/*        2*t5i */				__asm	addpd	xmm5,xmm5	/*        2*t5r */				\
			__asm	addpd	xmm4,xmm0	/* B3r =t2r+t5i */				__asm	addpd	xmm5,xmm3	/* B4i =t2i+t5r */				\
			__asm	movaps	[edx     ],xmm0	/* <-B4r */					__asm	movaps	[ecx+0x10],xmm3	/* <-B3i */		\
			__asm	movaps	[ecx     ],xmm4	/* <-B3r*/					__asm	movaps	[edx+0x10],xmm5	/* <-B4i */		\
		}

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix28_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix28_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

  #ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int tid;
		int ndivr;
		int _pad0;	// Pads to make sizeof this struct a multiple of 16 bytes
		int _pad1;
		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;
		int wts_idx_inc2;
		int icycle0;
		int icycle1;
		int icycle2;
		int icycle3;
		int icycle4;
		int icycle5;
		int icycle6;
		int jcycle0;
		int jcycle1;
		int jcycle2;
		int jcycle3;
		int jcycle4;
		int jcycle5;
		int jcycle6;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *rn0;
		struct complex *rn1;
		struct complex *s1p00r;

		int bjmodn00;
		int bjmodn01;
		int bjmodn02;
		int bjmodn03;
		int bjmodn04;
		int bjmodn05;
		int bjmodn06;
		int bjmodn07;
		int bjmodn08;
		int bjmodn09;
		int bjmodn10;
		int bjmodn11;
		int bjmodn12;
		int bjmodn13;
		int bjmodn14;
		int bjmodn15;
		int bjmodn16;
		int bjmodn17;
		int bjmodn18;
		int bjmodn19;
		int bjmodn20;
		int bjmodn21;
		int bjmodn22;
		int bjmodn23;
		int bjmodn24;
		int bjmodn25;
		int bjmodn26;
		int bjmodn27;

		// Pad to make size a multiple of 64 bytes:
		double dpad0;
		double dpad1;
		double dpad2;

		/* carries: */
		double cy_r00;
		double cy_r01;
		double cy_r02;
		double cy_r03;
		double cy_r04;
		double cy_r05;
		double cy_r06;
		double cy_r07;
		double cy_r08;
		double cy_r09;
		double cy_r10;
		double cy_r11;
		double cy_r12;
		double cy_r13;
		double cy_r14;
		double cy_r15;
		double cy_r16;
		double cy_r17;
		double cy_r18;
		double cy_r19;
		double cy_r20;
		double cy_r21;
		double cy_r22;
		double cy_r23;
		double cy_r24;
		double cy_r25;
		double cy_r26;
		double cy_r27;

		double cy_i00;
		double cy_i01;
		double cy_i02;
		double cy_i03;
		double cy_i04;
		double cy_i05;
		double cy_i06;
		double cy_i07;
		double cy_i08;
		double cy_i09;
		double cy_i10;
		double cy_i11;
		double cy_i12;
		double cy_i13;
		double cy_i14;
		double cy_i15;
		double cy_i16;
		double cy_i17;
		double cy_i18;
		double cy_i19;
		double cy_i20;
		double cy_i21;
		double cy_i22;
		double cy_i23;
		double cy_i24;
		double cy_i25;
		double cy_i26;
		double cy_i27;
	};

  #endif

#endif

#define RADIX_07_DFT_NUSS(\
	__A0r,__A0i,\
	__A1r,__A1i,\
	__A2r,__A2i,\
	__A3r,__A3i,\
	__A4r,__A4i,\
	__A5r,__A5i,\
	__A6r,__A6i,\
	__t0r,__t0i,\
	__t1r,__t1i,\
	__t2r,__t2i,\
	__t3r,__t3i,\
	__t4r,__t4i,\
	__t5r,__t5i,\
	__t6r,__t6i,\
	__B0r,__B0i,\
	__B1r,__B1i,\
	__B2r,__B2i,\
	__B3r,__B3i,\
	__B4r,__B4i,\
	__B5r,__B5i,\
	__B6r,__B6i,\
	__cx0,__sx0,\
	__cx1,__sx1,\
	__cx2,__sx2,\
	__cx3,__sx3,\
	__rt,__it)\
{\
	__t0r = __A0r;						__t0i = __A0i;				\
	__t6r = __A1r - __A6r;				__t6i = __A1i - __A6i;	/* x1 - x6 */	\
	__t1r = __A1r + __A6r;				__t1i = __A1i + __A6i;	/* x1 + x6 */	\
	\
	__t5r = __A2r - __A5r;				__t5i = __A2i - __A5i;	/* x2 - x5 */	\
	__t2r = __A2r + __A5r;				__t2i = __A2i + __A5i;	/* x2 + x5 */	\
	\
	__t4r = __A3r - __A4r;				__t4i = __A3i - __A4i;	/* x3 - x4 */	\
	__t3r = __A3r + __A4r;				__t3i = __A3i + __A4i;	/* x3 + x4 */	\
	\
	__rt = __t1r+__t2r+__t3r;			__it = __t1i+__t2i+__t3i;		\
	__B0r= __rt+__t0r;					__B0i= __it+__t0i;				\
	__t0r= __rt*__cx0+__t0r;			__t0i= __it*__cx0+__t0i;		\
	__t1r= __t1r-__t2r;					__t1i= __t1i-__t2i;				\
	__t2r= __t3r-__t2r;					__t2i= __t3i-__t2i;				\
	__t3r=(__t1r+__t2r)*__cx3;			__t3i=(__t1i+__t2i)*__cx3;		\
	__t1r= __t1r*__cx1;					__t1i= __t1i*__cx1;				\
	__t2r= __t2r*__cx2;					__t2i= __t2i*__cx2;				\
	__rt = __t1r-__t3r;					__it = __t1i-__t3i;				\
	__t2r= __t2r-__t3r;					__t2i= __t2i-__t3i;				\
																		\
	__t1r= __t0r-__rt-__t2r;			__t1i= __t0i-__it-__t2i;		\
	__t2r= __t0r+__t2r;					__t2i= __t0i+__t2i;				\
	__t0r= __t0r+__rt;					__t0i= __t0i+__it;				\
																		\
	__t3r=(__t6r-__t4r+__t5r)*__sx0;	__t3i=(__t6i-__t4i+__t5i)*__sx0;\
	__t6r= __t6r-__t5r;					__t6i= __t6i-__t5i;				\
	__t5r= __t4r+__t5r;					__t5i= __t4i+__t5i;				\
	__t4r=(__t5r-__t6r)*__sx3;			__t4i=(__t5i-__t6i)*__sx3;		\
	__t6r= __t6r*__sx1;					__t6i= __t6i*__sx1;				\
	__t5r= __t5r*__sx2;					__t5i= __t5i*__sx2;				\
	__t6r= __t4r+__t6r;					__t6i= __t4i+__t6i;				\
	__t5r= __t4r-__t5r;					__t5i= __t4i-__t5i;				\
																		\
	__t4r= __t3r-__t6r-__t5r;			__t4i= __t3i-__t6i-__t5i;		\
	__t5r= __t3r+__t5r;					__t5i= __t3i+__t5i;				\
	__t3r= __t3r+__t6r;					__t3i= __t3i+__t6i;				\
																		\
	__B1r =__t0r-__t3i;					__B1i =__t0i+__t3r;				\
	__B2r =__t1r-__t4i;					__B2i =__t1i+__t4r;				\
	__B3r =__t2r+__t5i;					__B3i =__t2i-__t5r;				\
	__B4r =__t2r-__t5i;					__B4i =__t2i+__t5r;				\
	__B5r =__t1r+__t4i;					__B5i =__t1i-__t4r;				\
	__B6r =__t0r+__t3i;					__B6i =__t0i-__t3r;				\
}

/**************/

int radix28_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-28 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-28 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const int RADIX = 28;
	const double crnd = 3.0*0x4000000*0x2000000;
	int NDIVR,i,j,j1,j2,jstart,jhi,full_pass,k,khi,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24;
	static double radix_inv, n2inv;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << 4;
	const int odd_radix = 7;
#endif

#if defined(USE_SSE2) || !defined (LO_ADD)
	/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
	static double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
				 	cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
				 	cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
				 	cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
				/* Switch the sign of ss3 in these: */
				 	sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
				 	sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
				 	sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
				 	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#else
	/* Non-SSE2 version assumes LO_ADD = 1 and uses the corresponding versions of the sincos constants: */
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
#endif

	double scale
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int l,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27;	/* indices into weights arrays (mod NWT) */
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD

	#ifdef USE_PTHREAD
		static struct complex *__r0;	// Base address for discrete per-thread local stores
		static struct cy_thread_data_t *tdat = 0x0;
		// Threadpool-based dispatch stuff:
		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)cy28_process_chunk, NULL, 0x0};
	#endif

  #else

//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */

  #endif

	static struct complex *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
  #ifndef COMPILER_TYPE_GCC
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i
  #endif
	, *tmp;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27;
	/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY #define below;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		, icycle0, icycle1, icycle2, icycle3, icycle4, icycle5, icycle6
		, jcycle0, jcycle1, jcycle2, jcycle3, jcycle4, jcycle5, jcycle6;

	static double wt_arr[7],wtinv_arr[7],bs_arr[7],bsinv_arr[7];

  #if(defined(USE_SCALAR_CARRY) || defined(DEBUG_SSE2))
	double temp,frac
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;
  #endif
  #ifdef USE_SCALAR_CARRY
	double cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
	      ,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27;
  #else
	static struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26;
	static struct complex *cy_i00,*cy_i02,*cy_i04,*cy_i06,*cy_i08,*cy_i10,*cy_i12,*cy_i14,*cy_i16,*cy_i18,*cy_i20,*cy_i22,*cy_i24,*cy_i26;
  #endif

  #ifdef DEBUG_SSE2
	int jt,jp;
  #endif

#else

	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it;
	int jt,jp,k1,k2,m,m2,ntmp;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	static double wt_re,wt_im;									/* Fermat-mod weights stuff */
  #if PFETCH
	double *addr, *addp;
  #endif
	double temp,frac;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27;
	double re,im
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,
	*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0;

#ifdef CTIME
	const double ICPS = 1.0/CLOCKS_PER_SEC;
	clock1 = clock();
	dt_fwd = dt_inv = dt_cy = dt_tot = 0.0;
#endif

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii00=ii01=ii02=ii03=ii04=ii05=ii06=ii07=ii08=ii09=ii10=ii11=ii12=ii13=ii14=ii15=ii16=ii17=ii18=ii19=ii20=ii21=ii22=ii23=ii24=ii25=ii26=ii27=-1;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/28 in radix28_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT(HERE, LO_ADD,"radix28_ditN_cy_dif1.c: LO_ADD");
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		if(0 != (j & 0xf)) {
			printf("sizeof(cy_thread_data_t) = %x\n",j);
			ASSERT(HERE, 0, "struct cy_thread_data_t not 16-byte size multiple!");
		}
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#ifdef OS_TYPE_MACOSX

			if(CY_THREADS > 1) {
				main_work_units = CY_THREADS/2;
				pool_work_units = CY_THREADS - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
			} else {
				main_work_units = 1;
				printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
			}

		#else

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifndef DEBUG_SSE2
	//	#error "USE_SCALAR_CARRY obsolete!"
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 128 dcomplex and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix28_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix28_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 56 16-byte slots of sc_arr for temporaries, next 8 for the nontrivial complex 16th roots,
	next 28 for the doubled carry pairs, next 2 for ROE and RND_CONST, next RADIX for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
	  #ifdef COMPILER_TYPE_GCC
		s1p00r = sc_ptr + 0x00;
		s1p01r = sc_ptr + 0x02;
		s1p02r = sc_ptr + 0x04;
		s1p03r = sc_ptr + 0x06;
		s1p04r = sc_ptr + 0x08;
		s1p05r = sc_ptr + 0x0a;
		s1p06r = sc_ptr + 0x0c;
		s1p07r = sc_ptr + 0x0e;
		s1p08r = sc_ptr + 0x10;
		s1p09r = sc_ptr + 0x12;
		s1p10r = sc_ptr + 0x14;
		s1p11r = sc_ptr + 0x16;
		s1p12r = sc_ptr + 0x18;
		s1p13r = sc_ptr + 0x1a;
		s1p14r = sc_ptr + 0x1c;
		s1p15r = sc_ptr + 0x1e;
		s1p16r = sc_ptr + 0x20;
		s1p17r = sc_ptr + 0x22;
		s1p18r = sc_ptr + 0x24;
		s1p19r = sc_ptr + 0x26;
		s1p20r = sc_ptr + 0x28;
		s1p21r = sc_ptr + 0x2a;
		s1p22r = sc_ptr + 0x2c;
		s1p23r = sc_ptr + 0x2e;
		s1p24r = sc_ptr + 0x30;
		s1p25r = sc_ptr + 0x32;
		s1p26r = sc_ptr + 0x34;
		s1p27r = sc_ptr + 0x36;
	  #else
		s1p00r = sc_ptr + 0x00;		s1p00i = sc_ptr + 0x01;
		s1p01r = sc_ptr + 0x02;		s1p01i = sc_ptr + 0x03;
		s1p02r = sc_ptr + 0x04;		s1p02i = sc_ptr + 0x05;
		s1p03r = sc_ptr + 0x06;		s1p03i = sc_ptr + 0x07;
		s1p04r = sc_ptr + 0x08;		s1p04i = sc_ptr + 0x09;
		s1p05r = sc_ptr + 0x0a;		s1p05i = sc_ptr + 0x0b;
		s1p06r = sc_ptr + 0x0c;		s1p06i = sc_ptr + 0x0d;
		s1p07r = sc_ptr + 0x0e;		s1p07i = sc_ptr + 0x0f;
		s1p08r = sc_ptr + 0x10;		s1p08i = sc_ptr + 0x11;
		s1p09r = sc_ptr + 0x12;		s1p09i = sc_ptr + 0x13;
		s1p10r = sc_ptr + 0x14;		s1p10i = sc_ptr + 0x15;
		s1p11r = sc_ptr + 0x16;		s1p11i = sc_ptr + 0x17;
		s1p12r = sc_ptr + 0x18;		s1p12i = sc_ptr + 0x19;
		s1p13r = sc_ptr + 0x1a;		s1p13i = sc_ptr + 0x1b;
		s1p14r = sc_ptr + 0x1c;		s1p14i = sc_ptr + 0x1d;
		s1p15r = sc_ptr + 0x1e;		s1p15i = sc_ptr + 0x1f;
		s1p16r = sc_ptr + 0x20;		s1p16i = sc_ptr + 0x21;
		s1p17r = sc_ptr + 0x22;		s1p17i = sc_ptr + 0x23;
		s1p18r = sc_ptr + 0x24;		s1p18i = sc_ptr + 0x25;
		s1p19r = sc_ptr + 0x26;		s1p19i = sc_ptr + 0x27;
		s1p20r = sc_ptr + 0x28;		s1p20i = sc_ptr + 0x29;
		s1p21r = sc_ptr + 0x2a;		s1p21i = sc_ptr + 0x2b;
		s1p22r = sc_ptr + 0x2c;		s1p22i = sc_ptr + 0x2d;
		s1p23r = sc_ptr + 0x2e;		s1p23i = sc_ptr + 0x2f;
		s1p24r = sc_ptr + 0x30;		s1p24i = sc_ptr + 0x31;
		s1p25r = sc_ptr + 0x32;		s1p25i = sc_ptr + 0x33;
		s1p26r = sc_ptr + 0x34;		s1p26i = sc_ptr + 0x35;
		s1p27r = sc_ptr + 0x36;		s1p27i = sc_ptr + 0x37;
	  #endif
		cc0		= sc_ptr + 0x38;
		ss0		= sc_ptr + 0x39;
		cc1		= sc_ptr + 0x3a;
		ss1		= sc_ptr + 0x3b;
		cc2		= sc_ptr + 0x3c;
		ss2		= sc_ptr + 0x3d;
		cc3  	= sc_ptr + 0x3e;
		ss3		= sc_ptr + 0x3f;	/* Extra 2 slots for scratch storage here */
	  #ifndef USE_SCALAR_CARRY
		cy_r00	= sc_ptr + 0x42;
		cy_r02	= sc_ptr + 0x43;
		cy_r04	= sc_ptr + 0x44;
		cy_r06	= sc_ptr + 0x45;
		cy_r08	= sc_ptr + 0x46;
		cy_r10	= sc_ptr + 0x47;
		cy_r12	= sc_ptr + 0x48;
		cy_r14	= sc_ptr + 0x49;
		cy_r16	= sc_ptr + 0x4a;
		cy_r18	= sc_ptr + 0x4b;
		cy_r20	= sc_ptr + 0x4c;
		cy_r22	= sc_ptr + 0x4d;
		cy_r24	= sc_ptr + 0x4e;
		cy_r26	= sc_ptr + 0x4f;
		cy_i00	= sc_ptr + 0x50;
		cy_i02	= sc_ptr + 0x51;
		cy_i04	= sc_ptr + 0x52;
		cy_i06	= sc_ptr + 0x53;
		cy_i08	= sc_ptr + 0x54;
		cy_i10	= sc_ptr + 0x55;
		cy_i12	= sc_ptr + 0x56;
		cy_i14	= sc_ptr + 0x57;
		cy_i16	= sc_ptr + 0x58;
		cy_i18	= sc_ptr + 0x59;
		cy_i20	= sc_ptr + 0x5a;
		cy_i22	= sc_ptr + 0x5b;
		cy_i24	= sc_ptr + 0x5c;
		cy_i26	= sc_ptr + 0x5d;
	  #endif
		max_err = sc_ptr + 0x5e;
		sse2_rnd= sc_ptr + 0x5f;
		half_arr= sc_ptr + 0x60;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */

		/* These remain fixed: */
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		cc0->re = cc0->im = cx0-1;	ss0->re = ss0->im = sx0;
		cc1->re = cc1->im = cx1;	ss1->re = ss1->im = sx1;
		cc2->re = cc2->im = cx2;	ss2->re = ss2->im = sx2;
		cc3->re = cc3->im = cx3;	ss3->re = ss3->im = sx3;

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = sse2_rnd->im = crnd;

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.
		*/
		tmp = half_arr;

	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		/* In Fermat-mod mode, init these in the same section where we compute the initial ii-values and icycle0-6 indices */
	}
	else
	{
		/* Forward-weight multipliers: */
		tmp->re = 1.0;	tmp->im = 1.0;	++tmp;
		tmp->re = .50;	tmp->im = 1.0;	++tmp;
		tmp->re = 1.0;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->re = .50;	tmp->im = .50;	++tmp;
		tmp->re = .25;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .25;	++tmp;
		tmp->re = .25;	tmp->im = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [0];	tmp->im = base   [1];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->re = baseinv[0];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[0];	tmp->im = baseinv[1];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[1];	++tmp;
	}

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		*sign_mask++ = (uint64)0x7FFFFFFFFFFFFFFFull;
		*sign_mask-- = (uint64)0x7FFFFFFFFFFFFFFFull;

	  #if 0
		// Set up the quadrupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + 2;
		__asm	mov	eax, bw
		__asm	mov	ebx, sse_bw
		__asm	movd	xmm0,eax	/* Move actual *value* of reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_sw  = sm_ptr + 4;
		__asm	lea	eax, sw
		__asm	mov	ebx, sse_sw
		__asm	movd	xmm0,[eax]	/* Variant 2: Move contents of address pointed to by reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_n   = sm_ptr + 6;
		__asm	lea	eax, n
		__asm	mov	ebx, sse_n
		__asm	movd	xmm0,[eax]
		__asm	pshufd	xmm0,xmm0,0	// Broadcast low 32 bits of xmm0 to all 4 slots of xmm0
		__asm	movaps	[ebx],xmm0
	  #else
		sse_bw  = sm_ptr + 2;
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_bw++ = tmp64;
		*sse_bw-- = tmp64;

		sse_sw  = sm_ptr + 4;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_sw++ = tmp64;
		*sse_sw-- = tmp64;

		sse_n   = sm_ptr + 6;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_n++ = tmp64;
		*sse_n-- = tmp64;
	  #endif

#ifdef USE_PTHREAD
	/* Populate the elements of the thread-specific data structs which don't change after init: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
	// int data:
		tdat[ithread].tid = ithread;
		tdat[ithread].ndivr = NDIVR;
	
		tdat[ithread].sw  = sw;
		tdat[ithread].nwt = nwt;
		tdat[ithread].wts_idx_inc2 = wts_idx_inc2;	// The real init of  this must wait until after we compute wts_idx_inc2 below

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		tdat[ithread].wt0 = wt0;
		tdat[ithread].wt1 = wt1;
		tdat[ithread].si  = si;
		tdat[ithread].rn0 = rn0;
		tdat[ithread].rn1 = rn1;
		tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
	}
#endif

		bjmodn00 = (uint32*)(sm_ptr + 8);
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;		free((void *)_cy_i00); _cy_i00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;		free((void *)_cy_i01); _cy_i01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;		free((void *)_cy_i02); _cy_i02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;		free((void *)_cy_i03); _cy_i03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;		free((void *)_cy_i04); _cy_i04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;		free((void *)_cy_i05); _cy_i05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;		free((void *)_cy_i06); _cy_i06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;		free((void *)_cy_i07); _cy_i07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;		free((void *)_cy_i08); _cy_i08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;		free((void *)_cy_i09); _cy_i09 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;		free((void *)_cy_i10); _cy_i10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;		free((void *)_cy_i11); _cy_i11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;		free((void *)_cy_i12); _cy_i12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;		free((void *)_cy_i13); _cy_i13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;		free((void *)_cy_i14); _cy_i14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;		free((void *)_cy_i15); _cy_i15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;		free((void *)_cy_i16); _cy_i16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;		free((void *)_cy_i17); _cy_i17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;		free((void *)_cy_i18); _cy_i18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;		free((void *)_cy_i19); _cy_i19 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;		free((void *)_cy_i20); _cy_i20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;		free((void *)_cy_i21); _cy_i21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;		free((void *)_cy_i22); _cy_i22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;		free((void *)_cy_i23); _cy_i23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;		free((void *)_cy_i24); _cy_i24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;		free((void *)_cy_i25); _cy_i25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;		free((void *)_cy_i26); _cy_i26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;		free((void *)_cy_i27); _cy_i27 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r27== 0x0);

		_cy_i00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i00== 0x0);
		_cy_i01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i01== 0x0);
		_cy_i02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i02== 0x0);
		_cy_i03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i03== 0x0);
		_cy_i04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i04== 0x0);
		_cy_i05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i05== 0x0);
		_cy_i06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i06== 0x0);
		_cy_i07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i07== 0x0);
		_cy_i08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i08== 0x0);
		_cy_i09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i09== 0x0);
		_cy_i10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i10== 0x0);
		_cy_i11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i11== 0x0);
		_cy_i12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i12== 0x0);
		_cy_i13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i13== 0x0);
		_cy_i14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i14== 0x0);
		_cy_i15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i15== 0x0);
		_cy_i16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i16== 0x0);
		_cy_i17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i17== 0x0);
		_cy_i18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i18== 0x0);
		_cy_i19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i19== 0x0);
		_cy_i20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i27== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix28_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix28_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = NDIVR/CY_THREADS;
		}
		else
		{
			jhi = NDIVR/CY_THREADS/2;
		}

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

#if FFT_DEBUG
	dbg_fname[0] += (char)(CY_THREADS - 1);	// 1-thread = "a.txt", 2-thread = "b.txt", etc.
	dbg_file = fopen(dbg_fname, "w");
	ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
	fprintf(dbg_file, "radix28_ditN_cy_dif1 DEBUG: fftlen = %d\n", n);
	fprintf(dbg_file,"CY_THREADS = %d\n", CY_THREADS);
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double rt = 1024.0*1024.0*1024.0*1024.0;
	for(j = 0; j < n; j++) {
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		a[j1  ] = rt*rng_isaac_rand_double_norm_pm1();
		a[j1+1] = rt*rng_isaac_rand_double_norm_pm1();
	}
#endif

/*...The radix-28 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;	_cy_i00[ithread] = 0;
		_cy_r01[ithread] = 0;	_cy_i01[ithread] = 0;
		_cy_r02[ithread] = 0;	_cy_i02[ithread] = 0;
		_cy_r03[ithread] = 0;	_cy_i03[ithread] = 0;
		_cy_r04[ithread] = 0;	_cy_i04[ithread] = 0;
		_cy_r05[ithread] = 0;	_cy_i05[ithread] = 0;
		_cy_r06[ithread] = 0;	_cy_i06[ithread] = 0;
		_cy_r07[ithread] = 0;	_cy_i07[ithread] = 0;
		_cy_r08[ithread] = 0;	_cy_i08[ithread] = 0;
		_cy_r09[ithread] = 0;	_cy_i09[ithread] = 0;
		_cy_r10[ithread] = 0;	_cy_i10[ithread] = 0;
		_cy_r11[ithread] = 0;	_cy_i11[ithread] = 0;
		_cy_r12[ithread] = 0;	_cy_i12[ithread] = 0;
		_cy_r13[ithread] = 0;	_cy_i13[ithread] = 0;
		_cy_r14[ithread] = 0;	_cy_i14[ithread] = 0;
		_cy_r15[ithread] = 0;	_cy_i15[ithread] = 0;
		_cy_r16[ithread] = 0;	_cy_i16[ithread] = 0;
		_cy_r17[ithread] = 0;	_cy_i17[ithread] = 0;
		_cy_r18[ithread] = 0;	_cy_i18[ithread] = 0;
		_cy_r19[ithread] = 0;	_cy_i19[ithread] = 0;
		_cy_r20[ithread] = 0;	_cy_i20[ithread] = 0;
		_cy_r21[ithread] = 0;	_cy_i21[ithread] = 0;
		_cy_r22[ithread] = 0;	_cy_i22[ithread] = 0;
		_cy_r23[ithread] = 0;	_cy_i23[ithread] = 0;
		_cy_r24[ithread] = 0;	_cy_i24[ithread] = 0;
		_cy_r25[ithread] = 0;	_cy_i25[ithread] = 0;
		_cy_r26[ithread] = 0;	_cy_i26[ithread] = 0;
		_cy_r27[ithread] = 0;	_cy_i27[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/

		khi = 1;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}

		/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
		so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
		*/
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*NDIVR/2) % nwt;
		MOD_ADD32(ii01,ii01,nwt,ii02);
		MOD_ADD32(ii02,ii01,nwt,ii03);
		MOD_ADD32(ii03,ii01,nwt,ii04);
		MOD_ADD32(ii04,ii01,nwt,ii05);
		MOD_ADD32(ii05,ii01,nwt,ii06);
		MOD_ADD32(ii06,ii01,nwt,ii07);
		MOD_ADD32(ii07,ii01,nwt,ii08);
		MOD_ADD32(ii08,ii01,nwt,ii09);
		MOD_ADD32(ii09,ii01,nwt,ii10);
		MOD_ADD32(ii10,ii01,nwt,ii11);
		MOD_ADD32(ii11,ii01,nwt,ii12);
		MOD_ADD32(ii12,ii01,nwt,ii13);
		MOD_ADD32(ii13,ii01,nwt,ii14);
		MOD_ADD32(ii14,ii01,nwt,ii15);
		MOD_ADD32(ii15,ii01,nwt,ii16);
		MOD_ADD32(ii16,ii01,nwt,ii17);
		MOD_ADD32(ii17,ii01,nwt,ii18);
		MOD_ADD32(ii18,ii01,nwt,ii19);
		MOD_ADD32(ii19,ii01,nwt,ii20);
		MOD_ADD32(ii20,ii01,nwt,ii21);
		MOD_ADD32(ii21,ii01,nwt,ii22);
		MOD_ADD32(ii22,ii01,nwt,ii23);
		MOD_ADD32(ii23,ii01,nwt,ii24);
		MOD_ADD32(ii24,ii01,nwt,ii25);
		MOD_ADD32(ii25,ii01,nwt,ii26);
		MOD_ADD32(ii26,ii01,nwt,ii27);
	}

	// In non-power-of-2-runlength case, both Mersenne and Fermat-mod share these next 2 loops:
	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	j = _bjmodnini[CY_THREADS];
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		MOD_ADD32(_bjmodn00[ithread], j, n, _bjmodn01[ithread]);
		MOD_ADD32(_bjmodn01[ithread], j, n, _bjmodn02[ithread]);
		MOD_ADD32(_bjmodn02[ithread], j, n, _bjmodn03[ithread]);
		MOD_ADD32(_bjmodn03[ithread], j, n, _bjmodn04[ithread]);
		MOD_ADD32(_bjmodn04[ithread], j, n, _bjmodn05[ithread]);
		MOD_ADD32(_bjmodn05[ithread], j, n, _bjmodn06[ithread]);
		MOD_ADD32(_bjmodn06[ithread], j, n, _bjmodn07[ithread]);
		MOD_ADD32(_bjmodn07[ithread], j, n, _bjmodn08[ithread]);
		MOD_ADD32(_bjmodn08[ithread], j, n, _bjmodn09[ithread]);
		MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn10[ithread]);
		MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);
		MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
		MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
		MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
		MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
		MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
		MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
		MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
		MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
		MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn20[ithread]);
		MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
		MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
		MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);
		MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
		MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
		MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
		MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);

		// Every (odd_radix)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn07[ithread] = n;
			_bjmodn14[ithread] = n;
			_bjmodn21[ithread] = n;
		}
	}

#ifdef USE_SSE2

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		/* Find the circular-index-shift described in the head-of-file comments of radix28_ditN_cy_dif1.c, by serching bjmodn01 ... bjmodn[nwt] for the one == bw: */
		if( _bjmodn01[0] == bw ) { wts_idx_incr = 1; };
		if( _bjmodn02[0] == bw ) { wts_idx_incr = 2; };
		if( _bjmodn03[0] == bw ) { wts_idx_incr = 3; };
		if( _bjmodn04[0] == bw ) { wts_idx_incr = 4; };
		if( _bjmodn05[0] == bw ) { wts_idx_incr = 5; };
		if( _bjmodn06[0] == bw ) { wts_idx_incr = 6; };

		ASSERT(HERE, wts_idx_incr != 0, "wts_idx_incr init failed!");

	  #ifndef USE_SCALAR_CARRY
		wts_idx_inc2 = (wts_idx_incr << 5);	/* In the SSE2 version, use icycle0-6 as actual address offsets, so wts_idx_incr includes a (<< 4) shift
											for the array-of-complex-doubles indexing, and another doubling to reflect the fact that the SSE2
											version of the loop is equivalent to 2 scalar loop executions, i.e. corresponds to 2 scalar-code
											incrementations of the icycle indices. */
		wts_idx_inc2 %= (nwt << 4);	/* Need an extra mod since 2*wts_idx_incr may be >= nwt */
	  #endif

		wts_idx_incr -= nwt;	/* Subtract nwt from the increments to ease fast-mod */
	  #ifndef USE_SCALAR_CARRY
		wts_idx_inc2 -= (nwt << 4);	/* Need an extra mod since 2*wts_idx_incr may be >= nwt */
	  #endif

		/* In this init section, both scalar and sse2 use the icycle values as array indices: */
		icycle0 = _bjmodn00[0] > sw;
		icycle1 = _bjmodn01[0] > sw;
		icycle2 = _bjmodn02[0] > sw;
		icycle3 = _bjmodn03[0] > sw;
		icycle4 = _bjmodn04[0] > sw;
		icycle5 = _bjmodn05[0] > sw;
		icycle6 = _bjmodn06[0] > sw;
	  /* Need this both in scalar mode and to ease the SSE2-array init */
		bs_arr[0] = base[icycle0];	bsinv_arr[0] = baseinv[icycle0];
		bs_arr[1] = base[icycle1];	bsinv_arr[1] = baseinv[icycle1];
		bs_arr[2] = base[icycle2];	bsinv_arr[2] = baseinv[icycle2];
		bs_arr[3] = base[icycle3];	bsinv_arr[3] = baseinv[icycle3];
		bs_arr[4] = base[icycle4];	bsinv_arr[4] = baseinv[icycle4];
		bs_arr[5] = base[icycle5];	bsinv_arr[5] = baseinv[icycle5];
		bs_arr[6] = base[icycle6];	bsinv_arr[6] = baseinv[icycle6];

		/* Now that are done using icycle indices as temps for the (bjmodn > sw) values, give them their proper starting values: */
		/* In this init section, both scalar and sse2 use the icycle values as array indices: */
		icycle0 = 0x0;
		icycle1 = 0x1;
		icycle2 = 0x2;
		icycle3 = 0x3;
		icycle4 = 0x4;
		icycle5 = 0x5;
		icycle6 = 0x6;

		/* Need this both in scalar mode and to ease the SSE2-array init */
		wt_arr[0] = wt0[ii00];	wtinv_arr[0] = scale*wt1[ii00];
		wt_arr[1] = wt0[ii01];	wtinv_arr[1] = scale*wt1[ii01];
		wt_arr[2] = wt0[ii02];	wtinv_arr[2] = scale*wt1[ii02];
		wt_arr[3] = wt0[ii03];	wtinv_arr[3] = scale*wt1[ii03];
		wt_arr[4] = wt0[ii04];	wtinv_arr[4] = scale*wt1[ii04];
		wt_arr[5] = wt0[ii05];	wtinv_arr[5] = scale*wt1[ii05];
		wt_arr[6] = wt0[ii06];	wtinv_arr[6] = scale*wt1[ii06];

	  #ifndef USE_SCALAR_CARRY
		tmp = half_arr;
		tmp->re = wt_arr[icycle0];	++tmp;
		tmp->re = wt_arr[icycle1];	++tmp;
		tmp->re = wt_arr[icycle2];	++tmp;
		tmp->re = wt_arr[icycle3];	++tmp;
		tmp->re = wt_arr[icycle4];	++tmp;
		tmp->re = wt_arr[icycle5];	++tmp;
		tmp->re = wt_arr[icycle6];	++tmp;

		tmp->re = wtinv_arr[icycle0];	++tmp;
		tmp->re = wtinv_arr[icycle1];	++tmp;
		tmp->re = wtinv_arr[icycle2];	++tmp;
		tmp->re = wtinv_arr[icycle3];	++tmp;
		tmp->re = wtinv_arr[icycle4];	++tmp;
		tmp->re = wtinv_arr[icycle5];	++tmp;
		tmp->re = wtinv_arr[icycle6];	++tmp;

		/* Now set the imaginary parts to the values corresponding to the 2nd of each pair of scalar-mode loop passes.
		Use this sequence for mod-add, as it is faster than general-mod '% nwt'. The reason we do not use the MOD_ADD32
		macor here is that wts_idx_incr is precomputed constant, so we have also pre-subtracted the modulus nwt from it:
		*/
	   #if 0
		jcycle0 = icycle0 + wts_idx_incr;		jcycle0 += ( (-(int)((uint32)jcycle0 >> 31)) & nwt);
		jcycle1 = icycle1 + wts_idx_incr;		jcycle1 += ( (-(int)((uint32)jcycle1 >> 31)) & nwt);
		jcycle2 = icycle2 + wts_idx_incr;		jcycle2 += ( (-(int)((uint32)jcycle2 >> 31)) & nwt);
		jcycle3 = icycle3 + wts_idx_incr;		jcycle3 += ( (-(int)((uint32)jcycle3 >> 31)) & nwt);
		jcycle4 = icycle4 + wts_idx_incr;		jcycle4 += ( (-(int)((uint32)jcycle4 >> 31)) & nwt);
		jcycle5 = icycle5 + wts_idx_incr;		jcycle5 += ( (-(int)((uint32)jcycle5 >> 31)) & nwt);
		jcycle6 = icycle6 + wts_idx_incr;		jcycle6 += ( (-(int)((uint32)jcycle6 >> 31)) & nwt);
	   #else
		jcycle0 = icycle0 + wts_idx_incr;		jcycle0 += ( (-(jcycle0 < 0)) & nwt);
		jcycle1 = icycle1 + wts_idx_incr;		jcycle1 += ( (-(jcycle1 < 0)) & nwt);
		jcycle2 = icycle2 + wts_idx_incr;		jcycle2 += ( (-(jcycle2 < 0)) & nwt);
		jcycle3 = icycle3 + wts_idx_incr;		jcycle3 += ( (-(jcycle3 < 0)) & nwt);
		jcycle4 = icycle4 + wts_idx_incr;		jcycle4 += ( (-(jcycle4 < 0)) & nwt);
		jcycle5 = icycle5 + wts_idx_incr;		jcycle5 += ( (-(jcycle5 < 0)) & nwt);
		jcycle6 = icycle6 + wts_idx_incr;		jcycle6 += ( (-(jcycle6 < 0)) & nwt);
	   #endif
		tmp = half_arr;
		tmp->im = wt_arr[jcycle0];	++tmp;
		tmp->im = wt_arr[jcycle1];	++tmp;
		tmp->im = wt_arr[jcycle2];	++tmp;
		tmp->im = wt_arr[jcycle3];	++tmp;
		tmp->im = wt_arr[jcycle4];	++tmp;
		tmp->im = wt_arr[jcycle5];	++tmp;
		tmp->im = wt_arr[jcycle6];	++tmp;

		tmp->im = wtinv_arr[jcycle0];	++tmp;
		tmp->im = wtinv_arr[jcycle1];	++tmp;
		tmp->im = wtinv_arr[jcycle2];	++tmp;
		tmp->im = wtinv_arr[jcycle3];	++tmp;
		tmp->im = wtinv_arr[jcycle4];	++tmp;
		tmp->im = wtinv_arr[jcycle5];	++tmp;
		tmp->im = wtinv_arr[jcycle6];	++tmp;

		tmp = half_arr + odd_radix*2;	/* Put the base-mini-arrays right after the weights */

		/* Because we apply doubled weights to data arranged as [a.re,b.re],[a.im,b.im] but apply doubled base
		multipliers to shuffled data [a.re,a.im],[b.re,b.im] (i.e. shuffled to yield same data layout as in the scalar
		case), the weights need to have disparate real and imag parts, whereas the base/baseinv terms do not: */
		tmp->re = tmp->im = bs_arr[0];	++tmp;
		tmp->re = tmp->im = bs_arr[1];	++tmp;
		tmp->re = tmp->im = bs_arr[2];	++tmp;
		tmp->re = tmp->im = bs_arr[3];	++tmp;
		tmp->re = tmp->im = bs_arr[4];	++tmp;
		tmp->re = tmp->im = bs_arr[5];	++tmp;
		tmp->re = tmp->im = bs_arr[6];	++tmp;

		tmp->re = tmp->im = bsinv_arr[0];	++tmp;
		tmp->re = tmp->im = bsinv_arr[1];	++tmp;
		tmp->re = tmp->im = bsinv_arr[2];	++tmp;
		tmp->re = tmp->im = bsinv_arr[3];	++tmp;
		tmp->re = tmp->im = bsinv_arr[4];	++tmp;
		tmp->re = tmp->im = bsinv_arr[5];	++tmp;
		tmp->re = tmp->im = bsinv_arr[6];	++tmp;

		icycle0 <<= 4;		jcycle0 <<= 4;
		icycle1 <<= 4;		jcycle1 <<= 4;
		icycle2 <<= 4;		jcycle2 <<= 4;
		icycle3 <<= 4;		jcycle3 <<= 4;
		icycle4 <<= 4;		jcycle4 <<= 4;
		icycle5 <<= 4;		jcycle5 <<= 4;
		icycle6 <<= 4;		jcycle6 <<= 4;
	#endif	// USE_SCALAR_CARRY
	}
/*	fprintf(stderr, "radix28_ditN_cy_dif1: wts_idx_incr = %d\n", wts_idx_incr);*/

  #ifdef USE_PTHREAD

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

	// This sentinel-value check ensures that these inits occur just once, and moreover in full-pass mode,
	// which is needed to get the array-index-offset values of the icycle/jcycle indices right:
	} else if(tdat[0].wts_idx_inc2 == 0) {	// Fermat-mod
//	printf("Init per-thread cyclic indices, full_pass = %d\n",full_pass);
		tdat[0].wts_idx_inc2 = wts_idx_inc2;
		tdat[0].icycle0 = icycle0;
		tdat[0].icycle1 = icycle1;
		tdat[0].icycle2 = icycle2;
		tdat[0].icycle3 = icycle3;
		tdat[0].icycle4 = icycle4;
		tdat[0].icycle5 = icycle5;
		tdat[0].icycle6 = icycle6;

		tdat[0].jcycle0 = jcycle0;
		tdat[0].jcycle1 = jcycle1;
		tdat[0].jcycle2 = jcycle2;
		tdat[0].jcycle3 = jcycle3;
		tdat[0].jcycle4 = jcycle4;
		tdat[0].jcycle5 = jcycle5;
		tdat[0].jcycle6 = jcycle6;
		// For remaining threads, simulate the loop-evolution of the above indices:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			jstart = _jstart[ithread];
			jhi    = _jhi[ithread];
			for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
			{
				for(j = jstart; j < jhi; j += 4)
				{
					icycle0 += wts_idx_inc2;		icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt16);
					icycle1 += wts_idx_inc2;		icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt16);
					icycle2 += wts_idx_inc2;		icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt16);
					icycle3 += wts_idx_inc2;		icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt16);
					icycle4 += wts_idx_inc2;		icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt16);
					icycle5 += wts_idx_inc2;		icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt16);
					icycle6 += wts_idx_inc2;		icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt16);

					jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(int)((uint32)jcycle0 >> 31)) & nwt16);
					jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(int)((uint32)jcycle1 >> 31)) & nwt16);
					jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(int)((uint32)jcycle2 >> 31)) & nwt16);
					jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(int)((uint32)jcycle3 >> 31)) & nwt16);
					jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(int)((uint32)jcycle4 >> 31)) & nwt16);
					jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(int)((uint32)jcycle5 >> 31)) & nwt16);
					jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(int)((uint32)jcycle6 >> 31)) & nwt16);
				}
			}
			tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
			tdat[ithread].icycle0 = icycle0;
			tdat[ithread].icycle1 = icycle1;
			tdat[ithread].icycle2 = icycle2;
			tdat[ithread].icycle3 = icycle3;
			tdat[ithread].icycle4 = icycle4;
			tdat[ithread].icycle5 = icycle5;
			tdat[ithread].icycle6 = icycle6;
	
			tdat[ithread].jcycle0 = jcycle0;
			tdat[ithread].jcycle1 = jcycle1;
			tdat[ithread].jcycle2 = jcycle2;
			tdat[ithread].jcycle3 = jcycle3;
			tdat[ithread].jcycle4 = jcycle4;
			tdat[ithread].jcycle5 = jcycle5;
			tdat[ithread].jcycle6 = jcycle6;
		}
	}

	/* Init thread 1-CY_THREADS's local stores and pointers: */
	s1p00r = __r0 + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
		memcpy(s1p00r, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
		s1p00r += cslots_in_local_store;
	}
  #endif

#endif	// USE_SSE2

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix28_ditN_cy_dif1: Cleanup Pass:\n");
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef USE_OMP
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,\
		n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,\
		rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,\
		a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,\
		a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,\
		bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,\
		cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,\
		cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27\
	) default(shared) schedule(static)
#endif

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");
	
		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

//	printf("cy28_process_chunk: thread %d tdat-init, khi = %d, jlo = %d, jhi = %d\n", ithread,tdat[ithread].khi,tdat[ithread].jstart,tdat[ithread].jhi);

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].s1p00r;
		ASSERT(HERE, ((tmp + 0x39)->re == sx0 && (tmp + 0x39)->im == sx0), "thread-local memcheck failed!");
		ASSERT(HERE, ((tmp + 0x5f)->re == crnd && (tmp + 0x5f)->im == crnd), "thread-local memcheck failed!");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			ASSERT(HERE, (tmp + 0x60+10)->re * (tmp + 0x60+14)->re == 1.0 && (tmp + 0x60+10)->im * (tmp + 0x60+14)->im == 1.0, "thread-local memcheck failed!");
			tdat[ithread].bjmodn00 = _bjmodn00[ithread];
			tdat[ithread].bjmodn01 = _bjmodn01[ithread];
			tdat[ithread].bjmodn02 = _bjmodn02[ithread];
			tdat[ithread].bjmodn03 = _bjmodn03[ithread];
			tdat[ithread].bjmodn04 = _bjmodn04[ithread];
			tdat[ithread].bjmodn05 = _bjmodn05[ithread];
			tdat[ithread].bjmodn06 = _bjmodn06[ithread];
			tdat[ithread].bjmodn07 = _bjmodn07[ithread];
			tdat[ithread].bjmodn08 = _bjmodn08[ithread];
			tdat[ithread].bjmodn09 = _bjmodn09[ithread];
			tdat[ithread].bjmodn10 = _bjmodn10[ithread];
			tdat[ithread].bjmodn11 = _bjmodn11[ithread];
			tdat[ithread].bjmodn12 = _bjmodn12[ithread];
			tdat[ithread].bjmodn13 = _bjmodn13[ithread];
			tdat[ithread].bjmodn14 = _bjmodn14[ithread];
			tdat[ithread].bjmodn15 = _bjmodn15[ithread];
			tdat[ithread].bjmodn16 = _bjmodn16[ithread];
			tdat[ithread].bjmodn17 = _bjmodn17[ithread];
			tdat[ithread].bjmodn18 = _bjmodn18[ithread];
			tdat[ithread].bjmodn19 = _bjmodn19[ithread];
			tdat[ithread].bjmodn20 = _bjmodn20[ithread];
			tdat[ithread].bjmodn21 = _bjmodn21[ithread];
			tdat[ithread].bjmodn22 = _bjmodn22[ithread];
			tdat[ithread].bjmodn23 = _bjmodn23[ithread];
			tdat[ithread].bjmodn24 = _bjmodn24[ithread];
			tdat[ithread].bjmodn25 = _bjmodn25[ithread];
			tdat[ithread].bjmodn26 = _bjmodn26[ithread];
			tdat[ithread].bjmodn27 = _bjmodn27[ithread];
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			ASSERT(HERE, (tmp + 0x60)->re * (tmp + 0x60+odd_radix)->re == scale && (tmp + 0x60)->im * (tmp + 0x60+odd_radix)->im == scale, "thread-local memcheck failed!");
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];	tdat[ithread].cy_i00 = _cy_i00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];	tdat[ithread].cy_i01 = _cy_i01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];	tdat[ithread].cy_i02 = _cy_i02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];	tdat[ithread].cy_i03 = _cy_i03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];	tdat[ithread].cy_i04 = _cy_i04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];	tdat[ithread].cy_i05 = _cy_i05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];	tdat[ithread].cy_i06 = _cy_i06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];	tdat[ithread].cy_i07 = _cy_i07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];	tdat[ithread].cy_i08 = _cy_i08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];	tdat[ithread].cy_i09 = _cy_i09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];	tdat[ithread].cy_i10 = _cy_i10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];	tdat[ithread].cy_i11 = _cy_i11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];	tdat[ithread].cy_i12 = _cy_i12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];	tdat[ithread].cy_i13 = _cy_i13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];	tdat[ithread].cy_i14 = _cy_i14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];	tdat[ithread].cy_i15 = _cy_i15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];	tdat[ithread].cy_i16 = _cy_i16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];	tdat[ithread].cy_i17 = _cy_i17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];	tdat[ithread].cy_i18 = _cy_i18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];	tdat[ithread].cy_i19 = _cy_i19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];	tdat[ithread].cy_i20 = _cy_i20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];	tdat[ithread].cy_i21 = _cy_i21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];	tdat[ithread].cy_i22 = _cy_i22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];	tdat[ithread].cy_i23 = _cy_i23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];	tdat[ithread].cy_i24 = _cy_i24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];	tdat[ithread].cy_i25 = _cy_i25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];	tdat[ithread].cy_i26 = _cy_i26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];	tdat[ithread].cy_i27 = _cy_i27[ithread];
		}
	}
#endif

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
	#ifdef USE_SSE2
		max_err->re = 0.0;	max_err->im = 0.0;
	#endif

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];
		*bjmodn01 = _bjmodn01[ithread];
		*bjmodn02 = _bjmodn02[ithread];
		*bjmodn03 = _bjmodn03[ithread];
		*bjmodn04 = _bjmodn04[ithread];
		*bjmodn05 = _bjmodn05[ithread];
		*bjmodn06 = _bjmodn06[ithread];
		*bjmodn07 = _bjmodn07[ithread];
		*bjmodn08 = _bjmodn08[ithread];
		*bjmodn09 = _bjmodn09[ithread];
		*bjmodn10 = _bjmodn10[ithread];
		*bjmodn11 = _bjmodn11[ithread];
		*bjmodn12 = _bjmodn12[ithread];
		*bjmodn13 = _bjmodn13[ithread];
		*bjmodn14 = _bjmodn14[ithread];
		*bjmodn15 = _bjmodn15[ithread];
		*bjmodn16 = _bjmodn16[ithread];
		*bjmodn17 = _bjmodn17[ithread];
		*bjmodn18 = _bjmodn18[ithread];
		*bjmodn19 = _bjmodn19[ithread];
		*bjmodn20 = _bjmodn20[ithread];
		*bjmodn21 = _bjmodn21[ithread];
		*bjmodn22 = _bjmodn22[ithread];
		*bjmodn23 = _bjmodn23[ithread];
		*bjmodn24 = _bjmodn24[ithread];
		*bjmodn25 = _bjmodn25[ithread];
		*bjmodn26 = _bjmodn26[ithread];
		*bjmodn27 = _bjmodn27[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];
		bjmodn01 = _bjmodn01[ithread];
		bjmodn02 = _bjmodn02[ithread];
		bjmodn03 = _bjmodn03[ithread];
		bjmodn04 = _bjmodn04[ithread];
		bjmodn05 = _bjmodn05[ithread];
		bjmodn06 = _bjmodn06[ithread];
		bjmodn07 = _bjmodn07[ithread];
		bjmodn08 = _bjmodn08[ithread];
		bjmodn09 = _bjmodn09[ithread];
		bjmodn10 = _bjmodn10[ithread];
		bjmodn11 = _bjmodn11[ithread];
		bjmodn12 = _bjmodn12[ithread];
		bjmodn13 = _bjmodn13[ithread];
		bjmodn14 = _bjmodn14[ithread];
		bjmodn15 = _bjmodn15[ithread];
		bjmodn16 = _bjmodn16[ithread];
		bjmodn17 = _bjmodn17[ithread];
		bjmodn18 = _bjmodn18[ithread];
		bjmodn19 = _bjmodn19[ithread];
		bjmodn20 = _bjmodn20[ithread];
		bjmodn21 = _bjmodn21[ithread];
		bjmodn22 = _bjmodn22[ithread];
		bjmodn23 = _bjmodn23[ithread];
		bjmodn24 = _bjmodn24[ithread];
		bjmodn25 = _bjmodn25[ithread];
		bjmodn26 = _bjmodn26[ithread];
		bjmodn27 = _bjmodn27[ithread];
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			col = _col[ithread];
			co2 = _co2[ithread];
			co3 = _co3[ithread];

			/* init carries	*/
		#if defined(USE_SSE2) && !defined(USE_SCALAR_CARRY)
			cy_r00->re = _cy_r00[ithread];
			cy_r00->im = _cy_r01[ithread];
			cy_r02->re = _cy_r02[ithread];
			cy_r02->im = _cy_r03[ithread];
			cy_r04->re = _cy_r04[ithread];
			cy_r04->im = _cy_r05[ithread];
			cy_r06->re = _cy_r06[ithread];
			cy_r06->im = _cy_r07[ithread];
			cy_r08->re = _cy_r08[ithread];
			cy_r08->im = _cy_r09[ithread];
			cy_r10->re = _cy_r10[ithread];
			cy_r10->im = _cy_r11[ithread];
			cy_r12->re = _cy_r12[ithread];
			cy_r12->im = _cy_r13[ithread];
			cy_r14->re = _cy_r14[ithread];
			cy_r14->im = _cy_r15[ithread];
			cy_r16->re = _cy_r16[ithread];
			cy_r16->im = _cy_r17[ithread];
			cy_r18->re = _cy_r18[ithread];
			cy_r18->im = _cy_r19[ithread];
			cy_r20->re = _cy_r20[ithread];
			cy_r20->im = _cy_r21[ithread];
			cy_r22->re = _cy_r22[ithread];
			cy_r22->im = _cy_r23[ithread];
			cy_r24->re = _cy_r24[ithread];
			cy_r24->im = _cy_r25[ithread];
			cy_r26->re = _cy_r26[ithread];
			cy_r26->im = _cy_r27[ithread];
		#else
			/* init carries	*/
			cy_r00 = _cy_r00[ithread];
			cy_r01 = _cy_r01[ithread];
			cy_r02 = _cy_r02[ithread];
			cy_r03 = _cy_r03[ithread];
			cy_r04 = _cy_r04[ithread];
			cy_r05 = _cy_r05[ithread];
			cy_r06 = _cy_r06[ithread];
			cy_r07 = _cy_r07[ithread];
			cy_r08 = _cy_r08[ithread];
			cy_r09 = _cy_r09[ithread];
			cy_r10 = _cy_r10[ithread];
			cy_r11 = _cy_r11[ithread];
			cy_r12 = _cy_r12[ithread];
			cy_r13 = _cy_r13[ithread];
			cy_r14 = _cy_r14[ithread];
			cy_r15 = _cy_r15[ithread];
			cy_r16 = _cy_r16[ithread];
			cy_r17 = _cy_r17[ithread];
			cy_r18 = _cy_r18[ithread];
			cy_r19 = _cy_r19[ithread];
			cy_r20 = _cy_r20[ithread];
			cy_r21 = _cy_r21[ithread];
			cy_r22 = _cy_r22[ithread];
			cy_r23 = _cy_r23[ithread];
			cy_r24 = _cy_r24[ithread];
			cy_r25 = _cy_r25[ithread];
			cy_r26 = _cy_r26[ithread];
			cy_r27 = _cy_r27[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#if defined(USE_SSE2) && !defined(USE_SCALAR_CARRY)
			cy_r00->re = _cy_r00[ithread];	cy_r00->im = _cy_i00[ithread];
			cy_r02->re = _cy_r01[ithread];	cy_r02->im = _cy_i01[ithread];
			cy_r04->re = _cy_r02[ithread];	cy_r04->im = _cy_i02[ithread];
			cy_r06->re = _cy_r03[ithread];	cy_r06->im = _cy_i03[ithread];
			cy_r08->re = _cy_r04[ithread];	cy_r08->im = _cy_i04[ithread];
			cy_r10->re = _cy_r05[ithread];	cy_r10->im = _cy_i05[ithread];
			cy_r12->re = _cy_r06[ithread];	cy_r12->im = _cy_i06[ithread];
			cy_r14->re = _cy_r07[ithread];	cy_r14->im = _cy_i07[ithread];
			cy_r16->re = _cy_r08[ithread];	cy_r16->im = _cy_i08[ithread];
			cy_r18->re = _cy_r09[ithread];	cy_r18->im = _cy_i09[ithread];
			cy_r20->re = _cy_r10[ithread];	cy_r20->im = _cy_i10[ithread];
			cy_r22->re = _cy_r11[ithread];	cy_r22->im = _cy_i11[ithread];
			cy_r24->re = _cy_r12[ithread];	cy_r24->im = _cy_i12[ithread];
			cy_r26->re = _cy_r13[ithread];	cy_r26->im = _cy_i13[ithread];
			cy_i00->re = _cy_r14[ithread];	cy_i00->im = _cy_i14[ithread];
			cy_i02->re = _cy_r15[ithread];	cy_i02->im = _cy_i15[ithread];
			cy_i04->re = _cy_r16[ithread];	cy_i04->im = _cy_i16[ithread];
			cy_i06->re = _cy_r17[ithread];	cy_i06->im = _cy_i17[ithread];
			cy_i08->re = _cy_r18[ithread];	cy_i08->im = _cy_i18[ithread];
			cy_i10->re = _cy_r19[ithread];	cy_i10->im = _cy_i19[ithread];
			cy_i12->re = _cy_r20[ithread];	cy_i12->im = _cy_i20[ithread];
			cy_i14->re = _cy_r21[ithread];	cy_i14->im = _cy_i21[ithread];
			cy_i16->re = _cy_r22[ithread];	cy_i16->im = _cy_i22[ithread];
			cy_i18->re = _cy_r23[ithread];	cy_i18->im = _cy_i23[ithread];
			cy_i20->re = _cy_r24[ithread];	cy_i20->im = _cy_i24[ithread];
			cy_i22->re = _cy_r25[ithread];	cy_i22->im = _cy_i25[ithread];
			cy_i24->re = _cy_r26[ithread];	cy_i24->im = _cy_i26[ithread];
			cy_i26->re = _cy_r27[ithread];	cy_i26->im = _cy_i27[ithread];
		#else
			/* init carries	*/
			cy_r00 = _cy_r00[ithread];	cy_i00 = _cy_i00[ithread];
			cy_r01 = _cy_r01[ithread];	cy_i01 = _cy_i01[ithread];
			cy_r02 = _cy_r02[ithread];	cy_i02 = _cy_i02[ithread];
			cy_r03 = _cy_r03[ithread];	cy_i03 = _cy_i03[ithread];
			cy_r04 = _cy_r04[ithread];	cy_i04 = _cy_i04[ithread];
			cy_r05 = _cy_r05[ithread];	cy_i05 = _cy_i05[ithread];
			cy_r06 = _cy_r06[ithread];	cy_i06 = _cy_i06[ithread];
			cy_r07 = _cy_r07[ithread];	cy_i07 = _cy_i07[ithread];
			cy_r08 = _cy_r08[ithread];	cy_i08 = _cy_i08[ithread];
			cy_r09 = _cy_r09[ithread];	cy_i09 = _cy_i09[ithread];
			cy_r10 = _cy_r10[ithread];	cy_i10 = _cy_i10[ithread];
			cy_r11 = _cy_r11[ithread];	cy_i11 = _cy_i11[ithread];
			cy_r12 = _cy_r12[ithread];	cy_i12 = _cy_i12[ithread];
			cy_r13 = _cy_r13[ithread];	cy_i13 = _cy_i13[ithread];
			cy_r14 = _cy_r14[ithread];	cy_i14 = _cy_i14[ithread];
			cy_r15 = _cy_r15[ithread];	cy_i15 = _cy_i15[ithread];
			cy_r16 = _cy_r16[ithread];	cy_i16 = _cy_i16[ithread];
			cy_r17 = _cy_r17[ithread];	cy_i17 = _cy_i17[ithread];
			cy_r18 = _cy_r18[ithread];	cy_i18 = _cy_i18[ithread];
			cy_r19 = _cy_r19[ithread];	cy_i19 = _cy_i19[ithread];
			cy_r20 = _cy_r20[ithread];	cy_i20 = _cy_i20[ithread];
			cy_r21 = _cy_r21[ithread];	cy_i21 = _cy_i21[ithread];
			cy_r22 = _cy_r22[ithread];	cy_i22 = _cy_i22[ithread];
			cy_r23 = _cy_r23[ithread];	cy_i23 = _cy_i23[ithread];
			cy_r24 = _cy_r24[ithread];	cy_i24 = _cy_i24[ithread];
			cy_r25 = _cy_r25[ithread];	cy_i25 = _cy_i25[ithread];
			cy_r26 = _cy_r26[ithread];	cy_i26 = _cy_i26[ithread];
			cy_r27 = _cy_r27[ithread];	cy_i27 = _cy_i27[ithread];
		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
		#ifdef USE_SSE2
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
		#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 = (j & mask01) + br4[j&3];
		#else
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 =  j;
		#endif
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

		#ifdef DEBUG_SSE2
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p04;	jp = j2 + p04;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p08;	jp = j2 + p08;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p12;	jp = j2 + p12;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p16;	jp = j2 + p16;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p20;	jp = j2 + p20;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p24;	jp = j2 + p24;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			jt = j1;		jp = j2;
			fprintf(stderr, "radix28_carry: A_in[00] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[01] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[02] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[03] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p04;	jp = j2 + p04;
			fprintf(stderr, "radix28_carry: A_in[04] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[05] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[06] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[07] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p08;	jp = j2 + p08;
			fprintf(stderr, "radix28_carry: A_in[08] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[09] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[10] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[11] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p12;	jp = j2 + p12;
			fprintf(stderr, "radix28_carry: A_in[12] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[13] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[14] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[15] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p16;	jp = j2 + p16;
			fprintf(stderr, "radix28_carry: A_in[16] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[17] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[18] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[19] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p20;	jp = j2 + p20;
			fprintf(stderr, "radix28_carry: A_in[20] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[21] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[22] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[23] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p24;	jp = j2 + p24;
			fprintf(stderr, "radix28_carry: A_in[24] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[25] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[26] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[27] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);
			fprintf(stderr, "\n");
		#endif

		/*...The radix-28 DIT pass is here:	*/
#ifdef CTIME
	clock2 = clock();
#endif
		/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
						  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
						  to properly re-use the ajp1 variables in the carry-pass version of this routine.
		*/
		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

			/* Since doing radix-7 in-place here, outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p00r,s1p03r,s1p02r,s1p01r)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p04r,s1p07r,s1p06r,s1p05r)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p08r,s1p11r,s1p10r,s1p09r)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p12r,s1p15r,s1p14r,s1p13r)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p16r,s1p19r,s1p18r,s1p17r)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p20r,s1p23r,s1p22r,s1p21r)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p24r,s1p27r,s1p26r,s1p25r)

			/*...and now do 4 radix-7 transforms...*/

			SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r);
			SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r);
			SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r);
			SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r);

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX28_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

		  #endif

			/* DEBUG: To do carries using the scalar [non-SSE2] carry macros, move appropriate half of the xmm temps into a1p's: */
			#ifdef USE_SCALAR_CARRY

			SSE_LOOP:

				if((j&3) == 0)
				{
					a1p00r = s1p00r->re;	a1p00i = s1p00i->re;
					a1p01r = s1p01r->re;	a1p01i = s1p01i->re;
					a1p02r = s1p02r->re;	a1p02i = s1p02i->re;
					a1p03r = s1p03r->re;	a1p03i = s1p03i->re;
					a1p04r = s1p04r->re;	a1p04i = s1p04i->re;
					a1p05r = s1p05r->re;	a1p05i = s1p05i->re;
					a1p06r = s1p06r->re;	a1p06i = s1p06i->re;
					a1p07r = s1p07r->re;	a1p07i = s1p07i->re;
					a1p08r = s1p08r->re;	a1p08i = s1p08i->re;
					a1p09r = s1p09r->re;	a1p09i = s1p09i->re;
					a1p10r = s1p10r->re;	a1p10i = s1p10i->re;
					a1p11r = s1p11r->re;	a1p11i = s1p11i->re;
					a1p12r = s1p12r->re;	a1p12i = s1p12i->re;
					a1p13r = s1p13r->re;	a1p13i = s1p13i->re;
					a1p14r = s1p14r->re;	a1p14i = s1p14i->re;
					a1p15r = s1p15r->re;	a1p15i = s1p15i->re;
					a1p16r = s1p16r->re;	a1p16i = s1p16i->re;
					a1p17r = s1p17r->re;	a1p17i = s1p17i->re;
					a1p18r = s1p18r->re;	a1p18i = s1p18i->re;
					a1p19r = s1p19r->re;	a1p19i = s1p19i->re;
					a1p20r = s1p20r->re;	a1p20i = s1p20i->re;
					a1p21r = s1p21r->re;	a1p21i = s1p21i->re;
					a1p22r = s1p22r->re;	a1p22i = s1p22i->re;
					a1p23r = s1p23r->re;	a1p23i = s1p23i->re;
					a1p24r = s1p24r->re;	a1p24i = s1p24i->re;
					a1p25r = s1p25r->re;	a1p25i = s1p25i->re;
					a1p26r = s1p26r->re;	a1p26i = s1p26i->re;
					a1p27r = s1p27r->re;	a1p27i = s1p27i->re;
				}
				else
				{
					a1p00r = s1p00r->im;	a1p00i = s1p00i->im;
					a1p01r = s1p01r->im;	a1p01i = s1p01i->im;
					a1p02r = s1p02r->im;	a1p02i = s1p02i->im;
					a1p03r = s1p03r->im;	a1p03i = s1p03i->im;
					a1p04r = s1p04r->im;	a1p04i = s1p04i->im;
					a1p05r = s1p05r->im;	a1p05i = s1p05i->im;
					a1p06r = s1p06r->im;	a1p06i = s1p06i->im;
					a1p07r = s1p07r->im;	a1p07i = s1p07i->im;
					a1p08r = s1p08r->im;	a1p08i = s1p08i->im;
					a1p09r = s1p09r->im;	a1p09i = s1p09i->im;
					a1p10r = s1p10r->im;	a1p10i = s1p10i->im;
					a1p11r = s1p11r->im;	a1p11i = s1p11i->im;
					a1p12r = s1p12r->im;	a1p12i = s1p12i->im;
					a1p13r = s1p13r->im;	a1p13i = s1p13i->im;
					a1p14r = s1p14r->im;	a1p14i = s1p14i->im;
					a1p15r = s1p15r->im;	a1p15i = s1p15i->im;
					a1p16r = s1p16r->im;	a1p16i = s1p16i->im;
					a1p17r = s1p17r->im;	a1p17i = s1p17i->im;
					a1p18r = s1p18r->im;	a1p18i = s1p18i->im;
					a1p19r = s1p19r->im;	a1p19i = s1p19i->im;
					a1p20r = s1p20r->im;	a1p20i = s1p20i->im;
					a1p21r = s1p21r->im;	a1p21i = s1p21i->im;
					a1p22r = s1p22r->im;	a1p22i = s1p22i->im;
					a1p23r = s1p23r->im;	a1p23i = s1p23i->im;
					a1p24r = s1p24r->im;	a1p24i = s1p24i->im;
					a1p25r = s1p25r->im;	a1p25i = s1p25i->im;
					a1p26r = s1p26r->im;	a1p26i = s1p26i->im;
					a1p27r = s1p27r->im;	a1p27i = s1p27i->im;
				}
			#endif

		#else	/* !USE_SSE2 */

		  /*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
							 /*                                      outputs                                      */ /*                          inputs                           */
			RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

		  /*...and now do 4 radix-7 transforms...*/
		  #if LO_ADD
							 /*                                                   inputs                                                  */ /*               intermediates              */ /*                                                                     outputs                                                                         */
			RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		  #else
			RADIX_07_DFT_NUSS(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
		  #endif

		#endif	/* USE_SSE2 */

#ifdef CTIME
	clock3 = clock();
	dt_fwd += (double)(clock3 - clock2);
	clock2 = clock3;
#endif
		/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

		#if defined(USE_SCALAR_CARRY) || !defined(USE_SSE2)

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			#if defined(USE_SCALAR_CARRY)

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,*bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,*bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,*bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,*bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,*bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,*bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,*bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,*bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,*bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,*bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,*bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,*bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,*bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,*bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,*bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,*bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,*bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,*bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,*bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,*bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,*bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,*bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,*bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,*bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,*bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,*bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,*bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,*bjmodn27,27);

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#else	/* SSE2 mode uses pointers for the bjmodn's, non-SSE2 uses scalars: */

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/

			#endif

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
			}
			else	/* MODULUS_TYPE_FERMAT */
			{
			#if defined(USE_SCALAR_CARRY)

			  #if 0
				if((iter == 1) && (j < 32) &&!full_pass)
				{
					fprintf(stderr, "J = %d:\n",j);
					fprintf(stderr, "radix28_carry_out: icycle0 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle0,wt_arr[icycle0],wtinv_arr[icycle0],bs_arr[icycle0]);
					fprintf(stderr, "radix28_carry_out: icycle1 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle1,wt_arr[icycle1],wtinv_arr[icycle1],bs_arr[icycle1]);
					fprintf(stderr, "radix28_carry_out: icycle2 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle2,wt_arr[icycle2],wtinv_arr[icycle2],bs_arr[icycle2]);
					fprintf(stderr, "radix28_carry_out: icycle3 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle3,wt_arr[icycle3],wtinv_arr[icycle3],bs_arr[icycle3]);
					fprintf(stderr, "radix28_carry_out: icycle4 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle4,wt_arr[icycle4],wtinv_arr[icycle4],bs_arr[icycle4]);
					fprintf(stderr, "radix28_carry_out: icycle5 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle5,wt_arr[icycle5],wtinv_arr[icycle5],bs_arr[icycle5]);
					fprintf(stderr, "radix28_carry_out: icycle6 = %1d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle6,wt_arr[icycle6],wtinv_arr[icycle6],bs_arr[icycle6]);
				}
			  #endif
				ntmp = 0;
				fermat_carry_norm_errcheckB(a1p00r,a1p00i,cy_r00,cy_i00,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p01r,a1p01i,cy_r01,cy_i01,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p02r,a1p02i,cy_r02,cy_i02,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p03r,a1p03i,cy_r03,cy_i03,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p04r,a1p04i,cy_r04,cy_i04,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p05r,a1p05i,cy_r05,cy_i05,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p06r,a1p06i,cy_r06,cy_i06,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p07r,a1p07i,cy_r07,cy_i07,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p08r,a1p08i,cy_r08,cy_i08,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p09r,a1p09i,cy_r09,cy_i09,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p10r,a1p10i,cy_r10,cy_i10,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p11r,a1p11i,cy_r11,cy_i11,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p12r,a1p12i,cy_r12,cy_i12,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p13r,a1p13i,cy_r13,cy_i13,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p14r,a1p14i,cy_r14,cy_i14,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p15r,a1p15i,cy_r15,cy_i15,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p16r,a1p16i,cy_r16,cy_i16,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p17r,a1p17i,cy_r17,cy_i17,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p18r,a1p18i,cy_r18,cy_i18,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p19r,a1p19i,cy_r19,cy_i19,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p20r,a1p20i,cy_r20,cy_i20,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p21r,a1p21i,cy_r21,cy_i21,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p22r,a1p22i,cy_r22,cy_i22,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p23r,a1p23i,cy_r23,cy_i23,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p24r,a1p24i,cy_r24,cy_i24,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p25r,a1p25i,cy_r25,cy_i25,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p26r,a1p26i,cy_r26,cy_i26,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p27r,a1p27i,cy_r27,cy_i27,icycle6,ntmp,NRTM1,NRT_BITS);

				icycle0 += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
				icycle1 += wts_idx_incr;
				icycle2 += wts_idx_incr;
				icycle3 += wts_idx_incr;
				icycle4 += wts_idx_incr;
				icycle5 += wts_idx_incr;
				icycle6 += wts_idx_incr;
				icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt);
				icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt);
				icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt);
				icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt);
				icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt);
				icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt);
				icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt);

			#else	/* SSE2 mode uses pointers for the bjmodn's, non-SSE2 uses scalars: */

			  #if 0
				if((iter == 1) && (j < 32) && full_pass)
				{
					fprintf(stderr, "J = %d:\n",j);
					fprintf(stderr, "radix28_carry_out: ii00 = %d, bjmodn00 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii00, bjmodn00, wt0[ii00],bjmodn00 > sw,base[bjmodn00 > sw]);
					fprintf(stderr, "radix28_carry_out: ii01 = %d, bjmodn01 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii01, bjmodn01, wt0[ii01],bjmodn01 > sw,base[bjmodn01 > sw]);
					fprintf(stderr, "radix28_carry_out: ii02 = %d, bjmodn02 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii02, bjmodn02, wt0[ii02],bjmodn02 > sw,base[bjmodn02 > sw]);
					fprintf(stderr, "radix28_carry_out: ii03 = %d, bjmodn03 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii03, bjmodn03, wt0[ii03],bjmodn03 > sw,base[bjmodn03 > sw]);
					fprintf(stderr, "radix28_carry_out: ii04 = %d, bjmodn04 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii04, bjmodn04, wt0[ii04],bjmodn04 > sw,base[bjmodn04 > sw]);
					fprintf(stderr, "radix28_carry_out: ii05 = %d, bjmodn05 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii05, bjmodn05, wt0[ii05],bjmodn05 > sw,base[bjmodn05 > sw]);
					fprintf(stderr, "radix28_carry_out: ii06 = %d, bjmodn06 = %d, wt = %6.4f, i = %d, base = %6.4f\n",ii06, bjmodn06, wt0[ii06],bjmodn06 > sw,base[bjmodn06 > sw]);
				}
			  #endif
				ntmp = 0;
				fermat_carry_norm_errcheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,ntmp,NRTM1,NRT_BITS);

			#endif
			}

		#elif defined(USE_SSE2)	/* !defined(USE_SCALAR_CARRY) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL				/* Updating i prior to the 2nd-7th macro calls allows use of the same 0_2B macro for all */
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);	// i =((uint32)(sw - *bjmodn04) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);	// i =((uint32)(sw - *bjmodn08) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);	// i =((uint32)(sw - *bjmodn12) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);	// i =((uint32)(sw - *bjmodn16) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);	// i =((uint32)(sw - *bjmodn20) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);	// i =((uint32)(sw - *bjmodn24) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				Dec 2011: Suspect this was a side effect of my gcc asm macros not including cc/memory in the clobber list, because
				the code now runs correctly without this hack ... but the code runs sign. faster with iy left in. So still "bizarre" but in a new way.
				*/
				if(j < 0)
				{
					fprintf(stderr, "Iter %3d\n",iter);
				}

		  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

		  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

		  #endif
				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			}
			else	/* Fermat-mod carry in SSE2 mode */
			{
			  #if 0
				if((iter == 1) && (j < 32) &&!full_pass)
				{
					tmp = half_arr;
					fprintf(stderr, "J = %d:\n",j);
					fprintf(stderr, "radix28_carry_out: icycle0 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle0>>4,(tmp+(icycle0>>4))->re,(tmp+(icycle0>>4)+7)->re,(tmp+(icycle0>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle1 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle1>>4,(tmp+(icycle1>>4))->re,(tmp+(icycle1>>4)+7)->re,(tmp+(icycle1>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle2 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle2>>4,(tmp+(icycle2>>4))->re,(tmp+(icycle2>>4)+7)->re,(tmp+(icycle2>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle3 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle3>>4,(tmp+(icycle3>>4))->re,(tmp+(icycle3>>4)+7)->re,(tmp+(icycle3>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle4 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle4>>4,(tmp+(icycle4>>4))->re,(tmp+(icycle4>>4)+7)->re,(tmp+(icycle4>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle5 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle5>>4,(tmp+(icycle5>>4))->re,(tmp+(icycle5>>4)+7)->re,(tmp+(icycle5>>4)+14)->re);
					fprintf(stderr, "radix28_carry_out: icycle6 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",icycle6>>4,(tmp+(icycle6>>4))->re,(tmp+(icycle6>>4)+7)->re,(tmp+(icycle6>>4)+14)->re);
					tmp = half_arr;
					fprintf(stderr, "J = %d:\n",j+2);
					fprintf(stderr, "radix28_carry_out: icycle0 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle0>>4,(tmp+(icycle0>>4))->im,(tmp+(icycle0>>4)+7)->im,(tmp+(jcycle0>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle1 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle1>>4,(tmp+(icycle1>>4))->im,(tmp+(icycle1>>4)+7)->im,(tmp+(jcycle1>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle2 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle2>>4,(tmp+(icycle2>>4))->im,(tmp+(icycle2>>4)+7)->im,(tmp+(jcycle2>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle3 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle3>>4,(tmp+(icycle3>>4))->im,(tmp+(icycle3>>4)+7)->im,(tmp+(jcycle3>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle4 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle4>>4,(tmp+(icycle4>>4))->im,(tmp+(icycle4>>4)+7)->im,(tmp+(jcycle4>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle5 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle5>>4,(tmp+(icycle5>>4))->im,(tmp+(icycle5>>4)+7)->im,(tmp+(jcycle5>>4)+14)->im);
					fprintf(stderr, "radix28_carry_out: icycle6 = %8d, wt = %6.4f, wi = %6.4f, base = %6.4f\n",jcycle6>>4,(tmp+(icycle6>>4))->im,(tmp+(icycle6>>4)+7)->im,(tmp+(jcycle6>>4)+14)->im);
				}
			  #endif
					/* Get the needed Nth root of -1: */
					add1 = &rn0[0];
					add2 = &rn1[0];

					idx_offset = j;
					idx_incr = NDIVR;

				#if defined(COMPILER_TYPE_MSVC)

				/* The cy_[r|i]_idx[A|B] names here are not meaningful, each simple stores one [re,im] carry pair,
				e.g. cy_r01 stores the carries our of [a0.re,a0.im], cy_r23 stores the carries our of [a1.re,a1.im], etc.
				Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
														2-vector				                                          Scalar
														--------	 		                                           -----------__ */
					SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r00,cy_i00 */
					SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r01,cy_i01 */
					SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r02,cy_i02 */
					SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r03,cy_i03 */
					SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r04,cy_i04 */
					SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r05,cy_i05 */
					SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r06,cy_i06 */
					SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r07,cy_i07 */
					SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r08,cy_i08 */
					SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r09,cy_i09 */
					SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r10,cy_i10 */
					SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r11,cy_i11 */
					SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r12,cy_i12 */
					SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r13,cy_i13 */
					SSE2_fermat_carry_norm_errcheck(s1p14r,cy_i00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r14,cy_i14 */
					SSE2_fermat_carry_norm_errcheck(s1p15r,cy_i02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r15,cy_i15 */
					SSE2_fermat_carry_norm_errcheck(s1p16r,cy_i04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r16,cy_i16 */
					SSE2_fermat_carry_norm_errcheck(s1p17r,cy_i06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r17,cy_i17 */
					SSE2_fermat_carry_norm_errcheck(s1p18r,cy_i08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r18,cy_i18 */
					SSE2_fermat_carry_norm_errcheck(s1p19r,cy_i10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r19,cy_i19 */
					SSE2_fermat_carry_norm_errcheck(s1p20r,cy_i12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r20,cy_i20 */
					SSE2_fermat_carry_norm_errcheck(s1p21r,cy_i14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r21,cy_i21 */
					SSE2_fermat_carry_norm_errcheck(s1p22r,cy_i16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r22,cy_i22 */
					SSE2_fermat_carry_norm_errcheck(s1p23r,cy_i18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r23,cy_i23 */
					SSE2_fermat_carry_norm_errcheck(s1p24r,cy_i20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r24,cy_i24 */
					SSE2_fermat_carry_norm_errcheck(s1p25r,cy_i22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r25,cy_i25 */
					SSE2_fermat_carry_norm_errcheck(s1p26r,cy_i24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r26,cy_i26 */
					SSE2_fermat_carry_norm_errcheck(s1p27r,cy_i26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r27,cy_i27 */
				#else
				  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit 
					SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p14r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p15r,cy_i02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p16r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p17r,cy_i06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p18r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p19r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p20r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p21r,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p22r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p23r,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p24r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p25r,cy_i22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p26r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p27r,cy_i26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				  #else
					SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck_X2(s1p14r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p16r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p18r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				  #endif
				#endif
				#if 1
					icycle0 += wts_idx_inc2;		icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt16);
					icycle1 += wts_idx_inc2;		icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt16);
					icycle2 += wts_idx_inc2;		icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt16);
					icycle3 += wts_idx_inc2;		icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt16);
					icycle4 += wts_idx_inc2;		icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt16);
					icycle5 += wts_idx_inc2;		icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt16);
					icycle6 += wts_idx_inc2;		icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt16);

					jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(int)((uint32)jcycle0 >> 31)) & nwt16);
					jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(int)((uint32)jcycle1 >> 31)) & nwt16);
					jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(int)((uint32)jcycle2 >> 31)) & nwt16);
					jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(int)((uint32)jcycle3 >> 31)) & nwt16);
					jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(int)((uint32)jcycle4 >> 31)) & nwt16);
					jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(int)((uint32)jcycle5 >> 31)) & nwt16);
					jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(int)((uint32)jcycle6 >> 31)) & nwt16);
				#else	// Above is actually faster in my MSVC builds
					icycle0 = icycle0 + wts_idx_inc2;		icycle0 += ( (-(icycle0 < 0)) & nwt16);
					icycle1 = icycle1 + wts_idx_inc2;		icycle1 += ( (-(icycle1 < 0)) & nwt16);
					icycle2 = icycle2 + wts_idx_inc2;		icycle2 += ( (-(icycle2 < 0)) & nwt16);
					icycle3 = icycle3 + wts_idx_inc2;		icycle3 += ( (-(icycle3 < 0)) & nwt16);
					icycle4 = icycle4 + wts_idx_inc2;		icycle4 += ( (-(icycle4 < 0)) & nwt16);
					icycle5 = icycle5 + wts_idx_inc2;		icycle5 += ( (-(icycle5 < 0)) & nwt16);
					icycle6 = icycle6 + wts_idx_inc2;		icycle6 += ( (-(icycle6 < 0)) & nwt16);

					jcycle0 = jcycle0 + wts_idx_inc2;		jcycle0 += ( (-(jcycle0 < 0)) & nwt16);
					jcycle1 = jcycle1 + wts_idx_inc2;		jcycle1 += ( (-(jcycle1 < 0)) & nwt16);
					jcycle2 = jcycle2 + wts_idx_inc2;		jcycle2 += ( (-(jcycle2 < 0)) & nwt16);
					jcycle3 = jcycle3 + wts_idx_inc2;		jcycle3 += ( (-(jcycle3 < 0)) & nwt16);
					jcycle4 = jcycle4 + wts_idx_inc2;		jcycle4 += ( (-(jcycle4 < 0)) & nwt16);
					jcycle5 = jcycle5 + wts_idx_inc2;		jcycle5 += ( (-(jcycle5 < 0)) & nwt16);
					jcycle6 = jcycle6 + wts_idx_inc2;		jcycle6 += ( (-(jcycle6 < 0)) & nwt16);
				  #endif
				}	/* if(MODULUS_TYPE == ...) */

		#endif	/* USE_SSE2 */

#if 0
	if(j < 2 && full_pass && iter <= 30)
	{
  #ifdef USE_SSE2
	#ifdef USE_SCALAR_CARRY
		fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n"
		,iter,s1p00r->re,s1p01r->re,s1p00r->im,s1p01r->im,cy_r00,cy_i00,cy_r01,cy_i01,maxerr);
	#else
		fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n"
		,iter,s1p00r->re,s1p01r->re,s1p00r->im,s1p01r->im,cy_r00->re,cy_r00->im,cy_r02->re,cy_r02->im,MAX(max_err->re,max_err->im));
	#endif
  #else
		fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f, maxerr = %10.8f\n"
		,iter,a1p00r,a1p00i,a1p01r,a1p01i,cy_r00,cy_i00,cy_r01,cy_i01,maxerr);
  #endif
	}
#endif
/*...The radix-28 DIF pass is here:	*/

		#ifdef USE_SSE2

		  #ifdef USE_SCALAR_CARRY
			/* Recover SSE2-packed data from scalar-double carry outputs: */
			if((j&3) == 0)
			{
				s1p00r->re = a1p00r;	s1p00i->re = a1p00i;
				s1p01r->re = a1p01r;	s1p01i->re = a1p01i;
				s1p02r->re = a1p02r;	s1p02i->re = a1p02i;
				s1p03r->re = a1p03r;	s1p03i->re = a1p03i;
				s1p04r->re = a1p04r;	s1p04i->re = a1p04i;
				s1p05r->re = a1p05r;	s1p05i->re = a1p05i;
				s1p06r->re = a1p06r;	s1p06i->re = a1p06i;
				s1p07r->re = a1p07r;	s1p07i->re = a1p07i;
				s1p08r->re = a1p08r;	s1p08i->re = a1p08i;
				s1p09r->re = a1p09r;	s1p09i->re = a1p09i;
				s1p10r->re = a1p10r;	s1p10i->re = a1p10i;
				s1p11r->re = a1p11r;	s1p11i->re = a1p11i;
				s1p12r->re = a1p12r;	s1p12i->re = a1p12i;
				s1p13r->re = a1p13r;	s1p13i->re = a1p13i;
				s1p14r->re = a1p14r;	s1p14i->re = a1p14i;
				s1p15r->re = a1p15r;	s1p15i->re = a1p15i;
				s1p16r->re = a1p16r;	s1p16i->re = a1p16i;
				s1p17r->re = a1p17r;	s1p17i->re = a1p17i;
				s1p18r->re = a1p18r;	s1p18i->re = a1p18i;
				s1p19r->re = a1p19r;	s1p19i->re = a1p19i;
				s1p20r->re = a1p20r;	s1p20i->re = a1p20i;
				s1p21r->re = a1p21r;	s1p21i->re = a1p21i;
				s1p22r->re = a1p22r;	s1p22i->re = a1p22i;
				s1p23r->re = a1p23r;	s1p23i->re = a1p23i;
				s1p24r->re = a1p24r;	s1p24i->re = a1p24i;
				s1p25r->re = a1p25r;	s1p25i->re = a1p25i;
				s1p26r->re = a1p26r;	s1p26i->re = a1p26i;
				s1p27r->re = a1p27r;	s1p27i->re = a1p27i;

			#if 0//DEBUG_SSE2
				fprintf(stderr, "J = %d:\n",j);
				fprintf(stderr, "radix28_carry_out: a00= %20.5f, %20.5f\n",a1p00r,a1p00i);
				fprintf(stderr, "radix28_carry_out: a01= %20.5f, %20.5f\n",a1p01r,a1p01i);
				fprintf(stderr, "radix28_carry_out: a02= %20.5f, %20.5f\n",a1p02r,a1p02i);
				fprintf(stderr, "radix28_carry_out: a03= %20.5f, %20.5f\n",a1p03r,a1p03i);
				fprintf(stderr, "radix28_carry_out: a04= %20.5f, %20.5f\n",a1p04r,a1p04i);
				fprintf(stderr, "radix28_carry_out: a05= %20.5f, %20.5f\n",a1p05r,a1p05i);
				fprintf(stderr, "radix28_carry_out: a06= %20.5f, %20.5f\n",a1p06r,a1p06i);
				fprintf(stderr, "radix28_carry_out: a07= %20.5f, %20.5f\n",a1p07r,a1p07i);
				fprintf(stderr, "radix28_carry_out: a08= %20.5f, %20.5f\n",a1p08r,a1p08i);
				fprintf(stderr, "radix28_carry_out: a09= %20.5f, %20.5f\n",a1p09r,a1p09i);
				fprintf(stderr, "radix28_carry_out: a10= %20.5f, %20.5f\n",a1p10r,a1p10i);
				fprintf(stderr, "radix28_carry_out: a11= %20.5f, %20.5f\n",a1p11r,a1p11i);
				fprintf(stderr, "radix28_carry_out: a12= %20.5f, %20.5f\n",a1p12r,a1p12i);
				fprintf(stderr, "radix28_carry_out: a13= %20.5f, %20.5f\n",a1p13r,a1p13i);
				fprintf(stderr, "radix28_carry_out: a14= %20.5f, %20.5f\n",a1p14r,a1p14i);
				fprintf(stderr, "radix28_carry_out: a15= %20.5f, %20.5f\n",a1p15r,a1p15i);
				fprintf(stderr, "radix28_carry_out: a16= %20.5f, %20.5f\n",a1p16r,a1p16i);
				fprintf(stderr, "radix28_carry_out: a17= %20.5f, %20.5f\n",a1p17r,a1p17i);
				fprintf(stderr, "radix28_carry_out: a18= %20.5f, %20.5f\n",a1p18r,a1p18i);
				fprintf(stderr, "radix28_carry_out: a19= %20.5f, %20.5f\n",a1p19r,a1p19i);
				fprintf(stderr, "radix28_carry_out: a20= %20.5f, %20.5f\n",a1p20r,a1p20i);
				fprintf(stderr, "radix28_carry_out: a21= %20.5f, %20.5f\n",a1p21r,a1p21i);
				fprintf(stderr, "radix28_carry_out: a22= %20.5f, %20.5f\n",a1p22r,a1p22i);
				fprintf(stderr, "radix28_carry_out: a23= %20.5f, %20.5f\n",a1p23r,a1p23i);
				fprintf(stderr, "radix28_carry_out: a24= %20.5f, %20.5f\n",a1p24r,a1p24i);
				fprintf(stderr, "radix28_carry_out: a25= %20.5f, %20.5f\n",a1p25r,a1p25i);
				fprintf(stderr, "radix28_carry_out: a26= %20.5f, %20.5f\n",a1p26r,a1p26i);
				fprintf(stderr, "radix28_carry_out: a27= %20.5f, %20.5f\n",a1p27r,a1p27i);
			#endif
				j += 2;
				goto SSE_LOOP;	/* Go back and do 2nd set of carries */
			}
			else
			{
				s1p00r->im = a1p00r;	s1p00i->im = a1p00i;
				s1p01r->im = a1p01r;	s1p01i->im = a1p01i;
				s1p02r->im = a1p02r;	s1p02i->im = a1p02i;
				s1p03r->im = a1p03r;	s1p03i->im = a1p03i;
				s1p04r->im = a1p04r;	s1p04i->im = a1p04i;
				s1p05r->im = a1p05r;	s1p05i->im = a1p05i;
				s1p06r->im = a1p06r;	s1p06i->im = a1p06i;
				s1p07r->im = a1p07r;	s1p07i->im = a1p07i;
				s1p08r->im = a1p08r;	s1p08i->im = a1p08i;
				s1p09r->im = a1p09r;	s1p09i->im = a1p09i;
				s1p10r->im = a1p10r;	s1p10i->im = a1p10i;
				s1p11r->im = a1p11r;	s1p11i->im = a1p11i;
				s1p12r->im = a1p12r;	s1p12i->im = a1p12i;
				s1p13r->im = a1p13r;	s1p13i->im = a1p13i;
				s1p14r->im = a1p14r;	s1p14i->im = a1p14i;
				s1p15r->im = a1p15r;	s1p15i->im = a1p15i;
				s1p16r->im = a1p16r;	s1p16i->im = a1p16i;
				s1p17r->im = a1p17r;	s1p17i->im = a1p17i;
				s1p18r->im = a1p18r;	s1p18i->im = a1p18i;
				s1p19r->im = a1p19r;	s1p19i->im = a1p19i;
				s1p20r->im = a1p20r;	s1p20i->im = a1p20i;
				s1p21r->im = a1p21r;	s1p21i->im = a1p21i;
				s1p22r->im = a1p22r;	s1p22i->im = a1p22i;
				s1p23r->im = a1p23r;	s1p23i->im = a1p23i;
				s1p24r->im = a1p24r;	s1p24i->im = a1p24i;
				s1p25r->im = a1p25r;	s1p25i->im = a1p25i;
				s1p26r->im = a1p26r;	s1p26i->im = a1p26i;
				s1p27r->im = a1p27r;	s1p27i->im = a1p27i;

			#if 0//DEBUG_SSE2
				fprintf(stderr, "J = %d:\n",j);
				fprintf(stderr, "radix28_carry_out: a00= %20.5f, %20.5f\n",a1p00r,a1p00i);
				fprintf(stderr, "radix28_carry_out: a01= %20.5f, %20.5f\n",a1p01r,a1p01i);
				fprintf(stderr, "radix28_carry_out: a02= %20.5f, %20.5f\n",a1p02r,a1p02i);
				fprintf(stderr, "radix28_carry_out: a03= %20.5f, %20.5f\n",a1p03r,a1p03i);
				fprintf(stderr, "radix28_carry_out: a04= %20.5f, %20.5f\n",a1p04r,a1p04i);
				fprintf(stderr, "radix28_carry_out: a05= %20.5f, %20.5f\n",a1p05r,a1p05i);
				fprintf(stderr, "radix28_carry_out: a06= %20.5f, %20.5f\n",a1p06r,a1p06i);
				fprintf(stderr, "radix28_carry_out: a07= %20.5f, %20.5f\n",a1p07r,a1p07i);
				fprintf(stderr, "radix28_carry_out: a08= %20.5f, %20.5f\n",a1p08r,a1p08i);
				fprintf(stderr, "radix28_carry_out: a09= %20.5f, %20.5f\n",a1p09r,a1p09i);
				fprintf(stderr, "radix28_carry_out: a10= %20.5f, %20.5f\n",a1p10r,a1p10i);
				fprintf(stderr, "radix28_carry_out: a11= %20.5f, %20.5f\n",a1p11r,a1p11i);
				fprintf(stderr, "radix28_carry_out: a12= %20.5f, %20.5f\n",a1p12r,a1p12i);
				fprintf(stderr, "radix28_carry_out: a13= %20.5f, %20.5f\n",a1p13r,a1p13i);
				fprintf(stderr, "radix28_carry_out: a14= %20.5f, %20.5f\n",a1p14r,a1p14i);
				fprintf(stderr, "radix28_carry_out: a15= %20.5f, %20.5f\n",a1p15r,a1p15i);
				fprintf(stderr, "radix28_carry_out: a16= %20.5f, %20.5f\n",a1p16r,a1p16i);
				fprintf(stderr, "radix28_carry_out: a17= %20.5f, %20.5f\n",a1p17r,a1p17i);
				fprintf(stderr, "radix28_carry_out: a18= %20.5f, %20.5f\n",a1p18r,a1p18i);
				fprintf(stderr, "radix28_carry_out: a19= %20.5f, %20.5f\n",a1p19r,a1p19i);
				fprintf(stderr, "radix28_carry_out: a20= %20.5f, %20.5f\n",a1p20r,a1p20i);
				fprintf(stderr, "radix28_carry_out: a21= %20.5f, %20.5f\n",a1p21r,a1p21i);
				fprintf(stderr, "radix28_carry_out: a22= %20.5f, %20.5f\n",a1p22r,a1p22i);
				fprintf(stderr, "radix28_carry_out: a23= %20.5f, %20.5f\n",a1p23r,a1p23i);
				fprintf(stderr, "radix28_carry_out: a24= %20.5f, %20.5f\n",a1p24r,a1p24i);
				fprintf(stderr, "radix28_carry_out: a25= %20.5f, %20.5f\n",a1p25r,a1p25i);
				fprintf(stderr, "radix28_carry_out: a26= %20.5f, %20.5f\n",a1p26r,a1p26i);
				fprintf(stderr, "radix28_carry_out: a27= %20.5f, %20.5f\n",a1p27r,a1p27i);
				exit(0);
			#endif
				j -= 2;
			}
		  #endif	/* USE_SCALAR_CARRY */

#ifdef CTIME
	clock3 = clock();
	dt_cy += (double)(clock3 - clock2);
	clock2 = clock3;
#endif
		  #if defined(COMPILER_TYPE_MSVC)

			SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r);
			SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r);
			SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r);
			SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r);

			/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

			add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add0,add1,add2,add3)
			add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add0,add1,add2,add3)
			add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add0,add1,add2,add3)
			add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add0,add1,add2,add3)
			add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add0,add1,add2,add3)
			add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add0,add1,add2,add3)
			add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add0,add1,add2,add3)

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX28_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

		  #endif

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "radix28_carry_out: a00= %20.5f, %20.5f\n",s1p00r->re,s1p00i->re);
			fprintf(stderr, "radix28_carry_out: a01= %20.5f, %20.5f\n",s1p01r->re,s1p01i->re);
			fprintf(stderr, "radix28_carry_out: a02= %20.5f, %20.5f\n",s1p02r->re,s1p02i->re);
			fprintf(stderr, "radix28_carry_out: a03= %20.5f, %20.5f\n",s1p03r->re,s1p03i->re);
			fprintf(stderr, "radix28_carry_out: a04= %20.5f, %20.5f\n",s1p04r->re,s1p04i->re);
			fprintf(stderr, "radix28_carry_out: a05= %20.5f, %20.5f\n",s1p05r->re,s1p05i->re);
			fprintf(stderr, "radix28_carry_out: a06= %20.5f, %20.5f\n",s1p06r->re,s1p06i->re);
			fprintf(stderr, "radix28_carry_out: a07= %20.5f, %20.5f\n",s1p07r->re,s1p07i->re);
			fprintf(stderr, "radix28_carry_out: a08= %20.5f, %20.5f\n",s1p08r->re,s1p08i->re);
			fprintf(stderr, "radix28_carry_out: a09= %20.5f, %20.5f\n",s1p09r->re,s1p09i->re);
			fprintf(stderr, "radix28_carry_out: a10= %20.5f, %20.5f\n",s1p10r->re,s1p10i->re);
			fprintf(stderr, "radix28_carry_out: a11= %20.5f, %20.5f\n",s1p11r->re,s1p11i->re);
			fprintf(stderr, "radix28_carry_out: a12= %20.5f, %20.5f\n",s1p12r->re,s1p12i->re);
			fprintf(stderr, "radix28_carry_out: a13= %20.5f, %20.5f\n",s1p13r->re,s1p13i->re);
			fprintf(stderr, "radix28_carry_out: a14= %20.5f, %20.5f\n",s1p14r->re,s1p14i->re);
			fprintf(stderr, "radix28_carry_out: a15= %20.5f, %20.5f\n",s1p15r->re,s1p15i->re);
			fprintf(stderr, "radix28_carry_out: a16= %20.5f, %20.5f\n",s1p16r->re,s1p16i->re);
			fprintf(stderr, "radix28_carry_out: a17= %20.5f, %20.5f\n",s1p17r->re,s1p17i->re);
			fprintf(stderr, "radix28_carry_out: a18= %20.5f, %20.5f\n",s1p18r->re,s1p18i->re);
			fprintf(stderr, "radix28_carry_out: a19= %20.5f, %20.5f\n",s1p19r->re,s1p19i->re);
			fprintf(stderr, "radix28_carry_out: a20= %20.5f, %20.5f\n",s1p20r->re,s1p20i->re);
			fprintf(stderr, "radix28_carry_out: a21= %20.5f, %20.5f\n",s1p21r->re,s1p21i->re);
			fprintf(stderr, "radix28_carry_out: a22= %20.5f, %20.5f\n",s1p22r->re,s1p22i->re);
			fprintf(stderr, "radix28_carry_out: a23= %20.5f, %20.5f\n",s1p23r->re,s1p23i->re);
			fprintf(stderr, "radix28_carry_out: a24= %20.5f, %20.5f\n",s1p24r->re,s1p24i->re);
			fprintf(stderr, "radix28_carry_out: a25= %20.5f, %20.5f\n",s1p25r->re,s1p25i->re);
			fprintf(stderr, "radix28_carry_out: a26= %20.5f, %20.5f\n",s1p26r->re,s1p26i->re);
			fprintf(stderr, "radix28_carry_out: a27= %20.5f, %20.5f\n",s1p27r->re,s1p27i->re);
		  #endif

		  #ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p04;	jp = j2+p04;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[04] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[05] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[06] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[07] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p08;	jp = j2+p08;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p12;	jp = j2+p12;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[12] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[13] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[14] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[15] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p16;	jp = j2+p16;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p20;	jp = j2+p20;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[20] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[21] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[22] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[23] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p24;	jp = j2+p24;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[24] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[25] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[26] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[27] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			exit(0);
		  #endif

		#else	/* !USE_SSE2 */

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "J = %d:\n",j);
			fprintf(stderr, "radix28_carry_out: a00= %20.5f, %20.5f\n",a1p00r,a1p00i);
			fprintf(stderr, "radix28_carry_out: a01= %20.5f, %20.5f\n",a1p01r,a1p01i);
			fprintf(stderr, "radix28_carry_out: a02= %20.5f, %20.5f\n",a1p02r,a1p02i);
			fprintf(stderr, "radix28_carry_out: a03= %20.5f, %20.5f\n",a1p03r,a1p03i);
			fprintf(stderr, "radix28_carry_out: a04= %20.5f, %20.5f\n",a1p04r,a1p04i);
			fprintf(stderr, "radix28_carry_out: a05= %20.5f, %20.5f\n",a1p05r,a1p05i);
			fprintf(stderr, "radix28_carry_out: a06= %20.5f, %20.5f\n",a1p06r,a1p06i);
			fprintf(stderr, "radix28_carry_out: a07= %20.5f, %20.5f\n",a1p07r,a1p07i);
			fprintf(stderr, "radix28_carry_out: a08= %20.5f, %20.5f\n",a1p08r,a1p08i);
			fprintf(stderr, "radix28_carry_out: a09= %20.5f, %20.5f\n",a1p09r,a1p09i);
			fprintf(stderr, "radix28_carry_out: a10= %20.5f, %20.5f\n",a1p10r,a1p10i);
			fprintf(stderr, "radix28_carry_out: a11= %20.5f, %20.5f\n",a1p11r,a1p11i);
			fprintf(stderr, "radix28_carry_out: a12= %20.5f, %20.5f\n",a1p12r,a1p12i);
			fprintf(stderr, "radix28_carry_out: a13= %20.5f, %20.5f\n",a1p13r,a1p13i);
			fprintf(stderr, "radix28_carry_out: a14= %20.5f, %20.5f\n",a1p14r,a1p14i);
			fprintf(stderr, "radix28_carry_out: a15= %20.5f, %20.5f\n",a1p15r,a1p15i);
			fprintf(stderr, "radix28_carry_out: a16= %20.5f, %20.5f\n",a1p16r,a1p16i);
			fprintf(stderr, "radix28_carry_out: a17= %20.5f, %20.5f\n",a1p17r,a1p17i);
			fprintf(stderr, "radix28_carry_out: a18= %20.5f, %20.5f\n",a1p18r,a1p18i);
			fprintf(stderr, "radix28_carry_out: a19= %20.5f, %20.5f\n",a1p19r,a1p19i);
			fprintf(stderr, "radix28_carry_out: a20= %20.5f, %20.5f\n",a1p20r,a1p20i);
			fprintf(stderr, "radix28_carry_out: a21= %20.5f, %20.5f\n",a1p21r,a1p21i);
			fprintf(stderr, "radix28_carry_out: a22= %20.5f, %20.5f\n",a1p22r,a1p22i);
			fprintf(stderr, "radix28_carry_out: a23= %20.5f, %20.5f\n",a1p23r,a1p23i);
			fprintf(stderr, "radix28_carry_out: a24= %20.5f, %20.5f\n",a1p24r,a1p24i);
			fprintf(stderr, "radix28_carry_out: a25= %20.5f, %20.5f\n",a1p25r,a1p25i);
			fprintf(stderr, "radix28_carry_out: a26= %20.5f, %20.5f\n",a1p26r,a1p26i);
			fprintf(stderr, "radix28_carry_out: a27= %20.5f, %20.5f\n",a1p27r,a1p27i);
			if(j==0) {
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					jstart += nwt;
					jhi    += nwt;

					col += RADIX;
					co3 -= RADIX;
				}
				continue;
			} else if(j==2) {
				exit(0);
			}
		  #endif

		  #if PFETCH
			addr = &a[j1];
			prefetch_p_doubles(addr);
		  #endif

		/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/
							 /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
		  #if PFETCH
			RADIX_07_DFT_PFETCH(a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p01,p02,p03);	jt=p04+p01;	jp=p04+p02;
			RADIX_07_DFT_PFETCH(a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p04, jt, jp);	jt=p08-p01;	jp=p08+p01;
			RADIX_07_DFT_PFETCH(a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp, jt,p08, jp);	jt=p08+p02;	jp=p08+p03;
			RADIX_07_DFT_PFETCH(a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp, jt, jp,p12);
		  #else
			RADIX_07_DFT       (a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		  #endif

		/*...and now do 7 radix-4 transforms...*/
							 /*                          inputs                           */ /*                                      outputs                                      */
		  #if PFETCH
			addp = addr+p12+p01;
			prefetch_p_doubles(addp);

			RADIX_04_DIF_PFETCH(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it,addr,addp,p12+p02,p12+p03);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIF_PFETCH(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p16    ,p16+p01);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIF_PFETCH(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p16+p02,p16+p03);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIF_PFETCH(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it,addr,addp,p20    ,p20+p01);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIF_PFETCH(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it,addr,addp,p20+p02,p20+p03);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIF_PFETCH(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p24    ,p20+p01);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIF_PFETCH(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p20+p02,p20+p03);
		  #else
			RADIX_04_DIF       (a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIF       (a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIF       (a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIF       (a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIF       (a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIF       (a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIF       (a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		  #endif

		#endif	/* !USE_SSE2 */

#ifdef CTIME
	clock3 = clock();
	dt_inv += (double)(clock3 - clock2);
	clock2 = clock3;
#endif
			}

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#if defined(USE_SSE2) && !defined(USE_SCALAR_CARRY)
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r00[ithread] = cy_r00->re;
			_cy_r01[ithread] = cy_r00->im;
			_cy_r02[ithread] = cy_r02->re;
			_cy_r03[ithread] = cy_r02->im;
			_cy_r04[ithread] = cy_r04->re;
			_cy_r05[ithread] = cy_r04->im;
			_cy_r06[ithread] = cy_r06->re;
			_cy_r07[ithread] = cy_r06->im;
			_cy_r08[ithread] = cy_r08->re;
			_cy_r09[ithread] = cy_r08->im;
			_cy_r10[ithread] = cy_r10->re;
			_cy_r11[ithread] = cy_r10->im;
			_cy_r12[ithread] = cy_r12->re;
			_cy_r13[ithread] = cy_r12->im;
			_cy_r14[ithread] = cy_r14->re;
			_cy_r15[ithread] = cy_r14->im;
			_cy_r16[ithread] = cy_r16->re;
			_cy_r17[ithread] = cy_r16->im;
			_cy_r18[ithread] = cy_r18->re;
			_cy_r19[ithread] = cy_r18->im;
			_cy_r20[ithread] = cy_r20->re;
			_cy_r21[ithread] = cy_r20->im;
			_cy_r22[ithread] = cy_r22->re;
			_cy_r23[ithread] = cy_r22->im;
			_cy_r24[ithread] = cy_r24->re;
			_cy_r25[ithread] = cy_r24->im;
			_cy_r26[ithread] = cy_r26->re;
			_cy_r27[ithread] = cy_r26->im;
		}
		else
		{
			_cy_r00[ithread] = cy_r00->re;	_cy_i00[ithread] = cy_r00->im;
			_cy_r01[ithread] = cy_r02->re;	_cy_i01[ithread] = cy_r02->im;
			_cy_r02[ithread] = cy_r04->re;	_cy_i02[ithread] = cy_r04->im;
			_cy_r03[ithread] = cy_r06->re;	_cy_i03[ithread] = cy_r06->im;
			_cy_r04[ithread] = cy_r08->re;	_cy_i04[ithread] = cy_r08->im;
			_cy_r05[ithread] = cy_r10->re;	_cy_i05[ithread] = cy_r10->im;
			_cy_r06[ithread] = cy_r12->re;	_cy_i06[ithread] = cy_r12->im;
			_cy_r07[ithread] = cy_r14->re;	_cy_i07[ithread] = cy_r14->im;
			_cy_r08[ithread] = cy_r16->re;	_cy_i08[ithread] = cy_r16->im;
			_cy_r09[ithread] = cy_r18->re;	_cy_i09[ithread] = cy_r18->im;
			_cy_r10[ithread] = cy_r20->re;	_cy_i10[ithread] = cy_r20->im;
			_cy_r11[ithread] = cy_r22->re;	_cy_i11[ithread] = cy_r22->im;
			_cy_r12[ithread] = cy_r24->re;	_cy_i12[ithread] = cy_r24->im;
			_cy_r13[ithread] = cy_r26->re;	_cy_i13[ithread] = cy_r26->im;
			_cy_r14[ithread] = cy_i00->re;	_cy_i14[ithread] = cy_i00->im;
			_cy_r15[ithread] = cy_i02->re;	_cy_i15[ithread] = cy_i02->im;
			_cy_r16[ithread] = cy_i04->re;	_cy_i16[ithread] = cy_i04->im;
			_cy_r17[ithread] = cy_i06->re;	_cy_i17[ithread] = cy_i06->im;
			_cy_r18[ithread] = cy_i08->re;	_cy_i18[ithread] = cy_i08->im;
			_cy_r19[ithread] = cy_i10->re;	_cy_i19[ithread] = cy_i10->im;
			_cy_r20[ithread] = cy_i12->re;	_cy_i20[ithread] = cy_i12->im;
			_cy_r21[ithread] = cy_i14->re;	_cy_i21[ithread] = cy_i14->im;
			_cy_r22[ithread] = cy_i16->re;	_cy_i22[ithread] = cy_i16->im;
			_cy_r23[ithread] = cy_i18->re;	_cy_i23[ithread] = cy_i18->im;
			_cy_r24[ithread] = cy_i20->re;	_cy_i24[ithread] = cy_i20->im;
			_cy_r25[ithread] = cy_i22->re;	_cy_i25[ithread] = cy_i22->im;
			_cy_r26[ithread] = cy_i24->re;	_cy_i26[ithread] = cy_i24->im;
			_cy_r27[ithread] = cy_i26->re;	_cy_i27[ithread] = cy_i26->im;
		}
		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r00[ithread] = cy_r00;
			_cy_r01[ithread] = cy_r01;
			_cy_r02[ithread] = cy_r02;
			_cy_r03[ithread] = cy_r03;
			_cy_r04[ithread] = cy_r04;
			_cy_r05[ithread] = cy_r05;
			_cy_r06[ithread] = cy_r06;
			_cy_r07[ithread] = cy_r07;
			_cy_r08[ithread] = cy_r08;
			_cy_r09[ithread] = cy_r09;
			_cy_r10[ithread] = cy_r10;
			_cy_r11[ithread] = cy_r11;
			_cy_r12[ithread] = cy_r12;
			_cy_r13[ithread] = cy_r13;
			_cy_r14[ithread] = cy_r14;
			_cy_r15[ithread] = cy_r15;
			_cy_r16[ithread] = cy_r16;
			_cy_r17[ithread] = cy_r17;
			_cy_r18[ithread] = cy_r18;
			_cy_r19[ithread] = cy_r19;
			_cy_r20[ithread] = cy_r20;
			_cy_r21[ithread] = cy_r21;
			_cy_r22[ithread] = cy_r22;
			_cy_r23[ithread] = cy_r23;
			_cy_r24[ithread] = cy_r24;
			_cy_r25[ithread] = cy_r25;
			_cy_r26[ithread] = cy_r26;
			_cy_r27[ithread] = cy_r27;
		}
		else
		{
			_cy_r00[ithread] = cy_r00;	_cy_i00[ithread] = cy_i00;
			_cy_r01[ithread] = cy_r01;	_cy_i01[ithread] = cy_i01;
			_cy_r02[ithread] = cy_r02;	_cy_i02[ithread] = cy_i02;
			_cy_r03[ithread] = cy_r03;	_cy_i03[ithread] = cy_i03;
			_cy_r04[ithread] = cy_r04;	_cy_i04[ithread] = cy_i04;
			_cy_r05[ithread] = cy_r05;	_cy_i05[ithread] = cy_i05;
			_cy_r06[ithread] = cy_r06;	_cy_i06[ithread] = cy_i06;
			_cy_r07[ithread] = cy_r07;	_cy_i07[ithread] = cy_i07;
			_cy_r08[ithread] = cy_r08;	_cy_i08[ithread] = cy_i08;
			_cy_r09[ithread] = cy_r09;	_cy_i09[ithread] = cy_i09;
			_cy_r10[ithread] = cy_r10;	_cy_i10[ithread] = cy_i10;
			_cy_r11[ithread] = cy_r11;	_cy_i11[ithread] = cy_i11;
			_cy_r12[ithread] = cy_r12;	_cy_i12[ithread] = cy_i12;
			_cy_r13[ithread] = cy_r13;	_cy_i13[ithread] = cy_i13;
			_cy_r14[ithread] = cy_r14;	_cy_i14[ithread] = cy_i14;
			_cy_r15[ithread] = cy_r15;	_cy_i15[ithread] = cy_i15;
			_cy_r16[ithread] = cy_r16;	_cy_i16[ithread] = cy_i16;
			_cy_r17[ithread] = cy_r17;	_cy_i17[ithread] = cy_i17;
			_cy_r18[ithread] = cy_r18;	_cy_i18[ithread] = cy_i18;
			_cy_r19[ithread] = cy_r19;	_cy_i19[ithread] = cy_i19;
			_cy_r20[ithread] = cy_r20;	_cy_i20[ithread] = cy_i20;
			_cy_r21[ithread] = cy_r21;	_cy_i21[ithread] = cy_i21;
			_cy_r22[ithread] = cy_r22;	_cy_i22[ithread] = cy_i22;
			_cy_r23[ithread] = cy_r23;	_cy_i23[ithread] = cy_i23;
			_cy_r24[ithread] = cy_r24;	_cy_i24[ithread] = cy_i24;
			_cy_r25[ithread] = cy_r25;	_cy_i25[ithread] = cy_i25;
			_cy_r26[ithread] = cy_r26;	_cy_i26[ithread] = cy_i26;
			_cy_r27[ithread] = cy_r27;	_cy_i27[ithread] = cy_i27;
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #ifdef OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy28_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that
	
	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix32_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;
		}
		else
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;	_cy_i00[ithread] = tdat[ithread].cy_i00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;	_cy_i01[ithread] = tdat[ithread].cy_i01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;	_cy_i02[ithread] = tdat[ithread].cy_i02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;	_cy_i03[ithread] = tdat[ithread].cy_i03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;	_cy_i04[ithread] = tdat[ithread].cy_i04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;	_cy_i05[ithread] = tdat[ithread].cy_i05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;	_cy_i06[ithread] = tdat[ithread].cy_i06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;	_cy_i07[ithread] = tdat[ithread].cy_i07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;	_cy_i08[ithread] = tdat[ithread].cy_i08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;	_cy_i09[ithread] = tdat[ithread].cy_i09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;	_cy_i10[ithread] = tdat[ithread].cy_i10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;	_cy_i11[ithread] = tdat[ithread].cy_i11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;	_cy_i12[ithread] = tdat[ithread].cy_i12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;	_cy_i13[ithread] = tdat[ithread].cy_i13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;	_cy_i14[ithread] = tdat[ithread].cy_i14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;	_cy_i15[ithread] = tdat[ithread].cy_i15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;	_cy_i16[ithread] = tdat[ithread].cy_i16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;	_cy_i17[ithread] = tdat[ithread].cy_i17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;	_cy_i18[ithread] = tdat[ithread].cy_i18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;	_cy_i19[ithread] = tdat[ithread].cy_i19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;	_cy_i20[ithread] = tdat[ithread].cy_i20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;	_cy_i21[ithread] = tdat[ithread].cy_i21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;	_cy_i22[ithread] = tdat[ithread].cy_i22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;	_cy_i23[ithread] = tdat[ithread].cy_i23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;	_cy_i24[ithread] = tdat[ithread].cy_i24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;	_cy_i25[ithread] = tdat[ithread].cy_i25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;	_cy_i26[ithread] = tdat[ithread].cy_i26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;	_cy_i27[ithread] = tdat[ithread].cy_i27;
		}
	}
#endif

#if FFT_DEBUG
	fprintf(dbg_file,"Iter = %d:\n\n",iter);
	for(j = 0; j < n; j++) {
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		fprintf(dbg_file,"a[%d] = %20.10e %20.10e\n", j,a[j1],a[j1+1]);
	}
	if(iter > 0 && !full_pass) {
		fclose(dbg_file);
		dbg_file = 0x0;
		sprintf(cbuf, "Wrote debug file %s", dbg_fname);
		fprintf(stderr, "%s\n", cbuf);
		exit(0);
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
	!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00= _cy_r00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];
		}

		_cy_r00[0] =+t54;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00;
		_cy_r02[0] = t02;
		_cy_r03[0] = t04;
		_cy_r04[0] = t06;
		_cy_r05[0] = t08;
		_cy_r06[0] = t10;
		_cy_r07[0] = t12;
		_cy_r08[0] = t14;
		_cy_r09[0] = t16;
		_cy_r10[0] = t18;
		_cy_r11[0] = t20;
		_cy_r12[0] = t22;
		_cy_r13[0] = t24;
		_cy_r14[0] = t26;
		_cy_r15[0] = t28;
		_cy_r16[0] = t30;
		_cy_r17[0] = t32;
		_cy_r18[0] = t34;
		_cy_r19[0] = t36;
		_cy_r20[0] = t38;
		_cy_r21[0] = t40;
		_cy_r22[0] = t42;
		_cy_r23[0] = t44;
		_cy_r24[0] = t46;
		_cy_r25[0] = t48;
		_cy_r26[0] = t50;
		_cy_r27[0] = t52;
	}
	else
	{
		t00= _cy_r00[CY_THREADS - 1];	t01= _cy_i00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t03= _cy_i01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t05= _cy_i02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t07= _cy_i03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t09= _cy_i04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];	t11= _cy_i05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];	t13= _cy_i06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];	t15= _cy_i07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];	t17= _cy_i08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];	t19= _cy_i09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];	t21= _cy_i10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];	t23= _cy_i11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];	t25= _cy_i12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];	t27= _cy_i13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];	t29= _cy_i14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];	t31= _cy_i15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];	t33= _cy_i16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];	t35= _cy_i17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];	t37= _cy_i18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];	t39= _cy_i19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];	t41= _cy_i20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];	t43= _cy_i21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];	t45= _cy_i22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];	t47= _cy_i23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];	t49= _cy_i24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];	t51= _cy_i25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];	t53= _cy_i26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];	t55= _cy_i27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];		_cy_i00[ithread] = _cy_i00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];		_cy_i01[ithread] = _cy_i01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];		_cy_i02[ithread] = _cy_i02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];		_cy_i03[ithread] = _cy_i03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];		_cy_i04[ithread] = _cy_i04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];		_cy_i05[ithread] = _cy_i05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];		_cy_i06[ithread] = _cy_i06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];		_cy_i07[ithread] = _cy_i07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];		_cy_i08[ithread] = _cy_i08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];		_cy_i09[ithread] = _cy_i09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];		_cy_i10[ithread] = _cy_i10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];		_cy_i11[ithread] = _cy_i11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];		_cy_i12[ithread] = _cy_i12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];		_cy_i13[ithread] = _cy_i13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];		_cy_i14[ithread] = _cy_i14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];		_cy_i15[ithread] = _cy_i15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];		_cy_i16[ithread] = _cy_i16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];		_cy_i17[ithread] = _cy_i17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];		_cy_i18[ithread] = _cy_i18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];		_cy_i19[ithread] = _cy_i19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];		_cy_i20[ithread] = _cy_i20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];		_cy_i21[ithread] = _cy_i21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];		_cy_i22[ithread] = _cy_i22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];		_cy_i23[ithread] = _cy_i23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];		_cy_i24[ithread] = _cy_i24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];		_cy_i25[ithread] = _cy_i25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];		_cy_i26[ithread] = _cy_i26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];		_cy_i27[ithread] = _cy_i27[ithread-1];
		}

		_cy_r00[0] =-t55;	_cy_i00[0] =+t54;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t10;	_cy_i06[0] = t11;
		_cy_r07[0] = t12;	_cy_i07[0] = t13;
		_cy_r08[0] = t14;	_cy_i08[0] = t15;
		_cy_r09[0] = t16;	_cy_i09[0] = t17;
		_cy_r10[0] = t18;	_cy_i10[0] = t19;
		_cy_r11[0] = t20;	_cy_i11[0] = t21;
		_cy_r12[0] = t22;	_cy_i12[0] = t23;
		_cy_r13[0] = t24;	_cy_i13[0] = t25;
		_cy_r14[0] = t26;	_cy_i14[0] = t27;
		_cy_r15[0] = t28;	_cy_i15[0] = t29;
		_cy_r16[0] = t30;	_cy_i16[0] = t31;
		_cy_r17[0] = t32;	_cy_i17[0] = t33;
		_cy_r18[0] = t34;	_cy_i18[0] = t35;
		_cy_r19[0] = t36;	_cy_i19[0] = t37;
		_cy_r20[0] = t38;	_cy_i20[0] = t39;
		_cy_r21[0] = t40;	_cy_i21[0] = t41;
		_cy_r22[0] = t42;	_cy_i22[0] = t43;
		_cy_r23[0] = t44;	_cy_i23[0] = t45;
		_cy_r24[0] = t46;	_cy_i24[0] = t47;
		_cy_r25[0] = t48;	_cy_i25[0] = t49;
		_cy_r26[0] = t50;	_cy_i26[0] = t51;
		_cy_r27[0] = t52;	_cy_i27[0] = t53;
	}

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			k = j;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
			k += p04;
			a[k    ] *= radix_inv;
			a[k+p01] *= radix_inv;
			a[k+p02] *= radix_inv;
			a[k+p03] *= radix_inv;
		}
	}
}	/* endfor(outer) */

#ifdef CTIME
	clock2 = clock();
	dt_tot = (double)(clock2 - clock1);
	printf("radix28_carry cycle times: total = %10.5f, fwd = %10.5f, inv = %10.5f, cy = %10.5f\n", dt_tot*ICPS, dt_fwd*ICPS, dt_inv*ICPS, dt_cy*ICPS);
#endif

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0])+fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0]);
		t00 += fabs(_cy_i00[0])+fabs(_cy_i01[0])+fabs(_cy_i02[0])+fabs(_cy_i03[0])+fabs(_cy_i04[0])+fabs(_cy_i05[0])+fabs(_cy_i06[0])+fabs(_cy_i07[0])+fabs(_cy_i08[0])+fabs(_cy_i09[0])+fabs(_cy_i10[0])+fabs(_cy_i11[0])+fabs(_cy_i12[0])+fabs(_cy_i13[0])+fabs(_cy_i14[0])+fabs(_cy_i15[0])+fabs(_cy_i16[0])+fabs(_cy_i17[0])+fabs(_cy_i18[0])+fabs(_cy_i19[0])+fabs(_cy_i20[0])+fabs(_cy_i21[0])+fabs(_cy_i22[0])+fabs(_cy_i23[0])+fabs(_cy_i24[0])+fabs(_cy_i25[0])+fabs(_cy_i26[0])+fabs(_cy_i27[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}

	if(t00 != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix28_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
			if(INTERACT)fprintf(stderr,"%s",cbuf);
			fp = fopen(   OFILE,"a");
			fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
			err=ERR_CARRY;
			return(err);
	}
	return(0);
}

/***************/

int radix28_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-28 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-28 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
/* For dual-use (Fermat / Mersenne-mod) carry routines, pack both the nwt/nrt and associated _bits params into a 32-bit int: */
	int n28, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k1,k2,k,khi,l,outer;
  #ifdef DEBUG_SSE2
	int jt,jp;
  #endif
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	/* The current versions of the macros in dft_macro.h doesn't allow a #define LO_ADD;
	it simply assumes LO_ADD = 1 (and this is explicitly checked at runtime),
	so must use the corresponding versions of the sincos constants :
	*/
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	static double radix_inv, n2inv;
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27
	,temp,scale;
#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27;	/* indices into weights arrays (mod NWT) */
	double wt_re,wt_im;									/* Fermat-mod weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii00=ii01=ii02=ii03=ii04=ii05=ii06=ii07=ii08=ii09=ii10=ii11=ii12=ii13=ii14=ii15=ii16=ii17=ii18=ii19=ii20=ii21=ii22=ii23=ii24=ii25=ii26=ii27=-1;

/*...change n28 and n_div_wt to non-static to work around a gcc compiler bug. */
	n28   = n/28;
	n_div_nwt = n28 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n28)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/28 in radix28_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT(HERE, LO_ADD,"radix28_ditN_cy_dif1.c: LO_ADD");
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)28));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < n28; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < n28/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-28 final DIT pass is here.	*/

	/* init carries	*/
	cy_r00= 0;	cy_i00= 0;
	cy_r01= 0;	cy_i01= 0;
	cy_r02= 0;	cy_i02= 0;
	cy_r03= 0;	cy_i03= 0;
	cy_r04= 0;	cy_i04= 0;
	cy_r05= 0;	cy_i05= 0;
	cy_r06= 0;	cy_i06= 0;
	cy_r07= 0;	cy_i07= 0;
	cy_r08= 0;	cy_i08= 0;
	cy_r09= 0;	cy_i09= 0;
	cy_r10= 0;	cy_i10= 0;
	cy_r11= 0;	cy_i11= 0;
	cy_r12= 0;	cy_i12= 0;
	cy_r13= 0;	cy_i13= 0;
	cy_r14= 0;	cy_i14= 0;
	cy_r15= 0;	cy_i15= 0;
	cy_r16= 0;	cy_i16= 0;
	cy_r17= 0;	cy_i17= 0;
	cy_r18= 0;	cy_i18= 0;
	cy_r19= 0;	cy_i19= 0;
	cy_r20= 0;	cy_i20= 0;
	cy_r21= 0;	cy_i21= 0;
	cy_r22= 0;	cy_i22= 0;
	cy_r23= 0;	cy_i23= 0;
	cy_r24= 0;	cy_i24= 0;
	cy_r25= 0;	cy_i25= 0;
	cy_r26= 0;	cy_i26= 0;
	cy_r27= 0;	cy_i27= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r00 = -2;
	}

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/

	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = 0;
		jhi = jstart+nwt-1;
		khi = n_div_nwt;
	}
	else
	{
		jstart = 0;
		jhi = n_div_nwt;
		khi = 1;
	}

for(outer=0; outer <= 1; outer++)
{
	i = 0;		/* Index into the BASE and BASEINV arrays. */
	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
	if(bw > 0)
		i = 1;

	bjmodn00= 0;
	bjmodn01= bjmodnini;
	bjmodn02= bjmodn01+bjmodnini-n; bjmodn02= bjmodn02+ ( (-(int)((uint32)bjmodn02>> 31)) & n);
	bjmodn03= bjmodn02+bjmodnini-n; bjmodn03= bjmodn03+ ( (-(int)((uint32)bjmodn03>> 31)) & n);
	bjmodn04= bjmodn03+bjmodnini-n; bjmodn04= bjmodn04+ ( (-(int)((uint32)bjmodn04>> 31)) & n);
	bjmodn05= bjmodn04+bjmodnini-n; bjmodn05= bjmodn05+ ( (-(int)((uint32)bjmodn05>> 31)) & n);
	bjmodn06= bjmodn05+bjmodnini-n; bjmodn06= bjmodn06+ ( (-(int)((uint32)bjmodn06>> 31)) & n);
	bjmodn07= bjmodn06+bjmodnini-n; bjmodn07= bjmodn07+ ( (-(int)((uint32)bjmodn07>> 31)) & n);
	bjmodn08= bjmodn07+bjmodnini-n; bjmodn08= bjmodn08+ ( (-(int)((uint32)bjmodn08>> 31)) & n);
	bjmodn09= bjmodn08+bjmodnini-n; bjmodn09= bjmodn09+ ( (-(int)((uint32)bjmodn09>> 31)) & n);
	bjmodn10= bjmodn09+bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+bjmodnini-n; bjmodn14= bjmodn14+ ( (-(int)((uint32)bjmodn14>> 31)) & n);
	bjmodn15= bjmodn14+bjmodnini-n; bjmodn15= bjmodn15+ ( (-(int)((uint32)bjmodn15>> 31)) & n);
	bjmodn16= bjmodn15+bjmodnini-n; bjmodn16= bjmodn16+ ( (-(int)((uint32)bjmodn16>> 31)) & n);
	bjmodn17= bjmodn16+bjmodnini-n; bjmodn17= bjmodn17+ ( (-(int)((uint32)bjmodn17>> 31)) & n);
	bjmodn18= bjmodn17+bjmodnini-n; bjmodn18= bjmodn18+ ( (-(int)((uint32)bjmodn18>> 31)) & n);
	bjmodn19= bjmodn18+bjmodnini-n; bjmodn19= bjmodn19+ ( (-(int)((uint32)bjmodn19>> 31)) & n);
	bjmodn20= bjmodn19+bjmodnini-n; bjmodn20= bjmodn20+ ( (-(int)((uint32)bjmodn20>> 31)) & n);
	bjmodn21= bjmodn20+bjmodnini-n; bjmodn21= bjmodn21+ ( (-(int)((uint32)bjmodn21>> 31)) & n);
	bjmodn22= bjmodn21+bjmodnini-n; bjmodn22= bjmodn22+ ( (-(int)((uint32)bjmodn22>> 31)) & n);
	bjmodn23= bjmodn22+bjmodnini-n; bjmodn23= bjmodn23+ ( (-(int)((uint32)bjmodn23>> 31)) & n);
	bjmodn24= bjmodn23+bjmodnini-n; bjmodn24= bjmodn24+ ( (-(int)((uint32)bjmodn24>> 31)) & n);
	bjmodn25= bjmodn24+bjmodnini-n; bjmodn25= bjmodn25+ ( (-(int)((uint32)bjmodn25>> 31)) & n);
	bjmodn26= bjmodn25+bjmodnini-n; bjmodn26= bjmodn26+ ( (-(int)((uint32)bjmodn26>> 31)) & n);
	bjmodn27= bjmodn26+bjmodnini-n; bjmodn27= bjmodn27+ ( (-(int)((uint32)bjmodn27>> 31)) & n);

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros that don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+28;
		co3=co2-28;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*n28/2) % nwt;
		ii02= (ii01+ ii01) % nwt;
		ii03= (ii02+ ii01) % nwt;
		ii04= (ii03+ ii01) % nwt;
		ii05= (ii04+ ii01) % nwt;
		ii06= (ii05+ ii01) % nwt;
		ii07= (ii06+ ii01) % nwt;
		ii08= (ii07+ ii01) % nwt;
		ii09= (ii08+ ii01) % nwt;
		ii10= (ii09+ ii01) % nwt;
		ii11= (ii10+ ii01) % nwt;
		ii12= (ii11+ ii01) % nwt;
		ii13= (ii12+ ii01) % nwt;
		ii14= (ii13+ ii01) % nwt;
		ii15= (ii14+ ii01) % nwt;
		ii16= (ii15+ ii01) % nwt;
		ii17= (ii16+ ii01) % nwt;
		ii18= (ii17+ ii01) % nwt;
		ii19= (ii18+ ii01) % nwt;
		ii20= (ii19+ ii01) % nwt;
		ii21= (ii20+ ii01) % nwt;
		ii22= (ii21+ ii01) % nwt;
		ii23= (ii22+ ii01) % nwt;
		ii24= (ii23+ ii01) % nwt;
		ii25= (ii24+ ii01) % nwt;
		ii26= (ii25+ ii01) % nwt;
		ii27= (ii26+ ii01) % nwt;

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn00= n;
		bjmodn07= n;
		bjmodn14= n;
		bjmodn21= n;
	}

	for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
	{
		for(j=jstart; j<jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
		{
		#ifdef USE_SSE2
			j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
		#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		#endif
			j2 = j1+RE_IM_STRIDE;

		#ifdef DEBUG_SSE2
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p04;	jp = j2 + p04;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p08;	jp = j2 + p08;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p12;	jp = j2 + p12;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p16;	jp = j2 + p16;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p20;	jp = j2 + p20;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p24;	jp = j2 + p24;
			a[jt    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    +1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jt+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03+1] = 1024.0*1024.0*1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			jt = j1;		jp = j2;
			fprintf(stderr, "radix28_carry: A_in[00] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[01] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[02] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[03] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p04;	jp = j2 + p04;
			fprintf(stderr, "radix28_carry: A_in[04] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[05] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[06] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[07] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p08;	jp = j2 + p08;
			fprintf(stderr, "radix28_carry: A_in[08] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[09] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[10] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[11] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p12;	jp = j2 + p12;
			fprintf(stderr, "radix28_carry: A_in[12] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[13] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[14] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[15] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p16;	jp = j2 + p16;
			fprintf(stderr, "radix28_carry: A_in[16] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[17] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[18] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[19] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p20;	jp = j2 + p20;
			fprintf(stderr, "radix28_carry: A_in[20] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[21] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[22] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[23] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);	jt = j1 + p24;	jp = j2 + p24;
			fprintf(stderr, "radix28_carry: A_in[24] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt    ],a[jt    +1],a[jp    ],a[jp    +1]);
			fprintf(stderr, "radix28_carry: A_in[25] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p01],a[jt+p01+1],a[jp+p01],a[jp+p01+1]);
			fprintf(stderr, "radix28_carry: A_in[26] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p02],a[jt+p02+1],a[jp+p02],a[jp+p02+1]);
			fprintf(stderr, "radix28_carry: A_in[27] = %20.5f, %20.5f, %20.5f, %20.5f\n",a[jt+p03],a[jt+p03+1],a[jp+p03],a[jp+p03+1]);
			fprintf(stderr, "\n");
		#endif

/*...The radix-28 DIT pass is here:	*/

/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
                  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
                  to properly re-use the ajp1 variables in the carry-pass version of this routine.
*/

/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
                     /*                                      outputs                                      */ /*                          inputs                           */
	RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);
	RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);
	RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);
	RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);
	RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);
	RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);
	RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

/*...and now do 4 radix-4 transforms...*/
                     /*                                                   inputs                                                  */ /*               intermediates              */ /*                                                                     outputs                                                                         */
	RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			l= j & (nwt-1);
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwt-l  ];
			sinwtm1 = si[nwt-l-1];

			wtl     =wt0[    l  ];
			wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
			wtlp1   =wt0[    l+1];
			wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			 cmplx_carry_norm_nocheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
			cmplx_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
			cmplx_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
			cmplx_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
			cmplx_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
			cmplx_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
			cmplx_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
			cmplx_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
			cmplx_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
			cmplx_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
			cmplx_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
			cmplx_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
			cmplx_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
			cmplx_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
			cmplx_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
			cmplx_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
			cmplx_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
			cmplx_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
			cmplx_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
			cmplx_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
			cmplx_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
			cmplx_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
			cmplx_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
			cmplx_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
			cmplx_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
			cmplx_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
			cmplx_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
			cmplx_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_nocheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,0 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,1 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,2 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,3 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,4 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,5 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,6 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,7 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,8 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,9 *n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,15*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,16*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,17*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,18*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,19*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,20*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,21*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,22*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,23*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,24*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,25*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,26*n28,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,27*n28,NRTM1,NRT_BITS);
		}

		  #ifdef DEBUG_SSE2
			fprintf(stderr, "J = %d:\n",j);
			fprintf(stderr, "radix28_carry_out: a00= %20.5f, %20.5f\n",a1p00r,a1p00i);
			fprintf(stderr, "radix28_carry_out: a01= %20.5f, %20.5f\n",a1p01r,a1p01i);
			fprintf(stderr, "radix28_carry_out: a02= %20.5f, %20.5f\n",a1p02r,a1p02i);
			fprintf(stderr, "radix28_carry_out: a03= %20.5f, %20.5f\n",a1p03r,a1p03i);
			fprintf(stderr, "radix28_carry_out: a04= %20.5f, %20.5f\n",a1p04r,a1p04i);
			fprintf(stderr, "radix28_carry_out: a05= %20.5f, %20.5f\n",a1p05r,a1p05i);
			fprintf(stderr, "radix28_carry_out: a06= %20.5f, %20.5f\n",a1p06r,a1p06i);
			fprintf(stderr, "radix28_carry_out: a07= %20.5f, %20.5f\n",a1p07r,a1p07i);
			fprintf(stderr, "radix28_carry_out: a08= %20.5f, %20.5f\n",a1p08r,a1p08i);
			fprintf(stderr, "radix28_carry_out: a09= %20.5f, %20.5f\n",a1p09r,a1p09i);
			fprintf(stderr, "radix28_carry_out: a10= %20.5f, %20.5f\n",a1p10r,a1p10i);
			fprintf(stderr, "radix28_carry_out: a11= %20.5f, %20.5f\n",a1p11r,a1p11i);
			fprintf(stderr, "radix28_carry_out: a12= %20.5f, %20.5f\n",a1p12r,a1p12i);
			fprintf(stderr, "radix28_carry_out: a13= %20.5f, %20.5f\n",a1p13r,a1p13i);
			fprintf(stderr, "radix28_carry_out: a14= %20.5f, %20.5f\n",a1p14r,a1p14i);
			fprintf(stderr, "radix28_carry_out: a15= %20.5f, %20.5f\n",a1p15r,a1p15i);
			fprintf(stderr, "radix28_carry_out: a16= %20.5f, %20.5f\n",a1p16r,a1p16i);
			fprintf(stderr, "radix28_carry_out: a17= %20.5f, %20.5f\n",a1p17r,a1p17i);
			fprintf(stderr, "radix28_carry_out: a18= %20.5f, %20.5f\n",a1p18r,a1p18i);
			fprintf(stderr, "radix28_carry_out: a19= %20.5f, %20.5f\n",a1p19r,a1p19i);
			fprintf(stderr, "radix28_carry_out: a20= %20.5f, %20.5f\n",a1p20r,a1p20i);
			fprintf(stderr, "radix28_carry_out: a21= %20.5f, %20.5f\n",a1p21r,a1p21i);
			fprintf(stderr, "radix28_carry_out: a22= %20.5f, %20.5f\n",a1p22r,a1p22i);
			fprintf(stderr, "radix28_carry_out: a23= %20.5f, %20.5f\n",a1p23r,a1p23i);
			fprintf(stderr, "radix28_carry_out: a24= %20.5f, %20.5f\n",a1p24r,a1p24i);
			fprintf(stderr, "radix28_carry_out: a25= %20.5f, %20.5f\n",a1p25r,a1p25i);
			fprintf(stderr, "radix28_carry_out: a26= %20.5f, %20.5f\n",a1p26r,a1p26i);
			fprintf(stderr, "radix28_carry_out: a27= %20.5f, %20.5f\n",a1p27r,a1p27i);
			if(j==0) {
				if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					jstart += nwt;
					jhi    += nwt;

					col += 28;
					co3 -= 28;
				}
				continue;
			} else if(j==2) {
				exit(0);
			}
		  #endif

#if 0
	if(j < 2 && (root_incr!=0) && iter <= 30)
	{
		fprintf(stderr, "Iter %3d: a0-3_in = %20.5f %20.5f %20.5f %20.5f, cy0-3 = %20.5f %20.5f %20.5f %20.5f\n"
		,iter,a1p00r,a1p00i,a1p01r,a1p01i,cy_r00,cy_i00,cy_r01,cy_i01);
	}
#endif

/*...The radix-28 DIF pass is here:	*/
#if PFETCH
addr = &a[j1];
prefetch_p_doubles(addr);
#endif

/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/
                     /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
#if PFETCH
	RADIX_07_DFT_PFETCH(a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p01,p02,p03);
	RADIX_07_DFT_PFETCH(a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p04,p05,p06);
	RADIX_07_DFT_PFETCH(a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p07,p08,p09);
	RADIX_07_DFT_PFETCH(a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p10,p11,p12);
#else
	RADIX_07_DFT       (a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
#endif

/*...and now do 7 radix-4 transforms...*/
                     /*                          inputs                           */ /*                                      outputs                                      */
#if PFETCH
	addp = addr+p13;
	prefetch_p_doubles(addp);

	RADIX_04_DIF_PFETCH(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it,addr,addp,p14,p15);
	RADIX_04_DIF_PFETCH(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it,addr,addp,p16,p17);
	RADIX_04_DIF_PFETCH(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it,addr,addp,p18,p19);
	RADIX_04_DIF_PFETCH(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it,addr,addp,p20,p21);
	RADIX_04_DIF_PFETCH(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p22,p23);
	RADIX_04_DIF_PFETCH(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it,addr,addp,p24,p25);
	RADIX_04_DIF_PFETCH(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p26,p27);
#else
	RADIX_04_DIF       (a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);
	RADIX_04_DIF       (a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it);
	RADIX_04_DIF       (a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it);
	RADIX_04_DIF       (a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it);
	RADIX_04_DIF       (a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
	RADIX_04_DIF       (a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it);
	RADIX_04_DIF       (a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
#endif

		iroot += root_incr;		/* increment sincos index.	*/

		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jstart += nwt;
			jhi    += nwt;

			col += 28;
			co3 -= 28;
		}
	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1    = cy_r27;
		cy_r27= cy_r26;
		cy_r26= cy_r25;
		cy_r25= cy_r24;
		cy_r24= cy_r23;
		cy_r23= cy_r22;
		cy_r22= cy_r21;
		cy_r21= cy_r20;
		cy_r20= cy_r19;
		cy_r19= cy_r18;
		cy_r18= cy_r17;
		cy_r17= cy_r16;
		cy_r16= cy_r15;
		cy_r15= cy_r14;
		cy_r14= cy_r13;
		cy_r13= cy_r12;
		cy_r12= cy_r11;
		cy_r11= cy_r10;
		cy_r10= cy_r09;
		cy_r09= cy_r08;
		cy_r08= cy_r07;
		cy_r07= cy_r06;
		cy_r06= cy_r05;
		cy_r05= cy_r04;
		cy_r04= cy_r03;
		cy_r03= cy_r02;
		cy_r02= cy_r01;
		cy_r01= cy_r00;
		cy_r00=    t1 ;
	}
	else
	{
		/* ...The 2 Mo"bius carries are here: */
		t1    = cy_r27;	t2    = cy_i27;
		cy_r27= cy_r26;	cy_i27= cy_i26;
		cy_r26= cy_r25;	cy_i26= cy_i25;
		cy_r25= cy_r24;	cy_i25= cy_i24;
		cy_r24= cy_r23;	cy_i24= cy_i23;
		cy_r23= cy_r22;	cy_i23= cy_i22;
		cy_r22= cy_r21;	cy_i22= cy_i21;
		cy_r21= cy_r20;	cy_i21= cy_i20;
		cy_r20= cy_r19;	cy_i20= cy_i19;
		cy_r19= cy_r18;	cy_i19= cy_i18;
		cy_r18= cy_r17;	cy_i18= cy_i17;
		cy_r17= cy_r16;	cy_i17= cy_i16;
		cy_r16= cy_r15;	cy_i16= cy_i15;
		cy_r15= cy_r14;	cy_i15= cy_i14;
		cy_r14= cy_r13;	cy_i14= cy_i13;
		cy_r13= cy_r12;	cy_i13= cy_i12;
		cy_r12= cy_r11;	cy_i12= cy_i11;
		cy_r11= cy_r10;	cy_i11= cy_i10;
		cy_r10= cy_r09;	cy_i10= cy_i09;
		cy_r09= cy_r08;	cy_i09= cy_i08;
		cy_r08= cy_r07;	cy_i08= cy_i07;
		cy_r07= cy_r06;	cy_i07= cy_i06;
		cy_r06= cy_r05;	cy_i06= cy_i05;
		cy_r05= cy_r04;	cy_i05= cy_i04;
		cy_r04= cy_r03;	cy_i04= cy_i03;
		cy_r03= cy_r02;	cy_i03= cy_i02;
		cy_r02= cy_r01;	cy_i02= cy_i01;
		cy_r01= cy_r00;	cy_i01= cy_i00;
		cy_r00=   -t2 ;	cy_i00=   +t1 ;
	}

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	}
	else
	{
		jhi = 7;
	}
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p01] *= radix_inv;
		a[j+p02] *= radix_inv;
		a[j+p03] *= radix_inv;
		a[j+p04] *= radix_inv;
		a[j+p05] *= radix_inv;
		a[j+p06] *= radix_inv;
		a[j+p07] *= radix_inv;
		a[j+p08] *= radix_inv;
		a[j+p09] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
		a[j+p15] *= radix_inv;
		a[j+p16] *= radix_inv;
		a[j+p17] *= radix_inv;
		a[j+p18] *= radix_inv;
		a[j+p19] *= radix_inv;
		a[j+p20] *= radix_inv;
		a[j+p21] *= radix_inv;
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
		a[j+p26] *= radix_inv;
		a[j+p27] *= radix_inv;
	}
}

	if(fabs(cy_r00)+fabs(cy_r01)+fabs(cy_r02)+fabs(cy_r03)+fabs(cy_r04)+fabs(cy_r05)+fabs(cy_r06)+fabs(cy_r07)+fabs(cy_r08)+fabs(cy_r09)+fabs(cy_r10)+fabs(cy_r11)+fabs(cy_r12)+fabs(cy_r13)+fabs(cy_r14)+fabs(cy_r15)+fabs(cy_r16)+fabs(cy_r17)+fabs(cy_r18)+fabs(cy_r19)+fabs(cy_r20)+fabs(cy_r21)+fabs(cy_r22)+fabs(cy_r23)+fabs(cy_r24)+fabs(cy_r25)+fabs(cy_r26)+fabs(cy_r27)
		+fabs(cy_i00)+fabs(cy_i01)+fabs(cy_i02)+fabs(cy_i03)+fabs(cy_i04)+fabs(cy_i05)+fabs(cy_i06)+fabs(cy_i07)+fabs(cy_i08)+fabs(cy_i09)+fabs(cy_i10)+fabs(cy_i11)+fabs(cy_i12)+fabs(cy_i13)+fabs(cy_i14)+fabs(cy_i15)+fabs(cy_i16)+fabs(cy_i17)+fabs(cy_i18)+fabs(cy_i19)+fabs(cy_i20)+fabs(cy_i21)+fabs(cy_i22)+fabs(cy_i23)+fabs(cy_i24)+fabs(cy_i25)+fabs(cy_i26)+fabs(cy_i27) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix28_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
			if(INTERACT)fprintf(stderr,"%s",cbuf);
			fp = fopen(   OFILE,"a");
			fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
			err=ERR_CARRY;
			return(err);
	}

	return(0);
}

/***************/

void radix28_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-28 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j=0; j < n28; j += 2)
	{
	#if defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,21,14, 7,24,17,10, 3,20,13, 6,27,16, 9, 2,23,12, 5,26,19, 8, 1,22,15, 4,25,18,11
		I.e. start out with first septet of indices {0,4,8,12,16,20,24}, permute those according to
		{0,4,8,12,16,20,24}*27%28 = {0,24,20,16,12,8,4}, then each is head of a length-4 list of indices with decrement 7.
	*/
						 /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
		RADIX_07_DFT(a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p21],a[j2+p21],a[j1+p17],a[j2+p17],a[j1+p13],a[j2+p13],a[j1+p09],a[j2+p09],a[j1+p05],a[j2+p05],a[j1+p01],a[j2+p01],a[j1+p25],a[j2+p25],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p14],a[j2+p14],a[j1+p10],a[j2+p10],a[j1+p06],a[j2+p06],a[j1+p02],a[j2+p02],a[j1+p26],a[j2+p26],a[j1+p22],a[j2+p22],a[j1+p18],a[j2+p18],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p07],a[j2+p07],a[j1+p03],a[j2+p03],a[j1+p27],a[j2+p27],a[j1+p23],a[j2+p23],a[j1+p19],a[j2+p19],a[j1+p15],a[j2+p15],a[j1+p11],a[j2+p11],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

	/*...and now do 7 radix-4 transforms...*/
						 /*                          inputs                           */ /*                                      outputs                                      */
		RADIX_04_DIF(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);
		RADIX_04_DIF(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it);
		RADIX_04_DIF(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it);
		RADIX_04_DIF(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it);
		RADIX_04_DIF(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it);
		RADIX_04_DIF(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
	}
}

/***************/

void radix28_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-28 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j=0; j < n28; j += 2)
	{
	#if defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,24,20,16,12, 8, 4,21,17,13, 9, 5, 1,25,14,10, 6, 2,26,22,18, 7, 3,27,23,19,15,11

		I.e. start out with first quartet of indices {0,7,14,21}, permute those according to
		  {0,7,14,21}*27%28 = {0,21,14,7}, then each is head of a length-7 list of indices with decrement 4

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] contain
		x[0,14, 7,21, 1,15, 8,22, 2,16, 9,23, 3,17,10,24, 4,18,11,25, 5,19,12,26, 6,20,13,27], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 24] ... in terms of a[] this translates to a[4 -> 15])
		x[0,14,21, 7,24,10,17, 3,20, 6,13,27,16, 2, 9,23,12,26, 5,19, 8,22, 1,15, 4,18,25,11], which means the a-indices get swapped as
		a[0, 1, 3, 2,15,14,13,12,25,24,26,27, 9, 8,10,11,22,23,20,21, 6, 7, 4, 5,16,17,19,18].
	*/
	/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
					  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
					  to properly re-use the ajp1 variables in the carry-pass version of this routine.
					  (NOTE: Didn't bother to do the swapping in the macro-less version of this routine - would only need to
					   do the swapping there if for some reason we wanted the macro-less DIT code in the carry routines.)
	*/
	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
					/*                                    inputs                                  */ /*                         outputs                   */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);
		RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);
		RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);
		RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);
		RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

	/*...and now do 4 radix-4 transforms...*/
					/*                                              inputs                                          */ /*               intermediates              */ /*                                                                     outputs                                                           */
		RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1    ],a[j2    ],a[j1+p08],a[j2+p08],a[j1+p16],a[j2+p16],a[j1+p24],a[j2+p24],a[j1+p04],a[j2+p04],a[j1+p12],a[j2+p12],a[j1+p20],a[j2+p20],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p07],a[j2+p07],a[j1+p15],a[j2+p15],a[j1+p23],a[j2+p23],a[j1+p03],a[j2+p03],a[j1+p11],a[j2+p11],a[j1+p19],a[j2+p19],a[j1+p27],a[j2+p27],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p14],a[j2+p14],a[j1+p22],a[j2+p22],a[j1+p02],a[j2+p02],a[j1+p10],a[j2+p10],a[j1+p18],a[j2+p18],a[j1+p26],a[j2+p26],a[j1+p06],a[j2+p06],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p21],a[j2+p21],a[j1+p01],a[j2+p01],a[j1+p09],a[j2+p09],a[j1+p17],a[j2+p17],a[j1+p25],a[j2+p25],a[j1+p05],a[j2+p05],a[j1+p13],a[j2+p13],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error pthreaded carry code requires SSE2-enabled build!
	#endif
	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy28_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 28;
		const double crnd = 3.0*0x4000000*0x2000000;
		const double sx0 = 0.44095855184409843174;
		const int odd_radix = 7;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24;
		double *add0, *add1, *add2, *add3;
		struct complex *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
		, *tmp;
		struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26;
		struct complex *cy_i00,*cy_i02,*cy_i04,*cy_i06,*cy_i08,*cy_i10,*cy_i12,*cy_i14,*cy_i16,*cy_i18,*cy_i20,*cy_i22,*cy_i24,*cy_i26;
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27;
		/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
	
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
#if FFT_DEBUG
	int ithread = thread_arg->tid;	/* unique thread index (use for debug) */
	fprintf(dbg_file,"cy28_process_chunk: thread %d, NDIVR = %d, NWT = %d, &rn0,1 = %llx %llx\n"\
		, ithread, thread_arg->ndivr, thread_arg->nwt, (uint64)thread_arg->rn0, (uint64)thread_arg->rn1);
#endif
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;	uint32 nwt16 = nwt << 4;
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int icycle0 = thread_arg->icycle0;
		int icycle1 = thread_arg->icycle1;
		int icycle2 = thread_arg->icycle2;
		int icycle3 = thread_arg->icycle3;
		int icycle4 = thread_arg->icycle4;
		int icycle5 = thread_arg->icycle5;
		int icycle6 = thread_arg->icycle6;
		int jcycle0 = thread_arg->jcycle0;
		int jcycle1 = thread_arg->jcycle1;
		int jcycle2 = thread_arg->jcycle2;
		int jcycle3 = thread_arg->jcycle3;
		int jcycle4 = thread_arg->jcycle4;
		int jcycle5 = thread_arg->jcycle5;
		int jcycle6 = thread_arg->jcycle6;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;	
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;	

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

		s1p00r = thread_arg->s1p00r;	cy_r00	= s1p00r + 0x42;
		s1p01r = s1p00r + 0x02;			cy_r02	= s1p00r + 0x43;
		s1p02r = s1p00r + 0x04;			cy_r04	= s1p00r + 0x44;
		s1p03r = s1p00r + 0x06;			cy_r06	= s1p00r + 0x45;
		s1p04r = s1p00r + 0x08;			cy_r08	= s1p00r + 0x46;
		s1p05r = s1p00r + 0x0a;			cy_r10	= s1p00r + 0x47;
		s1p06r = s1p00r + 0x0c;			cy_r12	= s1p00r + 0x48;
		s1p07r = s1p00r + 0x0e;			cy_r14	= s1p00r + 0x49;
		s1p08r = s1p00r + 0x10;			cy_r16	= s1p00r + 0x4a;
		s1p09r = s1p00r + 0x12;			cy_r18	= s1p00r + 0x4b;
		s1p10r = s1p00r + 0x14;			cy_r20	= s1p00r + 0x4c;
		s1p11r = s1p00r + 0x16;			cy_r22	= s1p00r + 0x4d;
		s1p12r = s1p00r + 0x18;			cy_r24	= s1p00r + 0x4e;
		s1p13r = s1p00r + 0x1a;			cy_r26	= s1p00r + 0x4f;
		s1p14r = s1p00r + 0x1c;			cy_i00	= s1p00r + 0x50;
		s1p15r = s1p00r + 0x1e;			cy_i02	= s1p00r + 0x51;
		s1p16r = s1p00r + 0x20;			cy_i04	= s1p00r + 0x52;
		s1p17r = s1p00r + 0x22;			cy_i06	= s1p00r + 0x53;
		s1p18r = s1p00r + 0x24;			cy_i08	= s1p00r + 0x54;
		s1p19r = s1p00r + 0x26;			cy_i10	= s1p00r + 0x55;
		s1p20r = s1p00r + 0x28;			cy_i12	= s1p00r + 0x56;
		s1p21r = s1p00r + 0x2a;			cy_i14	= s1p00r + 0x57;
		s1p22r = s1p00r + 0x2c;			cy_i16	= s1p00r + 0x58;
		s1p23r = s1p00r + 0x2e;			cy_i18	= s1p00r + 0x59;
		s1p24r = s1p00r + 0x30;			cy_i20	= s1p00r + 0x5a;
		s1p25r = s1p00r + 0x32;			cy_i22	= s1p00r + 0x5b;
		s1p26r = s1p00r + 0x34;			cy_i24	= s1p00r + 0x5c;
		s1p27r = s1p00r + 0x36;			cy_i26	= s1p00r + 0x5d;

		cc0		= s1p00r + 0x38;		max_err = s1p00r + 0x5e;
		ss0		= s1p00r + 0x39;		sse2_rnd= s1p00r + 0x5f;
		cc1		= s1p00r + 0x3a;		half_arr= s1p00r + 0x60;
		ss1		= s1p00r + 0x3b;
		cc2		= s1p00r + 0x3c;
		ss2		= s1p00r + 0x3d;
		cc3  	= s1p00r + 0x3e;
		ss3		= s1p00r + 0x3f;	

		ASSERT(HERE, (ss0->re == sx0 && ss0->im == sx0), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->re == crnd && sse2_rnd->im == crnd), "thread-local memcheck failed!");
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		ASSERT(HERE, (half_arr+10)->re * (half_arr+14)->re == 1.0 && (half_arr+10)->im * (half_arr+14)->im == 1.0, "thread-local memcheck failed!");
	} else {
		ASSERT(HERE, (half_arr)->re * (half_arr+odd_radix)->re == scale && (half_arr)->im * (half_arr+odd_radix)->im == scale, "thread-local memcheck failed!");
	}

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(s1p00r + radix28_creals_in_local_store);
		sse_bw  = sign_mask + 2;
		sse_sw  = sign_mask + 4;
		sse_n   = sign_mask + 6;
		bjmodn00 = (int*)(sign_mask + 8);
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */	/* init carries	*/
			*bjmodn00 = thread_arg->bjmodn00;	cy_r00->re = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;	cy_r00->im = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;	cy_r02->re = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;	cy_r02->im = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;	cy_r04->re = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;	cy_r04->im = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;	cy_r06->re = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;	cy_r06->im = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;	cy_r08->re = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;	cy_r08->im = thread_arg->cy_r09;
			*bjmodn10 = thread_arg->bjmodn10;	cy_r10->re = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;	cy_r10->im = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;	cy_r12->re = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;	cy_r12->im = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;	cy_r14->re = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;	cy_r14->im = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;	cy_r16->re = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;	cy_r16->im = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;	cy_r18->re = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;	cy_r18->im = thread_arg->cy_r19;
			*bjmodn20 = thread_arg->bjmodn20;	cy_r20->re = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;	cy_r20->im = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;	cy_r22->re = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;	cy_r22->im = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;	cy_r24->re = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;	cy_r24->im = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;	cy_r26->re = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;	cy_r26->im = thread_arg->cy_r27;
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
			cy_r00->re = thread_arg->cy_r00;	cy_r00->im = thread_arg->cy_i00;
			cy_r02->re = thread_arg->cy_r01;	cy_r02->im = thread_arg->cy_i01;
			cy_r04->re = thread_arg->cy_r02;	cy_r04->im = thread_arg->cy_i02;
			cy_r06->re = thread_arg->cy_r03;	cy_r06->im = thread_arg->cy_i03;
			cy_r08->re = thread_arg->cy_r04;	cy_r08->im = thread_arg->cy_i04;
			cy_r10->re = thread_arg->cy_r05;	cy_r10->im = thread_arg->cy_i05;
			cy_r12->re = thread_arg->cy_r06;	cy_r12->im = thread_arg->cy_i06;
			cy_r14->re = thread_arg->cy_r07;	cy_r14->im = thread_arg->cy_i07;
			cy_r16->re = thread_arg->cy_r08;	cy_r16->im = thread_arg->cy_i08;
			cy_r18->re = thread_arg->cy_r09;	cy_r18->im = thread_arg->cy_i09;
			cy_r20->re = thread_arg->cy_r10;	cy_r20->im = thread_arg->cy_i10;
			cy_r22->re = thread_arg->cy_r11;	cy_r22->im = thread_arg->cy_i11;
			cy_r24->re = thread_arg->cy_r12;	cy_r24->im = thread_arg->cy_i12;
			cy_r26->re = thread_arg->cy_r13;	cy_r26->im = thread_arg->cy_i13;
			cy_i00->re = thread_arg->cy_r14;	cy_i00->im = thread_arg->cy_i14;
			cy_i02->re = thread_arg->cy_r15;	cy_i02->im = thread_arg->cy_i15;
			cy_i04->re = thread_arg->cy_r16;	cy_i04->im = thread_arg->cy_i16;
			cy_i06->re = thread_arg->cy_r17;	cy_i06->im = thread_arg->cy_i17;
			cy_i08->re = thread_arg->cy_r18;	cy_i08->im = thread_arg->cy_i18;
			cy_i10->re = thread_arg->cy_r19;	cy_i10->im = thread_arg->cy_i19;
			cy_i12->re = thread_arg->cy_r20;	cy_i12->im = thread_arg->cy_i20;
			cy_i14->re = thread_arg->cy_r21;	cy_i14->im = thread_arg->cy_i21;
			cy_i16->re = thread_arg->cy_r22;	cy_i16->im = thread_arg->cy_i22;
			cy_i18->re = thread_arg->cy_r23;	cy_i18->im = thread_arg->cy_i23;
			cy_i20->re = thread_arg->cy_r24;	cy_i20->im = thread_arg->cy_i24;
			cy_i22->re = thread_arg->cy_r25;	cy_i22->im = thread_arg->cy_i25;
			cy_i24->re = thread_arg->cy_r26;	cy_i24->im = thread_arg->cy_i26;
			cy_i26->re = thread_arg->cy_r27;	cy_i26->im = thread_arg->cy_i27;
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

				add0 = &a[j1    ];
				SSE2_RADIX28_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			}
			else	/* Fermat-mod carry in SSE2 mode */
			{
				/* Get the needed Nth root of -1: */
				add1 = &rn0[0];
				add2 = &rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

			  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit 
				SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p14r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p15r,cy_i02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p16r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p17r,cy_i06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p18r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p19r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p20r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p21r,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p22r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p23r,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p24r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p25r,cy_i22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p26r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p27r,cy_i26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
			  #else
				SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck_X2(s1p14r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p16r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p18r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
			  #endif

				icycle0 += wts_idx_inc2;		icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt16);
				icycle1 += wts_idx_inc2;		icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt16);
				icycle2 += wts_idx_inc2;		icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt16);
				icycle3 += wts_idx_inc2;		icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt16);
				icycle4 += wts_idx_inc2;		icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt16);
				icycle5 += wts_idx_inc2;		icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt16);
				icycle6 += wts_idx_inc2;		icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt16);

				jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(int)((uint32)jcycle0 >> 31)) & nwt16);
				jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(int)((uint32)jcycle1 >> 31)) & nwt16);
				jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(int)((uint32)jcycle2 >> 31)) & nwt16);
				jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(int)((uint32)jcycle3 >> 31)) & nwt16);
				jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(int)((uint32)jcycle4 >> 31)) & nwt16);
				jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(int)((uint32)jcycle5 >> 31)) & nwt16);
				jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(int)((uint32)jcycle6 >> 31)) & nwt16);
			}	/* if(MODULUS_TYPE == ...) */


			/*...The radix-28 DIF pass is here:	*/

				add0 = &a[j1    ];
				SSE2_RADIX28_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			thread_arg->cy_r00 = cy_r00->re;
			thread_arg->cy_r01 = cy_r00->im;
			thread_arg->cy_r02 = cy_r02->re;
			thread_arg->cy_r03 = cy_r02->im;
			thread_arg->cy_r04 = cy_r04->re;
			thread_arg->cy_r05 = cy_r04->im;
			thread_arg->cy_r06 = cy_r06->re;
			thread_arg->cy_r07 = cy_r06->im;
			thread_arg->cy_r08 = cy_r08->re;
			thread_arg->cy_r09 = cy_r08->im;
			thread_arg->cy_r10 = cy_r10->re;
			thread_arg->cy_r11 = cy_r10->im;
			thread_arg->cy_r12 = cy_r12->re;
			thread_arg->cy_r13 = cy_r12->im;
			thread_arg->cy_r14 = cy_r14->re;
			thread_arg->cy_r15 = cy_r14->im;
			thread_arg->cy_r16 = cy_r16->re;
			thread_arg->cy_r17 = cy_r16->im;
			thread_arg->cy_r18 = cy_r18->re;
			thread_arg->cy_r19 = cy_r18->im;
			thread_arg->cy_r20 = cy_r20->re;
			thread_arg->cy_r21 = cy_r20->im;
			thread_arg->cy_r22 = cy_r22->re;
			thread_arg->cy_r23 = cy_r22->im;
			thread_arg->cy_r24 = cy_r24->re;
			thread_arg->cy_r25 = cy_r24->im;
			thread_arg->cy_r26 = cy_r26->re;
			thread_arg->cy_r27 = cy_r26->im;
		}
		else
		{
			thread_arg->cy_r00 = cy_r00->re;	thread_arg->cy_i00 = cy_r00->im;
			thread_arg->cy_r01 = cy_r02->re;	thread_arg->cy_i01 = cy_r02->im;
			thread_arg->cy_r02 = cy_r04->re;	thread_arg->cy_i02 = cy_r04->im;
			thread_arg->cy_r03 = cy_r06->re;	thread_arg->cy_i03 = cy_r06->im;
			thread_arg->cy_r04 = cy_r08->re;	thread_arg->cy_i04 = cy_r08->im;
			thread_arg->cy_r05 = cy_r10->re;	thread_arg->cy_i05 = cy_r10->im;
			thread_arg->cy_r06 = cy_r12->re;	thread_arg->cy_i06 = cy_r12->im;
			thread_arg->cy_r07 = cy_r14->re;	thread_arg->cy_i07 = cy_r14->im;
			thread_arg->cy_r08 = cy_r16->re;	thread_arg->cy_i08 = cy_r16->im;
			thread_arg->cy_r09 = cy_r18->re;	thread_arg->cy_i09 = cy_r18->im;
			thread_arg->cy_r10 = cy_r20->re;	thread_arg->cy_i10 = cy_r20->im;
			thread_arg->cy_r11 = cy_r22->re;	thread_arg->cy_i11 = cy_r22->im;
			thread_arg->cy_r12 = cy_r24->re;	thread_arg->cy_i12 = cy_r24->im;
			thread_arg->cy_r13 = cy_r26->re;	thread_arg->cy_i13 = cy_r26->im;
			thread_arg->cy_r14 = cy_i00->re;	thread_arg->cy_i14 = cy_i00->im;
			thread_arg->cy_r15 = cy_i02->re;	thread_arg->cy_i15 = cy_i02->im;
			thread_arg->cy_r16 = cy_i04->re;	thread_arg->cy_i16 = cy_i04->im;
			thread_arg->cy_r17 = cy_i06->re;	thread_arg->cy_i17 = cy_i06->im;
			thread_arg->cy_r18 = cy_i08->re;	thread_arg->cy_i18 = cy_i08->im;
			thread_arg->cy_r19 = cy_i10->re;	thread_arg->cy_i19 = cy_i10->im;
			thread_arg->cy_r20 = cy_i12->re;	thread_arg->cy_i20 = cy_i12->im;
			thread_arg->cy_r21 = cy_i14->re;	thread_arg->cy_i21 = cy_i14->im;
			thread_arg->cy_r22 = cy_i16->re;	thread_arg->cy_i22 = cy_i16->im;
			thread_arg->cy_r23 = cy_i18->re;	thread_arg->cy_i23 = cy_i18->im;
			thread_arg->cy_r24 = cy_i20->re;	thread_arg->cy_i24 = cy_i20->im;
			thread_arg->cy_r25 = cy_i22->re;	thread_arg->cy_i25 = cy_i22->im;
			thread_arg->cy_r26 = cy_i24->re;	thread_arg->cy_i26 = cy_i24->im;
			thread_arg->cy_r27 = cy_i26->re;	thread_arg->cy_i27 = cy_i26->im;
		}
		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}
/*
if(scale < 1.0) {
	printf("cy28_process_chunk: thread %d, maxerr = %20.15f\n", ithread, maxerr);
}
*/
#if FFT_DEBUG
	double cy_sum = 0;
	cy_sum += fabs(thread_arg->cy_r00)+fabs(thread_arg->cy_r01)+fabs(thread_arg->cy_r02)+fabs(thread_arg->cy_r03)+fabs(thread_arg->cy_r04)+fabs(thread_arg->cy_r05)+fabs(thread_arg->cy_r06)+fabs(thread_arg->cy_r07)+fabs(thread_arg->cy_r08)+fabs(thread_arg->cy_r09)+fabs(thread_arg->cy_r10)+fabs(thread_arg->cy_r11)+fabs(thread_arg->cy_r12)+fabs(thread_arg->cy_r13)+fabs(thread_arg->cy_r14)+fabs(thread_arg->cy_r15)+fabs(thread_arg->cy_r16)+fabs(thread_arg->cy_r17)+fabs(thread_arg->cy_r18)+fabs(thread_arg->cy_r19)+fabs(thread_arg->cy_r20)+fabs(thread_arg->cy_r21)+fabs(thread_arg->cy_r22)+fabs(thread_arg->cy_r23)+fabs(thread_arg->cy_r24)+fabs(thread_arg->cy_r25)+fabs(thread_arg->cy_r26)+fabs(thread_arg->cy_r27));
	cy_sum += fabs(thread_arg->cy_i00)+fabs(thread_arg->cy_i01)+fabs(thread_arg->cy_i02)+fabs(thread_arg->cy_i03)+fabs(thread_arg->cy_i04)+fabs(thread_arg->cy_i05)+fabs(thread_arg->cy_i06)+fabs(thread_arg->cy_i07)+fabs(thread_arg->cy_i08)+fabs(thread_arg->cy_i09)+fabs(thread_arg->cy_i10)+fabs(thread_arg->cy_i11)+fabs(thread_arg->cy_i12)+fabs(thread_arg->cy_i13)+fabs(thread_arg->cy_i14)+fabs(thread_arg->cy_i15)+fabs(thread_arg->cy_i16)+fabs(thread_arg->cy_i17)+fabs(thread_arg->cy_i18)+fabs(thread_arg->cy_i19)+fabs(thread_arg->cy_i20)+fabs(thread_arg->cy_i21)+fabs(thread_arg->cy_i22)+fabs(thread_arg->cy_i23)+fabs(thread_arg->cy_i24)+fabs(thread_arg->cy_i25)+fabs(thread_arg->cy_i26)+fabs(thread_arg->cy_i27));
	if(cy_sum != 0.0)
	{
		printf("cy_r00 = %20.10e %20.10e\n", cy_r00->re,cy_r00->im);
		printf("cy_r02 = %20.10e %20.10e\n", cy_r02->re,cy_r02->im);
		printf("cy_r04 = %20.10e %20.10e\n", cy_r04->re,cy_r04->im);
		printf("cy_r06 = %20.10e %20.10e\n", cy_r06->re,cy_r06->im);
		printf("cy_r08 = %20.10e %20.10e\n", cy_r08->re,cy_r08->im);
		printf("cy_r0A = %20.10e %20.10e\n", cy_r10->re,cy_r10->im);
		printf("cy_r0C = %20.10e %20.10e\n", cy_r12->re,cy_r12->im);
		printf("cy_r0E = %20.10e %20.10e\n", cy_r14->re,cy_r14->im);
		printf("cy_r10 = %20.10e %20.10e\n", cy_r16->re,cy_r16->im);
		printf("cy_r12 = %20.10e %20.10e\n", cy_r18->re,cy_r18->im);
		printf("cy_r14 = %20.10e %20.10e\n", cy_r20->re,cy_r20->im);
		printf("cy_r16 = %20.10e %20.10e\n", cy_r22->re,cy_r22->im);
		printf("cy_r18 = %20.10e %20.10e\n", cy_r24->re,cy_r24->im);
		printf("cy_r1A = %20.10e %20.10e\n", cy_r26->re,cy_r26->im);
		printf("cy_r1C = %20.10e %20.10e\n", cy_r00->re,cy_r00->im);
		printf("cy_r1E = %20.10e %20.10e\n", cy_r02->re,cy_r02->im);
		printf("cy_i00 = %20.10e %20.10e\n", cy_i04->re,cy_i04->im);
		printf("cy_i02 = %20.10e %20.10e\n", cy_i06->re,cy_i06->im);
		printf("cy_i04 = %20.10e %20.10e\n", cy_i08->re,cy_i08->im);
		printf("cy_i06 = %20.10e %20.10e\n", cy_i10->re,cy_i10->im);
		printf("cy_i08 = %20.10e %20.10e\n", cy_i12->re,cy_i12->im);
		printf("cy_i0A = %20.10e %20.10e\n", cy_i14->re,cy_i14->im);
		printf("cy_i0C = %20.10e %20.10e\n", cy_i16->re,cy_i16->im);
		printf("cy_i0E = %20.10e %20.10e\n", cy_i18->re,cy_i18->im);
		printf("cy_i10 = %20.10e %20.10e\n", cy_i20->re,cy_i20->im);
		printf("cy_i12 = %20.10e %20.10e\n", cy_i22->re,cy_i22->im);
		printf("cy_i14 = %20.10e %20.10e\n", cy_i24->re,cy_i24->im);
		printf("cy_i16 = %20.10e %20.10e\n", cy_i26->re,cy_i26->im);
		ASSERT(HERE, 0, "Nonzero exit carry in thread-function!");
	}
#endif
		return 0x0;
	}
#endif

