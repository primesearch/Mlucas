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

#include "Mlucas.h"

/***************/

int radix8_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[],double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-8 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-8 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	static int n8;
	int i,j,j1,j2,jstart,jhi,full_pass,k1,k2,k,khi,l,outer;
	// For some reason GCC treats these as possibly-uninited in cmplx_carry_norm_pow2_errcheck(), so init = 0;
	int bjmodn0 = 0,bjmodn1 = 0,bjmodn2 = 0,bjmodn3 = 0,bjmodn4 = 0,bjmodn5 = 0,bjmodn6 = 0,bjmodn7 = 0;
	static uint32 bjmodnini, nsave = 0;
	static uint64 psave = 0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv,n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7
		,temp,frac,scale;
	double maxerr = 0.0;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */

// Stuff for the multithreaded implementation is here: For this small radix we merely *simulate*
// key data-sharing aspects of the true multithreaded approach, which is deployed for larger radices such as 16:
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	static int *_bjmodn0 = 0x0, *_bjmodn1 = 0x0, *_bjmodn2 = 0x0, *_bjmodn3 = 0x0, *_bjmodn4 = 0x0, *_bjmodn5 = 0x0, *_bjmodn6 = 0x0, *_bjmodn7 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0
	, *_cy_r0 = 0x0, *_cy_r1 = 0x0, *_cy_r2 = 0x0, *_cy_r3 = 0x0, *_cy_r4 = 0x0, *_cy_r5 = 0x0, *_cy_r6 = 0x0, *_cy_r7 = 0x0
	, *_cy_i0 = 0x0, *_cy_i1 = 0x0, *_cy_i2 = 0x0, *_cy_i3 = 0x0, *_cy_i4 = 0x0, *_cy_i5 = 0x0, *_cy_i6 = 0x0, *_cy_i7 = 0x0;

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;

	if(RES_SHIFT) {
		WARN(HERE, "CY routines with radix < 16 do not support shifted residues!", "", 1);
		return(ERR_ASSERT);
	}

	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change n8 and n_div_wt to non-static to work around a gcc compiler bug. */
	n8   = n/8;
	n_div_nwt = n8 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n8)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/8 in radix8_ditN_cy_dif1.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave || n != nsave) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
		/* To-do: Support #thread change here! */
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
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
			CY_THREADS = MAX_THREADS;

		ASSERT(CY_THREADS >= NTHREADS,"radix8_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(isPow2(CY_THREADS)    ,"radix8_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(n8       %CY_THREADS == 0,"radix8_ditN_cy_dif1.c: n8      %CY_THREADS != 0 ... likely more threads than this leading radix can handle.");
			ASSERT(n_div_nwt%CY_THREADS == 0,"radix8_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0 ... likely more threads than this leading radix can handle.");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		psave = p;	nsave = n;
		first_entry=FALSE;
		radix_inv = 0.125;
		n2inv = 1.0/(n/2);

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		nm1   = n-1;

		/*   constant index offsets for load/stores are here.	*/
		pini = n8/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p1 = n8 + ( (n8 >> DAT_BITS) << PAD_BITS );
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;

		if(_cy_r0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn0); _bjmodn0 = 0x0;
			free((void *)_bjmodn1); _bjmodn1 = 0x0;
			free((void *)_bjmodn2); _bjmodn2 = 0x0;
			free((void *)_bjmodn3); _bjmodn3 = 0x0;
			free((void *)_bjmodn4); _bjmodn4 = 0x0;
			free((void *)_bjmodn5); _bjmodn5 = 0x0;
			free((void *)_bjmodn6); _bjmodn6 = 0x0;
			free((void *)_bjmodn7); _bjmodn7 = 0x0;

			free((void *)_cy_r0); _cy_r0 = 0x0;	free((void *)_cy_i0); _cy_i0 = 0x0;
			free((void *)_cy_r1); _cy_r1 = 0x0;	free((void *)_cy_i1); _cy_i1 = 0x0;
			free((void *)_cy_r2); _cy_r2 = 0x0;	free((void *)_cy_i2); _cy_i2 = 0x0;
			free((void *)_cy_r3); _cy_r3 = 0x0;	free((void *)_cy_i3); _cy_i3 = 0x0;
			free((void *)_cy_r4); _cy_r4 = 0x0;	free((void *)_cy_i4); _cy_i4 = 0x0;
			free((void *)_cy_r5); _cy_r5 = 0x0;	free((void *)_cy_i5); _cy_i5 = 0x0;
			free((void *)_cy_r6); _cy_r6 = 0x0;	free((void *)_cy_i6); _cy_i6 = 0x0;
			free((void *)_cy_r7); _cy_r7 = 0x0;	free((void *)_cy_i7); _cy_i7 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		_i       = (int *)malloc(CY_THREADS*sizeof(int)); if(!_i     ) { sprintf(cbuf,"ERROR: unable to allocate array _i       in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn0 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn0){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn0 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn1 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn1){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn1 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn2 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn2){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn2 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn3 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn3){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn3 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn4 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn4){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn4 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn5 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn5){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn5 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn6 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn6){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn6 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodn7 = (int *)malloc(CY_THREADS*sizeof(int)); if(!_bjmodn7){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodn7 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_jstart  = (int *)malloc(CY_THREADS*sizeof(int)); if(!_jstart ){ sprintf(cbuf,"ERROR: unable to allocate array _jstart  in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_jhi     = (int *)malloc(CY_THREADS*sizeof(int)); if(!_jhi    ){ sprintf(cbuf,"ERROR: unable to allocate array _jhi     in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_col     = (int *)malloc(CY_THREADS*sizeof(int)); if(!_col    ){ sprintf(cbuf,"ERROR: unable to allocate array _col    in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_co2     = (int *)malloc(CY_THREADS*sizeof(int)); if(!_co2    ){ sprintf(cbuf,"ERROR: unable to allocate array _co2    in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_co3     = (int *)malloc(CY_THREADS*sizeof(int)); if(!_co3    ){ sprintf(cbuf,"ERROR: unable to allocate array _co3    in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		_cy_r0  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r0){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r0 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r1  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r1){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r1 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r2  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r2){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r2 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r3  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r3){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r3 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r4  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r4){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r4 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r5  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r5){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r5 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r6  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r6){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r6 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_r7  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_r7){ sprintf(cbuf,"ERROR: unable to allocate array _cy_r7 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		_cy_i0  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i0){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i0 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i1  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i1){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i1 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i2  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i2){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i2 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i3  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i3){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i3 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i4  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i4){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i4 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i5  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i5){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i5 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i6  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i6){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i6 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_cy_i7  = (double *)malloc(CY_THREADS*sizeof(double)); if(!_cy_i7){ sprintf(cbuf,"ERROR: unable to allocate array _cy_i7 in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		_maxerr = (double *)malloc(CY_THREADS*sizeof(double)); if(!_maxerr){ sprintf(cbuf,"ERROR: unable to allocate array _maxerr in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/8-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int)); if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in radix8_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			_bjmodnini[0] = 0;
			_bjmodnini[1] = 0;
			for(j=0; j < n8/CY_THREADS; j++)
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
			for(j=0; j < n8; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
			ASSERT(_bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");
		}
	}	/* endif(first_entry) */

/*...The radix-8 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r0[ithread] = 0;	_cy_i0[ithread] = 0;
		_cy_r1[ithread] = 0;	_cy_i1[ithread] = 0;
		_cy_r2[ithread] = 0;	_cy_i2[ithread] = 0;
		_cy_r3[ithread] = 0;	_cy_i3[ithread] = 0;
		_cy_r4[ithread] = 0;	_cy_i4[ithread] = 0;
		_cy_r5[ithread] = 0;	_cy_i5[ithread] = 0;
		_cy_r6[ithread] = 0;	_cy_i6[ithread] = 0;
		_cy_r7[ithread] = 0;	_cy_i7[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r0[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		_maxerr[ithread] = 0.0;
    }

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		if(CY_THREADS > 1)
		{
			for(ithread = 1; ithread < CY_THREADS; ithread++)
			{
				_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
			}
		}

		/*
		Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
		then simply overwrite it with 1 prior to starting the k-loop.
		*/
		khi = n_div_nwt/CY_THREADS;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_bjmodn0[ithread] = _bjmodnini[ithread];
			_bjmodn1[ithread] = _bjmodn0[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn1[ithread] = _bjmodn1[ithread] + ( (-(int)((uint32)_bjmodn1[ithread] >> 31)) & n);
			_bjmodn2[ithread] = _bjmodn1[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn2[ithread] = _bjmodn2[ithread] + ( (-(int)((uint32)_bjmodn2[ithread] >> 31)) & n);
			_bjmodn3[ithread] = _bjmodn2[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn3[ithread] = _bjmodn3[ithread] + ( (-(int)((uint32)_bjmodn3[ithread] >> 31)) & n);
			_bjmodn4[ithread] = _bjmodn3[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn4[ithread] = _bjmodn4[ithread] + ( (-(int)((uint32)_bjmodn4[ithread] >> 31)) & n);
			_bjmodn5[ithread] = _bjmodn4[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn5[ithread] = _bjmodn5[ithread] + ( (-(int)((uint32)_bjmodn5[ithread] >> 31)) & n);
			_bjmodn6[ithread] = _bjmodn5[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn6[ithread] = _bjmodn6[ithread] + ( (-(int)((uint32)_bjmodn6[ithread] >> 31)) & n);
			_bjmodn7[ithread] = _bjmodn6[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn7[ithread] = _bjmodn7[ithread] + ( (-(int)((uint32)_bjmodn7[ithread] >> 31)) & n);

			_jstart[ithread] = ithread*n8/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*8);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+8 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-8;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		khi = 1;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*n8/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	}

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodn0 = _bjmodn0[ithread];
			bjmodn1 = _bjmodn1[ithread];
			bjmodn2 = _bjmodn2[ithread];
			bjmodn3 = _bjmodn3[ithread];
			bjmodn4 = _bjmodn4[ithread];
			bjmodn5 = _bjmodn5[ithread];
			bjmodn6 = _bjmodn6[ithread];
			bjmodn7 = _bjmodn7[ithread];
		}

		/* init carries	*/
		cy_r0 = _cy_r0[ithread];	cy_i0 = _cy_i0[ithread];
		cy_r1 = _cy_r1[ithread];	cy_i1 = _cy_i1[ithread];
		cy_r2 = _cy_r2[ithread];	cy_i2 = _cy_i2[ithread];
		cy_r3 = _cy_r3[ithread];	cy_i3 = _cy_i3[ithread];
		cy_r4 = _cy_r4[ithread];	cy_i4 = _cy_i4[ithread];
		cy_r5 = _cy_r5[ithread];	cy_i5 = _cy_i5[ithread];
		cy_r6 = _cy_r6[ithread];	cy_i6 = _cy_i6[ithread];
		cy_r7 = _cy_r7[ithread];	cy_i7 = _cy_i7[ithread];

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
			#ifdef USE_AVX
				j1 = (j & mask02) + br8[j&7];
			#elif defined(USE_SSE2)
				j1 = (j & mask01) + br4[j&3];
			#else
				j1 = j;
			#endif
				j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;

				t1 =a[j1   ];		t2 =a[j2   ];
				rt =a[j1+p1];		it =a[j2+p1];
				t3 =t1 -rt;			t4 =t2 -it;
				t1 =t1 +rt;			t2 =t2 +it;

				t5 =a[j1+p2];		t6 =a[j2+p2];
				rt =a[j1+p3];		it =a[j2+p3];
				t7 =t5 -rt;			t8 =t6 -it;
				t5 =t5 +rt;			t6 =t6 +it;

				t9 =a[j1+p4];		t10=a[j2+p4];
				rt =a[j1+p5];		it =a[j2+p5];
				t11=t9 -rt;			t12=t10-it;
				t9 =t9 +rt;			t10=t10+it;

				t13=a[j1+p6];		t14=a[j2+p6];
				rt =a[j1+p7];		it =a[j2+p7];
				t15=t13-rt;			t16=t14-it;
				t13=t13+rt;			t14=t14+it;

			/*       combine to get the 2 length-4 transforms...	*/
				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
						t8 =t4 +rt;	t4 =t4 -rt;

				rt =t13;t13=t9 -rt;	t9 =t9 +rt;
				it =t14;t14=t10-it;	t10=t10+it;

				rt =t15;t15=t11-t16;t11=t11+t16;
						t16=t12+rt;	t12=t12-rt;
			/*       now combine the two half-transforms	*/
				a1p0r=t1+t9;		a1p0i=t2+t10;
				a1p4r=t1-t9;		a1p4i=t2-t10;

				rt=(t11+t12)*ISRT2;	it=(t11-t12)*ISRT2;
				a1p1r=t3+rt;		a1p1i=t4-it;
				a1p5r=t3-rt;		a1p5i=t4+it;

				a1p2r=t5+t14;		a1p2i=t6-t13;
				a1p6r=t5-t14;		a1p6i=t6+t13;

				rt=(t15-t16)*ISRT2;	it=(t15+t16)*ISRT2;
				a1p3r=t7-rt;		a1p3i=t8-it;
				a1p7r=t7+rt;		a1p7i=t8+it;

				/*...and combine those to complete the radix-8 transform and do the carries. Since the outputs would
				normally be getting dispatched to 8 separate blocks of the A-array, we need 8 separate carries.	*/

				if(MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL)
				{	//			Indices in rightmost col are debug-usage only: vvv
					genfftmul_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,0x0);
					genfftmul_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,0x1);
					genfftmul_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,0x2);
					genfftmul_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,0x3);
					genfftmul_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,0x4);
					genfftmul_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,0x5);
					genfftmul_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,0x6);
					genfftmul_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,0x7);
				}
				else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
				{
					/* We can use the identity WT(J)*WT(N-J) = 2 here to obviate the
					need for an array of inverse weights - if we need the inverse
					of WT(J), just use WT(N-J)/2. Together with a reduced-array weights scheme
					this needs L* = (N - J) & (NWT - 1) = NWT - [J & (NWT - 1)] = NWT - L;
					also needs 	*/

					l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 8 carries, so precompute	*/
					n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
					n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
					sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/
	//	if(!j && full_pass)printf("Iter = %d: A01 = %18.10e,%18.10e\n",iter,a1p0r,a1p0i);
					/*...set0 is slightly different from others:	*/
					cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0,0x0,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p1r,a1p1i,cy_r1,bjmodn1,0x1,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p2r,a1p2i,cy_r2,bjmodn2,0x2,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p3r,a1p3i,cy_r3,bjmodn3,0x3,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p4r,a1p4i,cy_r4,bjmodn4,0x4,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p5r,a1p5i,cy_r5,bjmodn5,0x5,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p6r,a1p6i,cy_r6,bjmodn6,0x6,prp_mult);
					cmplx_carry_norm_pow2_errcheck (a1p7r,a1p7i,cy_r7,bjmodn7,0x7,prp_mult);

					i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
				}
				else
				{
					fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,0x0*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,0x1*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,0x2*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,0x3*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,0x4*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,0x5*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,0x6*n8,NRTM1,NRT_BITS,prp_mult);
					fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,0x7*n8,NRTM1,NRT_BITS,prp_mult);
				}

				/*...The radix-8 DIF pass is here:	*/
			#if PFETCH
				add0 = &a[j1];
				prefetch_p_doubles(add0);
			#endif
				/*       gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and combine to get the 4 length-2 transforms...	*/
				t3 =a1p0r -a1p4r;	t4 =a1p0i -a1p4i;
				t1 =a1p0r +a1p4r;	t2 =a1p0i +a1p4i;

				t7 =a1p2r -a1p6r;	t8 =a1p2i -a1p6i;
				t5 =a1p2r +a1p6r;	t6 =a1p2i +a1p6i;
			#if PFETCH
				addr = add0+p1;
				prefetch_p_doubles(addr);
			#endif
				t11=a1p1r -a1p5r;	t12=a1p1i -a1p5i;
				t9 =a1p1r +a1p5r;	t10=a1p1i +a1p5i;

				t15=a1p3r -a1p7r;	t16=a1p3i -a1p7i;
				t13=a1p3r +a1p7r;	t14=a1p3i +a1p7i;
			#if PFETCH
				addr = add0+p2;
				prefetch_p_doubles(addr);
			#endif
				/*       combine to get the 2 length-4 transforms...	*/
				rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
				it =t6;	t6 =t2 -it;	t2 =t2 +it;

				rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
						t8 =t4 -rt;	t4 =t4 +rt;
			#if PFETCH
				addr = add0+p3;
				prefetch_p_doubles(addr);
			#endif
				rt =t13;t13=t9 -rt;	t9 =t9 +rt;
				it =t14;t14=t10-it;	t10=t10+it;

				rt =t15;t15=t11+t16;t11=t11-t16;
						t16=t12-rt;	t12=t12+rt;
				/*       now combine the two half-transforms	*/
			#if PFETCH
				addr = add0+p4;
				prefetch_p_doubles(addr);
			#endif
				a[j1   ]=t1+t9;	a[j2   ]=t2+t10;
				a[j1+p1]=t1-t9;	a[j2+p1]=t2-t10;
			#if PFETCH
				addr = add0+p5;
				prefetch_p_doubles(addr);
			#endif
				a[j1+p2]=t5-t14;a[j2+p2]=t6+t13;
				a[j1+p3]=t5+t14;a[j2+p3]=t6-t13;
			#if PFETCH
				addr = add0+p6;
				prefetch_p_doubles(addr);
			#endif
				rt =(t11-t12)*ISRT2;	it =(t11+t12)*ISRT2;
				a[j1+p4]=t3+rt;	a[j2+p4]=t4+it;
				a[j1+p5]=t3-rt;	a[j2+p5]=t4-it;
			#if PFETCH
				addr = add0+p7;
				prefetch_p_doubles(addr);
			#endif
				rt =(t15+t16)*ISRT2;	it =(t16-t15)*ISRT2;
				a[j1+p6]=t7-rt;	a[j2+p6]=t8-it;
				a[j1+p7]=t7+rt;	a[j2+p7]=t8+it;

			}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += 8;
				co3 -= 8;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		_cy_r0[ithread] = cy_r0;	_cy_i0[ithread] = cy_i0;
		_cy_r1[ithread] = cy_r1;	_cy_i1[ithread] = cy_i1;
		_cy_r2[ithread] = cy_r2;	_cy_i2[ithread] = cy_i2;
		_cy_r3[ithread] = cy_r3;	_cy_i3[ithread] = cy_i3;
		_cy_r4[ithread] = cy_r4;	_cy_i4[ithread] = cy_i4;
		_cy_r5[ithread] = cy_r5;	_cy_i5[ithread] = cy_i5;
		_cy_r6[ithread] = cy_r6;	_cy_i6[ithread] = cy_i6;
		_cy_r7[ithread] = cy_r7;	_cy_i7[ithread] = cy_i7;

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/8 block into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the radix-8 forward DIF FFT of the first block of 8 complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the 8 outputs of (1);
	(3) Reweight and perform a radix-8 forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next 8 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1 = _cy_r0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(CY_THREADS > 1,"radix8_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];
		}

		_cy_r0[0] =+t15;	/* ...The wraparound carry is here: */
		_cy_r1[0] = t1 ;
		_cy_r2[0] = t3 ;
		_cy_r3[0] = t5 ;
		_cy_r4[0] = t7 ;
		_cy_r5[0] = t9 ;
		_cy_r6[0] = t11;
		_cy_r7[0] = t13;
	}
	else
	{
		t1 = _cy_r0[CY_THREADS - 1];	t2 = _cy_i0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];	t4 = _cy_i1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];	t6 = _cy_i2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];	t8 = _cy_i3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];	t10= _cy_i4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];	t12= _cy_i5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];	t14= _cy_i6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];	t16= _cy_i7[CY_THREADS - 1];

		// Handle valid case of high Re or Im-word < 0 and corr. cyout = 1 here, by positivizing the word:
		if(MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL && (t15 != 0.0 || t16 != 0.0)) {
			// Can't re-use t1-16 for DFT-temps, since those hold the carries:
			double _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16;
			// Must use n8 instead of p1 here since p1 may have pads which are not applied to element-2-slots-before
			j1 = n8-2;	j1 += ( (j1 >> DAT_BITS) << PAD_BITS );
			j2 = j1+RE_IM_STRIDE;
			ASSERT(t15 <= 1.0 && t16 <= 1.0, "genFFTmul expects carryouts = 0 or 1 at top!");
			// Undo the initial dif pass just for the 16 complex terms in question:
			RADIX_08_DIT(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7]
						,_t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16
						,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7]
						,rt,it)
			a[j1    ] *= radix_inv;	a[j2    ] *= radix_inv;
			a[j1+p1 ] *= radix_inv;	a[j2+p1 ] *= radix_inv;
			a[j1+p2 ] *= radix_inv;	a[j2+p2 ] *= radix_inv;
			a[j1+p3 ] *= radix_inv;	a[j2+p3 ] *= radix_inv;
			a[j1+p4 ] *= radix_inv;	a[j2+p4 ] *= radix_inv;
			a[j1+p5 ] *= radix_inv;	a[j2+p5 ] *= radix_inv;
			a[j1+p6 ] *= radix_inv;	a[j2+p6 ] *= radix_inv;
			a[j1+p7 ] *= radix_inv;	a[j2+p7 ] *= radix_inv;
			printf("CYHI.re = %10.3f, im = %10.3f, High words: Re = %10.3f, Im = %10.3f\n",t15,t16,a[j1+p7],a[j2+p7]);
			// Verify that any cyout = 1 has the corresponding high word < 0,
			// then absorb cyout back into the high word and zero the carry:
			if(t15 == 1.0) {
				ASSERT(a[j1+p7] < 0.0, "genFFTmul: Legal Re-cyout = 1 must have the corresponding high word < 0!");
				a[j1+p7] += FFT_MUL_BASE;	t15 = 0.0;
			}
			if(t16 == 1.0) {
				ASSERT(a[j2+p7] < 0.0, "genFFTmul: Legal Im-cyout = 1 must have the corresponding high word < 0!");
				a[j2+p7] += FFT_MUL_BASE;	t16 = 0.0;
			}
			// Redo the initial dif pass just for the 16 complex terms in question:
			RADIX_08_DIF(a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7]
						,_t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11,_t12,_t13,_t14,_t15,_t16
						,a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7]
						,rt,it)
		}

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(CY_THREADS > 1,"radix8_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];		_cy_i0[ithread] = _cy_i0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];		_cy_i1[ithread] = _cy_i1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];		_cy_i2[ithread] = _cy_i2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];		_cy_i3[ithread] = _cy_i3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];		_cy_i4[ithread] = _cy_i4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];		_cy_i5[ithread] = _cy_i5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];		_cy_i6[ithread] = _cy_i6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];		_cy_i7[ithread] = _cy_i7[ithread-1];
		}

		_cy_r0[0] =-t16;	_cy_i0[0] =+t15;	/* ...The 2 Mo"bius carries are here: */
		_cy_r1[0] = t1 ;	_cy_i1[0] = t2 ;
		_cy_r2[0] = t3 ;	_cy_i2[0] = t4 ;
		_cy_r3[0] = t5 ;	_cy_i3[0] = t6 ;
		_cy_r4[0] = t7 ;	_cy_i4[0] = t8 ;
		_cy_r5[0] = t9 ;	_cy_i5[0] = t10;
		_cy_r6[0] = t11;	_cy_i6[0] = t12;
		_cy_r7[0] = t13;	_cy_i7[0] = t14;
	}

	full_pass = 0;
	scale = prp_mult = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if((MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL) || (TRANSFORM_TYPE == RIGHT_ANGLE))
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
			a[j    ] *= radix_inv;
			a[j+p1 ] *= radix_inv;
			a[j+p2 ] *= radix_inv;
			a[j+p3 ] *= radix_inv;
			a[j+p4 ] *= radix_inv;
			a[j+p5 ] *= radix_inv;
			a[j+p6 ] *= radix_inv;
			a[j+p7 ] *= radix_inv;
		}
    }
}	/* endfor(outer) */

    t1 = 0;
    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		t1 += fabs(_cy_r0[ithread])+fabs(_cy_r1[ithread])+fabs(_cy_r2[ithread])+fabs(_cy_r3[ithread])+fabs(_cy_r4[ithread])+fabs(_cy_r5[ithread])+fabs(_cy_r6[ithread])+fabs(_cy_r7[ithread]);
		t1 += fabs(_cy_i0[ithread])+fabs(_cy_i1[ithread])+fabs(_cy_i2[ithread])+fabs(_cy_i3[ithread])+fabs(_cy_i4[ithread])+fabs(_cy_i5[ithread])+fabs(_cy_i6[ithread])+fabs(_cy_i7[ithread]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(t1 != 0.0)
	{
	    sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix8_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
	    mlucas_fprint(cbuf,INTERACT);
	    err=ERR_CARRY;
	    return(err);
	}

	return(0);
}

/***************/

void radix8_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-8 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
	int j,j1,j2;
	static int first_entry=TRUE;
	static uint32 n8,p1,p2,p3,p4,p5,p6,p7;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;

	if(!first_entry && (n >> 3) != n8)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n8 = n >> 3;

		p1 = n8;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-8 pass is here.	*/

    for(j = 0; j < n8; j += 2)
    {
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

/*       gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and combine to get the 4 length-2 transforms...	*/
		t1 =a[j1   ];		t2 =a[j2   ];
		rt =a[j1+p4];		it =a[j2+p4];
		t3 =t1 -rt;			t4 =t2 -it;
		t1 =t1 +rt;			t2 =t2 +it;

		t5 =a[j1+p2];		t6 =a[j2+p2];
		rt =a[j1+p6];		it =a[j2+p6];
		t7 =t5 -rt;			t8 =t6 -it;
		t5 =t5 +rt;			t6 =t6 +it;

		t9 =a[j1+p1];		t10=a[j2+p1];
		rt =a[j1+p5];		it =a[j2+p5];
		t11=t9 -rt;			t12=t10-it;
		t9 =t9 +rt;			t10=t10+it;

		t13=a[j1+p3];		t14=a[j2+p3];
		rt =a[j1+p7];		it =a[j2+p7];
		t15=t13-rt;			t16=t14-it;
		t13=t13+rt;			t14=t14+it;
/*       combine to get the 2 length-4 transforms...	*/
		rt =t5;				it =t6;
		t5 =t1 -rt;			t6 =t2 -it;
		t1 =t1 +rt;			t2 =t2 +it;

		rt =t7;				it =t8;
		t7 =t3 +it;			t8 =t4 -rt;
		t3 =t3 -it;			t4 =t4 +rt;

		rt =t13;			it =t14;
		t13=t9 -rt;			t14=t10-it;
		t9 =t9 +rt;			t10=t10+it;

		rt =t15;			it =t16;
		t15=t11+it;			t16=t12-rt;
		t11=t11-it;			t12=t12+rt;
/*       now combine the two half-transforms	*/
		a[j1   ]=t1+t9;		a[j2   ]=t2+t10;
		a[j1+p1]=t1-t9;		a[j2+p1]=t2-t10;

		a[j1+p2]=t5-t14;	a[j2+p2]=t6+t13;
		a[j1+p3]=t5+t14;	a[j2+p3]=t6-t13;

		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;
		a[j1+p4]=t3+rt;		a[j2+p4]=t4+it;
		a[j1+p5]=t3-rt;		a[j2+p5]=t4-it;

		rt =(t15+t16)*ISRT2;it =(t16-t15)*ISRT2;
		a[j1+p6]=t7-rt;		a[j2+p6]=t8-it;
		a[j1+p7]=t7+rt;		a[j2+p7]=t8+it;
	}

}

/***************/

void radix8_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-8 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix8_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix8_dif_pass1 for details on the algorithm.
*/
	int j,j1,j2;
	static int first_entry=TRUE;
	static uint32 n8,p1,p2,p3,p4,p5,p6,p7;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16;

	if(!first_entry && (n >> 3) != n8)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n8 = n/8;

		p1 = n8;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-8 pass is here.	*/

    for(j=0; j < n8; j += 2)
    {
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

/*       gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and combine to get the 4 length-2 transforms...	*/
		t1 =a[j1   ];		t2 =a[j2   ];
	    rt =a[j1+p1];		it =a[j2+p1];
	    t3 =t1 -rt;			t4 =t2 -it;
	    t1 =t1 +rt;			t2 =t2 +it;

	    t5 =a[j1+p2];		t6 =a[j2+p2];
	    rt =a[j1+p3];		it =a[j2+p3];
	    t7 =t5 -rt;			t8 =t6 -it;
	    t5 =t5 +rt;			t6 =t6 +it;

	    t9 =a[j1+p4];		t10=a[j2+p4];
	    rt =a[j1+p5];		it =a[j2+p5];
	    t11=t9 -rt;			t12=t10-it;
	    t9 =t9 +rt;			t10=t10+it;

	    t13=a[j1+p6];		t14=a[j2+p6];
	    rt =a[j1+p7];		it =a[j2+p7];
	    t15=t13-rt;			t16=t14-it;
	    t13=t13+rt;			t14=t14+it;
/*       combine to get the 2 length-4 transform...	*/
	    rt =t5;				it =t6;
	    t5 =t1 -rt;			t6 =t2 -it;
	    t1 =t1 +rt;			t2 =t2 +it;

	    rt =t7;				it =t8;
	    t7 =t3 -it;			t8 =t4 +rt;
	    t3 =t3 +it;			t4 =t4 -rt;

	    rt =t13;			it =t14;
	    t13=t9 -rt;			t14=t10-it;
	    t9 =t9 +rt;			t10=t10+it;

	    rt =t15;			it =t16;
	    t15=t11-it;			t16=t12+rt;
	    t11=t11+it;			t12=t12-rt;
/*       now combine the two half-transforms	*/
	    a[j1   ]=t1+t9;		a[j2   ]=t2+t10;
	    a[j1+p4]=t1-t9;		a[j2+p4]=t2-t10;

	    rt=(t11+t12)*ISRT2;	it=(t11-t12)*ISRT2;
	    a[j1+p1]=t3+rt;		a[j2+p1]=t4-it;
	    a[j1+p5]=t3-rt;		a[j2+p5]=t4+it;

	    a[j1+p2]=t5+t14;	a[j2+p2]=t6-t13;
	    a[j1+p6]=t5-t14;	a[j2+p6]=t6+t13;

	    rt=(t15-t16)*ISRT2;	it=(t15+t16)*ISRT2;
	    a[j1+p3]=t7-rt;		a[j2+p3]=t8-it;
	    a[j1+p7]=t7+rt;		a[j2+p7]=t8+it;
	}
}

