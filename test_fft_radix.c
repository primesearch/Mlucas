/* To test a particular DFT routine:
1. Edit this file to set the desired DFT params [RADIX and TEST_TYPE] and compile it via 'gcc -c test_fft_radix.c';
2. Compile util.c via 'gcc -c -DTEST_FFT_RADIX util.c' [Note: If -DUSE_THREADS was used for other sources, need to use it here, too];
3. Link into an Mlucas binary and run-as-if-timing-test-at-any-desired-length to test the DFT.
*/

#include "Mlucas.h"

#define RADIX	28

#define TEST_TYPE	2	// 0 = DIF, 1 = DIT-but-just-show-input-index-scramblings-needed, 2 = DIT,
						// 3 = DIF+DIT (which should return the original inputs after dividing by N)

void matmul_double (double **, double *, double *, int, int);
void matmul_complex(struct complex **, struct complex *, struct complex *, int, int);

double ISRT2 = .70710678118654752440;

#define CABS(a,b)	sqrt((a)*(a) + (b)*(b))

void test_fft_radix(void)
{
	double iradix = 1.0/RADIX;
	int rmul = 2*RE_IM_STRIDE;
	int i,j,k,l,nradices, *index = 0x0, nerr = 0, pow2, podd;
	int radix_prim[10];
	int *dit_scramble = 0x0;	/* This holds the input-perm for the DIF of the given length ... compute this using the primitive radices */
	double err_r, err_i, abserr, maxerr, avgerr;
	/* "Random" inputs are just decimal digits of Pi ... make this as big as needed, currently support up to complex length = 1024.
	For the curious, the number of occurrences of each decimal digit in these first 2048 is as follows:

	0	185
	1	213
	2	215
	3	191
	4	200
	5	212
	6	205
	7	200
	8	208
	9	219

	The latest occurrence of any digit 0-9 is 0, which does not occur until the 33rd digit.
	*/
	const char ref[2048] = {
	3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5,0,2,8,8,4,1,9,7,1,6,9,3,9,9,3,7,5,1,0,5,8,2,0,9,7,4,9,4,4,5,9,2,
	3,0,7,8,1,6,4,0,6,2,8,6,2,0,8,9,9,8,6,2,8,0,3,4,8,2,5,3,4,2,1,1,7,0,6,7,9,8,2,1,4,8,0,8,6,5,1,3,2,8,2,3,0,6,6,4,7,0,9,3,8,4,4,6,
	0,9,5,5,0,5,8,2,2,3,1,7,2,5,3,5,9,4,0,8,1,2,8,4,8,1,1,1,7,4,5,0,2,8,4,1,0,2,7,0,1,9,3,8,5,2,1,1,0,5,5,5,9,6,4,4,6,2,2,9,4,8,9,5,
	4,9,3,0,3,8,1,9,6,4,4,2,8,8,1,0,9,7,5,6,6,5,9,3,3,4,4,6,1,2,8,4,7,5,6,4,8,2,3,3,7,8,6,7,8,3,1,6,5,2,7,1,2,0,1,9,0,9,1,4,5,6,4,8,
	5,6,6,9,2,3,4,6,0,3,4,8,6,1,0,4,5,4,3,2,6,6,4,8,2,1,3,3,9,3,6,0,7,2,6,0,2,4,9,1,4,1,2,7,3,7,2,4,5,8,7,0,0,6,6,0,6,3,1,5,5,8,8,1,
	7,4,8,8,1,5,2,0,9,2,0,9,6,2,8,2,9,2,5,4,0,9,1,7,1,5,3,6,4,3,6,7,8,9,2,5,9,0,3,6,0,0,1,1,3,3,0,5,3,0,5,4,8,8,2,0,4,6,6,5,2,1,3,8,
	4,1,4,6,9,5,1,9,4,1,5,1,1,6,0,9,4,3,3,0,5,7,2,7,0,3,6,5,7,5,9,5,9,1,9,5,3,0,9,2,1,8,6,1,1,7,3,8,1,9,3,2,6,1,1,7,9,3,1,0,5,1,1,8,
	5,4,8,0,7,4,4,6,2,3,7,9,9,6,2,7,4,9,5,6,7,3,5,1,8,8,5,7,5,2,7,2,4,8,9,1,2,2,7,9,3,8,1,8,3,0,1,1,9,4,9,1,2,9,8,3,3,6,7,3,3,6,2,4,
	4,0,6,5,6,6,4,3,0,8,6,0,2,1,3,9,4,9,4,6,3,9,5,2,2,4,7,3,7,1,9,0,7,0,2,1,7,9,8,6,0,9,4,3,7,0,2,7,7,0,5,3,9,2,1,7,1,7,6,2,9,3,1,7,
	6,7,5,2,3,8,4,6,7,4,8,1,8,4,6,7,6,6,9,4,0,5,1,3,2,0,0,0,5,6,8,1,2,7,1,4,5,2,6,3,5,6,0,8,2,7,7,8,5,7,7,1,3,4,2,7,5,7,7,8,9,6,0,9,
	1,7,3,6,3,7,1,7,8,7,2,1,4,6,8,4,4,0,9,0,1,2,2,4,9,5,3,4,3,0,1,4,6,5,4,9,5,8,5,3,7,1,0,5,0,7,9,2,2,7,9,6,8,9,2,5,8,9,2,3,5,4,2,0,
	1,9,9,5,6,1,1,2,1,2,9,0,2,1,9,6,0,8,6,4,0,3,4,4,1,8,1,5,9,8,1,3,6,2,9,7,7,4,7,7,1,3,0,9,9,6,0,5,1,8,7,0,7,2,1,1,3,4,9,9,9,9,9,9,
	8,3,7,2,9,7,8,0,4,9,9,5,1,0,5,9,7,3,1,7,3,2,8,1,6,0,9,6,3,1,8,5,9,5,0,2,4,4,5,9,4,5,5,3,4,6,9,0,8,3,0,2,6,4,2,5,2,2,3,0,8,2,5,3,
	3,4,4,6,8,5,0,3,5,2,6,1,9,3,1,1,8,8,1,7,1,0,1,0,0,0,3,1,3,7,8,3,8,7,5,2,8,8,6,5,8,7,5,3,3,2,0,8,3,8,1,4,2,0,6,1,7,1,7,7,6,6,9,1,
	4,7,3,0,3,5,9,8,2,5,3,4,9,0,4,2,8,7,5,5,4,6,8,7,3,1,1,5,9,5,6,2,8,6,3,8,8,2,3,5,3,7,8,7,5,9,3,7,5,1,9,5,7,7,8,1,8,5,7,7,8,0,5,3,
	2,1,7,1,2,2,6,8,0,6,6,1,3,0,0,1,9,2,7,8,7,6,6,1,1,1,9,5,9,0,9,2,1,6,4,2,0,1,9,8,9,3,8,0,9,5,2,5,7,2,0,1,0,6,5,4,8,5,8,6,3,2,7,8,
	8,6,5,9,3,6,1,5,3,3,8,1,8,2,7,9,6,8,2,3,0,3,0,1,9,5,2,0,3,5,3,0,1,8,5,2,9,6,8,9,9,5,7,7,3,6,2,2,5,9,9,4,1,3,8,9,1,2,4,9,7,2,1,7,
	7,5,2,8,3,4,7,9,1,3,1,5,1,5,5,7,4,8,5,7,2,4,2,4,5,4,1,5,0,6,9,5,9,5,0,8,2,9,5,3,3,1,1,6,8,6,1,7,2,7,8,5,5,8,8,9,0,7,5,0,9,8,3,8,
	1,7,5,4,6,3,7,4,6,4,9,3,9,3,1,9,2,5,5,0,6,0,4,0,0,9,2,7,7,0,1,6,7,1,1,3,9,0,0,9,8,4,8,8,2,4,0,1,2,8,5,8,3,6,1,6,0,3,5,6,3,7,0,7,
	6,6,0,1,0,4,7,1,0,1,8,1,9,4,2,9,5,5,5,9,6,1,9,8,9,4,6,7,6,7,8,3,7,4,4,9,4,4,8,2,5,5,3,7,9,7,7,4,7,2,6,8,4,7,1,0,4,0,4,7,5,3,4,6,
	4,6,2,0,8,0,4,6,6,8,4,2,5,9,0,6,9,4,9,1,2,9,3,3,1,3,6,7,7,0,2,8,9,8,9,1,5,2,1,0,4,7,5,2,1,6,2,0,5,6,9,6,6,0,2,4,0,5,8,0,3,8,1,5,
	0,1,9,3,5,1,1,2,5,3,3,8,2,4,3,0,0,3,5,5,8,7,6,4,0,2,4,7,4,9,6,4,7,3,2,6,3,9,1,4,1,9,9,2,7,2,6,0,4,2,6,9,9,2,2,7,9,6,7,8,2,3,5,4,
	7,8,1,6,3,6,0,0,9,3,4,1,7,2,1,6,4,1,2,1,9,9,2,4,5,8,6,3,1,5,0,3,0,2,8,6,1,8,2,9,7,4,5,5,5,7,0,6,7,4,9,8,3,8,5,0,5,4,9,4,5,8,8,5,
	8,6,9,2,6,9,9,5,6,9,0,9,2,7,2,1,0,7,9,7,5,0,9,3,0,2,9,5,5,3,2,1,1,6,5,3,4,4,9,8,7,2,0,2,7,5,5,9,6,0,2,3,6,4,8,0,6,6,5,4,9,9,1,1,
	9,8,8,1,8,3,4,7,9,7,7,5,3,5,6,6,3,6,9,8,0,7,4,2,6,5,4,2,5,2,7,8,6,2,5,5,1,8,1,8,4,1,7,5,7,4,6,7,2,8,9,0,9,7,7,7,7,2,7,9,3,8,0,0,
	0,8,1,6,4,7,0,6,0,0,1,6,1,4,5,2,4,9,1,9,2,1,7,3,2,1,7,2,1,4,7,7,2,3,5,0,1,4,1,4,4,1,9,7,3,5,6,8,5,4,8,1,6,1,3,6,1,1,5,7,3,5,2,5,
	5,2,1,3,3,4,7,5,7,4,1,8,4,9,4,6,8,4,3,8,5,2,3,3,2,3,9,0,7,3,9,4,1,4,3,3,3,4,5,4,7,7,6,2,4,1,6,8,6,2,5,1,8,9,8,3,5,6,9,4,8,5,5,6,
	2,0,9,9,2,1,9,2,2,2,1,8,4,2,7,2,5,5,0,2,5,4,2,5,6,8,8,7,6,7,1,7,9,0,4,9,4,6,0,1,6,5,3,4,6,6,8,0,4,9,8,8,6,2,7,2,3,2,7,9,1,7,8,6,
	0,8,5,7,8,4,3,8,3,8,2,7,9,6,7,9,7,6,6,8,1,4,5,4,1,0,0,9,5,3,8,8,3,7,8,6,3,6,0,9,5,0,6,8,0,0,6,4,2,2,5,1,2,5,2,0,5,1,1,7,3,9,2,9,
	8,4,8,9,6,0,8,4,1,2,8,4,8,8,6,2,6,9,4,5,6,0,4,2,4,1,9,6,5,2,8,5,0,2,2,2,1,0,6,6,1,1,8,6,3,0,6,7,4,4,2,7,8,6,2,2,0,3,9,1,9,4,9,4,
	5,0,4,7,1,2,3,7,1,3,7,8,6,9,6,0,9,5,6,3,6,4,3,7,1,9,1,7,2,8,7,4,6,7,7,6,4,6,5,7,5,7,3,9,6,2,4,1,3,8,9,0,8,6,5,8,3,2,6,4,5,9,9,5,
	8,1,3,3,9,0,4,7,8,0,2,7,5,9,0,0,9,9,4,6,5,7,6,4,0,7,8,9,5,1,2,6,9,4,6,8,3,9,8,3,5,2,5,9,5,7,0,9,8,2,5,8,2,2,6,2,0,5,2,2,4,8,9,4};

	const char* test_info_str[] = {"DIF","DIT (but only show index-scramblings)","DIT","Combined DIF+DIT"};
	double *a = 0x0, *b = 0x0, *arrtmp = 0x0, *ptmp = 0x0;
	struct complex *ac, *bc;
	struct complex **mat = 0x0, **matp = 0x0, **ctmpp = 0x0, *ctmp = 0x0;
	double t0,t1,t2,t3;
	double theta, twopi = 6.2831853071795864769;

	/********* allocate all radix-dependent arrays dynamically: ********/
	index        = ALLOC_INT(index       , RADIX);
	dit_scramble = ALLOC_INT(dit_scramble, RADIX);
	/* double a[rmul*RADIX], b[rmul*RADIX], arrtmp[rmul*RADIX]: */
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT(HERE, (ptmp != 0x0), "FATAL: unable to allocate array A in test_fft_radix.\n");
	a    = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ac = (struct complex *)a;
	ASSERT(HERE, ((long)((void *)a) & 63) == 0x0,"test_fft_radix: A[] not aligned on 64-byte boundary!");
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT(HERE, (ptmp != 0x0), "FATAL: unable to allocate array B in test_fft_radix.\n");
	b    = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ASSERT(HERE, ((long)((void *)b) & 63) == 0x0,"test_fft_radix: B[] not aligned on 64-byte boundary!");
	bc = (struct complex *)b;
	ptmp = ALLOC_DOUBLE(ptmp, rmul*RADIX);	ASSERT(HERE, (ptmp != 0x0), "FATAL: unable to allocate array A_ptmp in test_fft_radix.\n");
	arrtmp = ALIGN_DOUBLE(ptmp);	ptmp = 0x0;
	ASSERT(HERE, ((long)((void *)arrtmp) & 63) == 0x0,"test_fft_radix: arrtmp[] not aligned on 64-byte boundary!");
	/* struct complex mat[radix][RADIX], *matp[RADIX]: */
	ctmpp = ALLOC_POINTER(ctmpp,struct complex*, RADIX);	ASSERT(HERE, (ctmpp != 0x0), "FATAL: unable to allocate array MATP in test_fft_radix.\n");
	matp  = ALIGN_POINTER(ctmpp,struct complex*);
	ctmpp = ALLOC_POINTER(ctmpp,struct complex*, RADIX);	ASSERT(HERE, (ctmpp != 0x0), "FATAL: unable to allocate array MAT[][] in test_fft_radix.\n");
	mat   = ALIGN_POINTER(ctmpp,struct complex*);
	for(i = 0; i < RADIX; ++i) {
		ctmp = ALLOC_COMPLEX(ctmp, RADIX);	ASSERT(HERE, (ctmp != 0x0), "FATAL: unable to allocate array Ctmp in test_fft_radix.\n");
		mat[i] = ALIGN_COMPLEX(ctmp);
		ctmp = 0x0;	/* Must re-init pointer so the realloc used by the ALLOC macro allocates new fresh memory for each row */
	}

	fprintf(stderr, "test_fft_radix: Testing radix-%d %s dft:\n", RADIX, test_info_str[TEST_TYPE]);

	/* Power-of-2 component of the DFT length: */
	pow2 = 1 << trailz32(RADIX);
	podd = RADIX >> trailz32(RADIX);
	ASSERT(HERE, RADIX == pow2*podd, "Radix decomposition failed!");
	ASSERT(HERE, (podd < 32), "test_fft_radix: Illegal radix; must be odd*2^n with odd < 16");
	/* These may not have been init'ed yet, so do it here: */
	DAT_BITS = DAT_BITS_DEF;
	PAD_BITS = PAD_BITS_DEF;

	/* Init dit_scramble array: This is the input permutation applied to DIT inputs, in top of the bit-reversal reordering: */
	k = 0;	/* Current perm index */
	l = 0;	/* Current perm value */
	for(i = 0; i < pow2; i++)
	{
		for(j = 0; j < podd; ++j)
		{
			dit_scramble[k++] = l;
			l = (l - pow2)%RADIX;	if(l < 0) { l += RADIX; }
		}
		l = (l - podd)%RADIX;	if(l < 0) { l += RADIX; }
	}

	/* Init data array in scalar-layout mode for reference matrix-multiply-DFT computation: */
	t0 = t1 = 0.;
	for(i = 0; i < RADIX ; i++)
	{
		a[2*i  ] = ref[2*i  ];	t0 += ref[2*i  ];
		a[2*i+1] = ref[2*i+1];	t1 += ref[2*i+1];
	}
	printf("DC signal components: sum[Re,Im] = %15.5f  %15.5f\n",t0,t1);

	/* Init DFT matrix */
	for(i = 0; i < RADIX; i++)
	{
		theta = i*twopi/RADIX;
		for(j = 0; j < RADIX; j++)
		{
			mat[i][j].re = cos(j*theta);
			mat[i][j].im = sin(j*theta);
	/*printf("mat[%4d][%4d] = %15.10f  %15.10f\n",i,j, mat[i][j].re, mat[i][j].im);*/
		}
	}
	matmul_complex(mat,ac,bc,RADIX,RADIX);

	/* In SSE2 mode re-Init data array, using [re,re,im,im] data layout: */
#ifdef USE_SSE2
	ASSERT(HERE, rmul == 4, "!");
	for(i = 0; i < RADIX ; i++)
	{
		a[ 2*i   *RE_IM_STRIDE  ] = ref[2*i  ];
		a[ 2*i   *RE_IM_STRIDE+1] = ref[2*i  ];
		a[(2*i+1)*RE_IM_STRIDE  ] = ref[2*i+1];
		a[(2*i+1)*RE_IM_STRIDE+1] = ref[2*i+1];
	}
#endif

	// If doing a DIF followed by a DIT, save a copy of the original data:
#if TEST_TYPE == 3
	for(i = 0; i < RADIX ; i++)
	{
		arrtmp[2*i  ] = ref[2*i  ];
		arrtmp[2*i+1] = ref[2*i+1];
	}
#endif

/*...Forward (DIF) FFT sincos data are in bit-reversed order.	*/

	l =0;

	switch(RADIX){
	case 2 :
		nradices = 1;
		radix_prim[l++] = 2; break;
	case 3 :
		nradices = 1;
		radix_prim[l++] = 3; break;
	case 4 :
		nradices = 2;
		radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 5 :
		nradices = 1;
		radix_prim[l++] = 5; break;
	case 6 :
		nradices = 2;
		radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 7 :
		nradices = 1;
		radix_prim[l++] = 7; break;
	case 8 :
		nradices = 3;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 9 :
		nradices = 2;
		radix_prim[l++] = 3; radix_prim[l++] = 3; break;
	case 10 :
		nradices = 2;
		radix_prim[l++] = 5; radix_prim[l++] = 2; break;
	case 11 :
		nradices = 1;
		radix_prim[l++] =11; break;
	case 12 :
		nradices = 3;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 13 :
		nradices = 1;
		radix_prim[l++] =13; break;
	case 14 :
		nradices = 2;
		radix_prim[l++] = 7; radix_prim[l++] = 2; break;
	case 15 :
		nradices = 2;
		radix_prim[l++] = 5; radix_prim[l++] = 3; break;
	case 16 :
		nradices = 4;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 18 :
		nradices = 3;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 20 :
		nradices = 3;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 22 :
		nradices = 2;
		radix_prim[l++] =11; radix_prim[l++] = 2; break;
	case 24 :
		nradices = 4;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 26 :
		nradices = 2;
		radix_prim[l++] =13; radix_prim[l++] = 2; break;
	case 28 :
		nradices = 3;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 30 :
		nradices = 3;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
	case 31 :
		nradices = 1;
		radix_prim[l++] =31; break;
	case 32 :
		nradices = 5;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 36 :
		nradices = 4;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 40:
		nradices = 4;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 44 :
		nradices = 3;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 48 :
		nradices = 5;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 52 :
		nradices = 3;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 56 :
		nradices = 4;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 60 :
		nradices = 4;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 64 :
		nradices = 6;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 72 :
		nradices = 5;
		radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 80 :
		nradices = 5;
		radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 88 :
		nradices = 4;
		radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 96 :
		nradices = 6;
		radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 104:
		nradices = 4;
		radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 112:
		nradices = 5;
		radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 120:
		nradices = 5;
		radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 128:
		nradices = 7;
		radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	case 992:
		nradices = 6;
		radix_prim[l++] =31; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
	default :
		printf("FATAL: radix %d not available. Halting...\n",RADIX); exit(EXIT_FAILURE);
	}

/*...Allocate and initialize an index array containing (RADIX) indices...	*/

	for(i=0; i < RADIX; i++)
	{
		index[i]=i;
	}

/*...then bit-reverse INDEX with respect to the [nradices] primitive radices.
		 The order of radices sent to bit_reverse_int is the reverse of that in which these radices are processed
		 in the forward (decimation in frequency) FFT. This is moot for a power-of-2 FFT (or any FFT whose length
		 is a prime power), but necessary for general vector lengths which are a product of 2 or more distinct primes.

		 If the current (Ith) radix is composite with distinct prime factors (e.g. 15 = 3*5), we must specify these
		 factors here in the opposite order from that which is used in the actual FFT-pass routine. For example,
		 if the radix-15 pass implementation does 5 radix-3 DFTs, followed by 3 radix-5 DFTs, then we send (3,5)
		 as the corresponding reverse-ordered prime radices to the bit-reversal routine, not (5,3).	*/

	bit_reverse_int(&index[0],RADIX,nradices,&radix_prim[nradices-1],-1,(int *)0x0);
/*
	printf("bit-reversal index array = [");
	for(i=0; i < RADIX; i++)
	{
		printf(" %d",index[i]);
	}
	printf("]\n\n");
*/
	/* DIT-specific index-permute diagnostics ... the final 'Combined' permutation is that needed for DIT inputs: */
	if(TEST_TYPE > 0)
	{
		printf("DIT input-scramble array = [");
		for(i=0; i < RADIX; i++)
		{
			printf("%3d,",dit_scramble[i]);
		}
		printf("]\n\n");

		printf("Bit-reversal array = [");
		for(i=0; i < RADIX; i++)
		{
			j = index[i];
			printf("%3d,",j);
		}
		printf("]\n\n");

		printf("DIT input-scramble + bit-reversal array = [");
		for(i=0; i < RADIX; i++)
		{
			j = dit_scramble[index[i]];
			printf("%3d,",j);
		}
		printf("]\n\n");

		/* Now find the location of each j-valued element above in the bit-reversal vector ... the *location* in that vector equals the final scrambled-DIT-input index: */
		printf("Combined DIT input-scramble array = [");
		for(i=0; i < RADIX; i++)
		{
			j = dit_scramble[index[i]];
			for(k = 0; k < RADIX; ++k)
			{
				if(j == index[k])
				{
					printf("%3d,",k);
					break;
				}
			}
		}
		printf("]\n\n");
	}
#if TEST_TYPE == 1
	exit(0);
#endif

#if TEST_TYPE == 0 || TEST_TYPE == 3

	#if RADIX == 2
		radix2_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 3
		radix3_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 4
		radix4_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 5
		radix5_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 6
		radix6_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 7
		radix7_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 8
		radix8_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 9
		radix9_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 10
		radix10_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 11
		radix11_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 12
		radix12_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 13
		radix13_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 14
		radix14_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 15
		radix15_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 16
		radix16_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 18
		radix18_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 20
		radix20_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 22
		radix22_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 24
		radix24_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 26
		radix26_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 28
		radix28_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 30
		radix30_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 31
		radix31_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 32
		radix32_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 36
		radix36_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 40
		radix40_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 44
		radix44_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 48
		radix48_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 52
		radix52_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 56
		radix56_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 60
		radix60_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 64
		radix64_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 72
		radix72_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 80
		radix80_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 88
		radix88_dif_pass1 (a,rmul*RADIX);
	#elif RADIX == 96
		radix96_dif_pass1 (a,rmul*RADIX);
	#elif RADIX ==104
		radix104_dif_pass1(a,rmul*RADIX);
	#elif RADIX ==112
		radix112_dif_pass1(a,rmul*RADIX);
	#elif RADIX ==120
		radix120_dif_pass1(a,rmul*RADIX);
	#elif RADIX ==128
		radix128_dif_pass1(a,rmul*RADIX);
	#elif RADIX ==992
		radix31x32_dif_pass1(a,rmul*RADIX);
	#endif

#endif

#if TEST_TYPE == 0

	nerr = 0;
	printf("To deduce the required output-idx ordering, sort Actual-outputs data by left col-of-real-parts, move\n");
	printf("resulting [re,im,i] block to file, repeat procedure for Expected-outputs, then compare the two files.\n");
	printf("If they differ only in their respective rcol indices, those two cols give the required index mapping.\n");
	printf("\n");
	printf("         Actual outputs:             i          Expected outputs:            i BR oidx:\n");
	printf(" -------------------------------   ---   -------------------------------   --- -------\n");
	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
	#ifdef USE_SSE2
		ASSERT(HERE, a[2*i*RE_IM_STRIDE] == a[2*i*RE_IM_STRIDE+1] && a[(2*i+1)*RE_IM_STRIDE] == a[(2*i+1)*RE_IM_STRIDE+1], "1/2 components of SSE2-pack mismatch!");
	#endif
		err_r = fabs(a[2*i*RE_IM_STRIDE]-b[2*j]);  err_i = fabs(a[(2*i+1)*RE_IM_STRIDE]-b[2*j+1]);
		abserr = CABS(err_r, err_i);
		if(abserr < 0.00001)
		{
			printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d  (%4d)\n",a[2*i*RE_IM_STRIDE],a[(2*i+1)*RE_IM_STRIDE],i,b[2*j  ],b[2*j+1],i,j);
		}
		else
		{
			++nerr;
			printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d  (%4d), ERR= %15.10e\n",a[2*i*RE_IM_STRIDE],a[(2*i+1)*RE_IM_STRIDE],i,b[2*j  ],b[2*j+1],i,j, abserr);
		}
		/*
		j=i;
		printf("%4d  %20.10f  %20.10f  %20.10f  %20.10f\n",i,bc[j].re,bc[j].im,a[2*i  ],a[2*i+1]);
		*/
	}
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIF transform!");

#endif

#if TEST_TYPE == 2

	// Bit-reverse the inputs to the transform...
  #ifdef USE_SSE2
	t0 = t1 = t2 = t3 = 0.;
	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
		arrtmp[ 2*i   *RE_IM_STRIDE  ] = a[2*j*RE_IM_STRIDE  ];
		arrtmp[ 2*i   *RE_IM_STRIDE+1] = a[2*j*RE_IM_STRIDE+1];
		arrtmp[(2*i+1)*RE_IM_STRIDE  ] = -a[(2*j+1)*RE_IM_STRIDE  ];
		arrtmp[(2*i+1)*RE_IM_STRIDE+1] = -a[(2*j+1)*RE_IM_STRIDE+1];
		t0 += arrtmp[ 2*i   *RE_IM_STRIDE  ];
		t1 += arrtmp[ 2*i   *RE_IM_STRIDE+1];
		t2 += arrtmp[(2*i+1)*RE_IM_STRIDE  ];
		t3 += arrtmp[(2*i+1)*RE_IM_STRIDE+1];
/*printf("J = [%3d]: add %6d, %6d\n",j,(int)a[2*j  ],(int)a[2*j+1]);*/
	}
/*printf("sum[Re,Im] = %15.5f  %15.5f\n",t0,t2);*/
	ASSERT(HERE, t0==t1 && t2==t3, "!");
	for(i = 0; i < rmul*RADIX ; i+=2)
	{
		a[i  ] = arrtmp[i  ];
		a[i+1] = arrtmp[i+1];
		ASSERT(HERE, a[i  ] == a[i+1], "!");
	}
  #else
	for(i = 0; i < RADIX ; i++)
	{
		j = index[i];
		arrtmp[2*i  ] = a[2*j  ];
		arrtmp[2*i+1] =-a[2*j+1];
	}
	for(i = 0; i < rmul*RADIX ; i++)
	{
		a[i] = arrtmp[i];
	}
  #endif

#endif

#if TEST_TYPE > 0

	#if RADIX == 2
		radix2_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 3
		radix3_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 4
		radix4_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 5
		radix5_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 6
		radix6_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 7
		radix7_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 8
		radix8_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 9
		radix9_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 10
		radix10_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 11
		radix11_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 12
		radix12_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 13
		radix13_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 14
		radix14_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 15
		radix15_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 16
		radix16_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 18
		radix18_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 20
		radix20_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 22
		radix22_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 24
		radix24_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 26
		radix26_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 28
		radix28_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 30
		radix30_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 31
		radix31_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 32
		radix32_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 36
		radix36_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 40
		radix40_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 44
		radix44_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 48
		radix48_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 52
		radix52_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 56
		radix56_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 60
		radix60_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 64
		radix64_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 72
		radix72_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 80
		radix80_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 88
		radix88_dit_pass1 (a,rmul*RADIX);
	#elif RADIX == 96
		radix96_dit_pass1 (a,rmul*RADIX);
	#elif RADIX ==104
		radix104_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==112
		radix112_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==120
		radix120_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==128
		radix128_dit_pass1(a,rmul*RADIX);
	#elif RADIX ==992
		radix31x32_dit_pass1(a,rmul*RADIX);
	#endif

#endif

#if TEST_TYPE == 2

	nerr = 0;
	printf("To deduce the required output-idx ordering, sort Actual-outputs data by left col-of-real-parts, move\n");
	printf("resulting [re,im,i] block to file, repeat procedure for Expected-outputs, then compare the two files.\n");
	printf("If they differ only in their respective rcol indices, those two cols give the required index mapping.\n");
	printf("\n");
	printf("         Actual outputs:             i          Expected outputs:            i\n");
	printf(" -------------------------------   ---   -------------------------------   ---\n");
	for(i = 0; i < RADIX ; i++)
	{
	#ifdef USE_SSE2
		t0 = a[2*i*RE_IM_STRIDE  ];
		t1 = a[2*i*RE_IM_STRIDE+1];
		t2 = a[(2*i+1)*RE_IM_STRIDE  ];
		t3 = a[(2*i+1)*RE_IM_STRIDE+1];
		ASSERT(HERE, t0 == t1 && t2 == t3, "1/2 components of SSE2-pack mismatch!");
	#else
		t0 = a[2*i*RE_IM_STRIDE  ];
		t2 = a[2*i*RE_IM_STRIDE+1];
	#endif
		/* Flip signs on imaginary parts of ref-outputs here: */
		err_r = fabs(t0-b[2*i]);  err_i = fabs(t2+b[2*i+1]);
		abserr = CABS(err_r, err_i);
		if(abserr < 0.00001)
		{
			printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d\n",t0,t2,i,b[2*i  ],-b[2*i+1],i);
		}
		else
		{
			++nerr;
			printf("%15.10f  %15.10f  %4d  %15.10f  %15.10f  %4d, ERR= %15.10e\n",t0,t2,i,b[2*i  ],-b[2*i+1],i, abserr);
		}
	}
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIT transform!");

#endif

#if TEST_TYPE == 3

	nerr = 0;
	avgerr = maxerr = 0.0;
	for(i = 0; i < RADIX ; i++)
	{
	#ifdef USE_SSE2
		ASSERT(HERE, a[2*i*RE_IM_STRIDE] == a[2*i*RE_IM_STRIDE+1] && a[(2*i+1)*RE_IM_STRIDE] == a[(2*i+1)*RE_IM_STRIDE+1], "1/2 components of SSE2-pack mismatch!");
	#endif
		a[2*i*RE_IM_STRIDE] *= iradix;	a[(2*i+1)*RE_IM_STRIDE] *= iradix;
		err_r = fabs(a[2*i*RE_IM_STRIDE]-arrtmp[2*i]);  err_i = fabs(a[(2*i+1)*RE_IM_STRIDE]-arrtmp[2*i+1]);
		abserr = CABS(err_r, err_i);
		avgerr += abserr;
		maxerr = MAX(maxerr, abserr);
		if(abserr < 0.00001)
		{
			printf("%4d  %25.15f  %25.15f\n",i,a[2*i*RE_IM_STRIDE], a[(2*i+1)*RE_IM_STRIDE]);
		}
		else
		{
			++nerr;
			printf("%4d  %25.15f  %25.15f, ERR= %15.10e\n",i,a[2*i*RE_IM_STRIDE], a[(2*i+1)*RE_IM_STRIDE], CABS(err_r, err_i));
		}
	}
	avgerr *= iradix;
	printf("test_fft_radix: %d Mismatches detected in DIF/DIT combo; maxerr = %15.10e, avgerr = %15.10e\n", nerr, maxerr, avgerr);
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIF/DIT combo!");

#endif
	printf("");
}

void matmul_double(double **mat, double vec_in[], double vec_out[], int nrow, int ncol)
{
	int i,j;
	double tmp;

	/*	Loop over [nrow] rows of matrix; inner product of each row with length-[ncol] vector is stored in a second vector of length [nrow]. */
	for(i = 0; i < nrow; i++)
	{
		tmp = 0.0;
		for(j = 0; j < ncol; j++)
		{
			tmp += mat[i][j]*vec_in[j];
		}
		vec_out[i] = tmp;
	}

	return;
}

void matmul_complex(struct complex **mat, struct complex vec_in[], struct complex vec_out[], int nrow, int ncol)
{
	int i,j;
	struct complex tmp;

	/*	Loop over [nrow] rows of matrix; inner product of each row with length-[ncol] vector is stored in a second vector of length [nrow]. */
	for(i = 0; i < nrow; i++)
	{
		vec_out[i].re = 0.0;
		vec_out[i].im = 0.0;
		for(j = 0; j < ncol; j++)
		{
			tmp = cmul(&mat[i][j], &vec_in[j]);
			vec_out[i].re += tmp.re;
			vec_out[i].im += tmp.im;
		/*
			printf("mat[%4d][%4d] = %15.10f  %15.10f\n",i,j, mat[i][j].re, mat[i][j].im);
			printf("cmul * vec_in = %15.10f  %15.10f\n",     vec_in[j].re, vec_in[j].im);
		*/
		}
	/*
		printf("--------------------------------\n");
		printf("vec_out[%4d] = %15.10f  %15.10f\n",i, vec_out[i].re,vec_out[i].im);
	*/
	}

	return;
}

#if 0
Here is the recipe for extracting the output perm from the above data.
NOTE: This assumes that the as-yet-unpermuted DFT outputs are in strictly ASCENDING-INDEX ORDER!

1. Copy the first column of ordered indices into the 2nd ()-bracketed column of indices
(which contains the index of the corresponding matrix-multiply-DFT output for reference);

2. Sort the 5th column (real parts of reference-DFT outputs) in numeric order, copy the
reordered 2nd column which that gives, then undo the sort. [If there are duplicate entries in col5, use col6 instead.
If col6 also has duplicates, do a lexical-sort-of-both-cols-5-and-6, or tweak your test inputs].

3. Sort the 3rd column (real parts of twiddleless-DFT outputs) in numeric order and paste the
reordered 2nd column you copied in step [2] into the resulting 2nd column. [If you used col6 in step 2, use col4 here.]

4. The resulting first 2 columns give the required from -> (to) output-index permutation for the twiddleless DFT.
Doing a final sort of the 1st column in ascending numeric order gives a convenient from/to lookup table.

radix-992 DIT:
   0   0 4585.0000000000 -4441.0000000000  4585.0000000000  4441.0000000000
   1  32 -175.7544491272   -83.0435496863    61.2766665287    -6.3736654598, ERR= 2.5333611695e+002
   2  64  -10.7035054535  -156.2475688113    71.4201695137   -97.2511435757, ERR= 2.6646931375e+002
   3  96   22.4376114922   -30.6040635425   100.2337280653    89.6547447405, ERR= 9.7668923941e+001
   4 128 -173.7103893583   176.6160059971     9.7351717174   -47.4252385631, ERR= 2.2437140698e+002
   5 160  -23.6252423946   124.8961100720   133.6578003338    52.9866766954, ERR= 2.3744523865e+002
   6 192 -221.5572951240    32.5920035235   -33.3816139138   -94.0252457237, ERR= 1.9794981749e+002
   7 224  -48.3990568341  -139.9377954914   -44.9722397106    61.6124332225, ERR= 7.8400289860e+001
   8 256  -11.5508698985  -101.4751476746   254.6244357784    93.9839904670, ERR= 2.6628069924e+002
   9 288   -5.5943597436   -13.7582451884   -33.4234078717   120.4449292356, ERR= 1.1025653936e+002
  10 320   80.4729029338   -22.1056173471   -31.3597970729     9.9740607928, ERR= 1.1248878813e+002
  11 352   45.4021962129   -95.4308787288   -60.9022072948   254.3855763868, ERR= 1.9122557913e+002
  12 384  110.8688320697   -67.7703925102    58.6356822268   -42.8125463845, ERR= 1.2229835779e+002
  13 416  177.0396872829   -13.3104721611    24.6482714243   196.7639848335, ERR= 2.3849179218e+002
  14 448  -70.6233558433    30.8690367099    46.6850711272   -29.3027420892, ERR= 1.1731888304e+002
  15 480 -152.7777427046    13.2622781849    70.5166962618   -53.4313247509, ERR= 2.2687873143e+002
  16 512  -87.9456575299    69.7446768582   -46.6344108680    74.2867139660, ERR= 1.4983878217e+002
  17 544  -63.3580368824    52.5655721701    90.3610859980  -192.9568218908, ERR= 2.0818086304e+002
  18 576 -159.6058876151     6.4738946531   -39.9372076688    27.3072582009, ERR= 1.2434532258e+002
  19 608  -52.1987816951    45.7849749219    41.5744645338   101.5787492886, ERR= 1.7466965656e+002
  20 640  -24.3460797650  -242.7649913886    52.7636212777  -108.0609990176, ERR= 3.5920019702e+002
  21 672   18.1137762912   -91.4202651393   104.7602117360    63.1748749694, ERR= 9.1134004856e+001
  22 704 -122.8456109311    82.9998403020   -77.2766206296   -87.6503708868, ERR= 4.5805679908e+001
  23 736  -38.5817626076    57.1572445255   -21.4067118012    34.8139965373, ERR= 9.3561164768e+001
  24 768   12.3778022375   -36.4916557750    98.3274702004  -168.8863372316, ERR= 2.2263752028e+002
  25 800  -56.5480564048   -78.7064263171  -110.6376061243  -171.8673385720, ERR= 2.5634525749e+002
  26 832  -32.1645223865   -82.1457451111   -84.4477093948    95.4656661980, ERR= 5.3953238471e+001
  27 864   72.2958357021   -20.9603178977   125.7489891039  -136.8840637361, ERR= 1.6664959773e+002
  28 896   84.1896309551    79.2352317856    53.7767140416   -85.6241883128, ERR= 3.1076748232e+001
  29 928  -18.3270693932   -24.2890218817   106.9672773912   102.0808460223, ERR= 1.4747962991e+002
  30 960   62.0194565151   -52.7347150515    19.9500223454   109.2984213506, ERR= 7.0493192307e+001
  31  31   62.7486777972   -66.4792772044    62.7486777972    66.4792772044
  32  63  -62.8438877266   -33.6988613374  -175.7544491272    83.0435496863, ERR= 1.2322212928e+002
  33  95   31.5519949540   163.9230736511   -14.3516643651   -62.8114343707, ERR= 1.1104372804e+002
  34 127   55.7845214542   112.6494619667   123.5512178722   114.0327996558, ERR= 2.3659495531e+002
  35 159  -15.1323694596    52.2566848342   -44.4586291356    29.6194452170, ERR= 8.6969708397e+001
  36 191   -8.1875214779   -31.5030425864    38.0634647854  -165.5296740870, ERR= 2.0238835236e+002
  37 223    0.8063542246    29.0397677353    49.3386067562   -15.5275726128, ERR= 5.0378159482e+001
  38 255   66.9668419310    97.9462217387    55.2495531296    65.1210919024, ERR= 1.6348774766e+002
  39 287  -85.5085313713   222.8801694807   -71.6476171782  -102.9963370061, ERR= 1.2068246862e+002
  40 319  128.6190428738    69.8655601852     1.8012169292   -68.7843781101, ERR= 1.2682243466e+002
  41 351    1.0098922110    -0.9667903071   -99.4508057263   -59.1940359984, ERR= 1.1709686952e+002
  42 383   27.4594938918  -132.7496009596  -150.4886197741  -133.3030866545, ERR= 3.2007743398e+002
  43 415   15.6525259180   199.6745957334   -48.2158384538    32.9630125018, ERR= 2.4124556935e+002
  44 447   83.1551850113  -122.7163304571    15.9527764301    99.0201058313, ERR= 7.1257805051e+001
  45 479   99.8540009867   -19.8724052885   -85.2238993264   134.7234102499, ERR= 2.1781777367e+002
  46 511  -16.5790371298   144.7888789940    77.4600241100   -20.0576160917, ERR= 1.5620894015e+002
  47 543  -21.1548487967    40.2059980251  -140.0162134317   -30.7315338950, ERR= 1.1923837249e+002
  48 575  139.3235123177   -49.7261937956    60.0432569270   -25.4598446636, ERR= 1.0926252456e+002
  49 607  -40.4873659003   -57.6284338338  -102.2090609818    70.4784903464, ERR= 6.3045155215e+001
  50 639    5.3616776642    51.4434309268    63.1710204308   -51.7197605574, ERR= 5.7810003193e+001
  51 671   91.8472831188  -131.7091257820     5.6436129984    33.5911112079, ERR= 1.3060711132e+002
  52 703   24.6247095236    65.2947475128   119.9776842004   214.8496237183, ERR= 2.9592745481e+002
  53 735 -107.2758095029   -62.2206468744   -49.6436635881    36.4320136379, ERR= 6.3138877460e+001
  54 767   48.0861070591   112.1276680037   135.8977024148   -30.9064986028, ERR= 1.1961502681e+002
  55 799   54.7641653963  -100.1553131297   -12.2400564335    62.2427766832, ERR= 7.6986532347e+001
  56 831   88.0407208249   133.5657471935   -18.5462047402   -94.3643629548, ERR= 1.1356725420e+002
  57 863  -46.2880040172   -65.0647586214   -95.2352553945    19.4414133832, ERR= 6.6912801825e+001
  58 895   55.3669560323  -110.4991619699  -165.4428823106    54.1459309757, ERR= 2.2788740938e+002
  59 927   14.4977342719    88.4168461351   -62.1934393657    62.8928500249, ERR= 1.6963537445e+002
  60 959  -94.7258363926   -11.3228204763   -72.3833877859   -30.6565515208, ERR= 4.7554733550e+001
  61 991  -90.0844492193   -34.4067375301   -37.8899251469  -146.4288670941, ERR= 1.8821738560e+002
  62  62  -45.7234336040    55.3409200715   -45.7234336040   -55.3409200715
  63  94    0.3072472892   108.7059975286   -62.8438877266    33.6988613374, ERR= 1.5577936218e+002
  64 126  -39.8534651439    94.0333265291   -10.7035054535   156.2475688113, ERR= 2.5197271027e+002
  65 158  130.6774126281    47.4255311359   -21.2281654747    28.9373064188, ERR= 1.7001937424e+002
  66 190   -8.2294332097   214.0906269146   -22.8911680430  -207.9751327320, ERR= 1.5886023336e+001
  67 222   51.5131463009    66.0982618769   214.1299125746   -75.8587566718, ERR= 1.6290942248e+002
  68 254   15.4795791094     8.1598723665   -14.9440435830    58.9122491562, ERR= 7.3649618487e+001
  69 286   -8.9205164376    90.9099539737   -40.0168582147    89.4404296065, ERR= 1.8301159343e+002
  70 318   43.0750444632    39.8971626494   -52.5932789888   -48.7415380689, ERR= 9.6076277450e+001
  71 350   42.7072193694    47.7601214284    14.7254134511    76.3525042643, ERR= 1.2722784805e+002
  72 382    8.2474112454   242.9778663334  -126.9095914770   -32.0724444885, ERR= 2.5049653161e+002
  73 414   73.9177719526   117.4252857816   -17.1500071430    15.5154828799, ERR= 1.6114151657e+002
  74 446  -32.1170349968   144.2502402357  -209.4105624303   119.6472123666, ERR= 3.1792272703e+002
  75 478   -6.9482002453   -72.3884178899    30.6578530397    -7.8848373255, ERR= 8.8645421464e+001
  76 510  -64.4234914096   -32.5498497100    18.1063867379    85.0879156445, ERR= 9.7833681108e+001
  77 542   10.1769913859  -163.6299785519   179.3410693348    -5.3501651201, ERR= 2.3910410750e+002
  78 574   -4.7405893507    45.7004698574   129.9959419133    65.8934927911, ERR= 1.7494897930e+002
  79 606  -69.0226212116    46.6010355615   -39.1120553458   105.7043863093, ERR= 1.5521463682e+002
  80 638    8.5189309390   -22.1847247178    96.2170850128   -50.0915248980, ERR= 1.1364340054e+002
  81 670  -95.2341253863   -57.2406696108   -68.6903471788  -127.3731858846, ERR= 1.8651232614e+002
  82 702   22.9699495157  -124.4615681833    -2.5518916808   -66.8791114744, ERR= 1.9303528193e+002
  83 734 -157.9457596062   100.5517435045   -24.3629728955   -73.3088577199, ERR= 1.3633244563e+002
  84 766  -43.1912433269   -63.2343071573   -98.2681735618   -94.0070657753, ERR= 1.6660827592e+002
  85 798 -106.5839436879   -44.6263916814   150.5444681530   -32.8299612683, ERR= 2.6854144333e+002
  86 830   -1.6459329191  -190.4025161390   -64.9369282751    91.6561004394, ERR= 1.1728855318e+002
  87 862   16.9471386866   -23.5766518275  -133.4397429767   -48.8048331426, ERR= 1.6689905195e+002
  88 894  -31.7960148099   -15.8856069170   -63.7530942716   -28.9590182072, ERR= 5.5066281246e+001
  89 926   61.4687192615   -47.7608519700    24.7716962374  -157.7268607865, ERR= 2.0873876399e+002
  90 958   19.6035368702   115.2050854185     7.2470788610  -137.7641825799, ERR= 2.5721487501e+001
  91 990   88.7702339249   -71.1158917353    25.9369938808  -148.0214855498, ERR= 2.2796755510e+002
  92  30   19.9500223454  -109.2984213506   -13.8445219401    68.8722193998, ERR= 5.2691071612e+001
  93  93   37.0022487954    -1.8311052237    37.0022487954     1.8311052237
  94 125  128.8136760636     4.8069347096     0.3072472892  -108.7059975286, ERR= 1.6525409977e+002
  95 157  -22.2488153457    -1.3792325073    31.5519949540  -163.9230736511, ERR= 1.7383722159e+002
  96 189  -27.9465783584   -71.8559581868    22.4376114922    30.6040635425, ERR= 6.5117473835e+001
  97 221  113.3199292424  -130.7505793247    63.8032397090   -14.6456690398, ERR= 1.5359678246e+002
  98 253  -16.7682406386    57.9555451515   -93.6489595925   -29.9758875926, ERR= 8.1813850808e+001
  99 285   34.4610243002   -44.4230388319   -69.8820090674   -75.6495720356, ERR= 1.5907514103e+002
 100 317  -24.8913521469   152.9737782522    64.2582359468  -124.8971682233, ERR= 9.3466277812e+001
 101 349   78.4561566198   185.6982557795    91.8760235141   122.3676141406, ERR= 3.0835802736e+002
 102 381 -229.6356238851   -56.9167923040  -196.6185172606  -105.2949927915, ERR= 1.6553788857e+002
 103 413  -32.5275484149  -166.5815897269   -24.5323590815   134.0720384126, ERR= 3.3478261292e+001
 104 445   66.6653265295     8.4026524512   -96.2668019293    44.8758828242, ERR= 1.7142193793e+002
 105 477  -90.4593459645  -147.3482939986    23.5135828294    35.0341028436, ERR= 1.6001345579e+002
 106 509   19.0124066614     1.2752339055    24.2848712248    76.5094619565, ERR= 7.7963182291e+001
 107 541   99.2643842286    54.6728099485    -1.0947296139   -24.5994721286, ERR= 1.0476811241e+002
 108 573  212.7724327188   159.8066752830   -69.3748381741    91.1477539047, ERR= 3.7760456565e+002
 109 605    6.8120130862   177.2204402938   -38.2928447017    74.8902345764, ERR= 2.5611372587e+002
 110 637   21.0682140378  -141.6987558174   -59.7688322050   -34.5226901581, ERR= 1.9387786379e+002
 111 669  174.8856910171    71.2695626488    36.4257775250    75.8071131267, ERR= 2.0199677275e+002
 112 701  -62.2063022701   151.6704176655    38.3171656177   124.8553989247, ERR= 2.9423034316e+002
 113 733  -18.5593042939  -185.3771951645    62.5805546206   -72.8824795943, ERR= 2.7070599608e+002
 114 765  -57.7934504337    33.3979214326    -9.8168847968   -92.7794192924, ERR= 7.6340769831e+001
 115 797  -45.7994691664   152.4830549595    88.4416395398    46.0191053593, ERR= 2.3963259986e+002
 116 829  -14.9988867292   -28.6053233795     2.4146048829  -174.9399871792, ERR= 2.0428882285e+002
 117 861   99.1212359417    -5.0407358560    22.1146828541   -70.6823045562, ERR= 1.0799994476e+002
 118 893    4.3487268956   -18.4763906375   114.7692121387   -33.0102397910, ERR= 1.2183413592e+002
 119 925   45.1852203248   139.6948011352    37.5218172603    76.5603160040, ERR= 2.1639085802e+002
 120 957  224.9993064788  -106.8318345369  -140.9046867460    63.6669613196, ERR= 3.6844122806e+002
 121 989  -65.5009613050     8.4262506129    39.9037799060    52.9236198450, ERR= 1.2195887042e+002
 122  29  106.9672773912  -102.0808460223    29.2706661440   -68.0052093855, ERR= 1.8699205770e+002
 123  61  -37.8899251469   146.4288670941    92.7432056879   -78.6155616885, ERR= 1.4718579844e+002
 124 124  191.4386001800    -2.4436508139   191.4386001800     2.4436508139
 125 156  115.2201454052   -12.1153973624   128.8136760636    -4.8069347096, ERR= 2.1705976101e+001
 126 188   10.1744291979   -79.3947729889   -39.8534651439   -94.0333265291, ERR= 1.8049957317e+002
 127 220   61.2828998653    63.1355603959    55.7845214542  -112.6494619667, ERR= 4.9818255830e+001
 128 252   84.2391905516  -108.2460232069  -173.7103893583  -176.6160059971, ERR= 3.8429723061e+002
 129 284  -68.8297297419   -34.1107099334   114.9982673715     1.2792834409, ERR= 1.8673680700e+002
 130 316   16.5337948499   -72.4974515632    -9.9089010834  -144.4883824085, ERR= 2.1859109843e+002
 131 348  -39.1730368880  -187.8932638024   -52.7177599948  -165.1522989155, ERR= 3.5330529133e+002
 132 380   43.4556007126   -44.4029007369   -20.5027043263   -50.9988151092, ERR= 1.1485709456e+002
 133 412 -210.8156717753   -16.4443762084    45.3504097312    69.8924099911, ERR= 2.6168254361e+002
 134 444   21.5060957968    35.7005759623   -31.0110700831  -124.0945579784, ERR= 1.0281803718e+002
 135 476  -40.4366227275  -131.5416702517   -28.4508097987   233.4959580635, ERR= 1.0265640026e+002
 136 508  -63.9730461227   -98.1364943816    39.2463778867   -44.8072978020, ERR= 1.7631556147e+002
 137 540   26.6116517219   -70.1732114750   162.9830102288   -90.8884805515, ERR= 2.1104031856e+002
 138 572   38.8873850867    54.6246113553   -56.3030862668    49.4826229170, ERR= 1.4106573668e+002
 139 604    3.1081495088   103.0683564912  -118.3278503909   -31.1953820402, ERR= 1.4111139758e+002
 140 636   97.6451862750   -32.5942626700    88.1557230046    37.2019719767, ERR= 1.0548976169e+001
 141 668    9.6527959791   -45.9338900510    -2.4746438748    32.3798895330, ERR= 1.8187515703e+001
 142 700  -60.0025056925    68.5806814543    46.0168875540   -44.4440340358, ERR= 1.0873219161e+002
 143 732   39.7759446097   -98.0248211208   106.6679417532    -0.0248916532, ERR= 1.1869408350e+002
 144 764  -21.5578132934   -33.9473891844   -53.4414051976     7.1651419666, ERR= 4.1639550895e+001
 145 796   34.6111546262    92.8253467318  -116.0751224361    32.8618955504, ERR= 1.9622343634e+002
 146 828  -96.4967721161    37.0666462075    27.5629867528    12.5893205521, ERR= 1.3362836078e+002
 147 860   57.8512675485   220.3998763073    55.0410679207   -19.5155322740, ERR= 2.0090399921e+002
 148 892   62.8144985249    27.1757934001   -31.1691791779   -71.1311736240, ERR= 1.0375455231e+002
 149 924  -94.6832745145    17.2837604583   110.6418347814   -89.8129783809, ERR= 2.1775878389e+002
 150 956  -77.5629540294    23.3171955267    -4.3131570739   149.6471676407, ERR= 1.8783557618e+002
 151 988  210.8987956141   -48.6429163723   -59.6896591910    73.1515479934, ERR= 2.7169612603e+002
 152  28   53.7767140416    85.6241883128  -140.7075024676    18.2575883488, ERR= 2.2048930585e+002
 153  60  -72.3833877859    30.6565515208   -12.8132075493    40.9942606391, ERR= 9.3179639711e+001
 154  92  -13.8445219401   -68.8722193998   -37.4466171491    71.8027214658, ERR= 2.3783329048e+001
 155 155  -36.4068175360   -80.0494397856   -36.4068175360    80.0494397856
 156 187  -11.3864913812    -7.7878482367   115.2201454052    12.1153973624, ERR= 1.2668057531e+002
 157 219   86.1341047453    40.1110141032   -22.2488153457     1.3792325073, ERR= 1.1605299622e+002
 158 251  -59.2042893745   -32.0550463443   130.6774126281   -47.4255311359, ERR= 2.0584514313e+002
 159 283  -73.0732862455   102.3611206421   -15.1323694596   -52.2566848342, ERR= 7.6600289331e+001
 160 315   96.4594532703   139.7396602257   -23.6252423946  -124.8961100720, ERR= 1.2099861617e+002
 161 347   97.7779770881    62.8943357111   -67.1998504328   -70.7799831625, ERR= 1.6516617998e+002
 162 379   50.8680088562   149.0553734156    11.8924343495     6.0568195512, ERR= 1.5993401081e+002
 163 411  170.6858365169    -4.7432363365   112.2704796655   -42.0455250043, ERR= 7.4843450641e+001
 164 443   71.3787448884    66.3111201100     4.0886620368  -113.0424629250, ERR= 8.1925415174e+001
 165 475  158.0833458108  -198.5408689069   -66.2805620128   -22.6635158453, ERR= 3.1507228213e+002
 166 507   90.3160767306    30.3136863772   161.5277766128  -124.5470263796, ERR= 1.1811447231e+002
 167 539  145.6052992258   -92.4832715927   -41.2065621632    14.3008611672, ERR= 2.0251212521e+002
 168 571  112.6003053642   -40.1034872700   113.2796226597   -68.7889302304, ERR= 1.0889453642e+002
 169 603  -61.6100439492    45.5912333239   -31.5374304549   210.1727176654, ERR= 2.5752584474e+002
 170 635  -79.3814277952  -146.4416909013    -0.0903255809    48.3039874686, ERR= 1.2616690424e+002
 171 667 -129.0707260860    41.2204889828   -37.4643362555    19.5371886108, ERR= 1.0992372830e+002
 172 699 -132.9691466068    62.5565314329    -1.1037634711    -5.1717548143, ERR= 1.4381061107e+002
 173 731   47.2833262939    15.0283821753   145.7388056449   -22.3564196161, ERR= 9.8727815467e+001
 174 763  -64.0850880389   194.0961577098   -14.8689426862   -11.3565003760, ERR= 1.8925118580e+002
 175 795  113.0101126586    39.8274926822   -41.2758220429   -28.2217934955, ERR= 1.5472182102e+002
 176 827   85.3616806521   137.4804207739   182.4548776218   -22.2757433040, ERR= 1.5066255875e+002
 177 859  -65.0714731736   148.8769960574    21.7136566903   -34.4783639544, ERR= 1.4359215087e+002
 178 891   92.5491754403   -65.5241203855    63.9772488334    43.0520848426, ERR= 3.6350342109e+001
 179 923   -8.8053533460    11.2879030657    30.4524993631    23.1344776878, ERR= 5.2211869302e+001
 180 955  -75.7092826998    16.0514678970   -10.3885363095   170.1887097935, ERR= 1.9736312648e+002
 181 987 -119.3877783951    60.0481801020    -8.7474703032  -100.5126663791, ERR= 1.1780769255e+002
 182  27  125.7489891039   136.8840637361    23.5491147899    96.5333710554, ERR= 2.5481073991e+002
 183  59  -62.1934393657   -62.8928500249     8.9734003684   -23.5815961775, ERR= 1.1199352179e+002
 184  91   25.9369938808   148.0214855498    87.9169485264  -116.8666409729, ERR= 6.9369583525e+001
 185 123   92.7432056879    78.6155616885    79.9185662219   185.6731554616, ERR= 2.6459969273e+002
 186 186   10.6660487136   214.3896580851    10.6660487136  -214.3896580851
 187 218  -92.8945635678    38.1210877076   -11.3864913812     7.7878482367, ERR= 9.3547828575e+001
 188 250 -119.8530856113    61.9558797237    10.1744291979    79.3947729889, ERR= 1.9206030727e+002
 189 282   21.8594678676  -107.9642522052   -27.9465783584    71.8559581868, ERR= 6.1517892825e+001
 190 314   58.8356473746   -44.2384107721    -8.2294332096  -214.0906269146, ERR= 2.6689251909e+002
 191 346  -69.3246399841   -64.7948426117    -8.1875214779    31.5030425864, ERR= 6.9613872239e+001
 192 378   48.5004819551   111.7172547295  -221.5572951240   -32.5920035235, ERR= 2.8141074667e+002
 193 410  -58.1176728410     0.5857860757   -48.0230007565  -241.1519767177, ERR= 2.4077789451e+002
 194 442  -70.3737558827   -37.5282186020   -14.6093213404   -85.1974432832, ERR= 1.3480081693e+002
 195 474  -76.7923662005   152.5809413882    23.5386308858   108.0500118802, ERR= 2.7927549620e+002
 196 506  133.1965103415    11.0559566058   116.2245334060   -89.9339861721, ERR= 8.0683279243e+001
 197 538   -9.4785717179   -71.8750527460   147.5876769751   -34.7004204429, ERR= 1.8981079517e+002
 198 570  -51.7080287052  -162.6118438085   111.3559855651    60.4301219893, ERR= 1.9243434471e+002
 199 602  -53.0143555036    42.3689102307    34.3731108822    -1.0768186151, ERR= 9.6651984518e+001
 200 634  -61.9316792040   172.9620384483     3.0319836131    54.5513428709, ERR= 2.3660645842e+002
 201 666   95.8660288772   194.0000692438   -40.9976180060    96.0937843605, ERR= 3.2075863470e+002
 202 698   25.0308406595   -64.7459048425     4.7623662932   -38.0179418107, ERR= 1.0474358802e+002
 203 730 -114.4094261505   -17.3455257390    24.2057843403   157.7855661774, ERR= 1.9732607921e+002
 204 762   40.3630677954    52.5563242375   -93.2363388606    -8.4691762416, ERR= 1.4068574227e+002
 205 794 -200.3409524506   -11.4237929294    62.7583671543   -45.1107873413, ERR= 2.6910483225e+002
 206 826 -103.1518468764   -22.1495834756   -58.5473123383  -228.7998422935, ERR= 2.5488267653e+002
 207 858  109.6300508979  -105.2940753317     7.1999588612    -9.4875938401, ERR= 1.5384003163e+002
 208 890   -4.1547077844   -36.9284545673    47.6360818164    24.0454461700, ERR= 5.3369071501e+001
 209 922  -33.5110793709   119.1626953969   153.4689767660    -0.6040438950, ERR= 2.2139940207e+002
 210 954  -59.2197321198   -60.4932902414  -121.3323607850   -52.4862790019, ERR= 1.2892773831e+002
 211 986   33.5762451482    52.5103750896   -85.7423411071    10.0955357990, ERR= 1.3474577954e+002
 212  26  -84.4477093948   -95.4656661980   -44.4454467820     1.2038880736, ERR= 1.0239855384e+002
 213  58 -165.4428823106   -54.1459309757    42.0771782221    26.5908646919, ERR= 2.0934148466e+002
 214  90    7.2470788610   137.7641825799    32.6965380690   -42.5081543144, ERR= 9.8597088674e+001
 215 122   29.2706661441    68.0052093855   -64.7594438588  -120.6823916247, ERR= 1.0778008682e+002
 216 154  -37.4466171491   -71.8027214658   -51.7325431673  -147.3305667636, ERR= 2.1959846469e+002
 217 217  -33.0198794986   -47.4311781386   -33.0198794985    47.4311781386
 218 249  -87.7895434713    33.2879999868   -92.8945635678   -38.1210877076, ERR= 7.0299336486e+000
 219 281  -21.7209194754   -95.9772905433    86.1341047453   -40.1110141032, ERR= 1.7364542295e+002
 220 313   36.9981175002    69.9223967026    61.2828998653   -63.1355603959, ERR= 2.5215308873e+001
 221 345  -79.4216376067    10.0314614709   113.3199292424   130.7505793247, ERR= 2.3868157575e+002
 222 377  117.0801505384   -78.3723594813    51.5131463008   -66.0982618769, ERR= 1.5865305695e+002
 223 409  -71.4231633150   -75.4114948798     0.8063542246   -29.0397677353, ERR= 1.2699279297e+002
 224 441  -52.7456905165    10.6129676860   -48.3990568341   139.9377954914, ERR= 1.5061349713e+002
 225 473  -65.0634549015   -59.5882957492    77.8883251143    21.2286242051, ERR= 1.4800903962e+002
 226 505   70.2730857046   -12.4132681536  -141.0482470272    -7.6127070611, ERR= 2.1226809782e+002
 227 537   12.9053758228  -114.5392707240  -118.7822927677     4.8130758447, ERR= 1.7141026778e+002
 228 569  -20.2485592225   -48.8757173194   -86.4109287920   -90.2110824494, ERR= 1.5402141740e+002
 229 601  -60.4017415115    40.7011431004   -32.3307867556    45.8057618158, ERR= 9.0947364443e+001
 230 633   -4.7666471712    45.9738369609    51.9713495007    46.7285521962, ERR= 1.0868731859e+002
 231 665  -14.5267510465  -121.3713924184   118.8540632655   106.3859631797, ERR= 1.3421998628e+002
 232 697  101.2940014546    65.0174512796  -230.0155290585   -79.6817255505, ERR= 3.3163390350e+002
 233 729  -68.6393654328   106.1006217125   -17.6909507710   176.6491389271, ERR= 2.8730326851e+002
 234 761   49.7056182453    24.0720192020   138.7551988875   -59.5495157096, ERR= 9.5856562483e+001
 235 793  -40.2313370517   -19.5215154981  -154.4820563748  -116.1006028429, ERR= 1.7733185233e+002
 236 825    6.1637495413   180.9373893612    69.2768835726   -38.4974717841, ERR= 1.5579601345e+002
 237 857   29.2666165354   -27.8242394537   -38.0602209681  -113.1514593024, ERR= 1.5622756059e+002
 238 889    0.0349074474    58.4851139068    15.2708904525  -136.2781356213, ERR= 7.9270987162e+001
 239 921   34.3498132862  -104.3487008951    99.5322578266   152.1840808675, ERR= 8.0851559375e+001
 240 953   14.1416583255    89.3255579892   -98.1358247058    -0.3001537040, ERR= 1.4328906380e+002
 241 985   22.0702149093    49.5442109582   -62.4597414604   -74.1666794594, ERR= 8.8043054689e+001
 242  25 -110.6376061243   171.8673385720    64.1308501926   -22.3199992776, ERR= 2.3001830365e+002
 243  57  -95.2352553945   -19.4414133832    44.0069866772   -36.5778069871, ERR= 1.5008849066e+002
 244  89   24.7716962374   157.7268607865    64.8561585645   -64.5325908071, ERR= 1.0144917978e+002
 245 121   39.9037799060   -52.9236198450  -143.6751990846    97.6338039062, ERR= 1.8894507690e+002
 246 153  -12.8132075493   -40.9942606391   -27.8377415037   -98.5745014149, ERR= 1.4037512587e+002
 247 185   79.9185662219  -185.6731554616   -49.6951002391  -162.9571603937, ERR= 3.7194461909e+002
 248 248 -101.0000000000  -147.0000000000  -101.0000000000   147.0000000000
 249 280  -54.4889439311    64.9459494389   -87.7895434713   -33.2879999868, ERR= 4.5947314320e+001
 250 312   -1.8146747530    98.7217758566  -119.8530856113   -61.9558797237, ERR= 1.2363170126e+002
 251 344   12.3439267356    94.7694523075   -59.2042893745    32.0550463443, ERR= 1.4561456207e+002
 252 376  -85.0634515569    56.3010888798    84.2391905516   108.2460232069, ERR= 2.3609137367e+002
 253 408  -13.9822662487    67.8442893614   -16.7682406386   -57.9555451515, ERR= 1.0273700178e+001
 254 440   20.6405670123    57.5392151584    15.4795791094    -8.1598723665, ERR= 4.9648316091e+001
 255 472  -94.6121182847   235.3489379391    66.9668419310   -97.9462217387, ERR= 2.1210201980e+002
 256 504  -55.5776924628   -34.8701945914   -11.5508698985   101.4751476746, ERR= 7.9840972441e+001
 257 536  -33.4302617907    60.5474342366   -72.2818730544   -82.6332267736, ERR= 4.4690378492e+001
 258 568  -56.8769516712    95.6807838185   -37.2972703930   155.0196777193, ERR= 2.5146388475e+002
 259 600  -33.8577692186  -105.5380548794   -56.4518512214   -83.3769464062, ERR= 1.9026132096e+002
 260 632  -93.5012996665   -58.4310900568   -38.7384701899   -24.6516877280, ERR= 9.9507363832e+001
 261 664   99.6028071174    27.5903205162   120.1408139276   -76.1220114177, ERR= 5.2698526977e+001
 262 696   61.0985383730   -30.7054482555   -21.0174084811   -58.5402515142, ERR= 1.2127581645e+002
 263 728   72.9592451509   -52.1707803413   155.8535416893   -78.3446984827, ERR= 1.5461485896e+002
 264 760  -67.7632481919    48.7120504859    30.3933692743   -56.3058641008, ERR= 9.8449924112e+001
 265 792   21.8984041628   -78.1456793130  -182.1866695066    65.2003129711, ERR= 2.0449523174e+002
 266 824  -59.0031436233    64.8601543539   191.9789867313    48.4629686901, ERR= 2.7538002828e+002
 267 856   18.8411536527    64.7331996290  -123.7097546392  -181.1062720416, ERR= 1.8402025279e+002
 268 888 -162.9650911576   -12.7415230818    72.9544975019     4.7688182434, ERR= 2.3605426566e+002
 269 920 -190.6283354803   -75.1370405484   -21.2761106528    58.5768332713, ERR= 1.7015997332e+002
 270 952  -46.8148933838   -87.1408882317   -12.3272397812    -9.2418279344, ERR= 1.0236711497e+002
 271 984  205.3951118824    37.8036708560    39.6717496511     5.2044716261, ERR= 1.7121312189e+002
 272  24   98.3274702004   168.8863372316    78.0865197071   -36.5951199570, ERR= 1.3383072235e+002
 273  56  -18.5462047402    94.3643629548    64.0082565387    12.4798638359, ERR= 1.3502195331e+002
 274  88  -63.7530942716    28.9590182072  -158.1277654214  -114.4209713527, ERR= 1.2731977062e+002
 275 120 -140.9046867460   -63.6669613196   107.3096636197   -93.0721068824, ERR= 2.9356004365e+002
 276 152 -140.7075024676   -18.2575883488    51.6353018912   124.6721888676, ERR= 2.1981770082e+002
 277 184   87.9169485264   116.8666409729     0.5540210372  -165.9706119657, ERR= 1.0021716952e+002
 278 216  -51.7325431673   147.3305667636   132.9106242781  -155.4167575364, ERR= 1.8482014437e+002
 279 279   12.3175369900   -78.1695392596    12.3175369900    78.1695392596
 280 311  101.6728098162   169.4273193432   -54.4889439311   -64.9459494389, ERR= 1.8789052661e+002
 281 343   78.2071257714    -2.8736731948   -21.7209194754    95.9772905433, ERR= 1.3657927292e+002
 282 375 -149.3277652403   161.6133605062    21.8594678676   107.9642522052, ERR= 3.1933862600e+002
 283 407   54.8812177312   -88.7260410597   -73.0732862455  -102.3611206421, ERR= 2.2997099481e+002
 284 439  -68.0429690321    -4.0783401975   -68.8297297419    34.1107099334, ERR= 3.0042673389e+001
 285 471  113.3150489569   173.4500546910    34.4610243002    44.4230388319, ERR= 2.3170378091e+002
 286 503   70.7185696008  -116.8991678346    -8.9205164376   -90.9099539737, ERR= 2.2254665833e+002
 287 535   52.1934258151   -17.9783820461   -85.5085313713  -222.8801694807, ERR= 2.7744309481e+002
 288 567   76.9306594948  -100.3648146631    -5.5943597436    13.7582451884, ERR= 1.1962891238e+002
 289 599 -121.5869897448   -88.2973667110   195.2825062950    56.4350043064, ERR= 3.1846740439e+002
 290 631    6.6354455823  -142.8578867304    85.3836468318   180.6613680339, ERR= 8.7352060071e+001
 291 663  103.9819663194    46.9373147011   -26.0572522312  -169.5310845415, ERR= 1.7871606185e+002
 292 695 -145.0804177658   -52.3371384397   105.5315534295  -127.6254719988, ERR= 3.0853346863e+002
 293 727  -83.4080182166  -150.1606952099   120.2520286753    60.8628593160, ERR= 2.2237697317e+002
 294 759  -53.2843234979    63.9520984695    29.4028481839  -104.9400010425, ERR= 9.2288550309e+001
 295 791  -68.2626686455   -11.1404889912   -77.8088455207   -99.7316074319, ERR= 1.1128230434e+002
 296 823  -34.6367661386   116.0651805398   -78.0897554554     7.5571075871, ERR= 1.3103675974e+002
 297 855 -135.9875143956    59.6500207493    53.0552405970   -72.7157472514, ERR= 1.8949373716e+002
 298 887   -5.3344019695    25.9738763351    83.2006174995   136.5171406413, ERR= 1.8504534652e+002
 299 919 -116.2484630952    84.5679207424    57.4346748619   -64.5876972577, ERR= 1.7482860676e+002
 300 951   -3.6012723596  -177.2717568937   -36.9109182955    -8.4523909362, ERR= 1.8868755020e+002
 301 983   58.7182794053  -168.1402197748   -75.4105703373    25.2734145043, ERR= 1.9596293624e+002
 302  23  -21.4067118012   -34.8139965373     4.2778969398   103.0285403241, ERR= 7.2889801140e+001
 303  55  -12.2400564335   -62.2427766832   106.1352505970  -124.7295377291, ERR= 2.2129473485e+002
 304  87 -133.4397429767    48.8048331426   -24.5373989062     9.2760485837, ERR= 1.2342248323e+002
 305 119   37.5218172603   -76.5603160040    30.4454002466     2.8023239143, ERR= 7.4096673845e+001
 306 151  -59.6896591910   -73.1515479934   -59.2696466743   -12.1969145195, ERR= 8.5349495979e+001
 307 183    8.9734003684    23.5815961775    87.6696979786   -19.0519901579, ERR= 7.8826547484e+001
 308 215  -64.7594438588   120.6823916247  -137.8021464804  -119.6394898668, ERR= 7.3050147504e+001
 309 247  -49.6951002391   162.9571603937   -84.0008207986    39.2327456559, ERR= 2.0507959570e+002
 310 310  -33.2611675057    54.6547470504   -33.2611675057   -54.6547470504
 311 342   68.6320957491   -38.0469395353   101.6728098162  -169.4273193432, ERR= 2.1008868814e+002
 312 374  -95.9260022695  -110.8607055820    -1.8146747530   -98.7217758566, ERR= 2.2974280945e+002
 313 406  -49.2577327894    86.8944148346    36.9981175002   -69.9223967026, ERR= 8.7909732730e+001
 314 438 -196.1328593744  -111.6652604342    58.8356473746    44.2384107721, ERR= 2.6373342505e+002
 315 470  103.5977636895   -91.6784937170    96.4594532703  -139.7396602257, ERR= 2.3152822171e+002
 316 502  -44.8454624344    29.6872566756    16.5337948499    72.4974515632, ERR= 1.1920204622e+002
 317 534  -93.1511863092    84.3196272733   -24.8913521469  -152.9737782522, ERR= 9.6813208843e+001
 318 566 -100.5806031729    12.3140696743    43.0750444632   -39.8971626494, ERR= 1.4627977343e+002
 319 598  -43.5742668663  -105.7625053624   128.6190428738   -69.8655601852, ERR= 2.4595884478e+002
 320 630 -114.4852501364   -79.2975697814    80.4729029338    22.1056173471, ERR= 2.0317381936e+002
 321 662   34.8660322862   -76.3710121138   -61.4761750087   -57.1611516884, ERR= 1.6465922287e+002
 322 694 -166.6816663189   -33.1355693522   -25.0402321749   112.6033207930, ERR= 1.6241126619e+002
 323 726  126.2321649472    11.5678673223   -80.6969757474   -52.3581302743, ERR= 2.1091115385e+002
 324 758  -85.1518931401    44.4591443131    -7.4531366408    16.8301493520, ERR= 9.8961984012e+001
 325 790 -113.8014385883    76.6935203523   -35.4389436091    99.6742301539, ERR= 1.9299291189e+002
 326 822  -27.4880795508    68.2372215474    42.3223432012   -92.0646270673, ERR= 7.3764763801e+001
 327 854  138.4848151749  -259.6573729437   -66.7630184058    28.7012598493, ERR= 3.0897799172e+002
 328 886 -211.8473998895   -46.6213564113   176.1818640942    61.3950991335, ERR= 3.8831040828e+002
 329 918  -51.8846138870   -63.5599100501   -73.5226968936  -120.7769306668, ERR= 1.8560247165e+002
 330 950   31.7591610290    46.1288445618   -24.3644992345   -60.6248137711, ERR= 5.7965492879e+001
 331 982   12.6264095494   -92.5280968119    -9.3605256118  -126.5305325050, ERR= 2.2015927960e+002
 332  22  -77.2766206296    87.6503708868   -94.9576871155   191.5352668958, ERR= 2.7974495609e+002
 333  54  135.8977024148    30.9064986028   -46.3271771162   -32.8012035946, ERR= 1.8223472948e+002
 334  86  -64.9369282751   -91.6561004394    19.6165238857   -17.1691522966, ERR= 1.3781227052e+002
 335 118  114.7692121387    33.0102397909   -79.0320454643    74.2789409704, ERR= 2.2151725837e+002
 336 150   -4.3131570739  -149.6471676407   -32.5070300009   -70.5127942141, ERR= 2.2195788626e+002
 337 182   23.5491147899   -96.5333710554   212.3227188907    78.7602558625, ERR= 1.8960843132e+002
 338 214   32.6965380690    42.5081543144   -35.7328333893   -52.5753973102, ERR= 6.9165947255e+001
 339 246  -27.8377415037    98.5745014149    -8.6324910088    26.3041538271, ERR= 1.2634682498e+002
 340 278  132.9106242781   155.4167575364   -75.8705986219    88.4347153133, ERR= 3.2101890886e+002
 341 341   -8.8258149745    35.2692910936    -8.8258149745   -35.2692910936
 342 373 -118.0838221775   -24.8674198674    68.6320957490    38.0469395353, ERR= 1.8718048441e+002
 343 405   -1.2509008019    -1.4004680284    78.2071257714     2.8736731948, ERR= 7.9471682506e+001
 344 437   14.0197197129   -47.5204398180    12.3439267356   -94.7694523075, ERR= 1.4229975995e+002
 345 469   92.4763133811    12.9422495624   -79.4216376067   -10.0314614709, ERR= 1.7192259375e+002
 346 501  -39.8806453645    39.3147448146   -69.3246399841    64.7948426117, ERR= 1.0819313755e+002
 347 533  120.4088387271    50.5356646671    97.7779770881   -62.8943357111, ERR= 2.5785512376e+001
 348 565  -78.0477216790  -185.6493875568   -39.1730368879   187.8932638024, ERR= 3.8939390060e+001
 349 597  -94.3025611897   -46.8417064536    78.4561566198  -185.6982557795, ERR= 2.8969019420e+002
 350 629  -36.5987071734    40.9660916665    42.7072193694   -47.7601214284, ERR= 7.9596412138e+001
 351 661   -9.2049854476  -107.3117324843     1.0098922110     0.9667903071, ERR= 1.0683440669e+002
 352 693  -65.3010118993  -146.8362030118    45.4021962129    95.4308787288, ERR= 1.2205616597e+002
 353 725   47.3921351909    -5.4744469718  -130.4266467311   -14.0225271928, ERR= 1.7888446329e+002
 354 757   42.8860168411  -221.2572528413    22.5083459188    37.3549091128, ERR= 1.8502789385e+002
 355 789   87.2791358539   -84.0230432699   -43.8463522062    30.3286246779, ERR= 1.4169327509e+002
 356 821  -51.4180868691    -3.0962044360   -78.3313334497    18.4126501263, ERR= 3.0966374507e+001
 357 853   -7.0870283665   119.2844002803   116.0296514701    59.9738702403, ERR= 2.1746550164e+002
 358 885    9.8449291512    55.8274542960    90.0512348089    73.5750120123, ERR= 1.5224338985e+002
 359 917  199.7365791657  -119.9072713535  -154.9605287449   211.3651556154, ERR= 3.6629848888e+002
 360 949   39.0642988736   -39.0991769115  -211.0829281403   123.2480280445, ERR= 2.6392170113e+002
 361 981  -72.1507395839   -93.6970217112   -92.1710797602    31.3192386034, ERR= 6.5511845083e+001
 362  21  104.7602117360   -63.1748749694   -26.2084693311   -32.4290211052, ERR= 1.6215085681e+002
 363  53  -49.6436635881   -36.4320136379  -163.5928784224   108.9487571498, ERR= 1.3506702651e+002
 364  85  150.5444681530    32.8299612683    51.9045752982   -24.5060973359, ERR= 9.8990480215e+001
 365 117   22.1146828541    70.6823045562    37.3180605370    -1.8605619989, ERR= 7.0481025401e+001
 366 149  110.6418347814    89.8129783809  -257.3838954820     1.6448131339, ERR= 3.7921954824e+002
 367 181   -8.7474703032   100.5126663791   -73.7579175120  -139.5522708037, ERR= 7.5831714737e+001
 368 213   42.0771782221   -26.5908646919    99.7194538368    28.2615538817, ERR= 5.7666481949e+001
 369 245 -143.6751990846   -97.6338039062   -55.9929830686   -47.4796214863, ERR= 1.6954668158e+002
 370 277    0.5540210372   165.9706119657    62.5511382455     4.1277219553, ERR= 1.8104443031e+002
 371 309  -84.0008207986   -39.2327456559   -34.9201071742   -60.4573421152, ERR= 1.1111719061e+002
 372 372  110.2842712475  -199.1126983722   110.2842712475   199.1126983722
 373 404  -47.9795950389    16.3139685521  -118.0838221775    24.8674198674, ERR= 8.1305039296e+001
 374 436  -52.0154327052    15.8485339717   -95.9260022695   110.8607055821, ERR= 1.3410208614e+002
 375 468  -62.3580155299  -115.5587917318  -149.3277652403  -161.6133605062, ERR= 2.9049636717e+002
 376 500 -129.3495365724   -16.2562196811   -85.0634515569   -56.3010888798, ERR= 8.5004825461e+001
 377 532  115.6773104023   -30.9244786585   117.0801505384    78.3723594813, ERR= 4.7468614421e+001
 378 564  109.9545941854   104.5381450919    48.5004819551  -111.7172547295, ERR= 6.1872025385e+001
 379 596  -15.4844937497    87.9543809582    50.8680088562  -149.0553734156, ERR= 9.0199700007e+001
 380 628   34.2360857936   161.3007713211    43.4556007126    44.4029007369, ERR= 2.0591017496e+002
 381 660   29.7538984457   -78.3217361299  -229.6356238851    56.9167923040, ERR= 2.6027119686e+002
 382 692  -72.4031984725   -28.2062928411     8.2474112454  -242.9778663334, ERR= 2.8292290299e+002
 383 724  150.3826454508   -78.8708603607    27.4594938918   132.7496009597, ERR= 1.3421259210e+002
 384 756   78.2524987757    84.7382943355   110.8688320697    67.7703925102, ERR= 1.5595744535e+002
 385 788   39.2468462200    55.2796825733    45.6716098939   226.9708981243, ERR= 2.8232369347e+002
 386 820  111.8927293123   -62.4389414084  -125.5097276676    11.7896964277, ERR= 2.4274528337e+002
 387 852   30.6947523285  -139.2346569361    66.5593079571   -88.9637986360, ERR= 2.3099957029e+002
 388 884 -134.4481019717   155.3166676613   117.5989824804   174.7488045845, ERR= 4.1529621808e+002
 389 916  -23.6878912726    71.7923261932   -19.7467221195   -11.9734199995, ERR= 5.9948597586e+001
 390 948   29.0062589615    58.8390070410   119.6104008403   -70.8905699738, ERR= 9.1402137255e+001
 391 980 -190.0401565606    54.6345135812   -86.5864577592    24.9306567037, ERR= 1.3051162446e+002
 392  20   52.7636212777   108.0609990176    61.9160284882    45.5467202595, ERR= 1.5388014160e+002
 393  52  119.9776842004  -214.8496237183    37.9214566651   -21.8716606215, ERR= 2.5053979911e+002
 394  84  -98.2681735618    94.0070657753   -22.1965850362    86.5289905510, ERR= 1.9590853533e+002
 395 116    2.4146048829   174.9399871792    51.7814435785   -66.7574155122, ERR= 1.1891405962e+002
 396 148  -31.1691791779    71.1311736240    30.9017430549    93.7490619187, ERR= 1.7617687550e+002
 397 180  -10.3885363095  -170.1887097935   -34.3646559848  -138.4138515934, ERR= 3.0953254305e+002
 398 212  -44.4454467820    -1.2038880736   168.0813773782    66.2920990714, ERR= 2.2227038984e+002
 399 244   64.8561585645    64.5325908071    75.8909078873  -261.8071087897, ERR= 1.9758289688e+002
 400 276   51.6353018912  -124.6721888676     2.2204695699    11.7754032527, ERR= 1.2323761542e+002
 401 308 -137.8021464804   119.6394898668   111.9883074242    37.0979219847, ERR= 2.9489300964e+002
 402 340  -75.8705986219   -88.4347153133  -153.8898569800   -42.6920241631, ERR= 1.5258186812e+002
 403 403    2.7097118752   -27.1447129921     2.7097118752    27.1447129921
 404 435 -226.6472008924    26.0820515135   -47.9795950389   -16.3139685521, ERR= 1.7893442605e+002
 405 467  -35.8192779546    97.8555263321    -1.2509008019     1.4004680284, ERR= 1.0510340202e+002
 406 499 -141.8493938577  -121.1367924361   -49.2577327894   -86.8944148346, ERR= 2.2770638748e+002
 407 531   80.4233209674   100.6389804265    54.8812177311    88.7260410597, ERR= 1.9107985347e+002
 408 563   66.0944735727   242.9096816132   -13.9822662487   -67.8442893613, ERR= 1.9251019668e+002
 409 595 -170.3458967297    10.8730802681   -71.4231633150    75.4114948798, ERR= 1.3126589464e+002
 410 627  174.1555850729    55.8029265043   -58.1176728410    -0.5857860756, ERR= 2.3874630665e+002
 411 659   64.3761437808   -34.5641642595   170.6858365169     4.7432363365, ERR= 1.1041303597e+002
 412 691   46.5261737243    22.4742806859  -210.8156717753    16.4443762084, ERR= 2.6026810657e+002
 413 723 -192.8421891065    71.9680433649   -32.5275484148   166.5815897269, ERR= 2.8741383312e+002
 414 755  107.5100349982   209.6346338488    73.9177719527  -117.4252857816, ERR= 9.8137678837e+001
 415 787   20.2454662176   -50.2782467495    15.6525259180  -199.6745957334, ERR= 2.4999503708e+002
 416 819 -101.1021420573    16.2685230825   177.0396872829    13.3104721611, ERR= 2.7971019679e+002
 417 851  -17.2461082482   -21.3125973294    14.0920160499    -1.3402017912, ERR= 3.8668169630e+001
 418 883  -32.6394473977   125.8076257732   -12.9318381412   -74.5564997193, ERR= 5.4909632893e+001
 419 915   54.3621135571   -94.5973873375   195.4674255197  -122.3772047245, ERR= 2.5882171985e+002
 420 947    1.9734867389    53.3842750262   -39.7896726581   243.7524783423, ERR= 3.0005734733e+002
 421 979   52.1539841648    37.8420824864    68.5919618926   181.8645181981, ERR= 2.2032066970e+002
 422  19   41.5744645338  -101.5787492886     8.4822501085   -77.7070101947, ERR= 1.8231422931e+002
 423  51    5.6436129984   -33.5911112079    24.1865450800   -40.2362748672, ERR= 7.6120452343e+001
 424  83  -24.3629728955    73.3088577199    61.5886938335   -66.8978203113, ERR= 8.6190431106e+001
 425 115   88.4416395398   -46.0191053593    97.9233990846   -84.9265555870, ERR= 1.3128849868e+002
 426 147   55.0410679207    19.5155322740  -167.6942282449   -60.7766760710, ERR= 2.2652482015e+002
 427 179   30.4524993631   -23.1344776878   -90.4533621824   -56.9940414523, ERR= 1.4504760231e+002
 428 211  -85.7423411072   -10.0955357991   -13.3930700417   -26.7952917372, ERR= 8.1211761340e+001
 429 243   44.0069866772    36.5778069871    36.0523084615   -17.7115943332, ERR= 2.0474640056e+001
 430 275  107.3096636197    93.0721068824   -38.5532737923   -11.7492980597, ERR= 1.6700118486e+002
 431 307   87.6696979786    19.0519901579   174.2149633615    50.7344018107, ERR= 1.1117654188e+002
 432 339   -8.6324910088   -26.3041538271   -78.3391252275     8.3575138504, ERR= 7.1979835652e+001
 433 371  -34.9201071742    60.4573421152   -70.9213404148     9.8578553256, ERR= 7.8995669413e+001
 434 434 -131.2295535871    80.2413807831  -131.2295535871   -80.2413807831
 435 466   36.2301589277    54.1211216919  -226.6472008924   -26.0820515135, ERR= 2.6436848481e+002
 436 498  -66.8936926857    10.4147341475   -52.0154327052   -15.8485339717, ERR= 1.5839469706e+001
 437 530  -35.1044804136   -24.4280796180    14.0197197129    47.5204398180, ERR= 5.4281158220e+001
 438 562 -100.9502871188   -71.7383604341  -196.1328593744   111.6652604342, ERR= 1.0321763127e+002
 439 594  -83.2705732123    -2.7259022057   -68.0429690321     4.0783401974, ERR= 1.5287544525e+001
 440 626  -31.9818065906    32.7498727015    20.6405670123   -57.5392151584, ERR= 5.8168941051e+001
 441 658  -11.7616851566   244.9165118792   -52.7456905165   -10.6129676860, ERR= 2.3786096678e+002
 442 690   52.5702202415   -55.1165912629   -70.3737558826    37.5282186020, ERR= 1.2419570088e+002
 443 722 -310.4385884374  -115.1000960037    71.3787448884   -66.3111201100, ERR= 4.2272272870e+002
 444 754 -105.9400508984     7.5302362529    21.5060957968   -35.7005759623, ERR= 1.3052236723e+002
 445 786 -108.9048894430    34.0156679550    66.6653265295    -8.4026524512, ERR= 1.7742865411e+002
 446 818  -65.7877990550  -136.0842831518   -32.1170349968  -144.2502402357, ERR= 2.8234936755e+002
 447 850   64.6387108854    23.6662567171    83.1551850113   122.7163304571, ERR= 1.4754904826e+002
 448 882  -44.9093634520    46.1353865186   -70.6233558433   -30.8690367099, ERR= 2.9904361575e+001
 449 914   57.3331592550   -99.9913869519   -84.0084616716   -15.5575942947, ERR= 1.8256237529e+002
 450 946  -44.6481798072  -204.7242705858  -127.6371998554   -48.1309546961, ERR= 2.6612580183e+002
 451 978  -21.9112008775    83.6464895876   111.3548605336   -71.4532118416, ERR= 1.3382271536e+002
 452  18  -39.9372076688   -27.3072582009    37.6523804561  -152.1617662614, ERR= 1.9552308029e+002
 453  50   63.1710204308    51.7197605574   111.2554690320   -75.0485051139, ERR= 5.3444780099e+001
 454  82   -2.5518916808    66.8791114744   182.9541453158  -150.7270288336, ERR= 2.0357544795e+002
 455 114   -9.8168847968    92.7794192924   -24.0810521189  -142.4141442460, ERR= 5.1643706205e+001
 456 146   27.5629867528   -12.5893205521    85.9255335574    35.2296287370, ERR= 6.2600083261e+001
 457 178   63.9772488334   -43.0520848426    -2.5505896441   162.9292053826, ERR= 1.3710024552e+002
 458 210 -121.3323607849    52.4862790019   120.1018326843    37.9617378331, ERR= 2.5782031248e+002
 459 242   64.1308501926    22.3199992776   -75.6301063350    40.5726689851, ERR= 1.5325995136e+002
 460 274 -158.1277654214   114.4209713527   -59.7929635031   181.2761948984, ERR= 3.1161923464e+002
 461 306  -59.2696466743    12.1969145195   206.0330976296    25.9535854173, ERR= 2.6803172719e+002
 462 338  -35.7328333893    52.5753973103   118.4287596945   109.9784273140, ERR= 2.2403022716e+002
 463 370   62.5511382455    -4.1277219553  -201.2746741105   -25.9899093180, ERR= 2.6553932097e+002
 464 402 -153.8898569800    42.6920241631    28.4793867760     8.9057286725, ERR= 1.8952801684e+002
 465 465  -18.4551119962   -30.1871091818   -18.4551119962    30.1871091818
 466 497 -164.4638281657    20.5811524184    36.2301589277   -54.1211216919, ERR= 2.0347728619e+002
 467 529  -94.5931467258    42.9671071629   -35.8192779546   -97.8555263321, ERR= 8.0418320109e+001
 468 561 -161.3960569361   118.3091991672   -62.3580155299   115.5587917318, ERR= 2.5397395696e+002
 469 593   15.3431197086   -94.6224068512    92.4763133811   -12.9422495624, ERR= 1.3236194648e+002
 470 625   68.2760671923  -168.9505914184   103.5977636895    91.6784937170, ERR= 8.4962340638e+001
 471 657   11.8052011661   121.6321331776   113.3150489569  -173.4500546909, ERR= 1.1397081288e+002
 472 689   39.6901890526    68.6801027007   -94.6121182848  -235.3489379391, ERR= 2.1404581378e+002
 473 721   48.7457348862   -74.9674494479   -65.0634549015    59.5882957492, ERR= 1.1484358950e+002
 474 753   19.4829854161   -35.9326207469   -76.7923662005  -152.5809413882, ERR= 2.1167500192e+002
 475 785   63.9226312488    -1.6829746964   158.0833458108   198.5408689069, ERR= 2.1821840133e+002
 476 817  -95.6958798394   126.8669693024   -40.4366227275   131.5416702517, ERR= 2.6425103688e+002
 477 849  -22.6239330287   -75.0422485516   -90.4593459645   147.3482939986, ERR= 9.9145385452e+001
 478 881 -138.7025005785    -5.2914696650    -6.9482002453    72.3884178899, ERR= 1.4785532157e+002
 479 913 -195.3131174961    11.0707948697    99.8540009867    19.8724052885, ERR= 2.9678461798e+002
 480 945   78.1938240255   -33.2966605488  -152.7777427046   -13.2622781849, ERR= 2.3561748537e+002
 481 977 -186.8448839835   -44.1554164722   -64.3275106574   139.7895968451, ERR= 1.5542330334e+002
 482  17   90.3610859980   192.9568218908   145.7183608596    97.0239136208, ERR= 2.9521730106e+002
 483  49 -102.2090609818   -70.4784903464   -85.0112773150   -23.5912161093, ERR= 9.5628831613e+001
 484  81  -68.6903471788   127.3731858846   -65.0012730816    84.8831942242, ERR= 2.1228843625e+002
 485 113   62.5805546206    72.8824795943  -137.5611615386   113.3709515545, ERR= 2.7339906211e+002
 486 145 -116.0751224361   -32.8618955504   -52.8657921780   -79.4495089363, ERR= 1.2887696074e+002
 487 177   21.7136566904    34.4783639544   -48.7377849716   174.7439733244, ERR= 2.2076546842e+002
 488 209  153.4689767660     0.6040438950   -68.7572330438    33.9119531089, ERR= 2.2489073430e+002
 489 241  -62.4597414604    74.1666794594     9.6798186706   -28.8216289633, ERR= 8.5207333842e+001
 490 273   64.0082565387   -12.4798638360   -57.7442262240   -14.2871857981, ERR= 1.2466010591e+002
 491 305   30.4454002466    -2.8023239143     2.8185902882   -92.1410285897, ERR= 9.8881144882e+001
 492 337  212.3227188907   -78.7602558625   -51.3944076454    68.6457233369, ERR= 2.6391102023e+002
 493 369  -55.9929830686    47.4796214863   115.4417899842    23.1676642116, ERR= 1.8542092759e+002
 494 401  111.9883074241   -37.0979219847   -14.8260786480  -217.6831647248, ERR= 2.8459671583e+002
 495 433  -70.9213404147    -9.8578553256   -39.2545502178  -144.6690535764, ERR= 1.5773823625e+002
 496 496  -61.0000000000    45.0000000000   -61.0000000000   -45.0000000000
 497 528  110.6299979494   -64.4326704517  -164.4638281657   -20.5811524184, ERR= 2.8793048336e+002
 498 560 -128.8563100312   145.5805538130   -66.8936926857   -10.4147341475, ERR= 1.4869150868e+002
 499 592 -114.4290670912   -29.2751396898  -141.8493938577   121.1367924361, ERR= 9.5866769870e+001
 500 624   66.1329047676    56.3933095165  -129.3495365724    16.2562196811, ERR= 2.0854577187e+002
 501 656  -57.5187439210    12.2699787461   -39.8806453645   -39.3147448146, ERR= 3.2288107600e+001
 502 688  167.3725240026   -38.9224412908   -44.8454624344   -29.6872566756, ERR= 2.2303310163e+002
 503 720 -125.3969899817   -18.1242033165    70.7185696008   116.8991678346, ERR= 2.1958553305e+002
 504 752   -8.5293997106   -61.6807383154   -55.5776924628    34.8701945914, ERR= 5.4151150548e+001
 505 784  -90.6729941639    60.9925799875    70.2730857046    12.4132681536, ERR= 1.7689561658e+002
 506 816  162.4122574508  -109.5126853604   133.1965103415   -11.0559566058, ERR= 1.2405787885e+002
 507 848  178.5820179915    83.2909101386    90.3160767306   -30.3136863771, ERR= 1.0294397809e+002
 508 880  -16.3148994801    55.7569685321   -63.9730461227    98.1364943816, ERR= 1.6110399396e+002
 509 912   49.9117815849   -35.7087342864    19.0124066614    -1.2752339056, ERR= 4.8193207757e+001
 510 944   78.5021375315    33.4385946500   -64.4234914096    32.5498497100, ERR= 1.5742366467e+002
 511 976 -138.4043188780    18.6274858548   -16.5790371298  -144.7888789940, ERR= 1.7537986313e+002
 512  16  -46.6344108680   -74.2867139660   -87.9456575299   -69.7446768582, ERR= 1.4983878217e+002
 513  48   60.0432569270    25.4598446636   187.3601307838    13.4031030645, ERR= 1.3311617135e+002
 514  80   96.2170850128    50.0915248980   130.2300428861   -71.4578660955, ERR= 4.0167173655e+001
 515 112   38.3171656177  -124.8553989247    48.7534245618   120.0622863913, ERR= 1.1484312278e+001
 516 144  -53.4414051976    -7.1651419666   -73.5293580120   -15.7890512611, ERR= 3.0502800445e+001
 517 176  182.4548776218    22.2757433040    73.2079639954   -78.4735525144, ERR= 1.2285390469e+002
 518 208   47.6360818164   -24.0454461700    69.6772866021     7.1474987504, ERR= 2.7773284563e+001
 519 240  -98.1358247058     0.3001537040   -85.0384869639   -40.8592698268, ERR= 4.2621381449e+001
 520 272   78.0865197071    36.5951199570     3.3669685084    24.6196574089, ERR= 9.6593272536e+001
 521 304  -24.5373989062    -9.2760485837   -32.7167936948   -63.0706971160, ERR= 7.2807651469e+001
 522 336  -32.5070300009    70.5127942141    54.1407065519    30.9809405343, ERR= 1.3344964759e+002
 523 368   99.7194538368   -28.2615538816   -40.0095004782   -15.3911385965, ERR= 1.4638899629e+002
 524 400    2.2204695699   -11.7754032527   -41.3236946627   -31.2362528742, ERR= 6.1205365790e+001
 525 432  -78.3391252275    -8.3575138504   -92.6268545952   -18.6350569749, ERR= 3.0540761127e+001
 526 464   28.4793867760    -8.9057286725    -8.3832179095    17.3961488241, ERR= 3.7827752491e+001
 527 527   68.4045442736    -2.0947947786    68.4045442736     2.0947947786
 528 559   14.0992544009    36.5735020698   110.6299979494    64.4326704517, ERR= 1.3971553721e+002
 529 591  -85.4957947251    56.0346808802   -94.5931467257   -42.9671071629, ERR= 1.5922414901e+001
 530 623  -61.8823980669    13.8763176663   -35.1044804136    24.4280796180, ERR= 4.6736321263e+001
 531 655 -101.4869837281   -11.4494554486    80.4233209674  -100.6389804265, ERR= 2.1367071959e+002
 532 687    6.1792688152    43.3215932599   115.6773104023    30.9244786585, ERR= 1.3229625961e+002
 533 719   -6.8222622430   -18.3260005855   120.4088387271   -50.5356646671, ERR= 1.4467094385e+002
 534 751  -27.8347393769   -35.7734510198   -93.1511863092   -84.3196272732, ERR= 1.3670620210e+002
 535 783  -76.5511452800    -5.7742505551    52.1934258151    17.9783820461, ERR= 1.2932171284e+002
 536 815  109.0387447618  -149.6481612637   -33.4302617907   -60.5474342366, ERR= 2.5392834855e+002
 537 847  -54.3806333455   -93.6391970664    12.9053758228   114.5392707240, ERR= 7.0457221835e+001
 538 879  182.5689356946   112.2600903473    -9.4785717179    71.8750527460, ERR= 2.6606013611e+002
 539 911   32.5352984457   214.5901593161   145.6052992258    92.4832715927, ERR= 3.2722915067e+002
 540 943  -53.9135947867   165.4822970907    26.6116517219    70.1732114750, ERR= 2.4903380100e+002
 541 975  -12.4885235845   -64.4906762335    99.2643842286   -54.6728099485, ERR= 1.6336660872e+002
 542  15   70.5166962618    53.4313247509    10.1769913859   163.6299785519, ERR= 2.2529200913e+002
 543  47 -140.0162134317    30.7315338950   -21.1548487967   -40.2059980252, ERR= 1.1923837249e+002
 544  79  -39.1120553458  -105.7043863093   -63.3580368824   -52.5655721701, ERR= 1.6011635575e+002
 545 111   36.4257775250   -75.8071131267   -36.5332342538    82.6188955315, ERR= 7.3276311174e+001
 546 143  106.6679417532     0.0248916532    -2.5718496017    72.9524681928, ERR= 1.3137361632e+002
 547 175  -41.2758220429    28.2217934955     7.4893734772   190.9041276345, ERR= 2.2448655551e+002
 548 207    7.1999588612     9.4875938401   -58.2753402406   151.1887455621, ERR= 1.7350475739e+002
 549 239   99.5322578266  -152.1840808675     6.7054265826   -84.6483055376, ERR= 2.5437452673e+002
 550 271   39.6717496511    -5.2044716261    20.6632197226   -13.6285480774, ERR= 2.6758304154e+001
 551 303  106.1352505970   124.7295377290    -3.0470134365   100.6656068147, ERR= 2.5044707617e+002
 552 335  -79.0320454643   -74.2789409704   -65.5102850901   -42.5114312766, ERR= 1.1757052800e+002
 553 367  -73.7579175120   139.5522708037   133.5456255789   177.0354765286, ERR= 3.7842114203e+002
 554 399   75.8909078873   261.8071087897    24.2119259873   -67.2806626945, ERR= 2.0127408030e+002
 555 431  174.2149633615   -50.7344018107   131.6864678551    85.0737486038, ERR= 5.4661354431e+001
 556 463 -201.2746741104    25.9899093180   -61.4505803857     7.1622369540, ERR= 1.4370052884e+002
 557 495  -39.2545502178   144.6690535764  -172.1291158796   -56.8505066930, ERR= 1.5927255688e+002
 558 558  -39.3060036475   -10.5124929467   -39.3060036475    10.5124929467
 559 590   75.2183029415    78.8939664866    14.0992544009   -36.5735020698, ERR= 7.4340835366e+001
 560 622  179.9455168424  -132.5289706784  -128.8563100312  -145.5805538130, ERR= 4.1557607714e+002
 561 654  -60.9669941647   130.6186502529  -161.3960569361  -118.3091991672, ERR= 1.0118062678e+002
 562 686   79.2474562525   -85.1694337539  -100.9502871188    71.7383604341, ERR= 1.8069759391e+002
 563 718   35.4999943784   157.7559448119    66.0944735727  -242.9096816132, ERR= 9.0483042878e+001
 564 750   47.0603237121    79.7661987334   109.9545941854  -104.5381450919, ERR= 6.7596882951e+001
 565 782   49.4096660408   -80.7792753135   -78.0477216790   185.6493875568, ERR= 1.6505491852e+002
 566 814  100.1393216656    84.1941976425  -100.5806031728   -12.3140696743, ERR= 2.1320234760e+002
 567 846  -31.7916143581   103.9072466153    76.9306594948   100.3648146631, ERR= 2.3140356058e+002
 568 878 -103.1751852076    39.7950945376   -56.8769516712   -95.6807838185, ERR= 7.2572285998e+001
 569 910    2.8137707232   -82.1285739159   -20.2485592225    48.8757173194, ERR= 4.0467561507e+001
 570 942   78.8931948279   -12.7206001722   -51.7080287052   162.6118438085, ERR= 1.9880660076e+002
 571 974  108.4487037408    63.3682182560   112.6003053642    40.1034872700, ERR= 1.0355495952e+002
 572  14   46.6850711272    29.3027420892    38.8873850867   -54.6246113554, ERR= 2.6495300918e+001
 573  46   77.4600241100    20.0576160918   212.7724327188  -159.8066752830, ERR= 1.9452312836e+002
 574  78  129.9959419133   -65.8934927911    -4.7405893507   -45.7004698574, ERR= 1.7494897930e+002
 575 110  -59.7688322050    34.5226901581   139.3235123177    49.7261937956, ERR= 2.1618426422e+002
 576 142   46.0168875540    44.4440340358  -159.6058876151    -6.4738946531, ERR= 2.0909915627e+002
 577 174  -14.8689426862    11.3565003760   -36.2344795589   -69.9457797698, ERR= 6.2363369262e+001
 578 206  -58.5473123383   228.7998422935  -100.7999988120    28.6295990354, ERR= 2.6087392890e+002
 579 238   15.2708904525   136.2781356213  -148.8902365930   141.6286709234, ERR= 3.2277092303e+002
 580 270  -12.3272397811     9.2418279344   138.8623878963    93.2869125777, ERR= 1.8267579519e+002
 581 302    4.2778969398  -103.0285403241    85.9417446380    12.1738122280, ERR= 1.2216204664e+002
 582 334   19.6165238857    17.1691522966    65.2270668339    89.9728444958, ERR= 1.1644624985e+002
 583 366 -257.3838954820    -1.6448131339   104.9029998474   105.3301377961, ERR= 3.7683211259e+002
 584 398  168.0813773782   -66.2920990714   -12.2348868969    82.6075515967, ERR= 1.8105289049e+002
 585 430  -38.5532737923    11.7492980597   105.1409361908    46.1734178366, ERR= 1.5492923223e+002
 586 462  118.4287596945  -109.9784273140    34.4811205655   166.1091077712, ERR= 1.0098445130e+002
 587 494  -14.8260786480   217.6831647248    48.1669080377  -120.1815524993, ERR= 1.1608049258e+002
 588 526   -8.3832179095   -17.3961488241    32.0863406463    36.7411364302, ERR= 4.4855475866e+001
 589 589  -54.1824085191    87.7972921430   -54.1824085192   -87.7972921430
 590 621  -24.8894963985   -16.5962251213    75.2183029415   -78.8939664866, ERR= 1.3834720157e+002
 591 653  102.4681164645    50.8874139590   -85.4957947251   -56.0346808802, ERR= 1.8803437523e+002
 592 685   89.3630046569   103.3916111706  -114.4290670912    29.2751396898, ERR= 2.4317005426e+002
 593 717  130.9466350838    63.5729235145    15.3431197086    94.6224068512, ERR= 1.9593349718e+002
 594 749 -133.4455264491    55.4607781249   -83.2705732123     2.7259022057, ERR= 7.6832387053e+001
 595 781   51.1043443841   -21.5967406713  -170.3458967297   -10.8730802680, ERR= 2.2381800321e+002
 596 813   36.6944614161  -111.7402926271   -15.4844937497   -87.9543809581, ERR= 2.0639914249e+002
 597 845   96.9781803813  -100.3163610977   -94.3025611897    46.8417064537, ERR= 1.9861485540e+002
 598 877  102.9038997649    83.7774934509   -43.5742668663   105.7625053624, ERR= 2.3954386749e+002
 599 909  -75.4363727332   -26.3775863324  -121.5869897448    88.2973667110, ERR= 7.7226541116e+001
 600 941   15.8696102086  -162.1231387255   -33.8577692186   105.5380548794, ERR= 7.5330498329e+001
 601 973  106.7355575654   -22.8974067795   -60.4017415115   -40.7011431004, ERR= 1.7882855558e+002
 602  13   24.6482714243  -196.7639848335   -53.0143555036   -42.3689102307, ERR= 2.5142797204e+002
 603  45  -85.2238993264  -134.7234102499   -61.6100439492   -45.5912333238, ERR= 1.8185429567e+002
 604  77  179.3410693348     5.3501651201     3.1081495088  -103.0683564912, ERR= 2.0151150576e+002
 605 109  -38.2928447017   -74.8902345764     6.8120130862  -177.2204402938, ERR= 2.5611372587e+002
 606 141   -2.4746438748   -32.3798895330   -69.0226212116   -46.6010355615, ERR= 1.0327932909e+002
 607 173  145.7388056449    22.3564196160   -40.4873659003    57.6284338338, ERR= 2.0267650024e+002
 608 205   62.7583671543    45.1107873413   -52.1987816951   -45.7849749219, ERR= 1.1495912578e+002
 609 237  -38.0602209681   113.1514593024  -102.4301065825     5.6427747292, ERR= 1.3511310896e+002
 610 269  -21.2761106528   -58.5768332713   119.1480736188   115.6689396802, ERR= 1.5158647744e+002
 611 301  -75.4105703373   -25.2734145043   -31.6784857462    79.2848038497, ERR= 6.9496225809e+001
 612 333  -46.3271771162    32.8012035946    24.7272782856  -120.2862470233, ERR= 1.1270478453e+002
 613 365   37.3180605370     1.8605619989   -52.8431024141    65.4961763973, ERR= 1.1254317177e+002
 614 397  -34.3646559848   138.4138515934   -60.2284365657   114.0239369107, ERR= 2.5375928005e+002
 615 429   36.0523084615    17.7115943332    63.0591160334    24.1669111201, ERR= 4.9831484768e+001
 616 461  206.0330976296   -25.9535854173   -34.8875741990    57.2082813269, ERR= 2.4293955242e+002
 617 493  115.4417899842   -23.1676642116    13.0950627294     2.7909976848, ERR= 1.0435545562e+002
 618 525  -92.6268545952    18.6350569749    36.7159023654    29.5195380412, ERR= 1.3801599110e+002
 619 557 -172.1291158796    56.8505066930   -93.8856129166   124.6714281497, ERR= 1.9766703970e+002
 620 620   24.5613998200   -33.5563491861    24.5613998200    33.5563491861
 621 652   15.3977875939    11.4345931107   -24.8894963985    16.5962251213, ERR= 4.9079446026e+001
 622 684 -135.9082500177   -32.1557911501   179.9455168424   132.5289706783, ERR= 3.3141873394e+002
 623 716  -90.8172462088   -86.1655582178   -61.8823980669   -13.8763176663, ERR= 1.0414222183e+002
 624 748 -128.7075623777    85.4229705474    66.1329047676   -56.3933095165, ERR= 1.9699118980e+002
 625 780  -24.1196327521    14.2567974556    68.2760671923   168.9505914184, ERR= 2.0518750621e+002
 626 812   94.6308582947  -189.3488476200   -31.9818065906   -32.7498727015, ERR= 2.5565329741e+002
 627 844   63.3439552773    40.8413880006   174.1555850728   -55.8029265043, ERR= 1.1181710483e+002
 628 876 -118.5189513120    68.2195845888    34.2360857936  -161.3007713211, ERR= 1.7888043125e+002
 629 908   12.5393551285   132.3563346580   -36.5987071734   -40.9660916665, ERR= 1.0376283381e+002
 630 940  -50.1088001830   -88.6439945509  -114.4852501364    79.2975697814, ERR= 6.5051387107e+001
 631 972   40.9457769509   -64.3096463109     6.6354455823   142.8578867304, ERR= 8.5714788173e+001
 632  12   58.6356822268    42.8125463845   -93.5012996665    58.4310900568, ERR= 1.8274554763e+002
 633  44   15.9527764301   -99.0201058313    -4.7666471712   -45.9738369610, ERR= 1.4646684936e+002
 634  76   18.1063867379   -85.0879156445   -61.9316792040  -172.9620384483, ERR= 2.7017748020e+002
 635 108  -69.3748381741   -91.1477539047   -79.3814277952   146.4416909013, ERR= 5.6192092900e+001
 636 140   88.1557230046   -37.2019719767    97.6451862749    32.5942626700, ERR= 1.0548976169e+001
 637 172   -1.1037634711     5.1717548143    21.0682140378   141.6987558174, ERR= 1.4853465414e+002
 638 204  -93.2363388606     8.4691762416     8.5189309390    22.1847247178, ERR= 1.0627227567e+002
 639 236   69.2768835726    38.4974717840     5.3616776642   -51.4434309268, ERR= 6.5213122947e+001
 640 268   72.9544975019    -4.7688182434   -24.3460797651   242.7649913886, ERR= 2.5711783440e+002
 641 300  -36.9109182956     8.4523909363    14.0717851313   -34.4801362382, ERR= 5.7242288338e+001
 642 332  -94.9576871156  -191.5352668957   -19.6675668206    19.1490604767, ERR= 1.8811062271e+002
 643 364   51.9045752981    24.5060973359  -206.6774951947    10.4177699433, ERR= 2.6092980605e+002
 644 396   30.9017430549   -93.7490619187   -65.7820855185   -78.1879124615, ERR= 1.9725639626e+002
 645 428  -13.3930700417    26.7952917372  -127.8824124384   100.8017687191, ERR= 1.7143167548e+002
 646 460  -59.7929635031  -181.2761948984    56.7758101827    58.5801965439, ERR= 1.6924120955e+002
 647 492  -51.3944076454   -68.6457233369   -50.0101874401  -151.9602370469, ERR= 2.2061030307e+002
 648 524  -41.3236946627    31.2362528742   123.2030658406    91.0450357280, ERR= 2.0499211805e+002
 649 556  -61.4505803857    -7.1622369539   -31.6543388962    97.3574492462, ERR= 9.4989432714e+001
 650 588   32.0863406463   -36.7411364302   -58.7362756745   199.4659790319, ERR= 1.8635482831e+002
 651 651   -1.4454608424   -62.9013973194    -1.4454608424    62.9013973194
 652 683  -19.1025388762   159.0453997564    15.3977875938   -11.4345931107, ERR= 1.5158899289e+002
 653 715  -29.1980749698    35.2498033982   102.4681164644   -50.8874139590, ERR= 1.3259155641e+002
 654 747  -13.7034296404   -24.0875362415   -60.9669941647  -130.6186502529, ERR= 1.6176479429e+002
 655 779   24.0972394464   109.0840066940  -101.4869837281    11.4494554486, ERR= 1.7406812634e+002
 656 811  103.6020691419   136.7319150679   -57.5187439210   -12.2699787461, ERR= 2.0359442526e+002
 657 843  -96.4401762399   141.6630024392    11.8052011661  -121.6321331776, ERR= 1.1008313882e+002
 658 875   77.6347746261    69.1680418373   -11.7616851566  -244.9165118792, ERR= 1.9717822330e+002
 659 907   89.0382740535    74.8198416071    64.3761437808    34.5641642595, ERR= 1.1212975256e+002
 660 939  -21.1396962853    85.5707975192    29.7538984457    78.3217361299, ERR= 1.7161270515e+002
 661 971  -52.3337677875   108.5182404232    -9.2049854476   107.3117324842, ERR= 2.2009695380e+002
 662  11  -60.9022072948  -254.3855763868    34.8660322862    76.3710121138, ERR= 2.0214039875e+002
 663  43  -48.2158384538   -32.9630125018   103.9819663194   -46.9373147011, ERR= 1.7189599782e+002
 664  75   30.6578530397     7.8848373255    99.6028071174   -27.5903205162, ERR= 7.1705737292e+001
 665 107   -1.0947296139    24.5994721286   -14.5267510465   121.3713924184, ERR= 1.4658755915e+002
 666 139 -118.3278503909    31.1953820402    95.8660288772  -194.0000692438, ERR= 2.6904346134e+002
 667 171  -37.4643362555   -19.5371886108  -129.0707260860   -41.2204889828, ERR= 1.0992372830e+002
 668 203   24.2057843403  -157.7855661774     9.6527959791    45.9338900510, ERR= 1.1279444544e+002
 669 235 -154.4820563748   116.1006028429   174.8856910171   -71.2695626488, ERR= 3.3240477612e+002
 670 267 -123.7097546392   181.1062720416   -95.2341253863    57.2406696108, ERR= 2.4004192562e+002
 671 299   57.4346748619    64.5876972577    91.8472831188   131.7091257820, ERR= 1.9929041709e+002
 672 331   -9.3605256118   126.5305325050    18.1137762912    91.4202651393, ERR= 2.1967564148e+002
 673 363 -163.5928784225  -108.9487571498   -80.7445213898   -32.1130062488, ERR= 1.6359178267e+002
 674 395   51.7814435785    66.7574155122    19.4336537231     7.8034027898, ERR= 8.1275427617e+001
 675 427  -90.4533621824    56.9940414523    70.3344387521  -225.8683103111, ERR= 2.3317640449e+002
 676 459  -75.6301063350   -40.5726689851   -79.1017539281   -81.7413508048, ERR= 1.2236327788e+002
 677 491    2.8185902882    92.1410285896   -16.3465860226   -50.2907194852, ERR= 4.6029907182e+001
 678 523  -40.0095004782    15.3911385965   -60.1098282069    60.4926471241, ERR= 7.8500777767e+001
 679 555  131.6864678551   -85.0737486039    24.0992496516   -69.1484981582, ERR= 1.8804124791e+002
 680 587   48.1669080377   120.1815524993    47.3424854408    14.1734178432, ERR= 1.3435749971e+002
 681 619  -93.8856129166  -124.6714281496    91.9794962647  -180.1271013670, ERR= 3.5699857480e+002
 682 682   63.1329411597   -52.5074070912    63.1329411597    52.5074070912
 683 714  177.4867039096   215.6687621854   -19.1025388762  -159.0453997564, ERR= 2.0458136658e+002
 684 746   78.1239315149  -132.0951253664  -135.9082500177    32.1557911501, ERR= 2.3621525195e+002
 685 778 -174.7637913531    -6.0543561335    89.3630046568  -103.3916111706, ERR= 2.8590450177e+002
 686 810   18.1881904740    31.1536076457    79.2474562525    85.1694337539, ERR= 1.3137459381e+002
 687 842   19.1641429604   -29.1917316499     6.1792688152   -43.3215932599, ERR= 7.3666744506e+001
 688 874   91.6750568682    46.4511511533   167.3725240026    38.9224412908, ERR= 1.1409976695e+002
 689 906  -69.9121676613   -55.7051778793    39.6901890525   -68.6801027008, ERR= 1.6578412054e+002
 690 938   44.6767990100   -66.9397406492    52.5702202415    55.1165912628, ERR= 1.4215940354e+001
 691 970   21.2847386873    53.6454399165    46.5261737243   -22.4742806860, ERR= 4.0109490280e+001
 692  10  -31.3597970729    -9.9740607928   -72.4031984725    28.2062928412, ERR= 4.4910745751e+001
 693  42 -150.4886197741   133.3030866545   -65.3010118993   146.8362030118, ERR= 2.9280531100e+002
 694  74 -209.4105624303  -119.6472123666  -166.6816663189    33.1355693522, ERR= 9.6488460139e+001
 695 106   24.2848712248   -76.5094619565  -145.0804177658    52.3371384397, ERR= 1.7108156633e+002
 696 138  -56.3030862668   -49.4826229170    61.0985383730    30.7054482554, ERR= 1.1889374986e+002
 697 170   -0.0903255809   -48.3039874686   101.2940014546   -65.0174512796, ERR= 1.5205436609e+002
 698 202    4.7623662932    38.0179418107    25.0308406595    64.7459048425, ERR= 1.0474358802e+002
 699 234  138.7551988875    59.5495157095  -132.9691466068   -62.5565314329, ERR= 2.7174098343e+002
 700 266  191.9789867313   -48.4629686901   -60.0025056925   -68.5806814543, ERR= 2.7783788180e+002
 701 298   83.2006174994  -136.5171406413   -62.2063022701  -151.6704176655, ERR= 3.2279287644e+002
 702 330  -24.3644992345    60.6248137711    22.9699495158   124.4615681834, ERR= 1.9104323810e+002
 703 362  -26.2084693311    32.4290211052    24.6247095236   -65.2947475128, ERR= 6.0532371875e+001
 704 394  -22.1965850362   -86.5289905510  -122.8456109311   -82.9998403020, ERR= 1.9715539786e+002
 705 426 -167.6942282449    60.7766760711    -9.7716154515     6.5823029156, ERR= 1.7168804175e+002
 706 458  120.1018326843   -37.9617378331   144.5183306635   -90.3924741532, ERR= 1.3065591876e+002
 707 490  -57.7442262240    14.2871857981    83.5610322290   -53.5019207154, ERR= 1.4664573468e+002
 708 522   54.1407065519   -30.9809405342   -24.2761615434   -81.6283622821, ERR= 1.3722266680e+002
 709 554   24.2119259873    67.2806626945   -75.3220827142   -50.7216504057, ERR= 1.0090203058e+002
 710 586   34.4811205655  -166.1091077712  -127.9497715281   -91.1016745809, ERR= 3.0420582056e+002
 711 618   36.7159023654   -29.5195380413    40.0762866267    27.0319451072, ERR= 4.1809449877e+000
 712 650  -58.7362756746  -199.4659790319   -49.0424522429    45.4622881274, ERR= 1.5430848008e+002
 713 713    2.3255876705   233.8411813954     2.3255876705  -233.8411813954
 714 745  -71.8376018135   111.3699805344   177.4867039096  -215.6687621854, ERR= 2.7026069873e+002
 715 777    0.6204457819   -40.5595660040   -29.1980749698   -35.2498033982, ERR= 8.1462903637e+001
 716 809   -1.7222557026   -17.2852587749   -90.8172462088    86.1655582178, ERR= 1.1261621990e+002
 717 841    6.6764861559   -54.4568263885   130.9466350837   -63.5729235146, ERR= 1.7138871543e+002
 718 873  -74.5629122658    38.2351344668    35.4999943784  -157.7559448119, ERR= 1.6247789857e+002
 719 905   83.0096889320   -51.1974310462    -6.8222622430    18.3260005855, ERR= 9.5657254782e+001
 720 937  -30.4764481307    21.8247685729  -125.3969899817    18.1242033165, ERR= 1.0298460866e+002
 721 969  -40.2559619407   -47.0178378923    48.7457348863    74.9674494479, ERR= 9.3287098916e+001
 722   9  -33.4234078717  -120.4449292356  -310.4385884374   115.1000960037, ERR= 2.7706673836e+002
 723  41  -99.4508057263    59.1940359984  -192.8421891064   -71.9680433649, ERR= 9.4260945008e+001
 724  73  -17.1500071430   -15.5154828799   150.3826454508    78.8708603606, ERR= 1.7911195812e+002
 725 105   23.5135828294   -35.0341028436    47.3921351909     5.4744469718, ERR= 3.7999454182e+001
 726 137  162.9830102288    90.8884805515   126.2321649472   -11.5678673223, ERR= 8.7420731591e+001
 727 169  -31.5374304549  -210.1727176654   -83.4080182165   150.1606952099, ERR= 7.9322132560e+001
 728 201  -40.9976180060   -96.0937843605    72.9592451509    52.1707803412, ERR= 1.2212860821e+002
 729 233  -17.6909507710  -176.6491389272   -68.6393654328  -106.1006217125, ERR= 2.8730326851e+002
 730 265 -182.1866695066   -65.2003129711  -114.4094261505    17.3455257391, ERR= 8.2968881986e+001
 731 297   53.0552405970    72.7157472514    47.2833262940   -15.0283821753, ERR= 5.7975400681e+001
 732 329  -73.5226968936   120.7769306668    39.7759446097    98.0248211209, ERR= 2.4639559402e+002
 733 361  -92.1710797602   -31.3192386034   -18.5593042939   185.3771951645, ERR= 1.7074117098e+002
 734 393   37.9214566651    21.8716606215  -157.9457596063  -100.5517435045, ERR= 2.1107942072e+002
 735 425   97.9233990846    84.9265555870  -107.2758095029    62.2206468744, ERR= 2.5250547399e+002
 736 457   -2.5505896441  -162.9292053826   -38.5817626076   -57.1572445255, ERR= 2.2301634662e+002
 737 489    9.6798186706    28.8216289633   -75.3083345899   -34.1406107461, ERR= 8.5154434775e+001
 738 521  -32.7167936948    63.0706971160   -40.9374327872  -162.6745394865, ERR= 9.9942505082e+001
 739 553  133.5456255789  -177.0354765287   -88.1474952141   -11.3029596307, ERR= 2.9089380595e+002
 740 585  105.1409361908   -46.1734178366   -10.7313355642   -19.8115890350, ERR= 1.3334318315e+002
 741 617   13.0950627294    -2.7909976848   -38.9758545708  -189.5185478438, ERR= 1.9923438892e+002
 742 649  -31.6543388962   -97.3574492462    -9.4712311637    90.5823895260, ERR= 2.3194648152e+001
 743 681   91.9794962648   180.1271013669   162.0205710583    36.5414829659, ERR= 2.2770820713e+002
 744 744  -83.0000000000   131.0000000000   -83.0000000000  -131.0000000000
 745 776  155.1137916720   -72.1492237628   -71.8376018135  -111.3699805344, ERR= 2.9186680755e+002
 746 808 -147.8363825331   -59.7288483239    78.1239315149   132.0951253664, ERR= 2.3726555076e+002
 747 840    7.8107946387  -195.8521254404   -13.7034296404    24.0875362416, ERR= 1.7310671838e+002
 748 872   20.7338997283   137.2104069575  -128.7075623777   -85.4229705473, ERR= 1.5816032741e+002
 749 904   37.0055254155   133.2823239064  -133.4455264491   -55.4607781248, ERR= 1.8737596983e+002
 750 936 -191.1581025531   -56.4782275171    47.0603237121   -79.7661987334, ERR= 2.7442769958e+002
 751 968   48.2258652606   -33.7182126258   -27.8347393769    35.7734510198, ERR= 7.6088366934e+001
 752   8  254.6244357784   -93.9839904670    -8.5293997105    61.6807383153, ERR= 2.6512910295e+002
 753  40    1.8012169292    68.7843781101    19.4829854161    35.9326207469, ERR= 1.0619931632e+002
 754  72 -126.9095914770    32.0724444885  -105.9400508984    -7.5302362529, ERR= 3.2280669404e+001
 755 104  -96.2668019293   -44.8758828242   107.5100349981  -209.6346338488, ERR= 3.2603773151e+002
 756 136   39.2463778867    44.8072978020    78.2524987757   -84.7382943355, ERR= 5.5820802135e+001
 757 168  113.2796226597    68.7889302304    42.8860168411   221.2572528414, ERR= 2.9846615898e+002
 758 200    3.0319836131   -54.5513428709   -85.1518931402   -44.4591443132, ERR= 1.3258760384e+002
 759 232 -230.0155290585    79.6817255505   -53.2843234979   -63.9520984695, ERR= 1.7742981764e+002
 760 264   30.3933692743    56.3058641008   -67.7632481918   -48.7120504859, ERR= 9.8449924112e+001
 761 296  -78.0897554554    -7.5571075871    49.7056182451   -24.0720192019, ERR= 1.3165127877e+002
 762 328  176.1818640942   -61.3950991335    40.3630677954   -52.5563242375, ERR= 1.7728979755e+002
 763 360 -211.0829281403  -123.2480280445   -64.0850880390  -194.0961577097, ERR= 3.4973661122e+002
 764 392   61.9160284882   -45.5467202595   -21.5578132934    33.9473891844, ERR= 8.4275896573e+001
 765 424   61.5886938336    66.8978203113   -57.7934504337   -33.3979214326, ERR= 1.2399330464e+002
 766 456   85.9255335573   -35.2296287370   -43.1912433269    63.2343071573, ERR= 1.3211890132e+002
 767 488  -68.7572330439   -33.9119531089    48.0861070592  -112.1276680037, ERR= 1.8702924119e+002
 768 520    3.3669685084   -24.6196574089    12.3778022376    36.4916557750, ERR= 1.4904343988e+001
 769 552  -65.5102850901    42.5114312766   164.4566926230    -5.1182315889, ERR= 2.3298725764e+002
 770 584  -12.2348868969   -82.6075515967   -51.3893275419    40.7650153161, ERR= 5.7305043972e+001
 771 616  -34.8875741990   -57.2082813269   143.4736406552    93.5875658127, ERR= 1.8203344556e+002
 772 648  123.2030658406   -91.0450357280    80.9913763918  -135.0682946327, ERR= 2.3001970544e+002
 773 680   47.3424854408   -14.1734178432  -154.1955070984    -5.0684887719, ERR= 2.0245447243e+002
 774 712  -49.0424522429   -45.4622881274    -4.7912827890    16.6064408934, ERR= 5.2828268168e+001
 775 775 -203.7062570730    92.1649674799  -203.7062570730   -92.1649674799
 776 807    7.9896222270   -32.0813005645   155.1137916720    72.1492237627, ERR= 1.5248265378e+002
 777 839  -45.3980005659    25.9373926662     0.6204457819    40.5595660040, ERR= 8.0867440399e+001
 778 871   27.9645753855    57.3620664398  -174.7637913530     6.0543561335, ERR= 2.1241570877e+002
 779 903    1.4143783056    23.0261437725    24.0972394464  -109.0840066940, ERR= 8.8997010962e+001
 780 935   76.9817150357    57.5948733348   -24.1196327521   -14.2567974556, ERR= 1.0999850611e+002
 781 967  -90.7074011746   -63.2451213174    51.1043443841    21.5967406713, ERR= 1.4780107844e+002
 782   7  -44.9722397106   -61.6124332225    49.4096660408    80.7792753136, ERR= 9.6308421070e+001
 783  39  -71.6476171782   102.9963370061   -76.5511452800     5.7742505551, ERR= 1.0888106037e+002
 784  71   14.7254134511   -76.3525042643   -90.6729941639   -60.9925799875, ERR= 1.7312566677e+002
 785 103  -24.5323590815  -134.0720384126    63.9226312488     1.6829746963, ERR= 1.5922044312e+002
 786 135  -28.4508097987  -233.4959580635  -108.9048894430   -34.0156679549, ERR= 2.7934804275e+002
 787 167  -41.2065621632   -14.3008611672    20.2454662176    50.2782467495, ERR= 7.1209016743e+001
 788 199   34.3731108822     1.0768186151    39.2468462199   -55.2796825733, ERR= 5.4421537624e+001
 789 231  118.8540632655  -106.3859631797    87.2791358539    84.0230432699, ERR= 3.8692069316e+001
 790 263  155.8535416893    78.3446984827  -113.8014385882   -76.6935203523, ERR= 2.6966003556e+002
 791 295  -77.8088455207    99.7316074319   -68.2626686455    11.1404889912, ERR= 1.1128230434e+002
 792 327  -66.7630184058   -28.7012598494    21.8984041628    78.1456793130, ERR= 1.0151649358e+002
 793 359 -154.9605287449  -211.3651556154   -40.2313370517    19.5215154980, ERR= 2.2353248015e+002
 794 391  -86.5864577592   -24.9306567037  -200.3409524506    11.4237929294, ERR= 1.1455357014e+002
 795 423   24.1865450800    40.2362748672   113.0101126586   -39.8274926823, ERR= 8.8824508219e+001
 796 455  -24.0810521189   142.4141442460    34.6111546262   -92.8253467318, ERR= 7.6836345381e+001
 797 487  -48.7377849716  -174.7439733244   -45.7994691664  -152.4830549594, ERR= 3.2724022023e+002
 798 519  -85.0384869639    40.8592698269  -106.5839436879    44.6263916815, ERR= 8.8158975884e+001
 799 551   -3.0470134365  -100.6656068147    54.7641653963   100.1553131297, ERR= 5.7813430945e+001
 800 583  104.9029998474  -105.3301377961   -56.5480564049    78.7064263171, ERR= 1.6363149323e+002
 801 615   63.0591160334   -24.1669111201    41.8542378569    73.6079470665, ERR= 5.3796495183e+001
 802 647  -50.0101874401   151.9602370468    13.1435495109   -99.8738911986, ERR= 8.1861968671e+001
 803 679   24.0992496516    69.1484981582   -47.9941906291     4.1195794900, ERR= 1.0278947093e+002
 804 711   40.0762866268   -27.0319451072   170.7060733563   -82.3322165479, ERR= 1.7036625557e+002
 805 743  162.0205710583   -36.5414829659    35.7740100494     7.2179958998, ERR= 1.2960733413e+002
 806 806  -85.7093952427   -15.4831741752   -85.7093952428    15.4831741752
 807 838  -28.6588419552  -137.4490529646     7.9896222269    32.0813005645, ERR= 1.1155928098e+002
 808 870   53.7926376050  -135.0049028668  -147.8363825331    59.7288483239, ERR= 2.1522255028e+002
 809 902 -107.6503743616   -39.0555071587    -1.7222557025    17.2852587749, ERR= 1.0814208264e+002
 810 934  107.3220120595    40.4864927733    18.1881904740   -31.1536076457, ERR= 8.9621096262e+001
 811 966   83.5177741428   180.4031841496   103.6020691419  -136.7319150680, ERR= 4.8068270708e+001
 812   6  -33.3816139138    94.0252457237    94.6308582947   189.3488476200, ERR= 3.1094705308e+002
 813  38   55.2495531296   -65.1210919024    36.6944614161   111.7402926270, ERR= 5.0176102925e+001
 814  70  -52.5932789888    48.7415380689   100.1393216655   -84.1941976425, ERR= 1.5679329824e+002
 815 102 -196.6185172606   105.2949927915   109.0387447618   149.6481612637, ERR= 3.9802308178e+002
 816 134  -31.0110700831   124.0945579784   162.4122574508   109.5126853603, ERR= 3.0329017092e+002
 817 166  161.5277766128   124.5470263796   -95.6958798394  -126.8669693023, ERR= 2.5723411821e+002
 818 198  111.3559855651   -60.4301219893   -65.7877990550   136.0842831518, ERR= 1.9262261687e+002
 819 230   51.9713495007   -46.7285521962  -101.1021420573   -16.2685230825, ERR= 1.6552983209e+002
 820 262  -21.0174084811    58.5402515142   111.8927293124    62.4389414083, ERR= 1.7972498393e+002
 821 294   29.4028481839   104.9400010426   -51.4180868692     3.0962044360, ERR= 1.3492162628e+002
 822 326   42.3223432013    92.0646270673   -27.4880795508   -68.2372215474, ERR= 7.3764763801e+001
 823 358   90.0512348089   -73.5750120123   -34.6367661386  -116.0651805398, ERR= 2.2695924791e+002
 824 390  119.6104008403    70.8905699738   -59.0031436233   -64.8601543539, ERR= 1.7871531602e+002
 825 422    8.4822501085    77.7070101947     6.1637495413  -180.9373893612, ERR= 1.0325641204e+002
 826 454  182.9541453158   150.7270288336  -103.1518468763    22.1495834757, ERR= 3.3427976584e+002
 827 486  -52.8657921780    79.4495089363    85.3616806521  -137.4804207739, ERR= 1.4991471233e+002
 828 518   69.6772866022    -7.1474987504   -96.4967721161   -37.0666462075, ERR= 1.7195554194e+002
 829 550   20.6632197226    13.6285480774   -14.9988867292    28.6053233796, ERR= 5.5276448283e+001
 830 582   65.2270668339   -89.9728444958    -1.6459329190   190.4025161390, ERR= 1.2065702235e+002
 831 614  -60.2284365657  -114.0239369107    88.0407208249  -133.5657471935, ERR= 2.8859035796e+002
 832 646   56.7758101827   -58.5801965439   -32.1645223865    82.1457451111, ERR= 9.2009335596e+001
 833 678  -60.1098282069   -60.4926471241    55.0919640052  -137.1106223950, ERR= 2.2873238742e+002
 834 710 -127.9497715281    91.1016745809  -121.1976627828   289.8490941120, ERR= 3.8101060240e+002
 835 742   -9.4712311637   -90.5823895259   170.2283690068    42.9272255942, ERR= 1.8591116414e+002
 836 774   -4.7912827889   -16.6064408934   -23.4133592605  -129.7919534385, ERR= 1.4757801867e+002
 837 837   63.5450799416    57.9604814800    63.5450799416   -57.9604814800
 838 869   31.4921743589   -85.7306938869   -28.6588419551   137.4490529646, ERR= 7.9328011631e+001
 839 901   -7.5872152116  -197.9992845809   -45.3980005659   -25.9373926661, ERR= 2.2710634272e+002
 840 933 -149.6716108073   -94.9770039111     7.8107946387   195.8521254404, ERR= 1.8702004750e+002
 841 965  -86.7763973799   101.1911985237     6.6764861559    54.4568263885, ERR= 1.8154820049e+002
 842   5  133.6578003338   -52.9866766954    19.1641429604    29.1917316499, ERR= 1.1694014276e+002
 843  37   49.3386067562    15.5275726128   -96.4401762399  -141.6630024392, ERR= 1.9277344275e+002
 844  69  -40.0168582147   -89.4404296065    63.3439552773   -40.8413880006, ERR= 1.6630336667e+002
 845 101   91.8760235141  -122.3676141406    96.9781803813   100.3163610977, ERR= 2.2633819065e+001
 846 133   45.3504097312   -69.8924099911   -31.7916143581  -103.9072466153, ERR= 1.9015049965e+002
 847 165  -66.2805620128    22.6635158453   -54.3806333455    93.6391970664, ERR= 1.1690991974e+002
 848 197  147.5876769751    34.7004204429   178.5820179915   -83.2909101387, ERR= 5.7634059929e+001
 849 229  -32.3307867555   -45.8057618158   -22.6239330287    75.0422485516, ERR= 3.0805765141e+001
 850 261  120.1408139277    76.1220114177    64.6387108854   -23.6662567171, ERR= 7.6368119286e+001
 851 293  120.2520286753   -60.8628593160   -17.2461082482    21.3125973294, ERR= 1.4307327102e+002
 852 325  -35.4389436091   -99.6742301539    30.6947523285   139.2346569361, ERR= 7.7062916539e+001
 853 357  116.0296514701   -59.9738702403    -7.0870283665  -119.2844002803, ERR= 2.1746550164e+002
 854 389  -19.7467221195    11.9734199995   138.4848151749   259.6573729436, ERR= 3.1435729206e+002
 855 421   68.5919618926  -181.8645181981  -135.9875143956   -59.6500207493, ERR= 3.1651545719e+002
 856 453  111.2554690320    75.0485051139    18.8411536527   -64.7331996290, ERR= 9.2988231590e+001
 857 485 -137.5611615386  -113.3709515545    29.2666165354    27.8242394537, ERR= 1.8748265917e+002
 858 517   73.2079639954    78.4735525144   109.6300508979   105.2940753317, ERR= 1.8734222551e+002
 859 549    6.7054265826    84.6483055376   -65.0714731737  -148.8769960575, ERR= 9.6318471876e+001
 860 581   85.9417446380   -12.1738122279    57.8512675485  -220.3998763073, ERR= 2.3426394409e+002
 861 613  -52.8431024141   -65.4961763973    99.1212359417     5.0407358560, ERR= 1.6354822048e+002
 862 645 -127.8824124384  -100.8017687191    16.9471386866    23.5766518275, ERR= 1.6413201259e+002
 863 677  -16.3465860226    50.2907194852   -46.2880040172    65.0647586214, ERR= 1.1917791255e+002
 864 709  -75.3220827142    50.7216504058    72.2958357020    20.9603178977, ERR= 1.6410165879e+002
 865 741  -38.9758545708   189.5185478438   -28.4937762126   -14.5539534471, ERR= 1.7527830231e+002
 866 773 -154.1955070984     5.0684887720  -115.5924082722   154.6365805119, ERR= 1.6430431642e+002
 867 805   35.7740100494    -7.2179958998   -43.3161555250  -158.8354146603, ERR= 1.8392658712e+002
 868 868   53.7157287525  -136.8873016278    53.7157287526   136.8873016278
 869 900 -113.2071590153    96.6028828896    31.4921743590    85.7306938869, ERR= 2.3277334534e+002
 870 932   78.2165438141    41.2895213835    53.7926376050   135.0049028668, ERR= 1.7797823242e+002
 871 964  -60.2030160842   137.3019546120    27.9645753855   -57.3620664398, ERR= 1.1901222587e+002
 872   4    9.7351717174    47.4252385631    20.7338997283  -137.2104069575, ERR= 9.0456334667e+001
 873  36   38.0634647855   165.5296740870   -74.5629122659   -38.2351344669, ERR= 1.6996646912e+002
 874  68  -14.9440435830   -58.9122491562    91.6750568682   -46.4511511533, ERR= 1.4989689358e+002
 875 100   64.2582359468   124.8971682233    77.6347746261   -69.1680418373, ERR= 5.7312017193e+001
 876 132  -20.5027043263    50.9988151092  -118.5189513118   -68.2195845888, ERR= 9.9517534006e+001
 877 164    4.0886620368   113.0424629250   102.9038997650   -83.7774934509, ERR= 1.0305770056e+002
 878 196  116.2245334060    89.9339861721  -103.1751852076   -39.7950945376, ERR= 2.2505587080e+002
 879 228  -86.4109287920    90.2110824494   182.5689356946  -112.2600903473, ERR= 2.6988205989e+002
 880 260  -38.7384701899    24.6516877280   -16.3148994802   -55.7569685321, ERR= 3.8345208531e+001
 881 292  105.5315534295   127.6254719988  -138.7025005785     5.2914696650, ERR= 2.7805968158e+002
 882 324   -7.4531366408   -16.8301493520   -44.9093634520   -46.1353865187, ERR= 7.3264095124e+001
 883 356  -78.3313334497   -18.4126501263   -32.6394473977  -125.8076257732, ERR= 1.5128528161e+002
 884 388  117.5989824804  -174.7488045845  -134.4481019718  -155.3166676613, ERR= 4.1529621808e+002
 885 420  -39.7896726581  -243.7524783423     9.8449291512   -55.8274542960, ERR= 3.0366384331e+002
 886 452   37.6523804561   152.1617662614  -211.8473998895    46.6213564113, ERR= 3.1900606617e+002
 887 484  -65.0012730816   -84.8831942243    -5.3344019695   -25.9738763351, ERR= 1.2589450187e+002
 888 516  -73.5293580120    15.7890512611  -162.9650911577    12.7415230818, ERR= 9.3876216560e+001
 889 548  -58.2753402405  -151.1887455622     0.0349074474   -58.4851139068, ERR= 2.1763090849e+002
 890 580  138.8623878962   -93.2869125778    -4.1547077844    36.9284545673, ERR= 1.5372106377e+002
 891 612   24.7272782856   120.2862470233    92.5491754403    65.5241203855, ERR= 1.9780116878e+002
 892 644  -65.7820855185    78.1879124615    62.8144985249   -27.1757934001, ERR= 1.3834492300e+002
 893 676  -79.1017539281    81.7413508048     4.3487268956    18.4763906375, ERR= 1.3041310689e+002
 894 708  -24.2761615434    81.6283622821   -31.7960148098    15.8856069170, ERR= 9.7803488599e+001
 895 740  -10.7313355642    19.8115890350    55.3669560322   110.4991619699, ERR= 1.4611596757e+002
 896 772   80.9913763918   135.0682946327    84.1896309551   -79.2352317856, ERR= 5.5924589754e+001
 897 804  170.7060733563    82.3322165479   122.1887600894    74.5936024394, ERR= 1.6425480922e+002
 898 836  -23.4133592605   129.7919534385   186.0707673679    77.9988569750, ERR= 2.9506036705e+002
 899 899 -296.6219354846  -167.3258061401  -296.6219354847   167.3258061401
 900 931   18.8418793995  -188.1531988978  -113.2071590153   -96.6028828896, ERR= 3.1388369607e+002
 901 963   -6.9843542264   152.6346885253    -7.5872152115   197.9992845810, ERR= 3.5063449137e+002
 902   3  100.2337280653   -89.6547447405  -107.6503743616    39.0555071587, ERR= 2.1395345963e+002
 903  35  -44.4586291356   -29.6194452170     1.4143783056   -23.0261437726, ERR= 6.9827579449e+001
 904  67  214.1299125746    75.8587566718    37.0055254155  -133.2823239064, ERR= 1.8620020032e+002
 905  99  -69.8820090674    75.6495720356    83.0096889320    51.1974310462, ERR= 1.9866059878e+002
 906 131  -52.7177599948   165.1522989155   -69.9121676613    55.7051778793, ERR= 2.2152578340e+002
 907 163  112.2704796655    42.0455250043    89.0382740535   -74.8198416071, ERR= 4.0173264821e+001
 908 195   23.5386308858  -108.0500118802    12.5393551284  -132.3563346580, ERR= 2.4065783911e+002
 909 227 -118.7822927677    -4.8130758446   -75.4363727332    26.3775863324, ERR= 4.8413808941e+001
 910 259  -56.4518512214    83.3769464062     2.8137707232    82.1285739159, ERR= 1.7579673262e+002
 911 291  -26.0572522312   169.5310845415    32.5352984457  -214.5901593161, ERR= 7.3914864637e+001
 912 323  -80.6969757474    52.3581302743    49.9117815849    35.7087342864, ERR= 1.5752593477e+002
 913 355  -43.8463522062   -30.3286246779  -195.3131174961   -11.0707948696, ERR= 1.5702258731e+002
 914 387   66.5593079571    88.9637986360    57.3331592550    99.9913869519, ERR= 1.8918029491e+002
 915 419  195.4674255197   122.3772047245    54.3621135572    94.5973873375, ERR= 2.5882171985e+002
 916 451  111.3548605336    71.4532118416   -23.6878912726   -71.7923261932, ERR= 1.3504317759e+002
 917 483  -85.0112773150    23.5912161093   199.7365791657   119.9072713534, ERR= 3.1886228638e+002
 918 515   48.7534245618  -120.0622863913   -51.8846138870    63.5599100501, ERR= 1.1541461482e+002
 919 547    7.4893734772  -190.9041276345  -116.2484630952   -84.5679207423, ERR= 3.0198659182e+002
 920 579 -148.8902365931  -141.6286709234  -190.6283354803    75.1370405485, ERR= 7.8506087717e+001
 921 611  -31.6784857462   -79.2848038497    34.3498132862   104.3487008951, ERR= 7.0625315633e+001
 922 643 -206.6774951947   -10.4177699433   -33.5110793709  -119.1626953969, ERR= 2.1628154005e+002
 923 675   70.3344387520   225.8683103111    -8.8053533460   -11.2879030656, ERR= 2.2870911190e+002
 924 707   83.5610322290    53.5019207154   -94.6832745145   -17.2837604583, ERR= 1.8188674503e+002
 925 739  -88.1474952141    11.3029596307    45.1852203248  -139.6948011352, ERR= 1.8510018368e+002
 926 771  143.4736406552   -93.5875658127    61.4687192615    47.7608519700, ERR= 9.3940911399e+001
 927 803  -47.9941906289    -4.1195794900    14.4977342719   -88.4168461351, ERR= 1.1166123206e+002
 928 835  170.2283690068   -42.9272255942   -18.3270693932    24.2890218817, ERR= 1.8947436763e+002
 929 867  -43.3161555249   158.8354146602   -58.0109731073   -61.7448404559, ERR= 9.8196320008e+001
 930 930  165.4305637139   -54.1236317770   165.4305637139    54.1236317769
 931 962  138.9331981736    47.8506727936    18.8418793996   188.1531988978, ERR= 2.6480134497e+002
 932   2   71.4201695137    97.2511435757    78.2165438141   -41.2895213835, ERR= 5.6372811372e+001
 933  34  123.5512178722  -114.0327996558  -149.6716108072    94.9770039111, ERR= 2.7388654122e+002
 934  66  -22.8911680430   207.9751327320   107.3220120595   -40.4864927733, ERR= 2.1215069358e+002
 935  98  -93.6489595925    29.9758875926    76.9817150357   -57.5948733348, ERR= 1.7285148393e+002
 936 130   -9.9089010834   144.4883824085  -191.1581025531    56.4782275172, ERR= 2.7062677498e+002
 937 162   11.8924343495    -6.0568195512   -30.4764481307   -21.8247685729, ERR= 5.0719869469e+001
 938 194  -14.6093213404    85.1974432832    44.6767990100    66.9397406493, ERR= 1.6328063817e+002
 939 226 -141.0482470272     7.6127070611   -21.1396962853   -85.5707975192, ERR= 1.4302281080e+002
 940 258  -37.2972703930  -155.0196777193   -50.1088001829    88.6439945509, ERR= 6.7600788543e+001
 941 290   85.3836468317  -180.6613680340    15.8696102086   162.1231387255, ERR= 7.1943500287e+001
 942 322  -25.0402321749  -112.6033207930    78.8931948279    12.7206001722, ERR= 1.4414824011e+002
 943 354   22.5083459189   -37.3549091128   -53.9135947867  -165.4822970907, ERR= 2.1675618847e+002
 944 386 -125.5097276676   -11.7896964277    78.5021375315   -33.4385946499, ERR= 2.0896516326e+002
 945 418  -12.9318381412    74.5564997194    78.1938240255    33.2966605488, ERR= 1.4119557530e+002
 946 450 -127.6371998554    48.1309546961   -44.6481798072   204.7242705858, ERR= 2.6612580183e+002
 947 482  145.7183608596   -97.0239136208     1.9734867389   -53.3842750262, ERR= 2.0805098425e+002
 948 514  130.2300428861    71.4578660955    29.0062589614   -58.8390070410, ERR= 1.0200730384e+002
 949 546   -2.5718496017   -72.9524681928    39.0642988735    39.0991769118, ERR= 5.3662036771e+001
 950 578 -100.7999988120   -28.6295990354    31.7591610291   -46.1288445617, ERR= 1.5218658202e+002
 951 610  119.1480736188  -115.6689396803    -3.6012723595   177.2717568937, ERR= 1.3734012169e+002
 952 642  -19.6675668206   -19.1490604767   -46.8148933838    87.1408882316, ERR= 7.3211105585e+001
 953 674   19.4336537230    -7.8034027899    14.1416583255   -89.3255579892, ERR= 9.7273019061e+001
 954 706  144.5183306635    90.3924741531   -59.2197321197    60.4932902414, ERR= 2.5352655112e+002
 955 738  -40.9374327872   162.6745394865   -75.7092826997   -16.0514678970, ERR= 1.5068976962e+002
 956 770  -51.3893275419   -40.7650153161   -77.5629540294   -23.3171955268, ERR= 6.9221300696e+001
 957 802   13.1435495109    99.8738911986   224.9993064788   106.8318345368, ERR= 2.9599006539e+002
 958 834 -121.1976627829  -289.8490941120    19.6035368702  -115.2050854185, ERR= 4.2882848107e+002
 959 866 -115.5924082722  -154.6365805119   -94.7258363926    11.3228204763, ERR= 1.4482488611e+002
 960 898  186.0707673680   -77.9988569751    62.0194565151    52.7347150515, ERR= 1.2659780642e+002
 961 961  -60.3897014192  -135.4188535678   -60.3897014191   135.4188535679
 962   1   61.2766665287     6.3736654598   138.9331981736   -47.8506727936, ERR= 8.8039076804e+001
 963  33  -14.3516643651    62.8114343707    -6.9843542264  -152.6346885253, ERR= 9.0124881390e+001
 964  65  -21.2281654747   -28.9373064188   -60.2030160843  -137.3019546119, ERR= 1.7074697915e+002
 965  97   63.8032397090    14.6456690398   -86.7763973799  -101.1911985237, ERR= 1.7367888697e+002
 966 129  114.9982673715    -1.2792834409    83.5177741428  -180.4031841496, ERR= 1.8438964310e+002
 967 161  -67.1998504328    70.7799831625   -90.7074011746    63.2451213174, ERR= 1.3607106075e+002
 968 193  -48.0230007565   241.1519767177    48.2258652606    33.7182126259, ERR= 2.9123438190e+002
 969 225   77.8883251143   -21.2286242051   -40.2559619407    47.0178378923, ERR= 1.2092624242e+002
 970 257  -72.2818730544    82.6332267737    21.2847386872   -53.6454399165, ERR= 9.7954084242e+001
 971 289  195.2825062950   -56.4350043064   -52.3337677874  -108.5182404232, ERR= 2.9752880892e+002
 972 321  -61.4761750087    57.1611516884    40.9457769510    64.3096463109, ERR= 1.5888804552e+002
 973 353 -130.4266467311    14.0225271929   106.7355575654    22.8974067795, ERR= 2.4001873400e+002
 974 385   45.6716098939  -226.9708981243   108.4487037407   -63.3682182560, ERR= 2.9704842368e+002
 975 417   14.0920160499     1.3402017912   -12.4885235845    64.4906762335, ERR= 7.0994574361e+001
 976 449  -84.0084616716    15.5575942947  -138.4043188780   -18.6274858547, ERR= 5.4482414736e+001
 977 481  -64.3275106574  -139.7895968451  -186.8448839835    44.1554164722, ERR= 1.5542330334e+002
 978 513  187.3601307838   -13.4031030645   -21.9112008775   -83.6464895876, ERR= 2.3067967767e+002
 979 545  -36.5332342538   -82.6188955315    52.1539841648   -37.8420824864, ERR= 1.4958699788e+002
 980 577  -36.2344795589    69.9457797698  -190.0401565606   -54.6345135812, ERR= 1.5456591199e+002
 981 609 -102.4301065825    -5.6427747292   -72.1507395839    93.6970217112, ERR= 9.3114931549e+001
 982 641   14.0717851313    34.4801362382    12.6264095495    92.5280968120, ERR= 1.2701645710e+002
 983 673  -80.7445213898    32.1130062488    58.7182794055   168.1402197748, ERR= 2.4403120157e+002
 984 705   -9.7716154515    -6.5823029156   205.3951118824   -37.8036708561, ERR= 2.1969714431e+002
 985 737  -75.3083345899    34.1406107461    22.0702149093   -49.5442109581, ERR= 9.8589313833e+001
 986 769  164.4566926230     5.1182315889    33.5762451482   -52.5103750896, ERR= 1.3919664794e+002
 987 801   41.8542378569   -73.6079470665  -119.3877783951   -60.0481801021, ERR= 2.0943483028e+002
 988 833   55.0919640053   137.1106223949   210.8987956141    48.6429163723, ERR= 2.4244617122e+002
 989 865  -28.4937762126    14.5539534471   -65.5009613049    -8.4262506128, ERR= 3.7511071572e+001
 990 897  122.1887600892   -74.5936024395    88.7702339252    71.1158917355, ERR= 3.3598993478e+001
 991 929  -58.0109731073    61.7448404560   -90.0844492193    34.4067375301, ERR= 1.0135992215e+002
#endif
