#include "Mlucas.h"

#define radix 44
#define TEST_TYPE 1	// 0 = DIF, 1 = DIT, 2 = BOTH

void matmul_double (double **, double *, double *, int, int);
void matmul_complex(struct complex **, struct complex *, struct complex *, int, int);

double ISRT2 = .70710678118654752440;

void test_fft_radix(void)
{
	int i,j,l,nradices, index[radix], nerr = 0;
	int radix_prim[10];
	double err_r, err_i;
	double a[2*radix] = {
	 3,1
#if radix >= 2
	,4,1
#endif
#if radix >= 3
	,5,9
#endif
#if radix >= 4
	,2,6
#endif
#if radix >= 5
	,5,3
#endif
#if radix >= 6
	,5,8
#endif
#if radix >= 7
	,9,7
#endif
#if radix >= 8
	,9,3
#endif
#if radix >= 9
	,2,3
#endif
#if radix >= 10
	,8,4
#endif
#if radix >= 11
	,6,2
#endif
#if radix >= 12
	,6,4
#endif
#if radix >= 13
	,3,3
#endif
#if radix >= 14
	,8,3
#endif
#if radix >= 15
	,2,7
#endif
#if radix >= 16
	,9,5
#endif
#if radix >= 17
	,0,2
#endif
#if radix >= 18
	,8,8
#endif
#if radix >= 19
	,4,1
#endif
#if radix >= 20
	,9,7
#endif
#if radix >= 21
	,1,6
#endif
#if radix >= 22
	,9,3
#endif
#if radix >= 23
	,9,9
#endif
#if radix >= 24
	,3,7
#endif
#if radix >= 25
	,5,1
#endif
#if radix >= 26
	,0,5
#endif
#if radix >= 27
	,8,2
#endif
#if radix >= 28
	,0,9
#endif
#if radix >= 29
	,7,4
#endif
#if radix >= 30
	,9,4
#endif
#if radix >= 31
	,4,5
#endif
#if radix >= 32
	,9,2
#endif
#if radix >= 33
	,3,0
#endif
#if radix >= 34
	,7,8
#endif
#if radix >= 35
	,1,6
#endif
#if radix >= 36
	,4,0
#endif
#if radix >= 37
	,6,2
#endif
#if radix >= 38
	,8,6
#endif
#if radix >= 39
	,2,0
#endif
#if radix >= 40
	,8,9
#endif
#if radix >= 41
	,9,8
#endif
#if radix >= 42
	,6,2
#endif
#if radix >= 43
	,8,0
#endif
#if radix >= 44
	,3,4
#endif
#if radix >= 45
	,8,2
#endif
#if radix >= 46
	,5,3
#endif
#if radix >= 47
	,4,2
#endif
#if radix >= 48
	,1,1
#endif
#if radix >= 49
	,7,0
#endif
#if radix >= 50
	,6,7
#endif
#if radix >= 51
	,9,8
#endif
#if radix >= 52
	,2,1
#endif
#if radix >= 53
	,4,8
#endif
#if radix >= 54
	,0,8
#endif
#if radix >= 55
	,6,5
#endif
#if radix >= 56
	,1,3
#endif
#if radix >= 57
	,2,8
#endif
#if radix >= 58
	,2,3
#endif
#if radix >= 59
	,0,6
#endif
#if radix >= 60
	,6,4
#endif
#if radix >= 61
	,7,0
#endif
#if radix >= 62
	,9,3
#endif
#if radix >= 63
	,8,4
#endif
#if radix >= 64
	,4,6
#endif
#if radix >= 65
	,6,2
#endif
#if radix >= 66
	,8,3
#endif
#if radix >= 67
	,1,8
#endif
#if radix >= 68
	,5,3
#endif
#if radix >= 69
	,0,7
#endif
#if radix >= 70
	,1,7
#endif
#if radix >= 71
	,9,5
#endif
#if radix >= 72
	,8,6
#endif
#if radix >= 73
	,4,7
#endif
#if radix >= 74
	,6,9
#endif
#if radix >= 75
	,7,8
#endif
#if radix >= 76
	,9,1
#endif
#if radix >= 77
	,5,1
#endif
#if radix >= 78
	,5,6
#endif
#if radix >= 79
	,3,0
#endif
	};
	struct complex *ac = (struct complex *)a;
	double b[2*radix], arrtmp[2*radix];
	struct complex *bc = (struct complex *)b;

	struct complex mat[radix][radix], *matp[radix];
	double theta;
	double twopi = 6.2831853071795864769;

	/* Make sure this utility is run only in non-SSE2 mode, since the latter's data layout hoses us here. */
#ifdef USE_SSE2
	ASSERT(HERE, 0, "test_fft_radix: USE_SSE2 must not be set!");
#endif
	ASSERT(HERE, ((radix >> trailz32(radix)) < 16), "test_fft_radix: Illegal radix; must be odd*2^n with odd < 16");
	/* These may not have been init'ed yet, so do it here: */
	DAT_BITS = DAT_BITS_DEF;
	PAD_BITS = PAD_BITS_DEF;

	/* Init DFT matrix */

	for(i = 0; i < radix; i++)
	{
		matp[i] = &mat[i][0];

		theta = i*twopi/radix;
		for(j = 0; j < radix; j++)
		{
			mat[i][j].re = cos(j*theta);
			mat[i][j].im = sin(j*theta);
		}
	}
	matmul_complex(matp,ac,bc,radix,radix);

	// If doing a DIF followed by a DIT, save a copy of the original data:
#if TEST_TYPE == 2
	for(i = 0; i < radix ; i++)
	{
		arrtmp[2*i	] = a[2*i	];
		arrtmp[2*i+1] = a[2*i+1];
	}
#endif

/*...Forward (DIF) FFT sincos data are in bit-reversed order.	*/

	l =0;

	switch(radix){
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
	default :
		printf("FATAL: radix %d not available. Halting...\n",radix); exit(EXIT_FAILURE);
	}

/*...Allocate and initialize an index array containing (radix) indices...	*/

	for(i=0; i < radix; i++)
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

	bit_reverse_int(&index[0],radix,nradices,&radix_prim[nradices-1],-1,(int *)0x0);
	printf("bit-reversal index array = [");
	for(i=0; i < radix; i++)
	{
		printf(" %d",index[i]);
	}
	printf("]\n\n");

#if TEST_TYPE == 0 || TEST_TYPE == 2

	#if radix == 2
		radix2_dif_pass1 (a,2*radix);
	#elif radix == 3
		radix3_dif_pass1 (a,2*radix);
	#elif radix == 4
		radix4_dif_pass1 (a,2*radix);
	#elif radix == 5
		radix5_dif_pass1 (a,2*radix);
	#elif radix == 6
		radix6_dif_pass1 (a,2*radix);
	#elif radix == 7
		radix7_dif_pass1 (a,2*radix);
	#elif radix == 8
		radix8_dif_pass1 (a,2*radix);
	#elif radix == 9
		radix9_dif_pass1 (a,2*radix);
	#elif radix == 10
		radix10_dif_pass1 (a,2*radix);
	#elif radix == 11
		radix11_dif_pass1 (a,2*radix);
	#elif radix == 12
		radix12_dif_pass1 (a,2*radix);
	#elif radix == 13
		radix13_dif_pass1 (a,2*radix);
	#elif radix == 14
		radix14_dif_pass1 (a,2*radix);
	#elif radix == 15
		radix15_dif_pass1 (a,2*radix);
	#elif radix == 16
		radix16_dif_pass1 (a,2*radix);
	#elif radix == 18
		radix18_dif_pass1 (a,2*radix);
	#elif radix == 20
		radix20_dif_pass1 (a,2*radix);
	#elif radix == 22
		radix22_dif_pass1 (a,2*radix);
	#elif radix == 24
		radix24_dif_pass1 (a,2*radix);
	#elif radix == 26
		radix26_dif_pass1 (a,2*radix);
	#elif radix == 28
		radix28_dif_pass1 (a,2*radix);
	#elif radix == 30
		radix30_dif_pass1 (a,2*radix);
	#elif radix == 32
		radix32_dif_pass1 (a,2*radix);
	#elif radix == 36
		radix36_dif_pass1 (a,2*radix);
	#elif radix == 40
		radix40_dif_pass1 (a,2*radix);
	#elif radix == 44
		radix44_dif_pass1 (a,2*radix);
	#elif radix == 48
		radix48_dif_pass1 (a,2*radix);
	#elif radix == 52
		radix52_dif_pass1 (a,2*radix);
	#elif radix == 56
		radix56_dif_pass1 (a,2*radix);
	#elif radix == 60
		radix60_dif_pass1 (a,2*radix);
	#elif radix == 64
		radix64_dif_pass1 (a,2*radix);
	#elif radix == 72
		radix72_dif_pass1 (a,2*radix);
	#elif radix == 80
		radix80_dif_pass1 (a,2*radix);
	#elif radix == 88
		radix88_dif_pass1 (a,2*radix);
	#elif radix == 96
		radix96_dif_pass1 (a,2*radix);
	#elif radix ==104
		radix104_dif_pass1(a,2*radix);
	#elif radix ==112
		radix112_dif_pass1(a,2*radix);
	#elif radix ==120
		radix120_dif_pass1(a,2*radix);
	#elif radix ==128
		radix128_dif_pass1(a,2*radix);
	#endif

#endif

#if TEST_TYPE == 0

	nerr = 0;
	for(i = 0; i < radix ; i++)
	{
		j = index[i];
    	printf("%4d (%4d) %15.10f  %15.10f  %15.10f  %15.10f\n",i,j,a[2*i  ],a[2*i+1],b[2*j  ],b[2*j+1]);

		err_r = fabs(a[2*i  ]-b[2*j]);  err_i = fabs(a[2*i+1]-b[2*j+1]);
		if(err_r > 0.00001 || err_i > 0.00001)
		{
			printf(" ERR= %15.10f  %15.10f\n", err_r, err_i);
			++nerr;
		}
		/*
		j=i;
		printf("%4d  %20.10f  %20.10f  %20.10f  %20.10f\n",i,bc[j].re,bc[j].im,a[2*i  ],a[2*i+1]);
		*/
	}
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIF transform!");

#endif

#if TEST_TYPE == 1

	// Bit-reverse the inputs to the transform...
	for(i = 0; i < radix ; i++)
	{
		j = index[i];
		arrtmp[2*i	] = a[2*j	];
		arrtmp[2*i+1] =-a[2*j+1];
	}
	for(i = 0; i < 2*radix ; i++)
	{
		a[i] = arrtmp[i];
	}

#endif

#if TEST_TYPE > 0

	#if radix == 2
		radix2_dit_pass1 (a,2*radix);
	#elif radix == 3
		radix3_dit_pass1 (a,2*radix);
	#elif radix == 4
		radix4_dit_pass1 (a,2*radix);
	#elif radix == 5
		radix5_dit_pass1 (a,2*radix);
	#elif radix == 6
		radix6_dit_pass1 (a,2*radix);
	#elif radix == 7
		radix7_dit_pass1 (a,2*radix);
	#elif radix == 8
		radix8_dit_pass1 (a,2*radix);
	#elif radix == 9
		radix9_dit_pass1 (a,2*radix);
	#elif radix == 10
		radix10_dit_pass1 (a,2*radix);
	#elif radix == 11
		radix11_dit_pass1 (a,2*radix);
	#elif radix == 12
		radix12_dit_pass1 (a,2*radix);
	#elif radix == 13
		radix13_dit_pass1 (a,2*radix);
	#elif radix == 14
		radix14_dit_pass1 (a,2*radix);
	#elif radix == 15
		radix15_dit_pass1 (a,2*radix);
	#elif radix == 16
		radix16_dit_pass1 (a,2*radix);
	#elif radix == 18
		radix18_dit_pass1 (a,2*radix);
	#elif radix == 20
		radix20_dit_pass1 (a,2*radix);
	#elif radix == 22
		radix22_dit_pass1 (a,2*radix);
	#elif radix == 24
		radix24_dit_pass1 (a,2*radix);
	#elif radix == 26
		radix26_dit_pass1 (a,2*radix);
	#elif radix == 28
		radix28_dit_pass1 (a,2*radix);
	#elif radix == 30
		radix30_dit_pass1 (a,2*radix);
	#elif radix == 32
		radix32_dit_pass1 (a,2*radix);
	#elif radix == 36
		radix36_dit_pass1 (a,2*radix);
	#elif radix == 40
		radix40_dit_pass1 (a,2*radix);
	#elif radix == 44
		radix44_dit_pass1 (a,2*radix);
	#elif radix == 48
		radix48_dit_pass1 (a,2*radix);
	#elif radix == 52
		radix52_dit_pass1 (a,2*radix);
	#elif radix == 56
		radix56_dit_pass1 (a,2*radix);
	#elif radix == 60
		radix60_dit_pass1 (a,2*radix);
	#elif radix == 64
		radix64_dit_pass1 (a,2*radix);
	#elif radix == 72
		radix72_dit_pass1 (a,2*radix);
	#elif radix == 80
		radix80_dit_pass1 (a,2*radix);
	#elif radix == 88
		radix88_dit_pass1 (a,2*radix);
	#elif radix == 96
		radix96_dit_pass1 (a,2*radix);
	#elif radix ==104
		radix104_dit_pass1(a,2*radix);
	#elif radix ==112
		radix112_dit_pass1(a,2*radix);
	#elif radix ==120
		radix120_dit_pass1(a,2*radix);
	#elif radix ==128
		radix128_dit_pass1(a,2*radix);
	#endif

#endif

#if TEST_TYPE == 1

	nerr = 0;
	for(i = 0; i < radix ; i++)
	{
		printf("%4d  %15.10f  %15.10f  %15.10f  %15.10f\n",i,a[2*i  ],a[2*i+1],b[2*i  ],-b[2*i+1]);

		err_r = fabs(a[2*i  ]-b[2*i]);  err_i = fabs(a[2*i+1]+b[2*i+1]);
		if(err_r > 0.00001 || err_i > 0.00001)
		{
			printf(" ERR= %15.10f  %15.10f\n", err_r, err_i);
			++nerr;
		}
	}
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIT transform!");

#endif

#if TEST_TYPE == 2

	nerr = 0;
	for(i = 0; i < radix ; i++)
	{
		a[2*i] /= radix;	a[2*i+1] /= radix;

		printf("%4d  %25.15f  %25.15f\n",i,a[2*i], a[2*i+1]);

		err_r = fabs(a[2*i  ]-arrtmp[2*i]);  err_i = fabs(a[2*i+1]-arrtmp[2*i+1]);
		if(err_r > 0.00001 || err_i > 0.00001)
		{
			printf(" ERR= %15.10f  %15.10f\n", err_r, err_i);
			++nerr;
		}
	}
	printf("\n");
	ASSERT(HERE, nerr == 0, "test_fft_radix: Mismatches detected in DIF/DIT combo!");

#endif
	printf("");
}

/* RADIX-24 test data:
normal output order:								twiddleless order:
	 0				 64.0000000000				 0				 64.0000000000	a0
	 1				 51.0000000000				 1				 51.0000000000
	 2				 -1.9019237886				 2				 -5.0000000000	a3
	 3				 -2.9019237886				 3				  8.0000000000
	 4				 -1.2679491924				 4				  6.0000000000	a6
	 5				 -4.3397459622				 5				-13.0000000000
	 6				 -5.0000000000				 6				 15.0000000000	a9
	 7				  8.0000000000				 7				  2.0000000000
	 8				 -4.0717967697				 8				 -4.0717967697	a4
	 9				 -5.1961524227				 9				 -5.1961524227
	10				 -7.0980762114				10				 -0.9019237886	a7
	11				 -8.0980762114				11				  8.2942286341
	12				  6.0000000000				12				 -4.7320508076	a10
	13				-13.0000000000				13				-21.6602540378
	14				 -0.9019237886				14				 -1.9019237886	a1
	15				  8.2942286341				15				 -2.9019237886
	16				-17.9282032303				16				-17.9282032303	a8
	17				  5.1961524227				17				  5.1961524227
	18				 15.0000000000				18				 -6.0980762114	a11
	19				  2.0000000000				19				 -7.2942286341
	20				 -4.7320508076				20				 -1.2679491924	a2
	21				-21.6602540378				21				 -4.3397459622
	22				 -6.0980762114				22				 -7.0980762114	a5
	23				 -7.2942286341				23				 -8.0980762114
*/

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
		}
	}

	return;
}
