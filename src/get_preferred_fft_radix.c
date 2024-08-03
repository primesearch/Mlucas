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

/*
	Given a target FFT length via the [kblocks] argument, looks for a runtime-optimized *.cfg
	preferred-FFT-radix file (assumed to have been created via a set of timing runs on the
	platform in question). If one is found containing an entry for the specified FFT length,
	returns the preferred FFT length as the function result and the preferred radix set
	(assumed to have been determined previously via self-test on the host platform)
	in the RADIX_VEC global.  We define "preferred FFT length" here as follows:

		Given an input FFT length NI, the preferred FFT length NP is the FFT length >= NI
		whose .cfg file entry contains the smallest per-iteration timing.

	This may or may not be the same as the input FFT length, depending on how well we've
	done our FFT coding (especially for the less-smooth leading FFT radices) and the quality
	of the platform (compiler + hardware) which the code is built and run on. It also implicitly
	assumes that if we have timing data for at least one FFT length NK > NI, then we also have timings
	for any supported FFT lengths NJ lying between NI and NK, so as to have a full range of options
	from which to choose. (E.g. if on some platform we had a faster timing for 2048K than for
	1792K, we'd want to know that the 2048K timing is also faster than that for 1920K.)

	Looks for the platform-specific *.cfg file for the number type stored in the MODULUS_TYPE global
	(e.g. mlucas.cfg for Mersennes and fermat.cfg for Fermats) and if the requisite file is found,
	look for an entry for the input FFT length (in units of K doubles).

	If an entry for the input FFT length is found in the .cfg file (i.e. a timing-test run has
	been done for that length), compares that timing to the best of any timings found for larger
	FFT lengths and returns the best of these as described above, but in the following manner:

	1) Stores the best-timing FFT radix set data for FFT length [kblocks] in the NRADICES and RADIX_VEC[]
	   globals, irrespective of whether a larger FFT length with a better timing was found in the .cfg file;

	2) If timing for the input FFT length [kblocks] is better than any larger ones, function returns value = kblocks;

	3) If a larger FFT length with a better timing was found in the .cfg file, returns the FFT radix set data
	for that FFT length in a compact bitwise form in the return value, as follows:

		- Bits <0:9> store (leading radix-1): We subtract the 1 so radices up to 1024 can be stored;
			<*** EWM Jan 2014: For F33 we want r0 ~= 4096, so will need to modify the code here.
				Once we are dealing with r0 > 1024, consider storing e.g. [odd part of r0, lg(pow2 part)]
				in some compact fashion, say 6 bits for any odd component (allowing these up to 63) and the
				rest for the lg2(pow2(r0)) term. Thus e.g. r0 = 4032 = 63*64 would map to a [63,6] pair,
				needing just 6+3 = 9 bits to store, as opposed to 11 bits for (r0-1) = 4095 = 0b11111111111 .
			***>
		- Bits <10:13> store (number of FFT radices);
		- Each successive pair of higher-order bits stores log2[(intermediate FFT radix)/8]: Since
		  our smallest permitted intermediate FFT radix is 8 and these must be powers of 2, this
		  again permits radices up to 64 to be stored using just 2 bits. Radix-8 of course maps to 0
		  under this scheme, but we know when to stop because bits <10:13> tell us the number of radices,
		  which can be as large as 10 under this scheme.

	In order to make it easy for the user to extract these bitwise FFT-radix data from the function
	return value, we define 2 handy utility functions in util.c:

		uint32	extractFFTlengthFrom32Bit (uint32 n) - returns the (real-vector) FFT length encoded by n according to the above scheme
		void	extractFFTradicesFrom32Bit(uint32 n) - extracts the FFT-radix data encoded by n and stores in the NRADICES and RADIX_VEC[] globals

	If the return FFT-length value differs from the input [kblocks] (which implies that a better timing
	datum was found for at least one larger FFT length in the .cfg file), caller must decide whether
	to reset the FFT length used for the computation to the larger value - if so, caller can use the
	2nd of the above-described functions to read the corr. FFT radix data into the NRADICES and RADIX_VEC globals.

	Returns 0 if no .cfg file is found or if the .cfg file contains no properly-formatted entry
	for the input FFT length. In this case the caller should do a timing test at the input
	FFT length so as to find the optimal radix set on-the-fly, or simply use the radix set index 0
	(guaranteed to be supported if the FFT length is) if maximal execution speed is not crucial.
*/
uint32	get_preferred_fft_radix(uint32 kblocks)
{
	uint32 i, j, k, kprod, found, retval = 0;
	double tbest = 0, tcurr;
	char *char_addr;

	/* FFT-radix configuration file is named mlucas.cfg or fermat.cfg,
	depending on whether Mersenne or Fermat number test is being done:
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		strcpy(CONFIGFILE,"fermat.cfg");
	else	/* For now, anything else gets done using the Mersenne-mod .cfg file */
		strcpy(CONFIGFILE,"mlucas.cfg");

	/*...Look for any FFT length >= [kblocks] and check the per-iteration timing: */
	found = 0;		/* Was an entry for the specified FFT length found in the .cfg file? */
	fp = mlucas_fopen(CONFIGFILE,"r");
	if(fp) {
		while(fgets(in_line, STR_MAX_LEN, fp)) {
		//	fprintf(stderr,"Current line: %s",in_line);
			/* Each FFT-length entry assumed to begin with an int followed by whitespace;
			any line not of that form is ignored, thus allowing pretty much any common comment-format: */
			if(sscanf(in_line, "%d", &i) == 1) {
				/* Consider any entry with an FFT length >= target, which further contains
				a per-iteration timing datum in the form 'msec/iter = [float arg]' in non-exponential form:
				*/
				if((i >= kblocks) && (char_addr = strstr(in_line, "msec/iter =")) != 0) {
					/* Stores whether we found an entry for the requested FFT length
					(whether that proves to have the best timing for lengths >= kblocks or not): */
					if(i == kblocks) {
						if(found) {
							sprintf(cbuf,"Multiple cfg-file entries for FFT length %uK encountered in %s - please delete or comment out all but one entry for this length, save the file and retry.",kblocks,CONFIGFILE);
							ASSERT(0,cbuf);
						} else
							found = TRUE;
					}
					if(sscanf(char_addr + 11, "%lf", &tcurr) == 1) {	// 11 chars in "msec/iter ="
						ASSERT(tcurr >= 0, "tcurr < 0!");
						if((tbest == 0.0) || ((tcurr > 0.0) && (tcurr < tbest))) {
							if((char_addr = strstr(in_line, "radices =")) == 0x0) {
								snprintf(cbuf,STR_MAX_LEN*2,"get_preferred_fft_radix: invalid format for %s file: 'radices =' not found in timing-data line %s", CONFIGFILE, in_line);
								ASSERT(0, cbuf);
							}
							char_addr += 9;	// 9 chars in "radices ="
							kprod = 1;	/* accumulate product of radices */
							for(j = 0; j < 10; j++) {	/* Read in the radices */
								if(sscanf(char_addr, "%d", &k) != 1) {
									snprintf(cbuf,STR_MAX_LEN*2,"get_preferred_fft_radix: invalid format for %s file: failed to read %dth element of radix set, offending input line %s", CONFIGFILE, j, in_line);
									ASSERT(0, cbuf);
								} else {
									// Advance to next WS char following the current numeric token - since sscanf skips leading WS,
									// Must do this in 2 steps. NOTE we *need* the trailing ; here to serve as executable-statement
									// loop bodies, otherwise the ensuing while or if() is treated so and each while() executes just once.
									// ***NOTE*** It is crucial to separate the loop-test from the ptr-incrementing here, because if e.g.
									// we have current radix k = 8, an opening while( isspace(*char_addr++)) increments char_addr to the WS
									// *following* the 8, and the loop continues, causing us to "lose the current radix",
									// leading to an eventual ASSERT in the kprod-based looping sanity checks.
									while( isspace(*char_addr)) char_addr++;	// 1. First skip any WS preceding current numeric token
									while(!isspace(*char_addr)) char_addr++;	// 2. Look for first WS char following current numeric token
									if(j == 0)
										ASSERT(k <= 1024, "get_preferred_fft_radix: Leading radix > 1024: out of range!");
									else if(k) {
										ASSERT(k <= 32  , "get_preferred_fft_radix: Intermediate radix > 32: out of range!");
										ASSERT(isPow2(k), "get_preferred_fft_radix: Intermediate FFT radix not a power of 2!");
									}
									/* If (i == kblocks), store the data directly into the NRADICES and RADIX_VEC[] globals: */
									if(i == kblocks) {
										if(k == 0) {
											ASSERT(!NRADICES, "Zero terminator of radix set found but NRADICES != 0 ... please check your mlucas.cfg file for duplicate FFT-length entries and remove the unwanted ones, or delete the file and rerun the self-test.");
											NRADICES = j;
											break;
										} else {
											kprod *= k;
											RADIX_VEC[j] = k;
										}
									} else {	/* Otherwise, store radix-set data into retval in above-described compact form */
										if(k == 0) {	/* Bits <10:13> store (number of FFT radices): */
											if(!((retval >> 10) & 0xf))	/* Set based only position of first zero in the list */
												retval += (j << 10);
										} else {
											kprod *= k;
										}
										/* Bits <0:9> store (leading radix-1): */
										if(j == 0)
											retval = k - 1;
										else if(k) {	/* Each successive pair of higher-order bits stores log2[(intermediate FFT radix)/8]: */
											k = trailz32(k) - 3;
											retval += (k << (12 + 2*j));
										}
									}
								}
							}
							/* Product of real-FFT radices (kblocks) must be divisible by 1K = 1024
							Since (kprod) here is product of complex radices, first multiply it by 2:
							*/
							kprod *= 2;
							if((kprod & 1023) != 0) {
								snprintf(cbuf,STR_MAX_LEN*2,"get_preferred_fft_radix: illegal data in %s file: product of complex radices (%d) not a multiple of 1K! Offending input line %s", CONFIGFILE, kprod, in_line);
								ASSERT(0, cbuf);
							}
							kprod >>= 10;
							tbest = tcurr;
							if(i == kblocks) {
								/* Product of radices must equal complex vector length (n/2): */
								if(kprod != kblocks) {
									snprintf(cbuf,STR_MAX_LEN*2,"get_preferred_fft_radix: mismatching data in %s file: (product of complex radices)/2^10 (%d) != kblocks/2 (%d), offending input line %s", CONFIGFILE, kprod, kblocks/2, in_line);
									ASSERT(0, cbuf);
								}
								retval = i;			/* Preferred FFT length */
							} else {
								ASSERT(i == extractFFTlengthFrom32Bit(retval), "get_preferred_fft_radix: i != extractFFTlengthFrom32Bit(retval)!");
							}
						}
					}
				}
			}
		}
		fclose(fp);	fp = 0x0;
	} else {
		sprintf(cbuf, "CONFIGFILE = %s: open failed -- please run the post-build self-tests as described in the README!", CONFIGFILE);
		ASSERT(0 , cbuf);
	}

	/* Only return nonzero if an entry for the specified FFT length was found.
	Otherwise clear RADIX_VEC and return 0:
	*/
	if(!found) {
		retval = 0;
		for(j=0; j<10; j++) { RADIX_VEC[j] = 0; }
		NRADICES = 0;
	}
	return retval;
}

/********* Functions related to FFT-radix-set compact 32-bit encoding ***********/

/* returns the (real-vector) FFT length encoded by n according to the above scheme */
uint32	extractFFTlengthFrom32Bit (uint32 n)
{
	uint32 i, nrad, retval;
	/* Bits <0:9> store (leading radix-1): We subtract the 1 so radices up to 1024 can be stored: */
	retval = (n & 0x3ff) + 1;	n >>= 10;
	ASSERT(retval > 4, "extractFFTlengthFrom32Bit: Leading radix must be 5 or larger!");
	/* Bits <10:13> store (number of FFT radices): */
	nrad   = (n & 0xf)    ;	n >>= 4;
	ASSERT(nrad >=  3, "extractFFTlengthFrom32Bit: Number of radices must be 3 or larger!");
	/* Each successive pair of higher-order bits stores log2[(intermediate FFT radix)/8]: */
	for(i = 1; i < nrad; i++)	/* Already done leading radix, so start at 1, not 0 */
	{
		retval *= ( 0x1 << ((n & 0x3)+3) );	n >>= 2;
	}
	/* return value is in units of K - combine div-by-2^10 with mul-by-2 (real-array length) here: */
	return (retval >> 9);
}

/* extracts the FFT-radix data encoded by n and stores in the NRADICES and RADIX_VEC[] globals */
void	extractFFTradicesFrom32Bit(uint32 n)
{
	uint32 i, nrad, retval;
	/* Bits <0:9> store (leading radix-1): We subtract the 1 so radices up to 1024 can be stored: */
	retval = (n & 0x3ff) + 1;	n >>= 10;
	ASSERT(retval > 4, "extractFFTradicesFrom32Bit: Leading radix must be 5 or larger!");
	RADIX_VEC[0] = retval;
	/* Bits <10:13> store (number of FFT radices): */
	nrad   = (n & 0xf)    ;	n >>= 4;
	ASSERT(nrad >=  3, "extractFFTradicesFrom32Bit: Number of radices must be 3 or larger!");
	ASSERT(nrad <= 10, "extractFFTradicesFrom32Bit: Number of radices must be 10 or smaller!");
	NRADICES = nrad;
	/* Each successive pair of higher-order bits stores log2[(intermediate FFT radix)/8]: */
	for(i = 1; i < 10; i++)	/* Already done leading radix, so start at 1, not 0 */
	{
		if(i < nrad)
			RADIX_VEC[i] = ( 0x1 << ((n & 0x3)+3) );
		else
			RADIX_VEC[i] = 0;

		n >>= 2;
	}
}

