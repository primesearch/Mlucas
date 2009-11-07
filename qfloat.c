/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

#include "qfloat.h"

/* Can turn on various qfloat-specific debugging features here: */
#ifndef QFDEBUG
	#define QFDEBUG		0
#endif

#if QFDEBUG
	#define QFONLY		0	/* Set to 1 to compile qfloat package in standalone mode - required if QFDEBUG is set! */
#else
	#define QFONLY		0	/* Set to 1 to compile qfloat package in standalone mode. */
#endif

#define DEBUG_QFMUL	0

/* Needed qfloat constants. (Remember, the high half of any qfloat has the same
	form as the IEEE double approximation to the qfloat, except for any differences
	due to rounding of the lowest bit in the latter.)
	Here's how to do this for any desired constant using the Unix 'bc' utility:

	0) Invoke 'bc -l' (i.e. bc in floating-point mode). Set (say) 75 decimal digits
	of precision and specify binary output, as follows:

	scale=75
	obase=2

	1) Generate the desired constant. We'll use 2*pi as out example, and we'll get
	this value as 8*arctan(1), which in bc syntax is:

	8*a(1)

	The output is

	110.010010000111111011010101000100010000101101000110000100011010011000\
	1001100011001100010100010111000000011011100000111001101000100101001000\
	0001001001110000010001000101001100111110011000111010000000010000010111\
	01111101010011000111011000100111001101011

	2) If the number in question is nonnegative, add the number of digits to the left of the binary point to the
	hex constant 0x3FE (= 1022 in decimal). In our example, we add 3, getting 0x401. If the number is negative,
	add to 0xBFE. The result is the 12-bit sign/exponent field of the qfloat. If the binary representation of the
	number has no nonzero digits to the left of the binary point, count the number of leading zeros to the right
	of the binary point and subtract this number from the appropriate one of the above two hex constants.
	(If the result is <= 0, your constant has underflowed, and special steps must be taken to represent it properly -
	see any reference on the IEEE floating-point standard for how to do this; most computer architecture texts have
	a section on this.)

	3) Delete the binary point and any leading zeros from the bc printout, and also delete the leading (leftmost)
	ones bit - this is the so-called "hidden bit" of the IEEE representation. Invoke another bc run, this time
	without the -l flag (i.e. now work in whole-number mode), and set the input base to be 2 by entering

	ibase=2

	Now cut the leading 52 bits left over after the above hidden-bit deletion and paste them into the bc shell.
	Then cut and paste the next 64 bits. If the leading bit of the remaining bits = 1, add one to the 64 bits
	that form the lower half of our qfloat. In our example, the leading 52 bits (sans hidden bit) are

	1001001000011111101101010100010001000010110100011000

	and the bc decimal output is

	2570638124657944 .

	The next 64 bits are

	0100011010011000100110001100110001010001011100000001101110000011

	or, in decimal:

	5086983782422027139 .

	The leftover bits are

	1001101000100101001000...,

	which have a leading ones bit, so we round by adding one to 5086983782422027139, getting 5086983782422027140 .


	4) This next step is separate from (3) since I've yet to find a bc implementation that allows one to switch
	both the input base and the output base from their defaults of 10 and correctly print things. That means that
	to convert binary to hex (which is what we are doing) we need to do 2 steps: (i) binary to decimal;
	(ii) decimal to hex. So, kill your bc run of step (3) and start another, again working in whole-number mode.
	Set the output format to hex via

	obase=16

	and paste in the two decimal outputs of step (3). In our example, the upper 52 bits of our significand are,
	in decimal:

	2570638124657944 .

	and the bc hex output is

	921FB54442D18 .

	Left-zero-pad (if necessary) the upper-52-bit result to be 13 hex characters long and concatenate
	the resulting string with the sign/exponent string you got in step (2). Thus, 0x401921FB54442D18
	is the high half of the qfloat approximation to 2*pi. The next 64 bits have decimal form

	5086983782422027140

	which in hex is

	469898CC51701B84 .

	This (left-zero-padded and appended to "Ox") is the hex representation of the lower half of our qfloat.
*/

const struct qfloat QZRO    = {0x0000000000000000ull, 0x0000000000000000ull};
const struct qfloat QONE    = {0x3FF0000000000000ull, 0x0000000000000000ull};
const struct qfloat QEPS    = {0x3E40000000000000ull, 0x0000000000000000ull};	/* Set QEPS to 2^-108 */
const struct qfloat QTWO    = {0x4000000000000000ull, 0x0000000000000000ull};
const struct qfloat QTHREE  = {0x4008000000000000ull, 0x0000000000000000ull};
const struct qfloat Q2PI    = {0x401921FB54442D18ull, 0x469898CC51701B84ull};
const struct qfloat QPI     = {0x400921FB54442D18ull, 0x469898CC51701B84ull};
const struct qfloat QPIHALF = {0x3FF921FB54442D18ull, 0x469898CC51701B84ull};
const struct qfloat QLN2    = {0x3FE62E42FEFA39EFull, 0x35793C7673007E5Full};
const struct qfloat QEXP    = {0x4005BF0A8B145769ull, 0x5355FB8AC404E7A8ull};
const struct qfloat QSQRT2  = {0x3FF6A09E667F3BCCull, 0x908B2FB1366EA958ull};


#if QFONLY	/* If set QFDEBUG to true, must compile in standalone mode. */
int main()
{
	return qtest();
}
#endif

int qtest(void)
{
	double c, d;
	int64 hidiff, lodiff;
	struct qfloat p, q ,r;
	double pi = 3.1415926535897932384;

#ifdef MUL_LOHI64_SUBROUTINE
	printf("INFO: qfloat routines using subroutine form of MUL_LOHI\n");
#endif

#if QFDEBUG
		printf("ABS(0x1234567890ABCDEF) gives %16llX\n",ABS((int64)0x1234567890ABCDEFull));
#endif
	if(!(ABS((int64)0x1234567890ABCDEFull) == 0x1234567890ABCDEFull)) ASSERT(HERE, 0,"ERROR 10 in qfloat.c");

	/*********** TEST THE TYPE CONVERSIONS **************/
	c = 0.0;	d = qfdbl(QZRO);
#if QFDEBUG
		printf("dble(0.0) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(hidiff == (int64)0)) ASSERT(HERE, 0,"ERROR 12 in qfloat.c");

	c = 1.0;	d = qfdbl(QONE);
#if QFDEBUG
		printf("dble(1.0) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(hidiff == (int64)0)) ASSERT(HERE, 0,"ERROR 14 in qfloat.c");

	c = 2.0;	d = qfdbl(QTWO);
#if QFDEBUG
		printf("dble(2.0) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(hidiff == (int64)0)) ASSERT(HERE, 0,"ERROR 16 in qfloat.c");

	c =-2.0;	d = qfdbl(qfneg(QTWO));
#if QFDEBUG
		printf("dble(-2.0)= %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(hidiff == (int64)0)) ASSERT(HERE, 0,"ERROR 18 in qfloat.c");

	c = 2*pi;	d = qfdbl(Q2PI);
#if QFDEBUG
		printf("dble(2pi) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(ABS(hidiff) < (int64)2)) ASSERT(HERE, 0,"ERROR 20 in qfloat.c");

	c =log(2.0);d = qfdbl(QLN2);
#if QFDEBUG
		printf("dble(ln2) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(ABS(hidiff) < (int64)2)) ASSERT(HERE, 0,"ERROR 22 in qfloat.c");

	c = exp(1.0);
	d = qfdbl(QEXP);
#if QFDEBUG
		printf("dble(exp) = %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(ABS(hidiff) < (int64)2)) ASSERT(HERE, 0,"ERROR 24 in qfloat.c");

	c = -c;		d = qfdbl(qfneg(QEXP));
#if QFDEBUG
		printf("dble(-exp)= %16llX  %16llX\n",*(int64 *)&c, *(int64 *)&d);
#endif
	hidiff = *(int64 *)&c - *(int64 *)&d;	if(!(ABS(hidiff) < (int64)2)) ASSERT(HERE, 0,"ERROR 26 in qfloat.c");

	/*********** TEST THE MULTIPLY ALGORITHM ************/
	/* e*e: 	0x401D8E64B8D4DDAD, 0xCC33A3BA206B68AC	*/
	q = qfmul(QEXP,QEXP);
#if QFDEBUG
		printf("      e*e  = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble( e*e) = %25.16e\n",qfdbl(q));
#endif
	hidiff = (int64)q.hi - (int64)0x401D8E64B8D4DDADull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 28 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0xCC33A3BA206B68ACull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 29 in qfloat.c");

	/* ln2*e:	0x3FFE258ECC242F82, 0x5DEC567E6A0E1111	*/
	q = qfmul(QLN2,QEXP);
#if QFDEBUG
		printf("     L2*e  = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble(L2*e) = %25.16e\n",qfdbl(q));
#endif
	hidiff = (int64)q.hi - (int64)0x3FFE258ECC242F82ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 30 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x5DEC567E6A0E1111ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 31 in qfloat.c");

	/* ln2^2:	0x3FDEBFBDFF82C58E, 0xA86F16B06EC97360	*/
	q = qfmul(QLN2,QLN2);
#if QFDEBUG
		printf("     L2^2  = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble(L2^2) = %25.16e\n",qfdbl(q));
#endif
	hidiff = (int64)q.hi - (int64)0x3FDEBFBDFF82C58Eull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 32 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0xA86F16B06EC97360ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 33 in qfloat.c");

	/* ln2*2pi:	0x40116BB24190A0B6, 0xE765BE0D06135E60	*/
	q = qfmul(QLN2,Q2PI);
#if QFDEBUG
		printf("     L2*pi = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble(L2*pi)= %25.16e\n",qfdbl(q));
#endif
	hidiff = (int64)q.hi - (int64)0x40116BB24190A0B6ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 34 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0xE765BE0D06135E60ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 35 in qfloat.c");

	/* 2pi*e:	0x403114580B45D474, 0x9E6108579A2D0CA7	*/
	q = qfmul(Q2PI,QEXP);
#if QFDEBUG
		printf("     pi*e  = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble(pi*e) = %25.16e\n",qfdbl(q));
#endif
	hidiff = (int64)q.hi - (int64)0x403114580B45D474ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 36 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x9E6108579A2D0CA7ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 37 in qfloat.c");

	/* 2pi*2pi:	0x4043BD3CC9BE45DE, 0x5A4ADC4D9B301183	*/
	q = qfmul(Q2PI,Q2PI);
#if QFDEBUG
		printf("  (2*pi)^2 = %16llX  %16llX\n",q.hi,q.lo);
		printf("dble(2pi^2)= %25.16e\n",qfdbl(q));
		printf("dble(2pi^2)= %25.16e\n",pi*pi);
#endif
	hidiff = (int64)q.hi - (int64)0x4043BD3CC9BE45DEull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 38 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x5A4ADC4D9B301183ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 39 in qfloat.c");

	/*********** TEST THE ADDITION ALGORITHM ************/
	/* 2*pi+e:	0x402200C04CE72C66, 0x7821CB48D9B947AC	*/
	q = qfadd(QEXP,Q2PI);
#if QFDEBUG
		printf("  2*pi + e = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x402200C04CE72C66ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 40 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x7821CB48D9B947ACull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 41 in qfloat.c");

	/********** TEST THE SUBTRACTION ALGORITHM **********/
	/* Both inputs   normalized, output   normalized, with just one significant bit. */
	q.hi = 0x3FEFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFFull;
	q = qfsub(q, q);
#if QFDEBUG
		printf("result1 = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 42 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0ull;	if(!(ABS(lodiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 43 in qfloat.c");

	p.hi = 0x3FEFFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x3FEFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result2 = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x38A0000000000000ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 44 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x0000000000000000ull;	if(!(ABS(lodiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 45 in qfloat.c");

	/* Both inputs   normalized, output denormalized, with just one significant bit. */
	p.hi = 0x00FFFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x00FFFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result3 = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x0000000000000000ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 46 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x0000000000004000ull;	if(!(ABS(lodiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 47 in qfloat.c");

	/* Both inputs denormalized, output denormalized, with just one significant bit. */
	q.hi = 0x000FFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFFull;
	q = qfsub(q, q);
#if QFDEBUG
		printf("result4 = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 48 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0ull;	if(!(ABS(lodiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 49 in qfloat.c");

	p.hi = 0x000FFFFFFFFFFFFFull;	p.lo = 0xFFFFFFFFFFFFFFFFull;
	q.hi = 0x000FFFFFFFFFFFFFull;	q.lo = 0xFFFFFFFFFFFFFFFEull;
	q = qfsub(p, q);
#if QFDEBUG
		printf("result5 = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 50 in qfloat.c");
	lodiff = (int64)q.lo - (int64)1ull;	if(!(ABS(lodiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 51 in qfloat.c");

	/* 2*pi-e:	0x400C84EC1D7402C7, 0x39DB360DDEDB4F60	*/
	q = qfsub(Q2PI,QEXP);
#if QFDEBUG
		printf("    2pi- e = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x400C84EC1D7402C7ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 52 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x39DB360DDEDB4F60ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 53 in qfloat.c");

	/* e-2*pi:	0xC00C84EC1D7402C7, 0x39DB360DDEDB4F60	*/
	r = qfsub(QEXP,Q2PI);
#if QFDEBUG
		printf("     e-2pi = %16llX  %16llX\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, qfneg(q)))) ASSERT(HERE, 0,"ERROR 54 in qfloat.c");

	/*********** TEST THE SQUARE ROOT ALGORITHM ************/
	/* sqrt(2):	0x3FF6A09E667F3BCC, 0x908B2FB1366EA958, qfsqrt gives ...956. */
	q = qfsqrt(QTWO);
#if QFDEBUG
		printf("sqrt(2) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FF6A09E667F3BCCull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 56 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x908B2FB1366EA958ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 57 in qfloat.c");

	/*********** TEST THE INVERSION AND DIVISION ALGORITHMS ************/
	/* 1/(2*pi):0x3FC45F306DC9C882, 0xA53F84EAFA3EA69B(B81B...), qfinv gives ...698. */
	q = qfinv(Q2PI);
#if QFDEBUG
		printf("1/(2*pi) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FC45F306DC9C882ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 58 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0xA53F84EAFA3EA69Bull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 59 in qfloat.c");

	/* 1/e:		0x3FD78B56362CEF37, 0xC6AEB7B1E0A4153E(4376...), qfinv gives ...53C. */
	q = qfinv(QEXP);
#if QFDEBUG
		printf("1/e      = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FD78B56362CEF37ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 60 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0xC6AEB7B1E0A4153Eull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 61 in qfloat.c");

	/* 1/ln2:	0x3FF71547652B82FE, 0x1777D0FFDA0D23A7(D11D...), qfinv gives ...3A6. */
	q = qfinv(QLN2);
#if QFDEBUG
		printf("1/ln(2)  = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FF71547652B82FEull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 62 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x1777D0FFDA0D23A7ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 63 in qfloat.c");

	/* 2*pi/ln2:0x40222123045B5DEB, 0x9C5398CE82C06E4B(80DB...), qfdiv gives ...E4A. */
	q = qfdiv(Q2PI, QLN2);
#if QFDEBUG
		printf("2*pi/ln(2)  = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x40222123045B5DEBull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 64 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x9C5398CE82C06E4Bull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 65 in qfloat.c");

	/*********** TEST THE EXPONENTIAL AND TRIG FUNCTIONS ************/
	/* exp(1):	0x4005BF0A8B145769, 0x5355FB8AC404E7A7(9E3B...), qfexp gives ...4E7A7, ~116 bits of accuracy. */
	q = qfexp(QONE);
#if QFDEBUG
		printf("qexp(1) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x4005BF0A8B145769ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 66 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x5355FB8AC404E7A7ull;	if(!(ABS(lodiff) < (int64)16)) ASSERT(HERE, 0,"ERROR 67 in qfloat.c");

	/* Sine and cosine are somewhat roundoff-error prone, so raise the error limit slightly. */
	/* cos(1):	0x3FE14A280FB5068B, 0x923848CDB2ED0E37(A534...), qfcs1 gives ...D0E38, ~116 bits of accuracy */
	q = qfcs1(QONE);
#if QFDEBUG
		printf("qcs1(1) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FE14A280FB5068Bull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 68 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x923848CDB2ED0E37ull;	if(!(ABS(lodiff) < (int64)64)) ASSERT(HERE, 0,"ERROR 69 in qfloat.c");
	r = qfcos(QONE);
#if QFDEBUG
		printf("qcos(1) = %16llX  %16llX\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, q))) ASSERT(HERE, 0,"ERROR 70 in qfloat.c");

	/* sin(1):	0x3FEAED548F090CEE, 0x0418DD3D2138A1E7(8651...), qfsn1 gives ...8A1E9, ~115 bits of accuracy */
	q = qfsn1(QONE);
	hidiff = (int64)q.hi - (int64)0x3FEAED548F090CEEull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 72 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x0418DD3D2138A1E7ull;	if(!(ABS(lodiff) < (int64)64)) ASSERT(HERE, 0,"ERROR 73 in qfloat.c");
#if QFDEBUG
		printf("qsn1(1) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	r = qfsin(QONE);
#if QFDEBUG
		printf("qsin(1) = %16llX  %16llX\n",r.hi,r.lo);
#endif
	if(!(qfcmpeq(r, q))) ASSERT(HERE, 0,"ERROR 74 in qfloat.c");

	/* cos(100):0x3FEB981DBF665FDF, 0x63F433736617A041(5D8A...), qfcos gives ...7A023, ~114 bits of accuracy */
	q = qfcos(i64_to_q((int64)100));
#if QFDEBUG
		printf("qcos(100) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0x3FEB981DBF665FDFull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 75 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x63F433736617A041ull;	if(!(ABS(lodiff) < (int64)64)) ASSERT(HERE, 0,"ERROR 76 in qfloat.c");

	/* sin(100):0xBFE03425B78C4DB8, 0x0708F6155D083EB2(1C6B...), qfsin gives ...83EE5, ~109 bits of accuracy */
	q = qfsin(i64_to_q((int64)100));
#if QFDEBUG
		printf("qsin(100) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	hidiff = (int64)q.hi - (int64)0xBFE03425B78C4DB8ull;	if(!(ABS(hidiff) ==(int64) 0)) ASSERT(HERE, 0,"ERROR 77 in qfloat.c");
	lodiff = (int64)q.lo - (int64)0x0708F6155D083EB2ull;	if(!(ABS(lodiff) < (int64)64)) ASSERT(HERE, 0,"ERROR 78 in qfloat.c");

	/*********** TEST THE ROUND-TOWARD-ZERO AND ROUND-TO-NEAREST FUNCTIONS ************/
	q = qfmul_pow2(QONE, -1);
	q = qfnint(q);
#if QFDEBUG
		printf("qfnint(0.5) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)1)) ASSERT(HERE, 0,"ERROR 80 in qfloat.c");

	q = qfmul_pow2(QONE, -1);
	q = qfint(q);
#if QFDEBUG
		printf("qfint(0.5) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)0)) ASSERT(HERE, 0,"ERROR 82 in qfloat.c");

	q = qfnint(QEXP);
#if QFDEBUG
		printf("qfnint(exp) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)3)) ASSERT(HERE, 0,"ERROR 84 in qfloat.c");

	q = qfint(QEXP);
#if QFDEBUG
		printf("qfint(exp) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)2)) ASSERT(HERE, 0,"ERROR 86 in qfloat.c");

	q = qfnint(Q2PI);
#if QFDEBUG
		printf("qfnint(2pi) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)6)) ASSERT(HERE, 0,"ERROR 88 in qfloat.c");

	q = qfmul_pow2(Q2PI, 20);
	q = qfnint(q);
#if QFDEBUG
		printf("qfnint(2pi*2^20) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(!q.hi && q.lo == (uint64)6588397)) ASSERT(HERE, 0,"ERROR 90 in qfloat.c");

	q = qfmul_pow2(Q2PI, 124);	/* This gives pi*2^125, which should still fit into a signed 128-bit int. */
	q = qfnint(q);
#if QFDEBUG
		printf("qfnint(2pi*2^125) = %16llX  %16llX\n",q.hi,q.lo);
#endif
	if(!(q.hi = (uint64)0x6487ED5110B4611Aull && q.lo == (uint64)0x62633145C06E1000ull)) ASSERT(HERE, 0,"ERROR 92 in qfloat.c");

	return 0;
}


/* qfloat negation. */
struct qfloat qfneg(struct qfloat q)
{
	q.hi ^= MASK_SIGN;
	return q;
}


/* qfloat absolute value. */
struct qfloat qfabs(struct qfloat q)
{
	q.hi &= ~MASK_SIGN;
	return q;
}


/* qfloat comparison: q1 == q2. 0 and -0 are treated as unequal, so one must be careful to distinguish between the two. */
uint32 qfcmpeq(struct qfloat q1, struct qfloat q2)
{
	if(q1.hi == q2.hi && q1.lo == q2.lo)
		return 1;
	else
		return 0;
}


/* qfloat comparison: q1 < q2 */
uint32 qfcmplt(struct qfloat q1, struct qfloat q2)
{
	uint32 both_signs = (uint32)(((q1.hi & MASK_SIGN) >> 62) + ((q2.hi & MASK_SIGN) >> 63));

	switch(both_signs)
	{
	case(0) :	/* Both q1 and q2 nonnegative */
		if(q1.hi < q2.hi || (q1.hi == q2.hi && q1.lo < q2.lo))
			return 1;
		else
			return 0;
	case(1) :	/* q1 >= 0, q2 < 0 */
		return 0;
	case(2) :	/* q1 < 0, q2 >= 0 */
		return 1;
	case(3) :	/* Both q1 and q2 negative, in which case a more-negative q1 looks larger w.r.to the unsigned compare */
		if(q1.hi > q2.hi || (q1.hi == q2.hi && q1.lo > q2.lo))
			return 1;
		else
			return 0;
	default:
		ASSERT(HERE, 0,"ERROR 98 in qfloat.c");
	}
	return 0;	/* Just to get the compiler to shut up ... this should never be reached. */
}


/* qfloat comparison: q1 <= q2 */
uint32 qfcmple(struct qfloat q1, struct qfloat q2)
{
	return (qfcmplt(q1, q2) | qfcmpeq(q1, q2));
}


/*
Return IEEE64-compliant floating double approximation to a qfloat.
Since the high part of a qfloat is already in IEEE double form, we only
need to round in the high bit of the low part and cast the result to a double.
*/
double qfdbl(struct qfloat q)
{
	uint64 hi;

	hi   = q.hi + (q.lo >> 63);
	/* If adding the carry rippled into the bottom bit of the exponent field, need to right-shift the significand one place.
	Note that this can only happen if bits <0:52> of hi were all flipped to 0, so don't need to worry about rounding the bit that gets shifted off.
	In fact, we don't even need to explicitly restore any hidden bit or do an actual right-shift while doing this, since:

	a) If number was normalized, the 1 that got carried out of the bottom 52 bits gets added to the exponent, which is what we wanted to do anyway;

	b) If number was denormalized, MSB of significand now gets treated as exponent field of 1, which again is what we wanted to do anyway.
	*/
	return *(double *)&hi;
}


/* Convert IEEE64 double to qfloat. */
struct qfloat dbl_to_q(double d)
{
	struct qfloat q;

	q.hi = *(uint64 *)&d;	/* Copy bit pattern of double into a uint64. */
	q.lo = (uint64)0;
	return q;
}


/* Convert 64-bit signed int to qfloat. */
struct qfloat i64_to_q(int64 i64)
{
	struct qfloat q;
	int32 lz, shift;
	uint64 sexp;

	if(!i64)
	{
		q.hi = (uint64)0;
		q.lo = (uint64)0;
		return q;
	}

	sexp = i64 & MASK_SIGN;
	if(sexp) i64 = -i64;
	lz = leadz64(i64);	/* Find # of leading zeros in |i64|. */

	/* Put leading ones bit of mantissa into bit position 52 of <0:63> and add to sign/exponent. */
	shift = lz - (int)11;	/* if lz = 11, int is already correctly positioned. If <> 11,
							i.e. shift <> 0, need to right/left-shift mantissa */

	sexp += ((uint64)1074 - shift) << 52;	/* Ex: i64 = 3, with lz = 62 and shift = 51, should yield exponent = 0x400 = 1024
											(prior to left-shift by 52 places); each bit less in lz means one more in exponent.
											Since the leftmost mantissa bit = 1 and gets added to the exponent, subtract an
											extra one, i.e. 1075 - 1 - shift = 1074 - shift. */

	if(shift < 0)
	{
		q.hi = sexp + (i64 >> (-shift));
		q.lo = i64 << (64+shift);
	}
	else	/* need to left-shift mantissa */
	{
		q.hi = sexp + (i64 << shift);
		q.lo = (uint64)0;
	}

	return q;
}


/* Multiply by an integer power of 2 is especially simple: */
struct qfloat qfmul_pow2(struct qfloat q, int32 pow)
{
	uint64 exp;
	struct qfloat return_val;

	/* Extract 11-bit exponent field and add sign-extended power-of-2 exponent: */

	exp = ((q.hi >> 52) & MASK_EXP);
if(!(exp)) ASSERT(HERE, 0,"ERROR 100 in qfloat.c");	/* Denormalized numbers not currently supported in this function. */
	exp += (int64)pow;	/* Cast pow to signed 64-bit int so it gets sign-extended properly. */

	/* Make sure the new exponent didn't over-or-underflow: */

	if(!((int64)exp > 0 && !(exp & ~MASK_EXP))) ASSERT(HERE, 0,"ERROR 110 in qfloat.c");

	/* Insert new exponent field into high part of q and return. */

	return_val.hi = (q.hi & ~(MASK_EXP << 52)) + (exp << 52);
	return_val.lo =  q.lo;

	return(return_val);
}


/* Nearest integer (round-to-nearest) of a qfloat. The result is stored in a qfloat,
but in this case the two halves should both be interpreted as 64-bit ints, according to:

	- if q.hi >= 0, result =  q.hi*2^64 + q.lo;
	- if q.hi <  0, result = -q.hi*2^64 - q.lo.

The allowable range of the result is [-2^127, +2^127-1].
*/
struct qfloat qfnint(struct qfloat q)
{
	int32 exp, sign, rshift, lshift;
	uint64 offword, carry;

	/* Separate upper part of the significand from the sign/exponent fields: */
	sign = (int32)(q.hi >> 63);
	exp  = (int32)(q.hi >> 52) & MASK_EXP;
	q.hi =         q.hi        & MASK_MANT;

	/* If number is normalized, restore hidden bit.
	Otherwise, number is so small that we can immediately return zero. */
	if(exp)
	{
		q.hi  += TWOE52;
	}
	else
	{
		return QZRO;
	}

	/* Get right-shift count. E.g. 1 has exp = 0x3FF and rshift = 116, 2 has exp = 0x400 and rshift = 115,
	so rshift = 0x400 + 115 - exp = 0x473 - exp.
	*/
	rshift = 0x473 - exp;
	carry = (uint64)0;

	if(rshift <= 0)		/* If rshift <= 0, require rshift strictly > -11, unless rshift = -11 and result = -2^127. */
	{
		if(rshift <= -11)
		{
			if(rshift == -11 && (!sign || (q.hi << -rshift) != MASK_SIGN || q.lo != (uint64)0))
			{
				ASSERT(HERE, 0,"ERROR: qfloat is too large to convert to 128-bit integer.");
			}
		}
		lshift =     - rshift;
		rshift = (64 + rshift) & 63;	/* shift count of 64 must get aliased to 0. */
		q.hi = (q.hi << lshift) + (q.lo >> rshift);
		q.lo = (q.lo << lshift);
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		offword = (q.lo << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
		q.lo = (q.lo >> rshift) + (q.hi << lshift);
		q.lo += offword;
		if(q.lo < offword)	/* Had a carry from the rounded bit being added into the right-shifted low part. */
		{
			carry = (uint64)1;
		}
		q.hi = (q.hi >> rshift) + carry;
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		offword = (q.lo << lshift) >> 63;
		q.lo = (q.lo >> rshift) + (q.hi << lshift) + offword;
		q.hi = (uint64)0;
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = q.lo >> 63;
		q.lo = q.hi + offword;
		q.hi = (uint64)0;
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift = 64 - rshift;
		offword = (q.hi << lshift) >> 63;
		q.lo = (q.hi >> rshift) + offword;
		q.hi = (uint64)0;
	}
	else			/* Entire significand shifted off, no need to round result. */
	{
		return QZRO;
	}

	/* If sign bit was set, negate high part of result. */
	if(sign)
	{
		q.hi = -q.hi;
	}

	return q;
}


/* Cast to integer (round-toward-zero) of a qfloat. The result is stored in a qfloat,
but in this case the two halves should both be interpreted as 64-bit ints, according to:

	- if q.hi >= 0, result =  q.hi*2^64 + q.lo;
	- if q.hi <  0, result = -q.hi*2^64 - q.lo.

The allowable range of the result is [-2^127, +2^127-1].
*/
struct qfloat qfint(struct qfloat q)
{
	int32 exp, sign, rshift, lshift;

	/* Separate upper part of the significand from the sign/exponent fields: */
	sign = (int32)(q.hi >> 63);
	exp  = (int32)(q.hi >> 52) & MASK_EXP;
	q.hi =         q.hi        & MASK_MANT;

	/* If number is normalized, restore hidden bit.
	Otherwise, number is so small that we can immediately return zero. */
	if(exp)
	{
		q.hi  += TWOE52;
	}
	else
	{
		return QZRO;
	}

	/* Get right-shift count. E.g. 1 has exp = 0x3FF and rshift = 116, 2 has exp = 0x400 and rshift = 115,
	so rshift = 0x400 + 115 - exp = 0x473 - exp.
	*/
	rshift = 0x473 - exp;

	if(rshift <= 0)		/* If rshift <= 0, require rshift strictly > -11, unless rshift = -11 and result = -2^127. */
	{
		if(rshift <= -11)
		{
			if(rshift == -11 && (!sign || (q.hi << -rshift) != MASK_SIGN || q.lo != (uint64)0))
			{
				ASSERT(HERE, 0,"ERROR: qfloat is too large to convert to 128-bit integer.");
			}
		}
		lshift =     - rshift;
		rshift = (64 + rshift) & 63;	/* shift count of 64 must get aliased to 0. */
		q.hi = (q.hi << lshift) + (q.lo >> rshift);
		q.lo = (q.lo << lshift);
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		q.lo = (q.lo >> rshift) + (q.hi << lshift);
		q.hi = (q.hi >> rshift);
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
	{
		lshift = 64 - rshift;
		q.lo = (q.lo >> rshift) + (q.hi << lshift);
		q.hi = (uint64)0;
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into q.hi, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		q.lo = q.hi;
		q.hi = (uint64)0;
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift = 64 - rshift;
		q.lo = (q.hi >> rshift);
		q.hi = (uint64)0;
	}
	else			/* Entire significand shifted off, no need to round result. */
	{
		return QZRO;
	}

	/* If sign bit was set, negate high part of result. */
	if(sign)
	{
		q.hi = -q.hi;
	}

	return q;
}


/*
Qfloat multiply algorithm. Needs only IEEE64-floating-compliant multiply and a macro to generate the
128-bit product of two unsigned 64-bit ints.

Rearrange the inputs as x = b + a*2^53, y = d + c*2^53, where a and c are 53 bits long,
and b and d are 64 bits long, in contrast with the obvious approach, which keeps 64 bits into the low part
of the 2-word representation and 53 bits into the high part. The multiply then looks like

        x*y = a*c       *2^106   (upper 128 bits)
            +(a*d + c*b)*2^53    (next-lower 117 bits, right-shifted 53 bits w.r.to a*c)
            + b*d                (lowest-order 106 bits, right-shifted 106 bits w.r.to a*c).

Notice what happens: superficially this looks like the same amount of work to calculate,
but hold on - since the middle (a*d + c*b) term of the above product is right-shifted 53 bits w.r.to
the leading order (a*c) term, only the upper 64 bits of a*d and c*b (considered separately) overlap
the a*c term, i.e. they align precisely with the lower half of a*c. Even better, since we only need
the upper 53 bits of these 3 64-bit aligned terms (i.e. MULL64(a*c), MULH64(a*d) and MULH64(c*b)) and since
MULH64(a*d) and MULH64(c*b) have only 53 bits to begin with, we can simply use an FMUL to emulate these
two product operations. Unless we absolutely need a result in which all 117 bits are correct, we don't
even care if the bottom bit of the FMUL result is wrong, i.e. need no expensive and messy error correction
operations to correct this bit. Our total cost is now just one 128-bit integer product (b*d), 2 FMUL,
and 2 adds with a potential carry into the upper 64 bits of a*c. The latter can be simplified by first
rounding MULL64(a*c) to its upper 53 bits, adding in FMUL(a*d) and FMUL(c*b) (properly scaled), and
combining the result with the lower 11 bits of MULH64(a*c). Sweet!

TOTAL COST (uisng Alpha multiply model): 6 FLTCVT, 4 FMUL, 2 IMUL + 50-70 IOPS (shift, add, etc.)
This is a lot of total operations, but the most-expensive operation, integer multiply, has been drastically
reduced in number from the all-integer model. Assuming a typical RISC multiple-instruction-issue model and
two 64-bit integer logic units which can operate in parallel on simple IOPS, we estimate 30-40 cycles per
qfloat multiply. That's a lot slower than the one cycle for multiplying two doubles, but still fast enough
to be useful. Also, were we to use an all-floating implementation, we'd need roughly a length-16 floating
convolution, followed by a carry step, and a bunch of floating-int conversions, which I estimate would be
~4x slower. (The latter approach is of course asymptotically better since we can use FFT, but that doesn't
begin to pay off until runlengths become much larger than the # of bits in a qfloat.)

IMPORTANT NOTE: since the input mantissas are in [2^116, 2^117)
(and this neglects the possibility of denormalized input operands, which requires additional special casing),
the product is in [2^232, 2^234), i.e. may have either 233 or 234 bits. We detect which by checking whether the
result of the MUL_LOHI64(b, d) has 127 or 128 bits, and then left-justify subsequent operations accordingly,
in order not to leave the MSB of the result mantissa vacant.


Example: in our qfloat format, the 117-bit approximation to e looks like

2^2 * 0.10101101111110000101010001011000101000101011101101001 0101001101010101111110111000101011000100000001001110011110101000

Rearranging the significand as b + a*2^53 gives

a = 1010110111111000010101000101100010100010101110110100101010011010 = 12535862302449814170
b = 10101111110111000101011000100000001001110011110101000            = 6187547923638184

Squaring, we get

a*a.hi = 8519001675203524399
a*a.lo = 16107666474326910116
a*b.hi = 4204874770886292	(rightmost digit should be 3 after rounding in low part).

Writing these in binary and aligning according to their relative power-of-2 multipliers:
																x
a*a*2^53 = 0111011000111001100100101110001101010011011101101011011100101111 1101111110001001111011010101000011111010101110010110010010100100
a*b      =                                                                  01110111100000100111110110011000010110010111010010100


*/
struct qfloat qfmul(struct qfloat q1, struct qfloat q2)
{
	const double fmul_scale = 1.0/((double)0x00400000*0x80000000);	/* 2^(-53) */
	struct qfloat return_val;
	uint64 sexp1, sexp2, hi1, hi2, lo1, lo2;
	uint64 a,b,c,d,lo,hi,bvac;
	double da,db,dc,dd,adhi,bchi;

	/* If either operand is zero, return zero. We mask off the sign bit so -0 is treated like 0. */
	if((!(q1.hi & ~MASK_SIGN) && !q1.lo) || (!(q2.hi & ~MASK_SIGN) && !q2.lo))
	{
		return QZRO;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */

	hi1 = q1.hi & MASK_MANT;
	hi2 = q2.hi & MASK_MANT;

	sexp1 = q1.hi - hi1;
	sexp2 = q2.hi - hi2;
	/* If number is normalized, restore hidden bit. */
	if(sexp1 & ~MASK_SIGN)
	{
		hi1 += TWOE52;
	}
	else
	{
		printf("Multiply by denormalized operand not supported!");
		ASSERT(HERE, 0,"ERROR in qfloat.c : qfmul");
	}

	if(sexp2 & ~MASK_SIGN)
	{
		hi2 += TWOE52;
	}
	else
	{
		printf("Multiply by denormalized operand not supported!");
		ASSERT(HERE, 0,"ERROR in qfloat.c : qfmul");
	}

	lo1 = q1.lo;
	lo2 = q2.lo;

	/* Need to subtract scaling factor prior to adding sexp1 and sexp2 together, so don't overflow into the sign bit. */
	sexp1 -= 0x3FF0000000000000ull;	/* Hex constant = ((uint64)1023 << 52) */

	a = (hi1<<11) + (lo1>>53);
	c = (hi2<<11) + (lo2>>53);

#if DEBUG_QFMUL
	printf("A,C = %20llu  %20llu\n", a, c);
#endif

#ifdef MUL_LOHI64_SUBROUTINE
	MUL_LOHI64(a,c,&lo,&hi);	/* start the integer product ASAP, since it will generally have a large latency. */
#else
	MUL_LOHI64(a,c,lo,hi);
#endif

#if DEBUG_QFMUL
	printf("L,H = %20llu  %20llu\n",lo,hi);
#endif

	b = (lo1 << 11) >> 11;
	d = (lo2 << 11) >> 11;

#if DEBUG_QFMUL
	printf("B,D = %20llu  %20llu\n", b, d);
#endif

	/* Make these as accurate as possible by rounding rather than truncating.
	If operand overflows into bit position 53 (of <0:63>) on rounding, the high part was = 2^53 - 1
	and is now = 2^53, which can be exactly represented as a double, i.e. we lose no precision if this occurs.
	*/
	da = (double)(hi1 + (lo1 >> 63)) * fmul_scale;
	dc = (double)(hi2 + (lo2 >> 63)) * fmul_scale;

#if DEBUG_QFMUL
	printf("da,dc = %20.15lf %20.15lf\n",da,dc);
#endif

	db = (double)(b);
	dd = (double)(d);

#if DEBUG_QFMUL
	printf("db,dd = %20.15lf %20.15lf\n",db,dd);
#endif

/* DEBUG: make sure we didn't lose any bits of b or d during the conversion to float. */
#if QFDEBUG
	if(!((uint64)db == b)) ASSERT(HERE, 0,"ERROR 120 in qfloat.c");
	if(!((uint64)dd == d)) ASSERT(HERE, 0,"ERROR 121 in qfloat.c");
#endif

	adhi = da*dd;
	bchi = db*dc;

#if DEBUG_QFMUL
	printf("ad,bc = %20.15lf %20.15lf\n",adhi,bchi);
#endif

	bvac = (uint64)leadz64(hi);
	if(!(bvac < 2)) ASSERT(HERE, 0,"ERROR 130 in qfloat.c");	/* Make sure at most the leftmost bit of high part is vacant. This check
						needs to remain in place until support for denormalized oprands is added. */
#if DEBUG_QFMUL
	printf("bvac = %20llu\n",bvac);
#endif

	/*
	printf("a*d_hi = %16llX = %20llu\n", (uint64)adhi, (uint64)adhi );
	printf("c*b_hi = %16llX = %20llu\n", (uint64)bchi, (uint64)bchi );
	printf("b*d_hi = %16llX = %20llu\n", hi, hi);
	printf("b*d_lo = %16llX = %20llu\n", lo, lo);
	*/
	/*
	Now need to right-shift MUL_LOHI result (12-bvac) places, FMUL results (1-bvac) place, and add together.
	(We could contemplate FADDing the 2 FMUL results together, which saves a float->int conversion
	but risks losing a bit of precision at the lower end. If FLTCVT is slow, it's probaby a good trade.)
	First add the 2 FMUL results to (LO >> (11-bvac)), which preserves an extra low bit, then right-shift
	one place, round and add to the lower (12-bvac) bits of HI.
	*/
	/*
	lo = (hi << 52) + (lo >> 12) + (lo & 1) + (uint64)(adhi + bchi);
	hi = (hi >> 12);
	printf("mul_hi = %16llX = %20llu\n", hi, hi);
	printf("mul_lo = %16llX = %20llu\n", lo, lo);
	*/

	return_val.hi = (hi >> (11-bvac));
if(!((return_val.hi >> 52) == (uint64)1)) ASSERT(HERE, 0,"ERROR 140 in qfloat.c");
	return_val.lo = (hi << (53+bvac)) + (lo >> (11-bvac)) + ((lo >> (10-bvac)) & (uint64)1) + (((uint64)adhi + (uint64)bchi) << bvac);
                                                            /* ^^^^rounding is here^^^^^ */   /* Maximize accuracy by converting to int prior to add. */
	if(return_val.lo < (hi << (53+bvac)))	/* had a carry out of lo part. */
	{
		return_val.hi++;

		/* If adding the carry rippled into the bottom bit of the exponent field, need to right-shift the significand one place.
		(This actually happens quite frequently, e.g. in qfinv, when we multiply x by its estimated inverse, and the product may
		be extremely close to unity, which may show up as return_val.hi = 0x001FFFFFFFFFFFFF prior to adding carrylo, and then
		we get 0x0020000000000000 after. Note that we don't actually need to right-shift the hi part in this case, since the
		bottom bits are all zero and the carry needs to be added to the exponent field to account for the right-shift, anyway.
		*/
		if(return_val.hi == 0x0020000000000000ull)	/* 	Note that this can only happen if bits <0:52> of hi were all flipped to 0... */
		{
			return_val.lo = (return_val.lo >> 1) + (return_val.lo & 1);	/* ...so don't need to worry about carrying in the lo bit of hi here. */
		}
		else
		{
if(!((return_val.hi >> 52) == (uint64)1)) ASSERT(HERE, 0,"ERROR 150 in qfloat.c");
		}
	}

	return_val.hi = (sexp1 + sexp2 - (bvac << 52)) + return_val.hi;

	return(return_val);
}


/*
Top-level add and subtract routines seen by the caller - these examine the
signs of the inputs, send the proper combination of +-q1 and +-q2 to the
low-level qfsum and qfdif functions, and properly sign the result.
*/
struct qfloat qfadd	(struct qfloat q1, struct qfloat q2)
{
	struct qfloat q;
	uint32 both_signs = (uint32)(((q1.hi & MASK_SIGN) >> 62) + ((q2.hi & MASK_SIGN) >> 63));

	if(both_signs == 3)		/* Both inputs negative */
	{
		q1.hi ^= MASK_SIGN;	/* Send -q1 and -q2 to low-level add routine... */
		q2.hi ^= MASK_SIGN;
		q = qfsum(q1, q2); q.hi ^= MASK_SIGN;	/* ...and return negated sum. */
	}
	else if(both_signs == 2)	/* q1 negative */
	{
		q1.hi ^= MASK_SIGN; /* Send q2 and -q1 to low-level sub routine. */
		q = qfdif(q2, q1);
	}
	else if(both_signs == 1)	/* q2 negative */
	{
		q2.hi ^= MASK_SIGN; /* Send q1 and -q2 to low-level sub routine. */
		q = qfdif(q1, q2);
	}
	else if(both_signs == 0)	/* Both inputs positive */
	{
		q = qfsum(q1, q2);
	}
	else
	{
		ASSERT(HERE, 0,"ERROR: unrecognized sign combination in QFADD");
	}

	return q;
}


struct qfloat qfsub	(struct qfloat q1, struct qfloat q2)
{
	struct qfloat q;
	uint32 both_signs = (uint32)(((q1.hi & MASK_SIGN) >> 62) + ((q2.hi & MASK_SIGN) >> 63));

	if(both_signs == 3)		/* Both inputs negative */
	{
		q1.hi ^=MASK_SIGN;	/* Send -q1 and -q2 to low-level sub routine... */
		q2.hi ^=MASK_SIGN;
		q = qfdif(q1, q2); q.hi ^= MASK_SIGN;	/* ...and return negated difference. */
	}
	else if(both_signs == 2)	/* q1 negative */
	{
		q1.hi ^= MASK_SIGN; /* Send q2 and -q1 to low-level add routine. */
		q = qfsum(q2, q2);
	}
	else if(both_signs == 1)	/* q2 negative */
	{
		q2.hi ^= MASK_SIGN; /* Send q1 and -q2 to low-level add routine. */
		q = qfsum(q1, q2);
	}
	else if(both_signs == 0)	/* Both inputs positive */
	{
		q = qfdif(q1, q2);
	}
	else
	{
		ASSERT(HERE, 0,"ERROR: unrecognized sign combination in QFSUB");
	}

	return q;
}


/*
Qfloat addition algorithm. Assumes both summands are nonnegative. A fair amount of code,
but should need only about 10 cycles in a nicely optimized machine implementation.

Example: in our qfloat format,

pi = 2^2 * 0.1100100100001111110110101010001000100001011010001100	0010001101001100010011000110011000101000101110000000110111000010,
 e = 2^2 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100,

and the exact sum pi + e looks like

     2^2 * 1.0111011100001000001011101111101011000100001001000000	1100110011110111010010100010101110001010101110101000000110010110
   = 2^3 * 0.1011101110000100000101110111110101100010000100100000	0110011001111011101001010001010111000101010111010100000011001011.

Example 2: the exact sum pi + 1024*e looks like

     2^12 * 0.0000000000110010010000111111011010101000100010000101	1010001100001000110100110001001100011001100010100010111000000011
   + 2^12 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100

   = 2^12 * 0.1010111000101010100110000100111101001011010000111010	0100110010110011110100001101100001111011100011001010000111010111.
*/
struct qfloat qfsum(struct qfloat q1, struct qfloat q2)
{
	struct qfloat return_val;
	struct qfloat *ptr0, *ptr1;
	uint32 rshift, lshift;
	uint64 exp0, exp1, hi0, hi1, lo0, lo1, offword, carry;

	/* Make sure both inputs are nonnegative. */
	if(!((int64)q1.hi >= 0 && (int64)q2.hi >= 0)) ASSERT(HERE, 0,"ERROR 160 in qfloat.c");

	/* Larger of the two operands gets index 0 in our local length-2 arrays: */
	if(qfcmple(q2, q1))
	{
		ptr0 = &q1;	ptr1 = &q2;
	}
	else
	{
		ptr0 = &q2;	ptr1 = &q1;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */

	hi0  = ptr0->hi & MASK_MANT;
	hi1  = ptr1->hi & MASK_MANT;

	/* Sign = zero, so get exponent by simply subtracting off low bits. */
	exp0 = ptr0->hi - hi0;
	exp1 = ptr1->hi - hi1;

	lo0  = ptr0->lo;
	lo1  = ptr1->lo;

	/* If number is normalized, restore hidden bit. Otherwise, add one to the exponent field to get the true exponent. */
	if(exp0)
	{
		hi0  += TWOE52;
	}
	else
	{
		exp0 += TWOE52;
	}

	if(exp1)
	{
		hi1  += TWOE52;
	}
	else
	{
		exp1 += TWOE52;
	}

	/* ...then right-shift the smaller summand and add. */

	rshift = (exp0 - exp1) >> 52;
	lshift = 64 - rshift;

	carry = 0;

	if(rshift == 0)		/* Operands perfectly aligned */
	{
		lo0 += lo1;
		if(lo0 < lo1)		/* Had a carry. */
		{
			carry = (uint64)1;
		}
		hi0 += hi1 + carry;
	}
	else if(rshift <= 53)	/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
		lo1 = (lo1 >> rshift) + (hi1 << lshift);
		lo1 += offword;
		if(lo1 < offword)	/* Had a carry from the rounded bit being added into the right-shifted low part. */
		{
			carry = (uint64)1;
		}
		hi1 = (hi1 >> rshift) + carry;

		lo0 += lo1;
		if(lo0 < lo1)		/* Had a carry. */
		{
			carry = (uint64)1;
		}
		hi0 += hi1 + carry;
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;
		lo1 = (lo1 >> rshift) + (hi1 << lshift) + offword;

		lo0 += lo1;
		if(lo0 < lo1)		/* Had a carry. */
		{
			hi0 += (uint64)1;
		}
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = lo1 >> 63;
		lo1 = hi1 + offword;

		lo0 += lo1;
		if(lo0 < lo1)		/* Had a carry. */
		{
			hi0 += (uint64)1;
		}
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes hi1 behave like lo1 does in the above segments. */

		offword = (hi1 << lshift) >> 63;
		lo1 = (hi1 >> rshift) + offword;

		lo0 += lo1;
		if(lo0 < lo1)		/* Had a carry. */
		{
			hi0 += (uint64)1;
		}
	}
	else			/* Entire summand completely shifted off,  no need to round result. */
	{
		return *ptr0;
	}

	/* If carry a one into 53-bit <53:53> of the high part, right-shift sum one place and add one to result exponent.
	Potentially need to do this step twice, if round bit from low word ripples all the way back into 53-bit, but we'll
	conveniently ignore that complication.
	*/
#if 1
	rshift = hi0 >> 52;	/* rshift = 3, 2, 1 or 0 */
	if(rshift > 1)
	{
		return_val.hi = exp0 + (hi0 >> 1);			/* Right-shift high word one place, then simply add result to exponent field,
													exp. gets incremented automatically if 53-bit was set. */
		return_val.lo = (lo0 >> 1) + (hi0 << 63);
	}
	else if(rshift == 1)
	{
		return_val.hi = exp0 + (hi0 & MASK_MANT);	/* Mask off MSB of result significand prior to adding to exponent field. */
		return_val.lo = lo0;
	}
	else	/* Both inputs were denormalized, result is also denormalized. Subtract one from exponent field prior to returning. */
	{
		return_val.hi = exp0 - TWOE52 + hi0;
		return_val.lo = lo0;
	}
#else
  #error Invalid code segment in qfloat.c::qfsum!
	/* Branchless version of the above sequence: */
  #if 0
	/*** 9/9/2005: WTF was this? How could this glaring of an error go unnoticed for so long?***/
	rshift = hi0 >> 52;
	lshift = (64 - rshift) & 63;
	return_val.hi = exp0 - TWOE52 + (hi0 >> rshift);
	/*******************************************************************************************/
  #endif
	rshift = hi0 >> 53;				/* rshift = 1 if there was a carry into 53-bit, 0 otherwise. */
	lshift = (64 - rshift) & 63;	/* Need to mask off lower 6 bits, since C doesn't specify what happens if shift count = 64. */
	return_val.hi = (exp0 + ((uint64)rshift << 52)) + ((hi0 >> rshift) - TWOE52);	/* exponent field, gets incremented by 1 if 53-bit was set,
													unchanged if 53-bit = 0 but 52-bit = 1, decremented by 1 if both = 0. */
	return_val.lo = (lo0 >> rshift) + (hi0 << lshift);
#endif
	return(return_val);
}


/*
Qfloat subtraction algorithm. Assumes both inputs are nonnegative.

Example: in our qfloat format,

pi = 2^2 * 0.1100100100001111110110101010001000100001011010001100	0010001101001100010011000110011000101000101110000000110111000010,
 e = 2^2 * 0.1010110111111000010101000101100010100010101110110100	1010100110101010111111011100010101100010000000100111001111010100,

and the difference pi - e looks like (with a borrow from the low-part subtract)

	 2^2 * 0.0001101100010111100001100100100101111110101011010111   0111100110100001010011101010000011000110101101011001100111110010
   = 2^-1* 0.1101100010111100001100100100101111110101011010111011   1100110100001010011101010000011000110101101011001100111101110000
In hex qfloat format, this is 0x399D8BC324BF56BB  0xCD0A750635ACCF70 .
*/
struct qfloat qfdif(struct qfloat q1, struct qfloat q2)
{
	struct qfloat return_val;
	struct qfloat *ptr0, *ptr1;
	uint32 rshift, lshift, lshift_max;
	uint64 exp0, exp1, hi0, hi1, lo0, lo1, offword, carry;

	/* Make sure both inputs are nonnegative. */
	if(!((int64)q1.hi >= 0 && (int64)q2.hi >= 0)) ASSERT(HERE, 0,"ERROR 170 in qfloat.c");

	/* Larger of the two operands gets index 0 in our local length-2 arrays: */
	if(qfcmple(q2, q1))
	{
		/* One more reason to make function arguments qfloat rather than qfloat * : if q1, q2 were ptrs
		and pointed to the same thing, ptr0 and ptr1 would also point to that one address, which would
		trigger the sign-flip conditional prior to returning, causing a return value of -0 here.
		*/
		ptr0 = &q1;	ptr1 = &q2;
	}
	else
	{
		ptr0 = &q2;	ptr1 = &q1;
	}

	/* Separate upper parts of the 2 significands from the sign/exponent fields: */

	hi0  = ptr0->hi & MASK_MANT;
	hi1  = ptr1->hi & MASK_MANT;

	/* Sign assumed zero, so get exponent by simply subtracting off low bits. */
	exp0 = ptr0->hi - hi0;
	exp1 = ptr1->hi - hi1;

	lo0  = ptr0->lo;
	lo1  = ptr1->lo;

	/* If number is normalized, restore hidden bit.
	Otherwise, add one to the exponent field to get the true exponent. */
	if(exp0)
	{
		hi0  += TWOE52;
	}
	else
	{
		exp0 += TWOE52;
	}

	if(exp1)
	{
		hi1  += TWOE52;
	}
	else
	{
		exp1 += TWOE52;
	}

	/* ...then right-shift the smaller summand and add. */

	rshift = (exp0 - exp1) >> 52;
	lshift = 64 - rshift;

	carry = 0;

	if(rshift == 0)		/* Operands perfectly aligned */
	{
		lo0 -= lo1;
		if(lo0 > ptr0->lo)		/* Had a borrow. */
		{
			carry = (uint64)1;
		}
		/* Since ptr0 >= ptr1, hi0 >= 0, no check for borrow-out needed here: */
		hi0 -= hi1 + carry;
	}
	else if(rshift <= 53)		/* Hi part partially shifted into lo, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;	/* Round by adding the MSB of off-shifted low word into next-higher word. */
		lo1 = (lo1 >> rshift) + (hi1 << lshift);
		lo1 += offword;
		if(lo1 < offword)	/* Had a carry from the rounded bit being added into the right-shifted low part. */
		{
			carry = (uint64)1;
		}
		hi1 = (hi1 >> rshift) + carry;

		lo0 -= lo1;
		if(lo0 > ptr0->lo)		/* Had a borrow. */
		{
			carry = (uint64)1;
		}
		hi0 -= hi1 + carry;
	}
	else if(rshift < 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{
		offword = (lo1 << lshift) >> 63;
		lo1 = (lo1 >> rshift) + (hi1 << lshift) + offword;

		lo0 -= lo1;
		if(lo0 > ptr0->lo)		/* Had a borrow. */
		{
			hi0 -= (uint64)1;
		}
	}
	else if(rshift == 64)	/* Hi part completely shifted into lo with no chance of carry-from-round into hi1, lo part partially shifted off. */
	{				/* Treat this case specially, since e.g. x86 does bizarre things when shift count== 64. */
		offword = lo1 >> 63;
		lo1 = hi1 + offword;

		lo0 -= lo1;
		if(lo0 > ptr0->lo)		/* Had a borrow. */
		{
			hi0 -= (uint64)1;
		}
	}
	else if(rshift <= 117)	/* Hi part completely shifted into lo and partially shifted off, lo part completely shifted off. */
	{
		rshift -= 64;
		lshift += 64;	/* lshift treated as unsigned, so put back into range. This makes hi1 behave like lo1 does in the above segments. */

		offword = (hi1 << lshift) >> 63;
		lo1 = (hi1 >> rshift) + offword;

		lo0 -= lo1;
		if(lo0 > ptr0->lo)		/* Had a borrow. */
		{
			hi0 -= (uint64)1;
		}
	}
	else			/* Entire summand completely shifted off, no need to round result. */
	{
		return *ptr0;
	}

	/* Now left-shift result so MSB is in bit <52:52> of high word. */
	if(hi0)
	{
		lshift = leadz64(hi0) - 11;
	}
	else
	{
		lshift = leadz64(lo0) + 53;
	}

	if(lshift && (exp0 > TWOE52))	/* Left-justify mantissa if subtract caused us to lose bits at the left end AND larger input was normalized. */
	{
		lshift_max = (exp0 - TWOE52) >> 52;	/* If result underflows, left-shift mantissa just enough to make exp field = 1. */
		lshift = MIN(lshift, lshift_max);
		rshift = 64 - lshift;				/* Don't have to worry about this being  = 64, since lshift > 0 here. */

		if(lshift < 53)	/* Hi part is nonzero, gets left-shifted <= 52 bits, upper (lshift) bits of lo shifted in to fill vacated bits. */
		{
			hi0 = (hi0 << lshift) + (lo0 >> rshift);
			lo0 = (lo0 << lshift);
			exp0 -= (uint64)lshift << 52;
		}
		else			/* Hi part is zero. Assuming lo part has lzlo lead zeros, right-shift lo (11-lzlo) bits and put that into hi. */
		{
#if QFDEBUG
	printf("WARNING: catastrophic loss of precision in subtract:\n %16llX %16llX -\n %16llX %16llX\n", ptr0->hi, ptr0->lo, ptr1->hi, ptr1->lo);
#endif
			if((int32)rshift > 0)	/* Lo part has > 53 SB, upper 53 of which get put into hi part. */
			{						/* That is, lo gets right-shifted <= 11 bits, result put into hi. Lshift now in [52, 63), i.e. in range. */
				hi0 = (lo0 >> rshift);
				lo0 = (lo0 << lshift);
				exp0 -= (uint64)lshift << 52;
			}
			else					/* Lo part has <= 53 SB, all of which get put into hi part. */
			{
				if(lshift < 117)	/* Lo part has at least one SB */
				{
					hi0 = (lo0 << -rshift);
					lo0 = (uint64)0;
					exp0 -= (uint64)lshift << 52;
				}
				else				/* Result = 0 */
				{
					exp0 = 	TWOE52;	/* Set exponent = 1, since we're going to subtract one prior to returning. */
				}
			}
		}
	}

	return_val.hi = exp0 - TWOE52 + hi0;	/* exponent field unchanged if 52-bit = 1, decremented by 1 if both = 0. */
	if(ptr0 == &q2)
	{
		return_val.hi ^= MASK_SIGN;	/* If flipped order of inputs, negate result. */
	}
	return_val.lo = lo0;

	return(return_val);
}


/*
Newtonian iterative inverse scheme: given nonzero x, whose multiplicative
inverse we wish to find, we begin with an initial guess y and iterate as follows:

		y = y*(2 - x*y) .

Assuming the initial guess is reasonably close to 1/x, the number of significant
bits in y should roughly double on each iteration. If x is a qfloat and our initial
guess is the double-precision-accurate inverse of x, we need just 2 such iterations.
*/
struct qfloat qfinv(struct qfloat q)
{
	double qinv_dble;
	struct qfloat qinv, xyprod;

/* Make sure q is properly normalized. This also catches potential divides-by-zero. */
if((q.hi & ~(MASK_SIGN + MASK_MANT)) == (uint64)0)
{
	ASSERT(HERE, 0,"ERROR: divide by denormalized input not supported.");
}

	/* Get double-precision approximation to qinv: */

	qinv_dble = (double)1.0/qfdbl(q);
	qinv = dbl_to_q(qinv_dble);

	/* Do 2 Newton iterations: */

	xyprod = qfmul(q, qinv);
	xyprod = qfsub(QTWO, xyprod);
	qinv   = qfmul(qinv, xyprod);

	xyprod = qfmul(q, qinv);
	xyprod = qfsub(QTWO, xyprod);
	qinv   = qfmul(qinv, xyprod);

	return qinv;
}


/*
Divide q1/q2, using iterative inverse and a multiply.
*/
struct qfloat qfdiv(struct qfloat q1, struct qfloat q2)
{
	struct qfloat qinv;

	qinv = qfinv(q2);
	qinv = qfmul(q1, qinv);

	return qinv;
}

/*
Return a quotient p/q of 64-bit ints in qfloat form.
*/
struct qfloat qf_rational_quotient(int64 p, int64 q)
{
	return qfdiv(i64_to_q(p),i64_to_q(q));
}

/*
Newtonian iterative square root scheme: given nonnegative x, whose square root
we wish to find, we begin with an initial guess to 1/sqrt(x) and iterate as follows:

		y = y*(3 - x*y^2)/2 .

Assuming the initial guess is reasonably close, the iteration should converge to 1/sqrt(x),
with the number of significant bits in y roughly doubling on each iteration. If x is a qfloat and
our initial guess is the double-precision-accurate approximation to sqrt(x), we need just 2 such
iterations.
*/
struct qfloat qfsqrt(struct qfloat q)
{
	double qisrt_dble;
	struct qfloat qisrt, xysqr;

/* Make sure q is nonnegative. This also catches potential divides-by-zero. */
if(q.hi >> 63)
{
	ASSERT(HERE, 0,"ERROR: sqrt of a negative number not supported.");
}
else if(qfcmpeq(q, QZRO))
{
	return QZRO;
}

	/* Get double-precision approximation to 1/sqrt(x): */

	qisrt_dble = (double)1.0/sqrt(qfdbl(q));
	qisrt = dbl_to_q(qisrt_dble);

	/* Do 2 Newton iterations: */

	xysqr = qfmul(q, qfmul(qisrt, qisrt));
	xysqr = qfsub(QTHREE, xysqr);
	qisrt = qfmul_pow2(qfmul(qisrt, xysqr), -1);

	xysqr = qfmul(q, qfmul(qisrt, qisrt));
	xysqr = qfsub(QTHREE, xysqr);
	qisrt = qfmul_pow2(qfmul(qisrt, xysqr), -1);

	qisrt = qfmul(q, qisrt);	/* Multiply 1/sqrt(x) by x to get sqrt(x). */
	return qisrt;
}


/*
Exponential. Algorithm is not sophisticated - given qfloat x, use the Taylor series about 0.
*/
struct qfloat qfexp(struct qfloat q)
{
	uint32 i, e_sum, e_new;
	struct qfloat sum, curr_term, mult, denom;

	/* Initialize partial sum to 1 + x... */

	curr_term = q;
	sum = qfadd(QONE, curr_term);

	/* get next term of series by multiplying current one by x/i */

	/* exp(1):	0x4005BF0A8B145769, 0x5355FB8AC404E7A7(9E3B...), qfexp gives ...4E7A7, ~116 bits of accuracy. */
	for(i = 2; i < 100; i++)
	{
		denom     = i64_to_q((int64)i);
		mult      = qfdiv(q, denom);
		curr_term = qfmul(curr_term, mult);
		sum       = qfadd(curr_term, sum);

		/* If the current term differs from the sum by more than a factor of 2^110, stop. */
		e_sum = (uint32)((sum      .hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		e_new = (uint32)((curr_term.hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		if((int32)(e_sum - e_new) > 115)
		{
#if QFDEBUG
			printf("End exp(x) summation after %d terms.\n", i);
#endif
			break;
		}
	}
	return sum;
}


/* Top-level sine and cosine routines seen by the caller - these examine the
   sign and magnitude of the input and map that to the proper call to either
   +-sin(arg) or +-cos(arg), where arg is in [0, pi/2).

   Use the following symmetries, all based on the general formulae

	cos(x+a) = cos(x)*cos(a) - sin(x)*sin(a);	sin(x+a) = sin(x)*cos(a) + cos(x)*sin(a);

   1) exp(I*(x + 2*k*pi)) = exp(I*x);						(full-rotation symmetry, a = 2*k*pi)

   2) cos(x + pi  ) = -cos(x);      sin(x + pi  ) = -sin(x);	(half-rotation symmetry, a = pi  )

   3) cos(x + pi/2) = -sin(x);      sin(x + pi/2) = +cos(x);	(symmetry about y-axis,  a = pi/2)

*/
struct qfloat qfcos	(struct qfloat q)
{
	struct qfloat qt;
	uint32 quad;

	qt = qfabs(q);				/* Need deal only with nonnegative numbers. */
	qt = qfdiv(qt, QPIHALF);	/* Get q/(pi/2)... */
	qt = qfint(qt);				/* ...And truncate that to the next-smaller integer. */

	/* For the quadrant, we only need the result modulo 4. */
	quad = qt.lo & (uint64)3;

if(!(       qt.hi == 0)) ASSERT(HERE, 0,"ERROR 180 in qfloat.c");
if(!((int64)qt.lo >= 0)) ASSERT(HERE, 0,"ERROR 181 in qfloat.c");
	qt = i64_to_q((int64)qt.lo);

	/* The calling argument is q mod pi/2 */
	q = qfsub(q, qfmul(qt, QPIHALF));

	if(quad == 0)
		q = qfcs1(q);
	else if(quad == 1)
		q = qfneg(qfsn1(q));
	else if(quad == 2)
		q = qfneg(qfcs1(q));
	else if(quad == 3)
		q = qfsn1(q);

	return q;
}

struct qfloat qfsin	(struct qfloat q)
{
	struct qfloat qt;
	uint32 quad, sign = 0;

	if(qfcmplt(q, QZRO))
	{
		sign = 1;
	}
	qt = qfabs(q);				/* Deal only with nonnegative numbers. */
	qt = qfdiv(qt, QPIHALF);	/* Get q/(pi/2)... */
	qt = qfint(qt);				/* ...And truncate that to the next-smaller integer. */

	/* For the quadrant, we only need the result modulo 4. */
	quad = qt.lo & (uint64)3;

if(!(       qt.hi == 0)) ASSERT(HERE, 0,"ERROR 190 in qfloat.c");
if(!((int64)qt.lo >= 0)) ASSERT(HERE, 0,"ERROR 191 in qfloat.c");
	qt = i64_to_q((int64)qt.lo);

	/* The calling argument is q mod pi/2 */
	q = qfsub(q, qfmul(qt, QPIHALF));

	if(quad == 0)
		q = qfsn1(q);
	else if(quad == 1)
		q = qfcs1(q);
	else if(quad == 2)
		q = qfneg(qfsn1(q));
	else if(quad == 3)
		q = qfneg(qfcs1(q));

	if(sign)
		return qfneg(q);
	else
		return q;
}


/* low-level sine and cosine routines. */
/*
Cosine function. Algorithm is not sophisticated - given qfloat x, use the Taylor series about 0.
Requires an argument in [0, pi/2), and then maps that to [0, pi/4) using the symmetry

		cos(pi/2 - x) = sin(x);

Actually, the Taylor series itself doesn't require any limitations whatsoever on the input
argument, but for the sake of fast convergence and high accuracy, we desire smallish arguments.
Since this routine is designed strictly to be called from the high-level qfcos routine, which
maps an arbitrary qfloat argument into [0, pi/2), the restriction on inputs is no problem.
*/
struct qfloat qfcs1(struct qfloat q)
{
	uint32 i, e_sum, e_new;
	struct qfloat qsqr, sum, curr_term, q1, mult, denom;

	/* Make sure argument is in range... */
	if(!(qfcmple(qfneg(QEPS), q) && qfcmplt(q, qfadd(QPIHALF, QEPS))))
	{
		qfcmple(qfneg(QEPS), q);
		qfcmplt(q, qfadd(QPIHALF, QEPS));
		ASSERT(HERE, 0,"ERROR 200 in qfloat.c");
	}

	if(qfcmpeq(QZRO, q))
	{
		return QONE;
	}
	else if(qfcmplt(qfmul_pow2(QPIHALF, -1), q))	/* Argument is > pi/4. Return cos (pi/2 - q). */
	{
		return qfsn1(qfsub(QPIHALF, q));
	}

	/* Initialize partial sum to 1 - x^2/2... */

	qsqr= qfmul(q, q);		curr_term = qfmul_pow2(qsqr, -1);
	sum = qfsub(QONE, curr_term);

	/* get next term of series by multiplying current one by x^2/(i*(i-1)) */

	for(i = 4; i < 100; i += 4)
	{
		/*
		denom     = i64_to_q((int64)i*(int64)(i-1));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		sum       = qfadd(sum, curr_term);

		denom     = i64_to_q((int64)(i+1)*(int64)(i+2));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		sum       = qfsub(sum, curr_term);
		*/
		denom     = i64_to_q((int64)i*(int64)(i-1));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		q1        = curr_term;

		denom     = i64_to_q((int64)(i+1)*(int64)(i+2));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		q1        = qfsub(q1, curr_term);	/* First calculate the difference of the current 2 terms... */
		sum       = qfadd(sum, q1);			/* ...and only then update the partial sum. */

		/* If the current term differs from the sum by more than a factor of 2^110, stop. */
		e_sum = (uint32)((sum      .hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		e_new = (uint32)((curr_term.hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		if((int32)(e_sum - e_new) > 115)
		{
#if QFDEBUG
			printf("End cos(x) summation after %d terms.\n", i/2);
#endif
			break;
		}
	}
	return sum;
}


/*
Sine function. Algorithm is not sophisticated - given qfloat x, use the Taylor series about 0.
Requires an argument in [0, pi/2), and then maps that to [0, pi/4) using the symmetry

	    sin(pi/2 - x) = cos(x).

Actually, the Taylor series itself doesn't require any limitations whatsoever on the input
argument, but for the sake of fast convergence and high accuracy, we desire smallish arguments.
Since this routine is designed strictly to be called from the high-level qfsin routine, which
maps an arbitrary qfloat argument into [0, pi/2), the restriction on inputs is no problem.
*/
struct qfloat qfsn1(struct qfloat q)
{
	uint32 i, e_sum, e_new;
	struct qfloat qsqr, sum, curr_term, q1, mult, denom;

	/* Make sure argument is in range... */
	if(!(qfcmple(qfneg(QEPS), q) && qfcmplt(q, qfadd(QPIHALF, QEPS))))
	{
		qfcmple(qfneg(QEPS), q);
		qfcmplt(q, qfadd(QPIHALF, QEPS));
		ASSERT(HERE, 0,"ERROR 210 in qfloat.c");
	}

	if(qfcmpeq(QZRO, q))
	{
		return QZRO;
	}
	else if(qfcmplt(qfmul_pow2(QPIHALF, -1), q))	/* Argument is > pi/4. Return cos (pi/2 - q). */
	{
		return qfcs1(qfsub(QPIHALF, q));
	}

	/* Initialize partial sum to x... */

	qsqr= qfmul(q, q);
	curr_term = q;
	sum       = q;

	/* get next term of series by multiplying current one by x^2/(i*(i-1)) */

	for(i = 3; i < 100; i += 4)
	{
		/*
		denom     = i64_to_q((int64)i*(int64)(i-1));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		sum       = qfsub(sum, curr_term);

		denom     = i64_to_q((int64)(i+1)*(int64)(i+2));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		sum       = qfadd(sum, curr_term);
		*/
		denom     = i64_to_q((int64)i*(int64)(i-1));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		q1        = curr_term;

		denom     = i64_to_q((int64)(i+1)*(int64)(i+2));
		mult      = qfdiv(qsqr, denom);
		curr_term = qfmul(curr_term, mult);
		q1        = qfsub(q1, curr_term);	/* First calculate the difference of the current 2 terms... */
		sum       = qfsub(sum, q1);			/* ...and only then update the partial sum. */

		/* If the current term differs from the sum by more than a factor of 2^110, stop. */
		e_sum = (uint32)((sum      .hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		e_new = (uint32)((curr_term.hi & ~(MASK_SIGN + MASK_MANT)) >> 52);
		if((int32)(e_sum - e_new) > 115)
		{
#if QFDEBUG
			printf("End sin(x) summation after %d terms.\n", i/2);
#endif
			break;
		}
	}
	return sum;
}

