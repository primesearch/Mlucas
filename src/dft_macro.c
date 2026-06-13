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

#include "masterdefs.h"
#include "Mdata.h"
#include "align.h"
#include "dft_macro.h"
#include "radix128_twiddles.h"
#include "radix256_twiddles.h"
// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC)
	#include "sse2_macro_gcc64.h"
	#include "radix09_sse_macro.h"
#endif

#ifndef COMPACT_OBJ	// Toggle for parametrized-loop-DFT compact-object code scheme
	#define COMPACT_OBJ	1
#endif

const uint8 reverse8[8] = {0,4,2,6,1,5,3,7};
const uint8 reverse16[16] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};

/************** RADIX-32 DIF/DIT: *****************************/

// Inlining these larger macros leads to very long compiles, so instead prototype here and define them as functions in dft_macro.c.

/* Totals: 376 ADD, 88 MUL	*/
/* Because of the way the original (non-macro-ized) code was written, it's convenient to index the A-inputs in terms
of 4 length-8 blocks with octal indices, and the B-outputs in terms of 2 length-16 blocks with hexadecimal indices.
MSVC allows a maximum of 'only' 127 macro args (which is probably a good thing), so unlike the smaller-radix DFT
macros which use actual array-indexed terms as args, here we use pointers to the real part of each complex arg:
*/
void RADIX_32_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 32 (index) offsets */
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);
		__rt =*(__A+__idx[0x10]);			__it =*(Aim+__idx[0x10]);
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(__A+__idx[0x08]);			__t05=*(Aim+__idx[0x08]);
		__rt =*(__A+__idx[0x18]);			__it =*(Aim+__idx[0x18]);
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02+__it;					__t07=__t03-__rt;
		__t02=__t02-__it;					__t03=__t03+__rt;

		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);
		__rt =*(__A+__idx[0x14]);			__it =*(Aim+__idx[0x14]);
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(__A+__idx[0x0C]);			__t0D=*(Aim+__idx[0x0C]);
		__rt =*(__A+__idx[0x1C]);			__it =*(Aim+__idx[0x1C]);
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A+__it;					__t0F=__t0B-__rt;
		__t0A=__t0A-__it;					__t0B=__t0B+__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04+__it;					__t0D=__t05-__rt;
		__t04=__t04-__it;					__t05=__t05+__rt;

		__rt =(__t0A-__t0B)*ISRT2;			__it =(__t0A+__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03-__it;
		__t02=__t02+__rt;					__t03=__t03+__it;

		__rt =(__t0E+__t0F)*ISRT2;			__it =(__t0F-__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x02]);			__t11=*(Aim+__idx[0x02]);
		__rt =*(__A+__idx[0x12]);			__it =*(Aim+__idx[0x12]);
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);
		__rt =*(__A+__idx[0x1A]);			__it =*(Aim+__idx[0x1A]);
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12+__it;					__t17=__t13-__rt;
		__t12=__t12-__it;					__t13=__t13+__rt;

		__t18=*(__A+__idx[0x06]);			__t19=*(Aim+__idx[0x06]);
		__rt =*(__A+__idx[0x16]);			__it =*(Aim+__idx[0x16]);
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);
		__rt =*(__A+__idx[0x1E]);			__it =*(Aim+__idx[0x1E]);
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A+__it;					__t1F=__t1B-__rt;
		__t1A=__t1A-__it;					__t1B=__t1B+__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14+__it;					__t1D=__t15-__rt;
		__t14=__t14-__it;					__t15=__t15+__rt;

		__rt =(__t1A-__t1B)*ISRT2;			__it =(__t1A+__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13-__it;
		__t12=__t12+__rt;					__t13=__t13+__it;

		__rt =(__t1E+__t1F)*ISRT2;			__it =(__t1F-__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x01]);			__t21=*(Aim+__idx[0x01]);
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(__A+__idx[0x09]);			__t25=*(Aim+__idx[0x09]);
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22+__it;					__t27=__t23-__rt;
		__t22=__t22-__it;					__t23=__t23+__rt;

		__t28=*(__A+__idx[0x05]);			__t29=*(Aim+__idx[0x05]);
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(__A+__idx[0x0D]);			__t2D=*(Aim+__idx[0x0D]);
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A+__it;					__t2F=__t2B-__rt;
		__t2A=__t2A-__it;					__t2B=__t2B+__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24+__it;					__t2D=__t25-__rt;
		__t24=__t24-__it;					__t25=__t25+__rt;

		__rt =(__t2A-__t2B)*ISRT2;			__it =(__t2A+__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23-__it;
		__t22=__t22+__rt;					__t23=__t23+__it;

		__rt =(__t2E+__t2F)*ISRT2;			__it =(__t2F-__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x03]);			__t31=*(Aim+__idx[0x03]);
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(__A+__idx[0x0B]);			__t35=*(Aim+__idx[0x0B]);
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32+__it;					__t37=__t33-__rt;
		__t32=__t32-__it;					__t33=__t33+__rt;

		__t38=*(__A+__idx[0x07]);			__t39=*(Aim+__idx[0x07]);
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(__A+__idx[0x0F]);			__t3D=*(Aim+__idx[0x0F]);
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A+__it;					__t3F=__t3B-__rt;
		__t3A=__t3A-__it;					__t3B=__t3B+__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34+__it;					__t3D=__t35-__rt;
		__t34=__t34-__it;					__t35=__t35+__rt;

		__rt =(__t3A-__t3B)*ISRT2;			__it =(__t3A+__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33-__it;
		__t32=__t32+__rt;					__t33=__t33+__it;

		__rt =(__t3E+__t3F)*ISRT2;			__it =(__t3F-__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(__B+__odx[0x00])=__t00+__t20;		*(Bim+__odx[0x00])=__t01+__t21;
		*(__B+__odx[0x01])=__t00-__t20;		*(Bim+__odx[0x01])=__t01-__t21;
		*(__B+__odx[0x02])=__t10-__t31;		*(Bim+__odx[0x02])=__t11+__t30;
		*(__B+__odx[0x03])=__t10+__t31;		*(Bim+__odx[0x03])=__t11-__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08+__t19;					__t08=__t08-__t19;
		__t19=__t09-__rt;					__t09=__t09+__rt;

		__rt =(__t28-__t29)*ISRT2;			__t29=(__t28+__t29)*ISRT2;	__t28=__rt;
		__rt =(__t39+__t38)*ISRT2;			__it =(__t39-__t38)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		*(__B+__odx[0x04])=__t08+__t28;		*(Bim+__odx[0x04])=__t09+__t29;
		*(__B+__odx[0x05])=__t08-__t28;		*(Bim+__odx[0x05])=__t09-__t29;
		*(__B+__odx[0x06])=__t18-__t39;		*(Bim+__odx[0x06])=__t19+__t38;
		*(__B+__odx[0x07])=__t18+__t39;		*(Bim+__odx[0x07])=__t19-__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t14-__t15)*ISRT2;			__it =(__t14+__t15)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 - __t25*s16;			__t25=__t25*c16 + __t24*s16;	__t24=__rt;
		__rt =__t34*s16 - __t35*c16;			__it =__t35*s16 + __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x08])=__t04+__t24;		*(Bim+__odx[0x08])=__t05+__t25;
		*(__B+__odx[0x09])=__t04-__t24;		*(Bim+__odx[0x09])=__t05-__t25;
		*(__B+__odx[0x0A])=__t14-__t35;		*(Bim+__odx[0x0A])=__t15+__t34;
		*(__B+__odx[0x0B])=__t14+__t35;		*(Bim+__odx[0x0B])=__t15-__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1D+__t1C)*ISRT2;			__it =(__t1D-__t1C)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 - __t2D*c16;			__t2D=__t2D*s16 + __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 - __t3D*s16;			__it =__t3D*c16 + __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x0C])=__t0C+__t2C;		*(Bim+__odx[0x0C])=__t0D+__t2D;
		*(__B+__odx[0x0D])=__t0C-__t2C;		*(Bim+__odx[0x0D])=__t0D-__t2D;
		*(__B+__odx[0x0E])=__t1C-__t3D;		*(Bim+__odx[0x0E])=__t1D+__t3C;
		*(__B+__odx[0x0F])=__t1C+__t3D;		*(Bim+__odx[0x0F])=__t1D-__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 - __t13*s16;			__it =__t13*c16 + __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 - __t23*s32_1;	__t23=__t23*c32_1 + __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 - __t33*s32_3;	__it =__t33*c32_3 + __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x10])=__t02+__t22;		*(Bim+__odx[0x10])=__t03+__t23;
		*(__B+__odx[0x11])=__t02-__t22;		*(Bim+__odx[0x11])=__t03-__t23;
		*(__B+__odx[0x12])=__t12-__t33;		*(Bim+__odx[0x12])=__t13+__t32;
		*(__B+__odx[0x13])=__t12+__t33;		*(Bim+__odx[0x13])=__t13-__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 + __t1B*c16;			__it =__t1B*s16 - __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 - __t2B*c32_3;	__t2B=__t2B*s32_3 + __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 + __t3B*s32_1;	__it =__t3B*c32_1 - __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x14])=__t0A+__t2A;		*(Bim+__odx[0x14])=__t0B+__t2B;
		*(__B+__odx[0x15])=__t0A-__t2A;		*(Bim+__odx[0x15])=__t0B-__t2B;
		*(__B+__odx[0x16])=__t1A-__t3B;		*(Bim+__odx[0x16])=__t1B+__t3A;
		*(__B+__odx[0x17])=__t1A+__t3B;		*(Bim+__odx[0x17])=__t1B-__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 - __t17*c16;			__it =__t17*s16 + __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 - __t27*s32_3;	__t27=__t27*c32_3 + __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 + __t37*c32_1;	__it =__t37*s32_1 - __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x18])=__t06+__t26;		*(Bim+__odx[0x18])=__t07+__t27;
		*(__B+__odx[0x19])=__t06-__t26;		*(Bim+__odx[0x19])=__t07-__t27;
		*(__B+__odx[0x1A])=__t16-__t37;		*(Bim+__odx[0x1A])=__t17+__t36;
		*(__B+__odx[0x1B])=__t16+__t37;		*(Bim+__odx[0x1B])=__t17-__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 + __t1F*s16;			__it =__t1F*c16 - __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 - __t2F*c32_1;	__t2F=__t2F*s32_1 + __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 - __t3F*c32_3;	__it =__t3F*s32_3 + __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x1C])=__t0E+__t2E;		*(Bim+__odx[0x1C])=__t0F+__t2F;
		*(__B+__odx[0x1D])=__t0E-__t2E;		*(Bim+__odx[0x1D])=__t0F-__t2F;
		*(__B+__odx[0x1E])=__t1E-__t3F;		*(Bim+__odx[0x1E])=__t1F+__t3E;
		*(__B+__odx[0x1F])=__t1E+__t3F;		*(Bim+__odx[0x1F])=__t1F-__t3E;
}

// With-twiddles out-of-place analog of above twiddleless DIF macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIF version of this macro processes the twiddles in bit-reversed order: 0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f.
void RADIX_32_DIF_TWIDDLE_OOP(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c10,double __s10,
	double __c08,double __s08,
	double __c18,double __s18,
	/**/
	double __c04,double __s04,
	double __c14,double __s14,
	double __c0C,double __s0C,
	double __c1C,double __s1C,
	/**/
	double __c02,double __s02,
	double __c12,double __s12,
	double __c0A,double __s0A,
	double __c1A,double __s1A,
	/**/
	double __c06,double __s06,
	double __c16,double __s16,
	double __c0E,double __s0E,
	double __c1E,double __s1E,
	/**/
	double __c01,double __s01,
	double __c11,double __s11,
	double __c09,double __s09,
	double __c19,double __s19,
	/**/
	double __c05,double __s05,
	double __c15,double __s15,
	double __c0D,double __s0D,
	double __c1D,double __s1D,
	/**/
	double __c03,double __s03,
	double __c13,double __s13,
	double __c0B,double __s0B,
	double __c1B,double __s1B,
	/**/
	double __c07,double __s07,
	double __c17,double __s17,
	double __c0F,double __s0F,
	double __c1F,double __s1F
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Aim = __A + RE_IM_STRIDE, *Bim = __B + RE_IM_STRIDE;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);	__t01=*(Aim+__idx[0x00]);
		__t02=*(__A+__idx[0x10]);	__t03=*(Aim+__idx[0x10]);	__rt = __t02*__c10 -__t03*__s10;	__it = __t03*__c10 +__t02*__s10;
		__t02=__t00-__rt;			__t03=__t01-__it;
		__t00=__t00+__rt;			__t01=__t01+__it;

		__t04=*(__A+__idx[0x08]);	__t05=*(Aim+__idx[0x08]);	__rt = __t04*__c08 -__t05*__s08;	__t05= __t05*__c08 +__t04*__s08;	__t04 = __rt;
		__t06=*(__A+__idx[0x18]);	__t07=*(Aim+__idx[0x18]);	__rt = __t06*__c18 -__t07*__s18;	__it = __t07*__c18 +__t06*__s18;
		__t06=__t04-__rt;			__t07=__t05-__it;
		__t04=__t04+__rt;			__t05=__t05+__it;

		__rt =__t04;				__it =__t05;
		__t04=__t00-__rt;			__t05=__t01-__it;
		__t00=__t00+__rt;			__t01=__t01+__it;

		__rt =__t06;				__it =__t07;
		__t06=__t02+__it;			__t07=__t03-__rt;
		__t02=__t02-__it;			__t03=__t03+__rt;

		__t08=*(__A+__idx[0x04]);	__t09=*(Aim+__idx[0x04]);	__rt = __t08*__c04 -__t09*__s04;	__t09= __t09*__c04 +__t08*__s04;	__t08 = __rt;
		__t0A=*(__A+__idx[0x14]);	__t0B=*(Aim+__idx[0x14]);	__rt = __t0A*__c14 -__t0B*__s14;	__it = __t0B*__c14 +__t0A*__s14;
		__t0A=__t08-__rt;			__t0B=__t09-__it;
		__t08=__t08+__rt;			__t09=__t09+__it;

		__t0C=*(__A+__idx[0x0C]);	__t0D=*(Aim+__idx[0x0C]);	__rt = __t0C*__c0C -__t0D*__s0C;	__t0D= __t0D*__c0C +__t0C*__s0C;	__t0C = __rt;
		__t0E=*(__A+__idx[0x1C]);	__t0F=*(Aim+__idx[0x1C]);	__rt = __t0E*__c1C -__t0F*__s1C;	__it = __t0F*__c1C +__t0E*__s1C;
		__t0E=__t0C-__rt;			__t0F=__t0D-__it;
		__t0C=__t0C+__rt;			__t0D=__t0D+__it;

		__rt =__t0C;				__it =__t0D;
		__t0C=__t08-__rt;			__t0D=__t09-__it;
		__t08=__t08+__rt;			__t09=__t09+__it;

		__rt =__t0E;				__it =__t0F;
		__t0E=__t0A+__it;			__t0F=__t0B-__rt;
		__t0A=__t0A-__it;			__t0B=__t0B+__rt;

		__rt =__t08;				__it =__t09;
		__t08=__t00-__rt;			__t09=__t01-__it;
		__t00=__t00+__rt;			__t01=__t01+__it;

		__rt =__t0C;				__it =__t0D;
		__t0C=__t04+__it;			__t0D=__t05-__rt;
		__t04=__t04-__it;			__t05=__t05+__rt;

		__rt =(__t0A-__t0B)*ISRT2;	__it =(__t0A+__t0B)*ISRT2;
		__t0A=__t02-__rt;			__t0B=__t03-__it;
		__t02=__t02+__rt;			__t03=__t03+__it;

		__rt =(__t0E+__t0F)*ISRT2;	__it =(__t0F-__t0E)*ISRT2;
		__t0E=__t06+__rt;			__t0F=__t07+__it;
		__t06=__t06-__rt;			__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x02]);	__t11=*(Aim+__idx[0x02]);	__rt = __t10*__c02 -__t11*__s02;	__t11= __t11*__c02 +__t10*__s02;	__t10 = __rt;
		__t12=*(__A+__idx[0x12]);	__t13=*(Aim+__idx[0x12]);	__rt = __t12*__c12 -__t13*__s12;	__it = __t13*__c12 +__t12*__s12;
		__t12=__t10-__rt;			__t13=__t11-__it;
		__t10=__t10+__rt;			__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);	__t15=*(Aim+__idx[0x0A]);	__rt = __t14*__c0A -__t15*__s0A;	__t15= __t15*__c0A +__t14*__s0A;	__t14 = __rt;
		__t16=*(__A+__idx[0x1A]);	__t17=*(Aim+__idx[0x1A]);	__rt = __t16*__c1A -__t17*__s1A;	__it = __t17*__c1A +__t16*__s1A;
		__t16=__t14-__rt;			__t17=__t15-__it;
		__t14=__t14+__rt;			__t15=__t15+__it;

		__rt =__t14;				__it =__t15;
		__t14=__t10-__rt;			__t15=__t11-__it;
		__t10=__t10+__rt;			__t11=__t11+__it;

		__rt =__t16;				__it =__t17;
		__t16=__t12+__it;			__t17=__t13-__rt;
		__t12=__t12-__it;			__t13=__t13+__rt;

		__t18=*(__A+__idx[0x06]);	__t19=*(Aim+__idx[0x06]);	__rt = __t18*__c06 -__t19*__s06;	__t19= __t19*__c06 +__t18*__s06;	__t18 = __rt;
		__t1A=*(__A+__idx[0x16]);	__t1B=*(Aim+__idx[0x16]);	__rt = __t1A*__c16 -__t1B*__s16;	__it = __t1B*__c16 +__t1A*__s16;
		__t1A=__t18-__rt;			__t1B=__t19-__it;
		__t18=__t18+__rt;			__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);	__t1D=*(Aim+__idx[0x0E]);	__rt = __t1C*__c0E -__t1D*__s0E;	__t1D= __t1D*__c0E +__t1C*__s0E;	__t1C = __rt;
		__t1E=*(__A+__idx[0x1E]);	__t1F=*(Aim+__idx[0x1E]);	__rt = __t1E*__c1E -__t1F*__s1E;	__it = __t1F*__c1E +__t1E*__s1E;
		__t1E=__t1C-__rt;			__t1F=__t1D-__it;
		__t1C=__t1C+__rt;			__t1D=__t1D+__it;

		__rt =__t1C;				__it =__t1D;
		__t1C=__t18-__rt;			__t1D=__t19-__it;
		__t18=__t18+__rt;			__t19=__t19+__it;

		__rt =__t1E;				__it =__t1F;
		__t1E=__t1A+__it;			__t1F=__t1B-__rt;
		__t1A=__t1A-__it;			__t1B=__t1B+__rt;

		__rt =__t18;				__it =__t19;
		__t18=__t10-__rt;			__t19=__t11-__it;
		__t10=__t10+__rt;			__t11=__t11+__it;

		__rt =__t1C;				__it =__t1D;
		__t1C=__t14+__it;			__t1D=__t15-__rt;
		__t14=__t14-__it;			__t15=__t15+__rt;

		__rt =(__t1A-__t1B)*ISRT2;	__it =(__t1A+__t1B)*ISRT2;
		__t1A=__t12-__rt;			__t1B=__t13-__it;
		__t12=__t12+__rt;			__t13=__t13+__it;

		__rt =(__t1E+__t1F)*ISRT2;	__it =(__t1F-__t1E)*ISRT2;
		__t1E=__t16+__rt;			__t1F=__t17+__it;
		__t16=__t16-__rt;			__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x01]);	__t21=*(Aim+__idx[0x01]);	__rt = __t20*__c01 -__t21*__s01;	__t21= __t21*__c01 +__t20*__s01;	__t20 = __rt;
		__t22=*(__A+__idx[0x11]);	__t23=*(Aim+__idx[0x11]);	__rt = __t22*__c11 -__t23*__s11;	__it = __t23*__c11 +__t22*__s11;
		__t22=__t20-__rt;			__t23=__t21-__it;
		__t20=__t20+__rt;			__t21=__t21+__it;

		__t24=*(__A+__idx[0x09]);	__t25=*(Aim+__idx[0x09]);	__rt = __t24*__c09 -__t25*__s09;	__t25= __t25*__c09 +__t24*__s09;	__t24 = __rt;
		__t26=*(__A+__idx[0x19]);	__t27=*(Aim+__idx[0x19]);	__rt = __t26*__c19 -__t27*__s19;	__it = __t27*__c19 +__t26*__s19;
		__t26=__t24-__rt;			__t27=__t25-__it;
		__t24=__t24+__rt;			__t25=__t25+__it;

		__rt =__t24;				__it =__t25;
		__t24=__t20-__rt;			__t25=__t21-__it;
		__t20=__t20+__rt;			__t21=__t21+__it;

		__rt =__t26;				__it =__t27;
		__t26=__t22+__it;			__t27=__t23-__rt;
		__t22=__t22-__it;			__t23=__t23+__rt;

		__t28=*(__A+__idx[0x05]);	__t29=*(Aim+__idx[0x05]);	__rt = __t28*__c05 -__t29*__s05;	__t29= __t29*__c05 +__t28*__s05;	__t28 = __rt;
		__t2A=*(__A+__idx[0x15]);	__t2B=*(Aim+__idx[0x15]);	__rt = __t2A*__c15 -__t2B*__s15;	__it = __t2B*__c15 +__t2A*__s15;
		__t2A=__t28-__rt;			__t2B=__t29-__it;
		__t28=__t28+__rt;			__t29=__t29+__it;

		__t2C=*(__A+__idx[0x0D]);	__t2D=*(Aim+__idx[0x0D]);	__rt = __t2C*__c0D -__t2D*__s0D;	__t2D= __t2D*__c0D +__t2C*__s0D;	__t2C = __rt;
		__t2E=*(__A+__idx[0x1D]);	__t2F=*(Aim+__idx[0x1D]);	__rt = __t2E*__c1D -__t2F*__s1D;	__it = __t2F*__c1D +__t2E*__s1D;
		__t2E=__t2C-__rt;			__t2F=__t2D-__it;
		__t2C=__t2C+__rt;			__t2D=__t2D+__it;

		__rt =__t2C;				__it =__t2D;
		__t2C=__t28-__rt;			__t2D=__t29-__it;
		__t28=__t28+__rt;			__t29=__t29+__it;

		__rt =__t2E;				__it =__t2F;
		__t2E=__t2A+__it;			__t2F=__t2B-__rt;
		__t2A=__t2A-__it;			__t2B=__t2B+__rt;

		__rt =__t28;				__it =__t29;
		__t28=__t20-__rt;			__t29=__t21-__it;
		__t20=__t20+__rt;			__t21=__t21+__it;

		__rt =__t2C;				__it =__t2D;
		__t2C=__t24+__it;			__t2D=__t25-__rt;
		__t24=__t24-__it;			__t25=__t25+__rt;

		__rt =(__t2A-__t2B)*ISRT2;	__it =(__t2A+__t2B)*ISRT2;
		__t2A=__t22-__rt;			__t2B=__t23-__it;
		__t22=__t22+__rt;			__t23=__t23+__it;

		__rt =(__t2E+__t2F)*ISRT2;	__it =(__t2F-__t2E)*ISRT2;
		__t2E=__t26+__rt;			__t2F=__t27+__it;
		__t26=__t26-__rt;			__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x03]);	__t31=*(Aim+__idx[0x03]);	__rt = __t30*__c03 -__t31*__s03;	__t31= __t31*__c03 +__t30*__s03;	__t30 = __rt;
		__t32=*(__A+__idx[0x13]);	__t33=*(Aim+__idx[0x13]);	__rt = __t32*__c13 -__t33*__s13;	__it = __t33*__c13 +__t32*__s13;
		__t32=__t30-__rt;			__t33=__t31-__it;
		__t30=__t30+__rt;			__t31=__t31+__it;

		__t34=*(__A+__idx[0x0B]);	__t35=*(Aim+__idx[0x0B]);	__rt = __t34*__c0B -__t35*__s0B;	__t35= __t35*__c0B +__t34*__s0B;	__t34 = __rt;
		__t36=*(__A+__idx[0x1B]);	__t37=*(Aim+__idx[0x1B]);	__rt = __t36*__c1B -__t37*__s1B;	__it = __t37*__c1B +__t36*__s1B;
		__t36=__t34-__rt;			__t37=__t35-__it;
		__t34=__t34+__rt;			__t35=__t35+__it;

		__rt =__t34;				__it =__t35;
		__t34=__t30-__rt;			__t35=__t31-__it;
		__t30=__t30+__rt;			__t31=__t31+__it;

		__rt =__t36;				__it =__t37;
		__t36=__t32+__it;			__t37=__t33-__rt;
		__t32=__t32-__it;			__t33=__t33+__rt;

		__t38=*(__A+__idx[0x07]);	__t39=*(Aim+__idx[0x07]);	__rt = __t38*__c07 -__t39*__s07;	__t39= __t39*__c07 +__t38*__s07;	__t38 = __rt;
		__t3A=*(__A+__idx[0x17]);	__t3B=*(Aim+__idx[0x17]);	__rt = __t3A*__c17 -__t3B*__s17;	__it = __t3B*__c17 +__t3A*__s17;
		__t3A=__t38-__rt;			__t3B=__t39-__it;
		__t38=__t38+__rt;			__t39=__t39+__it;

		__t3C=*(__A+__idx[0x0F]);	__t3D=*(Aim+__idx[0x0F]);	__rt = __t3C*__c0F -__t3D*__s0F;	__t3D= __t3D*__c0F +__t3C*__s0F;	__t3C = __rt;
		__t3E=*(__A+__idx[0x1F]);	__t3F=*(Aim+__idx[0x1F]);	__rt = __t3E*__c1F -__t3F*__s1F;	__it = __t3F*__c1F +__t3E*__s1F;
		__t3E=__t3C-__rt;			__t3F=__t3D-__it;
		__t3C=__t3C+__rt;			__t3D=__t3D+__it;

		__rt =__t3C;				__it =__t3D;
		__t3C=__t38-__rt;			__t3D=__t39-__it;
		__t38=__t38+__rt;			__t39=__t39+__it;

		__rt =__t3E;				__it =__t3F;
		__t3E=__t3A+__it;			__t3F=__t3B-__rt;
		__t3A=__t3A-__it;			__t3B=__t3B+__rt;

		__rt =__t38;				__it =__t39;
		__t38=__t30-__rt;			__t39=__t31-__it;
		__t30=__t30+__rt;			__t31=__t31+__it;

		__rt =__t3C;				__it =__t3D;
		__t3C=__t34+__it;			__t3D=__t35-__rt;
		__t34=__t34-__it;			__t35=__t35+__rt;

		__rt =(__t3A-__t3B)*ISRT2;	__it =(__t3A+__t3B)*ISRT2;
		__t3A=__t32-__rt;			__t3B=__t33-__it;
		__t32=__t32+__rt;			__t33=__t33+__it;

		__rt =(__t3E+__t3F)*ISRT2;	__it =(__t3F-__t3E)*ISRT2;
		__t3E=__t36+__rt;			__t3F=__t37+__it;
		__t36=__t36-__rt;			__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(__B+__odx[0x00])=__t00+__t20;		*(Bim+__odx[0x00])=__t01+__t21;
		*(__B+__odx[0x01])=__t00-__t20;		*(Bim+__odx[0x01])=__t01-__t21;
		*(__B+__odx[0x02])=__t10-__t31;		*(Bim+__odx[0x02])=__t11+__t30;
		*(__B+__odx[0x03])=__t10+__t31;		*(Bim+__odx[0x03])=__t11-__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08+__t19;					__t08=__t08-__t19;
		__t19=__t09-__rt;					__t09=__t09+__rt;

		__rt =(__t28-__t29)*ISRT2;			__t29=(__t28+__t29)*ISRT2;	__t28=__rt;
		__rt =(__t39+__t38)*ISRT2;			__it =(__t39-__t38)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		*(__B+__odx[0x04])=__t08+__t28;		*(Bim+__odx[0x04])=__t09+__t29;
		*(__B+__odx[0x05])=__t08-__t28;		*(Bim+__odx[0x05])=__t09-__t29;
		*(__B+__odx[0x06])=__t18-__t39;		*(Bim+__odx[0x06])=__t19+__t38;
		*(__B+__odx[0x07])=__t18+__t39;		*(Bim+__odx[0x07])=__t19-__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t14-__t15)*ISRT2;			__it =(__t14+__t15)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 - __t25*s16;__t25=__t25*c16 + __t24*s16;	__t24=__rt;
		__rt =__t34*s16 - __t35*c16;__it =__t35*s16 + __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x08])=__t04+__t24;		*(Bim+__odx[0x08])=__t05+__t25;
		*(__B+__odx[0x09])=__t04-__t24;		*(Bim+__odx[0x09])=__t05-__t25;
		*(__B+__odx[0x0A])=__t14-__t35;		*(Bim+__odx[0x0A])=__t15+__t34;
		*(__B+__odx[0x0B])=__t14+__t35;		*(Bim+__odx[0x0B])=__t15-__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1D+__t1C)*ISRT2;			__it =(__t1D-__t1C)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 - __t2D*c16;__t2D=__t2D*s16 + __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 - __t3D*s16;__it =__t3D*c16 + __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x0C])=__t0C+__t2C;		*(Bim+__odx[0x0C])=__t0D+__t2D;
		*(__B+__odx[0x0D])=__t0C-__t2C;		*(Bim+__odx[0x0D])=__t0D-__t2D;
		*(__B+__odx[0x0E])=__t1C-__t3D;		*(Bim+__odx[0x0E])=__t1D+__t3C;
		*(__B+__odx[0x0F])=__t1C+__t3D;		*(Bim+__odx[0x0F])=__t1D-__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 - __t13*s16;__it =__t13*c16 + __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 - __t23*s32_1;__t23=__t23*c32_1 + __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 - __t33*s32_3;__it =__t33*c32_3 + __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x10])=__t02+__t22;		*(Bim+__odx[0x10])=__t03+__t23;
		*(__B+__odx[0x11])=__t02-__t22;		*(Bim+__odx[0x11])=__t03-__t23;
		*(__B+__odx[0x12])=__t12-__t33;		*(Bim+__odx[0x12])=__t13+__t32;
		*(__B+__odx[0x13])=__t12+__t33;		*(Bim+__odx[0x13])=__t13-__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 + __t1B*c16;__it =__t1B*s16 - __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 - __t2B*c32_3;__t2B=__t2B*s32_3 + __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 + __t3B*s32_1;__it =__t3B*c32_1 - __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x14])=__t0A+__t2A;		*(Bim+__odx[0x14])=__t0B+__t2B;
		*(__B+__odx[0x15])=__t0A-__t2A;		*(Bim+__odx[0x15])=__t0B-__t2B;
		*(__B+__odx[0x16])=__t1A-__t3B;		*(Bim+__odx[0x16])=__t1B+__t3A;
		*(__B+__odx[0x17])=__t1A+__t3B;		*(Bim+__odx[0x17])=__t1B-__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 - __t17*c16;__it =__t17*s16 + __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 - __t27*s32_3;__t27=__t27*c32_3 + __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 + __t37*c32_1;__it =__t37*s32_1 - __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x18])=__t06+__t26;		*(Bim+__odx[0x18])=__t07+__t27;
		*(__B+__odx[0x19])=__t06-__t26;		*(Bim+__odx[0x19])=__t07-__t27;
		*(__B+__odx[0x1A])=__t16-__t37;		*(Bim+__odx[0x1A])=__t17+__t36;
		*(__B+__odx[0x1B])=__t16+__t37;		*(Bim+__odx[0x1B])=__t17-__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 + __t1F*s16;__it =__t1F*c16 - __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 - __t2F*c32_1;__t2F=__t2F*s32_1 + __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 - __t3F*c32_3;__it =__t3F*s32_3 + __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x1C])=__t0E+__t2E;		*(Bim+__odx[0x1C])=__t0F+__t2F;
		*(__B+__odx[0x1D])=__t0E-__t2E;		*(Bim+__odx[0x1D])=__t0F-__t2F;
		*(__B+__odx[0x1E])=__t1E-__t3F;		*(Bim+__odx[0x1E])=__t1F+__t3E;
		*(__B+__odx[0x1F])=__t1E+__t3F;		*(Bim+__odx[0x1F])=__t1F-__t3E;
}


void RADIX_32_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 32 (index) offsets */
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);
		__rt =*(__A+__idx[0x01]);			__it =*(Aim+__idx[0x01]);
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(__A+__idx[0x02]);			__t05=*(Aim+__idx[0x02]);
		__rt =*(__A+__idx[0x03]);			__it =*(Aim+__idx[0x03]);
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02-__it;					__t07=__t03+__rt;
		__t02=__t02+__it;					__t03=__t03-__rt;

		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);
		__rt =*(__A+__idx[0x05]);			__it =*(Aim+__idx[0x05]);
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(__A+__idx[0x06]);			__t0D=*(Aim+__idx[0x06]);
		__rt =*(__A+__idx[0x07]);			__it =*(Aim+__idx[0x07]);
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A-__it;					__t0F=__t0B+__rt;
		__t0A=__t0A+__it;					__t0B=__t0B-__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04-__it;					__t0D=__t05+__rt;
		__t04=__t04+__it;					__t05=__t05-__rt;

		__rt =(__t0A+__t0B)*ISRT2;			__it =(__t0A-__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03+__it;
		__t02=__t02+__rt;					__t03=__t03-__it;

		__rt =(__t0E-__t0F)*ISRT2;			__it =(__t0F+__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x08]);			__t11=*(Aim+__idx[0x08]);
		__rt =*(__A+__idx[0x09]);			__it =*(Aim+__idx[0x09]);
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);
		__rt =*(__A+__idx[0x0B]);			__it =*(Aim+__idx[0x0B]);
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12-__it;					__t17=__t13+__rt;
		__t12=__t12+__it;					__t13=__t13-__rt;

		__t18=*(__A+__idx[0x0C]);			__t19=*(Aim+__idx[0x0C]);
		__rt =*(__A+__idx[0x0D]);			__it =*(Aim+__idx[0x0D]);
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);
		__rt =*(__A+__idx[0x0F]);			__it =*(Aim+__idx[0x0F]);
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A-__it;					__t1F=__t1B+__rt;
		__t1A=__t1A+__it;					__t1B=__t1B-__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14-__it;					__t1D=__t15+__rt;
		__t14=__t14+__it;					__t15=__t15-__rt;

		__rt =(__t1A+__t1B)*ISRT2;			__it =(__t1A-__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13+__it;
		__t12=__t12+__rt;					__t13=__t13-__it;

		__rt =(__t1E-__t1F)*ISRT2;			__it =(__t1F+__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x10]);			__t21=*(Aim+__idx[0x10]);
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(__A+__idx[0x12]);			__t25=*(Aim+__idx[0x12]);
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22-__it;					__t27=__t23+__rt;
		__t22=__t22+__it;					__t23=__t23-__rt;

		__t28=*(__A+__idx[0x14]);			__t29=*(Aim+__idx[0x14]);
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(__A+__idx[0x16]);			__t2D=*(Aim+__idx[0x16]);
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A-__it;					__t2F=__t2B+__rt;
		__t2A=__t2A+__it;					__t2B=__t2B-__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24-__it;					__t2D=__t25+__rt;
		__t24=__t24+__it;					__t25=__t25-__rt;

		__rt =(__t2A+__t2B)*ISRT2;			__it =(__t2A-__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23+__it;
		__t22=__t22+__rt;					__t23=__t23-__it;

		__rt =(__t2E-__t2F)*ISRT2;			__it =(__t2F+__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x18]);			__t31=*(Aim+__idx[0x18]);
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(__A+__idx[0x1A]);			__t35=*(Aim+__idx[0x1A]);
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32-__it;					__t37=__t33+__rt;
		__t32=__t32+__it;					__t33=__t33-__rt;

		__t38=*(__A+__idx[0x1C]);			__t39=*(Aim+__idx[0x1C]);
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(__A+__idx[0x1E]);			__t3D=*(Aim+__idx[0x1E]);
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A-__it;					__t3F=__t3B+__rt;
		__t3A=__t3A+__it;					__t3B=__t3B-__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34-__it;					__t3D=__t35+__rt;
		__t34=__t34+__it;					__t35=__t35-__rt;

		__rt =(__t3A+__t3B)*ISRT2;			__it =(__t3A-__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33+__it;
		__t32=__t32+__rt;					__t33=__t33-__it;

		__rt =(__t3E-__t3F)*ISRT2;			__it =(__t3F+__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(__B+__odx[0x00])=__t00+__t20;			*(Bim+__odx[0x00])=__t01+__t21;
		*(__B+__odx[0x10])=__t00-__t20;			*(Bim+__odx[0x10])=__t01-__t21;
		*(__B+__odx[0x08])=__t10+__t31;			*(Bim+__odx[0x08])=__t11-__t30;
		*(__B+__odx[0x18])=__t10-__t31;			*(Bim+__odx[0x18])=__t11+__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08-__t19;					__t08=__t08+__t19;
		__t19=__t09+__rt;					__t09=__t09-__rt;

		__rt =(__t29+__t28)*ISRT2;			__t29=(__t29-__t28)*ISRT2;	__t28=__rt;
		__rt =(__t38-__t39)*ISRT2;			__it =(__t38+__t39)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		*(__B+__odx[0x04])=__t08+__t28;			*(Bim+__odx[0x04])=__t09+__t29;
		*(__B+__odx[0x14])=__t08-__t28;			*(Bim+__odx[0x14])=__t09-__t29;
		*(__B+__odx[0x0C])=__t18+__t39;			*(Bim+__odx[0x0C])=__t19-__t38;
		*(__B+__odx[0x1C])=__t18-__t39;			*(Bim+__odx[0x1C])=__t19+__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t15+__t14)*ISRT2;			__it =(__t15-__t14)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 + __t25*s16;			__t25=__t25*c16 - __t24*s16;	__t24=__rt;
		__rt =__t34*s16 + __t35*c16;			__it =__t35*s16 - __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x02])=__t04+__t24;			*(Bim+__odx[0x02])=__t05+__t25;
		*(__B+__odx[0x12])=__t04-__t24;			*(Bim+__odx[0x12])=__t05-__t25;
		*(__B+__odx[0x0A])=__t14+__t35;			*(Bim+__odx[0x0A])=__t15-__t34;
		*(__B+__odx[0x1A])=__t14-__t35;			*(Bim+__odx[0x1A])=__t15+__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1C-__t1D)*ISRT2;			__it =(__t1C+__t1D)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 + __t2D*c16;			__t2D=__t2D*s16 - __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 + __t3D*s16;			__it =__t3D*c16 - __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x06])=__t0C+__t2C;			*(Bim+__odx[0x06])=__t0D+__t2D;
		*(__B+__odx[0x16])=__t0C-__t2C;			*(Bim+__odx[0x16])=__t0D-__t2D;
		*(__B+__odx[0x0E])=__t1C+__t3D;			*(Bim+__odx[0x0E])=__t1D-__t3C;
		*(__B+__odx[0x1E])=__t1C-__t3D;			*(Bim+__odx[0x1E])=__t1D+__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 + __t13*s16;			__it =__t13*c16 - __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 + __t23*s32_1;	__t23=__t23*c32_1 - __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 + __t33*s32_3;	__it =__t33*c32_3 - __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x01])=__t02+__t22;			*(Bim+__odx[0x01])=__t03+__t23;
		*(__B+__odx[0x11])=__t02-__t22;			*(Bim+__odx[0x11])=__t03-__t23;
		*(__B+__odx[0x09])=__t12+__t33;			*(Bim+__odx[0x09])=__t13-__t32;
		*(__B+__odx[0x19])=__t12-__t33;			*(Bim+__odx[0x19])=__t13+__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 - __t1B*c16;			__it =__t1B*s16 + __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 + __t2B*c32_3;	__t2B=__t2B*s32_3 - __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 - __t3B*s32_1;	__it =__t3B*c32_1 + __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x05])=__t0A+__t2A;			*(Bim+__odx[0x05])=__t0B+__t2B;
		*(__B+__odx[0x15])=__t0A-__t2A;			*(Bim+__odx[0x15])=__t0B-__t2B;
		*(__B+__odx[0x0D])=__t1A+__t3B;			*(Bim+__odx[0x0D])=__t1B-__t3A;
		*(__B+__odx[0x1D])=__t1A-__t3B;			*(Bim+__odx[0x1D])=__t1B+__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 + __t17*c16;			__it =__t17*s16 - __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 + __t27*s32_3;	__t27=__t27*c32_3 - __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 - __t37*c32_1;	__it =__t37*s32_1 + __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x03])=__t06+__t26;			*(Bim+__odx[0x03])=__t07+__t27;
		*(__B+__odx[0x13])=__t06-__t26;			*(Bim+__odx[0x13])=__t07-__t27;
		*(__B+__odx[0x0B])=__t16+__t37;			*(Bim+__odx[0x0B])=__t17-__t36;
		*(__B+__odx[0x1B])=__t16-__t37;			*(Bim+__odx[0x1B])=__t17+__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 - __t1F*s16;			__it =__t1F*c16 + __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 + __t2F*c32_1;	__t2F=__t2F*s32_1 - __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 + __t3F*c32_3;	__it =__t3F*s32_3 - __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x07])=__t0E+__t2E;			*(Bim+__odx[0x07])=__t0F+__t2F;
		*(__B+__odx[0x17])=__t0E-__t2E;			*(Bim+__odx[0x17])=__t0F-__t2F;
		*(__B+__odx[0x0F])=__t1E+__t3F;			*(Bim+__odx[0x0F])=__t1F-__t3E;
		*(__B+__odx[0x1F])=__t1E-__t3F;			*(Bim+__odx[0x1F])=__t1F+__t3E;
}

// Post-twiddles out-of-place-permitting analog of above twiddleless DIT macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIT version of this macro processes the twiddles in bit-reversed order: 0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f.
void RADIX_32_DIT_TWIDDLE(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c10,double __s10,
	double __c08,double __s08,
	double __c18,double __s18,
	/**/
	double __c04,double __s04,
	double __c14,double __s14,
	double __c0C,double __s0C,
	double __c1C,double __s1C,
	/**/
	double __c02,double __s02,
	double __c12,double __s12,
	double __c0A,double __s0A,
	double __c1A,double __s1A,
	/**/
	double __c06,double __s06,
	double __c16,double __s16,
	double __c0E,double __s0E,
	double __c1E,double __s1E,
	/**/
	double __c01,double __s01,
	double __c11,double __s11,
	double __c09,double __s09,
	double __c19,double __s19,
	/**/
	double __c05,double __s05,
	double __c15,double __s15,
	double __c0D,double __s0D,
	double __c1D,double __s1D,
	/**/
	double __c03,double __s03,
	double __c13,double __s13,
	double __c0B,double __s0B,
	double __c1B,double __s1B,
	/**/
	double __c07,double __s07,
	double __c17,double __s17,
	double __c0F,double __s0F,
	double __c1F,double __s1F
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Aim = __A + RE_IM_STRIDE, *Bim = __B + RE_IM_STRIDE;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);
		__rt =*(__A+__idx[0x01]);			__it =*(Aim+__idx[0x01]);
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(__A+__idx[0x02]);			__t05=*(Aim+__idx[0x02]);
		__rt =*(__A+__idx[0x03]);			__it =*(Aim+__idx[0x03]);
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02-__it;					__t07=__t03+__rt;
		__t02=__t02+__it;					__t03=__t03-__rt;

		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);
		__rt =*(__A+__idx[0x05]);			__it =*(Aim+__idx[0x05]);
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(__A+__idx[0x06]);			__t0D=*(Aim+__idx[0x06]);
		__rt =*(__A+__idx[0x07]);			__it =*(Aim+__idx[0x07]);
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A-__it;					__t0F=__t0B+__rt;
		__t0A=__t0A+__it;					__t0B=__t0B-__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04-__it;					__t0D=__t05+__rt;
		__t04=__t04+__it;					__t05=__t05-__rt;

		__rt =(__t0A+__t0B)*ISRT2;			__it =(__t0A-__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03+__it;
		__t02=__t02+__rt;					__t03=__t03-__it;

		__rt =(__t0E-__t0F)*ISRT2;			__it =(__t0F+__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x08]);			__t11=*(Aim+__idx[0x08]);
		__rt =*(__A+__idx[0x09]);			__it =*(Aim+__idx[0x09]);
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);
		__rt =*(__A+__idx[0x0B]);			__it =*(Aim+__idx[0x0B]);
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12-__it;					__t17=__t13+__rt;
		__t12=__t12+__it;					__t13=__t13-__rt;

		__t18=*(__A+__idx[0x0C]);			__t19=*(Aim+__idx[0x0C]);
		__rt =*(__A+__idx[0x0D]);			__it =*(Aim+__idx[0x0D]);
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);
		__rt =*(__A+__idx[0x0F]);			__it =*(Aim+__idx[0x0F]);
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A-__it;					__t1F=__t1B+__rt;
		__t1A=__t1A+__it;					__t1B=__t1B-__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14-__it;					__t1D=__t15+__rt;
		__t14=__t14+__it;					__t15=__t15-__rt;

		__rt =(__t1A+__t1B)*ISRT2;			__it =(__t1A-__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13+__it;
		__t12=__t12+__rt;					__t13=__t13-__it;

		__rt =(__t1E-__t1F)*ISRT2;			__it =(__t1F+__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x10]);			__t21=*(Aim+__idx[0x10]);
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(__A+__idx[0x12]);			__t25=*(Aim+__idx[0x12]);
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22-__it;					__t27=__t23+__rt;
		__t22=__t22+__it;					__t23=__t23-__rt;

		__t28=*(__A+__idx[0x14]);			__t29=*(Aim+__idx[0x14]);
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(__A+__idx[0x16]);			__t2D=*(Aim+__idx[0x16]);
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A-__it;					__t2F=__t2B+__rt;
		__t2A=__t2A+__it;					__t2B=__t2B-__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24-__it;					__t2D=__t25+__rt;
		__t24=__t24+__it;					__t25=__t25-__rt;

		__rt =(__t2A+__t2B)*ISRT2;			__it =(__t2A-__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23+__it;
		__t22=__t22+__rt;					__t23=__t23-__it;

		__rt =(__t2E-__t2F)*ISRT2;			__it =(__t2F+__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x18]);			__t31=*(Aim+__idx[0x18]);
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(__A+__idx[0x1A]);			__t35=*(Aim+__idx[0x1A]);
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32-__it;					__t37=__t33+__rt;
		__t32=__t32+__it;					__t33=__t33-__rt;

		__t38=*(__A+__idx[0x1C]);			__t39=*(Aim+__idx[0x1C]);
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(__A+__idx[0x1E]);			__t3D=*(Aim+__idx[0x1E]);
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A-__it;					__t3F=__t3B+__rt;
		__t3A=__t3A+__it;					__t3B=__t3B-__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34-__it;					__t3D=__t35+__rt;
		__t34=__t34+__it;					__t35=__t35-__rt;

		__rt =(__t3A+__t3B)*ISRT2;			__it =(__t3A-__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33+__it;
		__t32=__t32+__rt;					__t33=__t33-__it;

		__rt =(__t3E-__t3F)*ISRT2;			__it =(__t3F+__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;
	/**********************************************************************************/
	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/**********************************************************************************/
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(__B+__odx[0x00])=__t00+__t20;				*(Bim+__odx[0x00])=__t01+__t21;
		__t00=__t00-__t20;		__t01  =__t01-__t21;
		*(__B+__odx[0x10])=__t00*__c10+__t01*__s10;	*(Bim+__odx[0x10])=__t01*__c10-__t00*__s10;

		__rt =__t10+__t31;		__it   =__t11-__t30;	/* mpy by E^-4 = -I is inlined here...	*/
		__t10=__t10-__t31;		__t11  =__t11+__t30;
		*(__B+__odx[0x08])=__rt *__c08+__it *__s08;	*(Bim+__odx[0x08])=__it *__c08-__rt *__s08;
		*(__B+__odx[0x18])=__t10*__c18+__t11*__s18;	*(Bim+__odx[0x18])=__t11*__c18-__t10*__s18;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;	__t18=__t08-__t19;	__t08=__t08+__t19;		/* twiddle mpy by E^8 =-I	*/
						__t19=__t09+__rt;	__t09=__t09-__rt;

		__rt =(__t29+__t28)*ISRT2;	__t29=(__t29-__t28)*ISRT2;		__t28=__rt;	/* twiddle mpy by E^-4	*/
		__rt =(__t38-__t39)*ISRT2;	__it =(__t38+__t39)*ISRT2;		/* twiddle mpy by E^4 = -E^-12 is here...	*/
		__t38=__t28+__rt;			__t28=__t28-__rt;				/* ...and get E^-12 by flipping signs here.	*/
		__t39=__t29+__it;			__t29=__t29-__it;

		__rt =__t08+__t28;		__it   =__t09+__t29;
		__t08=__t08-__t28;		__t09  =__t09-__t29;
		*(__B+__odx[0x04])=__rt *__c04+__it *__s04;	*(Bim+__odx[0x04])=__it *__c04-__rt *__s04;
		*(__B+__odx[0x14])=__t08*__c14+__t09*__s14;	*(Bim+__odx[0x14])=__t09*__c14-__t08*__s14;

		__rt =__t18+__t39;		__it   =__t19-__t38;	/* mpy by E^-4 = -I is inlined here...	*/
		__t18=__t18-__t39;		__t19  =__t19+__t38;
		*(__B+__odx[0x0C])=__rt *__c0C+__it *__s0C;	*(Bim+__odx[0x0C])=__it *__c0C-__rt *__s0C;
		*(__B+__odx[0x1C])=__t18*__c1C+__t19*__s1C;	*(Bim+__odx[0x1C])=__t19*__c1C-__t18*__s1C;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t15+__t14)*ISRT2;	__it =(__t15-__t14)*ISRT2;			/* twiddle mpy by E^-4	*/
		__t14=__t04-__rt;			__t04=__t04+__rt;
		__t15=__t05-__it;			__t05=__t05+__it;

		__rt =__t24*c16 + __t25*s16;		__t25=__t25*c16 - __t24*s16;		__t24=__rt;	/* twiddle mpy by E^-2	*/
		__rt =__t34*s16 + __t35*c16;		__it =__t35*s16 - __t34*c16;			/* twiddle mpy by E^-6	*/
		__t34=__t24-__rt;			__t24=__t24+__rt;
		__t35=__t25-__it;			__t25=__t25+__it;

		__rt =__t04+__t24;		__it   =__t05+__t25;
		__t04=__t04-__t24;		__t05  =__t05-__t25;
		*(__B+__odx[0x02])=__rt *__c02+__it *__s02;	*(Bim+__odx[0x02])=__it *__c02-__rt *__s02;
		*(__B+__odx[0x12])=__t04*__c12+__t05*__s12;	*(Bim+__odx[0x12])=__t05*__c12-__t04*__s12;

		__rt =__t14+__t35;		__it   =__t15-__t34;	/* mpy by E^-4 = -I is inlined here...	*/
		__t14=__t14-__t35;		__t15  =__t15+__t34;
		*(__B+__odx[0x0A])=__rt *__c0A+__it *__s0A;	*(Bim+__odx[0x0A])=__it *__c0A-__rt *__s0A;
		*(__B+__odx[0x1A])=__t14*__c1A+__t15*__s1A;	*(Bim+__odx[0x1A])=__t15*__c1A-__t14*__s1A;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1C-__t1D)*ISRT2;	__it =(__t1C+__t1D)*ISRT2;			/* twiddle mpy by E^4 = -E^-12 is here...	*/
		__t1C=__t0C+__rt;			__t0C=__t0C-__rt;				/* ...and get E^-12 by flipping signs here.	*/
		__t1D=__t0D+__it;			__t0D=__t0D-__it;

		__rt =__t2C*s16 + __t2D*c16;		__t2D=__t2D*s16 - __t2C*c16;		__t2C=__rt;	/* twiddle mpy by E^-6	*/
		__rt =__t3C*c16 + __t3D*s16;		__it =__t3D*c16 - __t3C*s16;			/* twiddle mpy by E^-18 is here...	*/
		__t3C=__t2C+__rt;			__t2C=__t2C-__rt;				/* ...and get E^-18 by flipping signs here.	*/
		__t3D=__t2D+__it;			__t2D=__t2D-__it;

		__rt =__t0C+__t2C;		__it   =__t0D+__t2D;
		__t0C=__t0C-__t2C;		__t0D  =__t0D-__t2D;
		*(__B+__odx[0x06])=__rt *__c06+__it *__s06;	*(Bim+__odx[0x06])=__it *__c06-__rt *__s06;
		*(__B+__odx[0x16])=__t0C*__c16+__t0D*__s16;	*(Bim+__odx[0x16])=__t0D*__c16-__t0C*__s16;

		__rt =__t1C+__t3D;		__it   =__t1D-__t3C;	/* mpy by E^-4 = -I is inlined here...	*/
		__t1C=__t1C-__t3D;		__t1D  =__t1D+__t3C;
		*(__B+__odx[0x0E])=__rt *__c0E+__it *__s0E;	*(Bim+__odx[0x0E])=__it *__c0E-__rt *__s0E;
		*(__B+__odx[0x1E])=__t1C*__c1E+__t1D*__s1E;	*(Bim+__odx[0x1E])=__t1D*__c1E-__t1C*__s1E;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 + __t13*s16;			__it =__t13*c16 - __t12*s16;
		__t12=__t02-__rt;			__t02=__t02+__rt;
		__t13=__t03-__it;			__t03=__t03+__it;

		__rt =__t22*c32_1 + __t23*s32_1;	__t23=__t23*c32_1 - __t22*s32_1;	__t22=__rt;	/* twiddle mpy by E^-1	*/
		__rt =__t32*c32_3 + __t33*s32_3;	__it =__t33*c32_3 - __t32*s32_3;		/* twiddle mpy by E^-3	*/
		__t32=__t22-__rt;			__t22=__t22+__rt;
		__t33=__t23-__it;			__t23=__t23+__it;

		__rt =__t02+__t22;		__it   =__t03+__t23;
		__t02=__t02-__t22;		__t03  =__t03-__t23;
		*(__B+__odx[0x01])=__rt *__c01+__it *__s01;	*(Bim+__odx[0x01])=__it *__c01-__rt *__s01;
		*(__B+__odx[0x11])=__t02*__c11+__t03*__s11;	*(Bim+__odx[0x11])=__t03*__c11-__t02*__s11;

		__rt =__t12+__t33;		__it   =__t13-__t32;	/* mpy by E^-4 = -I is inlined here...	*/
		__t12=__t12-__t33;		__t13  =__t13+__t32;
		*(__B+__odx[0x09])=__rt *__c09+__it *__s09;	*(Bim+__odx[0x09])=__it *__c09-__rt *__s09;
		*(__B+__odx[0x19])=__t12*__c19+__t13*__s19;	*(Bim+__odx[0x19])=__t13*__c19-__t12*__s19;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 - __t1B*c16;			__it =__t1B*s16 + __t1A*c16;
		__t1A=__t0A+__rt;			__t0A =__t0A-__rt;				/* ...and get E^-10 by flipping signs here.	*/
		__t1B=__t0B+__it;			__t0B =__t0B-__it;

		__rt =__t2A*s32_3 + __t2B*c32_3;	__t2B=__t2B*s32_3 - __t2A*c32_3;	__t2A=__rt;	/* twiddle mpy by E^-5	*/
		__rt =__t3A*c32_1 - __t3B*s32_1;	__it =__t3B*c32_1 + __t3A*s32_1;		/* twiddle mpy by -E^-15 is here...	*/
		__t3A=__t2A+__rt;			__t2A=__t2A-__rt;				/* ...and get E^-15 by flipping signs here.	*/
		__t3B=__t2B+__it;			__t2B=__t2B-__it;

		__rt =__t0A+__t2A;		__it   =__t0B+__t2B;
		__t0A=__t0A-__t2A;		__t0B  =__t0B-__t2B;
		*(__B+__odx[0x05])=__rt *__c05+__it *__s05;	*(Bim+__odx[0x05])=__it *__c05-__rt *__s05;
		*(__B+__odx[0x15])=__t0A*__c15+__t0B*__s15;	*(Bim+__odx[0x15])=__t0B*__c15-__t0A*__s15;

		__rt =__t1A+__t3B;		__it   =__t1B-__t3A;	/* mpy by E^-4 = -I is inlined here...	*/
		__t1A=__t1A-__t3B;		__t1B  =__t1B+__t3A;
		*(__B+__odx[0x0D])=__rt *__c0D+__it *__s0D;	*(Bim+__odx[0x0D])=__it *__c0D-__rt *__s0D;
		*(__B+__odx[0x1D])=__t1A*__c1D+__t1B*__s1D;	*(Bim+__odx[0x1D])=__t1B*__c1D-__t1A*__s1D;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 + __t17*c16;			__it =__t17*s16 - __t16*c16;
		__t16=__t06-__rt;			__t06 =__t06+__rt;
		__t17=__t07-__it;			__t07 =__t07+__it;

		__rt =__t26*c32_3 + __t27*s32_3;	__t27=__t27*c32_3 - __t26*s32_3;	__t26=__rt;	/* twiddle mpy by E^-3	*/
		__rt =__t36*s32_1 - __t37*c32_1;	__it =__t37*s32_1 + __t36*c32_1;		/* twiddle mpy by -E^-9 is here...	*/
		__t36=__t26+__rt;			__t26=__t26-__rt;				/* ...and get E^-9 by flipping signs here.	*/
		__t37=__t27+__it;			__t27=__t27-__it;

		__rt =__t06+__t26;		__it   =__t07+__t27;
		__t06=__t06-__t26;		__t07  =__t07-__t27;
		*(__B+__odx[0x03])=__rt *__c03+__it *__s03;	*(Bim+__odx[0x03])=__it *__c03-__rt *__s03;
		*(__B+__odx[0x13])=__t06*__c13+__t07*__s13;	*(Bim+__odx[0x13])=__t07*__c13-__t06*__s13;

		__rt =__t16+__t37;		__it   =__t17-__t36;	/* mpy by E^-4 = -I is inlined here...	*/
		__t16=__t16-__t37;		__t17  =__t17+__t36;
		*(__B+__odx[0x0B])=__rt *__c0B+__it *__s0B;	*(Bim+__odx[0x0B])=__it *__c0B-__rt *__s0B;
		*(__B+__odx[0x1B])=__t16*__c1B+__t17*__s1B;	*(Bim+__odx[0x1B])=__t17*__c1B-__t16*__s1B;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 - __t1F*s16;			__it =__t1F*c16 + __t1E*s16;
		__t1E=__t0E+__rt;			__t0E =__t0E-__rt;				/* ...and get E^-14 by flipping signs here.	*/
		__t1F=__t0F+__it;			__t0F =__t0F-__it;

		__rt =__t2E*s32_1 + __t2F*c32_1;	__t2F=__t2F*s32_1 - __t2E*c32_1;	__t2E=__rt;	/* twiddle mpy by E^-7	*/
		__rt =__t3E*s32_3 + __t3F*c32_3;	__it =__t3F*s32_3 - __t3E*c32_3;		/* twiddle mpy by -E^-21 is here...	*/
		__t3E=__t2E+__rt;			__t2E=__t2E-__rt;				/* ...and get E^-21 by flipping signs here.	*/
		__t3F=__t2F+__it;			__t2F=__t2F-__it;

		__rt =__t0E+__t2E;		__it   =__t0F+__t2F;
		__t0E=__t0E-__t2E;		__t0F  =__t0F-__t2F;
		*(__B+__odx[0x07])=__rt *__c07+__it *__s07;	*(Bim+__odx[0x07])=__it *__c07-__rt *__s07;
		*(__B+__odx[0x17])=__t0E*__c17+__t0F*__s17;	*(Bim+__odx[0x17])=__t0F*__c17-__t0E*__s17;

		__rt =__t1E+__t3F;		__it   =__t1F-__t3E;	/* mpy by E^-4 = -I is inlined here...	*/
		__t1E=__t1E-__t3F;		__t1F  =__t1F+__t3E;
		*(__B+__odx[0x0F])=__rt *__c0F+__it *__s0F;	*(Bim+__odx[0x0F])=__it *__c0F-__rt *__s0F;
		*(__B+__odx[0x1F])=__t1E*__c1F+__t1F*__s1F;	*(Bim+__odx[0x1F])=__t1F*__c1F-__t1E*__s1F;
}

// With-twiddles out-of-place analog of above twiddleless DIT macro: 31 nontrivial complex input twiddles E01-1f [E0 assumed = 1],
// The DIT version of this macro processes the twiddles in order.
// *** NOTE: This macro is a pre-twiddles DIT, designed as the 2nd level in a large-power-of-2 DIT DFT, e.g. radix-1024. **********
// *** It is NOT suitable as standalone radix-32 pass, since my FFT does such passes on the iFFT (DIT) side with post-twiddles. ***
void RADIX_32_DIT_TWIDDLE_OOP(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx,	/* Outputs: Base address plus 32 (index) offsets */
	/* Twiddles: */
	double __c01,double __s01,
	double __c02,double __s02,
	double __c03,double __s03,
	double __c04,double __s04,
	double __c05,double __s05,
	double __c06,double __s06,
	double __c07,double __s07,
	double __c08,double __s08,
	double __c09,double __s09,
	double __c0A,double __s0A,
	double __c0B,double __s0B,
	double __c0C,double __s0C,
	double __c0D,double __s0D,
	double __c0E,double __s0E,
	double __c0F,double __s0F,
	double __c10,double __s10,
	double __c11,double __s11,
	double __c12,double __s12,
	double __c13,double __s13,
	double __c14,double __s14,
	double __c15,double __s15,
	double __c16,double __s16,
	double __c17,double __s17,
	double __c18,double __s18,
	double __c19,double __s19,
	double __c1A,double __s1A,
	double __c1B,double __s1B,
	double __c1C,double __s1C,
	double __c1D,double __s1D,
	double __c1E,double __s1E,
	double __c1F,double __s1F
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Aim = __A + RE_IM_STRIDE, *Bim = __B + RE_IM_STRIDE;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);
		__t02=*(__A+__idx[0x01]);			__t03=*(Aim+__idx[0x01]);	__rt = __t02*__c01 +__t03*__s01;	__it = __t03*__c01 -__t02*__s01;
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(__A+__idx[0x02]);			__t05=*(Aim+__idx[0x02]);	__rt = __t04*__c02 +__t05*__s02;	__t05= __t05*__c02 -__t04*__s02;	__t04 = __rt;
		__t06=*(__A+__idx[0x03]);			__t07=*(Aim+__idx[0x03]);	__rt = __t06*__c03 +__t07*__s03;	__it = __t07*__c03 -__t06*__s03;
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02-__it;					__t07=__t03+__rt;
		__t02=__t02+__it;					__t03=__t03-__rt;

		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);	__rt = __t08*__c04 +__t09*__s04;	__t09= __t09*__c04 -__t08*__s04;	__t08 = __rt;
		__t0A=*(__A+__idx[0x05]);			__t0B=*(Aim+__idx[0x05]);	__rt = __t0A*__c05 +__t0B*__s05;	__it = __t0B*__c05 -__t0A*__s05;
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(__A+__idx[0x06]);			__t0D=*(Aim+__idx[0x06]);	__rt = __t0C*__c06 +__t0D*__s06;	__t0D= __t0D*__c06 -__t0C*__s06;	__t0C = __rt;
		__t0E=*(__A+__idx[0x07]);			__t0F=*(Aim+__idx[0x07]);	__rt = __t0E*__c07 +__t0F*__s07;	__it = __t0F*__c07 -__t0E*__s07;
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A-__it;					__t0F=__t0B+__rt;
		__t0A=__t0A+__it;					__t0B=__t0B-__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04-__it;					__t0D=__t05+__rt;
		__t04=__t04+__it;					__t05=__t05-__rt;

		__rt =(__t0A+__t0B)*ISRT2;			__it =(__t0A-__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03+__it;
		__t02=__t02+__rt;					__t03=__t03-__it;

		__rt =(__t0E-__t0F)*ISRT2;			__it =(__t0F+__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x08]);			__t11=*(Aim+__idx[0x08]);	__rt = __t10*__c08 +__t11*__s08;	__t11= __t11*__c08 -__t10*__s08;	__t10 = __rt;
		__t12=*(__A+__idx[0x09]);			__t13=*(Aim+__idx[0x09]);	__rt = __t12*__c09 +__t13*__s09;	__it = __t13*__c09 -__t12*__s09;
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);	__rt = __t14*__c0A +__t15*__s0A;	__t15= __t15*__c0A -__t14*__s0A;	__t14 = __rt;
		__t16=*(__A+__idx[0x0B]);			__t17=*(Aim+__idx[0x0B]);	__rt = __t16*__c0B +__t17*__s0B;	__it = __t17*__c0B -__t16*__s0B;
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12-__it;					__t17=__t13+__rt;
		__t12=__t12+__it;					__t13=__t13-__rt;

		__t18=*(__A+__idx[0x0C]);			__t19=*(Aim+__idx[0x0C]);	__rt = __t18*__c0C +__t19*__s0C;	__t19= __t19*__c0C -__t18*__s0C;	__t18 = __rt;
		__t1A=*(__A+__idx[0x0D]);			__t1B=*(Aim+__idx[0x0D]);	__rt = __t1A*__c0D +__t1B*__s0D;	__it = __t1B*__c0D -__t1A*__s0D;
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);	__rt = __t1C*__c0E +__t1D*__s0E;	__t1D= __t1D*__c0E -__t1C*__s0E;	__t1C = __rt;
		__t1E=*(__A+__idx[0x0F]);			__t1F=*(Aim+__idx[0x0F]);	__rt = __t1E*__c0F +__t1F*__s0F;	__it = __t1F*__c0F -__t1E*__s0F;
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A-__it;					__t1F=__t1B+__rt;
		__t1A=__t1A+__it;					__t1B=__t1B-__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14-__it;					__t1D=__t15+__rt;
		__t14=__t14+__it;					__t15=__t15-__rt;

		__rt =(__t1A+__t1B)*ISRT2;			__it =(__t1A-__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13+__it;
		__t12=__t12+__rt;					__t13=__t13-__it;

		__rt =(__t1E-__t1F)*ISRT2;			__it =(__t1F+__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x10]);			__t21=*(Aim+__idx[0x10]);	__rt = __t20*__c10 +__t21*__s10;	__t21= __t21*__c10 -__t20*__s10;	__t20 = __rt;
		__t22=*(__A+__idx[0x11]);			__t23=*(Aim+__idx[0x11]);	__rt = __t22*__c11 +__t23*__s11;	__it = __t23*__c11 -__t22*__s11;
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(__A+__idx[0x12]);			__t25=*(Aim+__idx[0x12]);	__rt = __t24*__c12 +__t25*__s12;	__t25= __t25*__c12 -__t24*__s12;	__t24 = __rt;
		__t26=*(__A+__idx[0x13]);			__t27=*(Aim+__idx[0x13]);	__rt = __t26*__c13 +__t27*__s13;	__it = __t27*__c13 -__t26*__s13;
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22-__it;					__t27=__t23+__rt;
		__t22=__t22+__it;					__t23=__t23-__rt;

		__t28=*(__A+__idx[0x14]);			__t29=*(Aim+__idx[0x14]);	__rt = __t28*__c14 +__t29*__s14;	__t29= __t29*__c14 -__t28*__s14;	__t28 = __rt;
		__t2A=*(__A+__idx[0x15]);			__t2B=*(Aim+__idx[0x15]);	__rt = __t2A*__c15 +__t2B*__s15;	__it = __t2B*__c15 -__t2A*__s15;
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(__A+__idx[0x16]);			__t2D=*(Aim+__idx[0x16]);	__rt = __t2C*__c16 +__t2D*__s16;	__t2D= __t2D*__c16 -__t2C*__s16;	__t2C = __rt;
		__t2E=*(__A+__idx[0x17]);			__t2F=*(Aim+__idx[0x17]);	__rt = __t2E*__c17 +__t2F*__s17;	__it = __t2F*__c17 -__t2E*__s17;
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A-__it;					__t2F=__t2B+__rt;
		__t2A=__t2A+__it;					__t2B=__t2B-__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24-__it;					__t2D=__t25+__rt;
		__t24=__t24+__it;					__t25=__t25-__rt;

		__rt =(__t2A+__t2B)*ISRT2;			__it =(__t2A-__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23+__it;
		__t22=__t22+__rt;					__t23=__t23-__it;

		__rt =(__t2E-__t2F)*ISRT2;			__it =(__t2F+__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x18]);			__t31=*(Aim+__idx[0x18]);	__rt = __t30*__c18 +__t31*__s18;	__t31= __t31*__c18 -__t30*__s18;	__t30 = __rt;
		__t32=*(__A+__idx[0x19]);			__t33=*(Aim+__idx[0x19]);	__rt = __t32*__c19 +__t33*__s19;	__it = __t33*__c19 -__t32*__s19;
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(__A+__idx[0x1A]);			__t35=*(Aim+__idx[0x1A]);	__rt = __t34*__c1A +__t35*__s1A;	__t35= __t35*__c1A -__t34*__s1A;	__t34 = __rt;
		__t36=*(__A+__idx[0x1B]);			__t37=*(Aim+__idx[0x1B]);	__rt = __t36*__c1B +__t37*__s1B;	__it = __t37*__c1B -__t36*__s1B;
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32-__it;					__t37=__t33+__rt;
		__t32=__t32+__it;					__t33=__t33-__rt;

		__t38=*(__A+__idx[0x1C]);			__t39=*(Aim+__idx[0x1C]);	__rt = __t38*__c1C +__t39*__s1C;	__t39= __t39*__c1C -__t38*__s1C;	__t38 = __rt;
		__t3A=*(__A+__idx[0x1D]);			__t3B=*(Aim+__idx[0x1D]);	__rt = __t3A*__c1D +__t3B*__s1D;	__it = __t3B*__c1D -__t3A*__s1D;
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(__A+__idx[0x1E]);			__t3D=*(Aim+__idx[0x1E]);	__rt = __t3C*__c1E +__t3D*__s1E;	__t3D= __t3D*__c1E -__t3C*__s1E;	__t3C = __rt;
		__t3E=*(__A+__idx[0x1F]);			__t3F=*(Aim+__idx[0x1F]);	__rt = __t3E*__c1F +__t3F*__s1F;	__it = __t3F*__c1F -__t3E*__s1F;
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A-__it;					__t3F=__t3B+__rt;
		__t3A=__t3A+__it;					__t3B=__t3B-__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34-__it;					__t3D=__t35+__rt;
		__t34=__t34+__it;					__t35=__t35-__rt;

		__rt =(__t3A+__t3B)*ISRT2;			__it =(__t3A-__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33+__it;
		__t32=__t32+__rt;					__t33=__t33-__it;

		__rt =(__t3E-__t3F)*ISRT2;			__it =(__t3F+__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(__B+__odx[0x00])=__t00+__t20;			*(Bim+__odx[0x00])=__t01+__t21;
		*(__B+__odx[0x10])=__t00-__t20;			*(Bim+__odx[0x10])=__t01-__t21;
		*(__B+__odx[0x08])=__t10+__t31;			*(Bim+__odx[0x08])=__t11-__t30;
		*(__B+__odx[0x18])=__t10-__t31;			*(Bim+__odx[0x18])=__t11+__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08-__t19;					__t08=__t08+__t19;
		__t19=__t09+__rt;					__t09=__t09-__rt;

		__rt =(__t29+__t28)*ISRT2;			__t29=(__t29-__t28)*ISRT2;	__t28=__rt;
		__rt =(__t38-__t39)*ISRT2;			__it =(__t38+__t39)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		*(__B+__odx[0x04])=__t08+__t28;			*(Bim+__odx[0x04])=__t09+__t29;
		*(__B+__odx[0x14])=__t08-__t28;			*(Bim+__odx[0x14])=__t09-__t29;
		*(__B+__odx[0x0C])=__t18+__t39;			*(Bim+__odx[0x0C])=__t19-__t38;
		*(__B+__odx[0x1C])=__t18-__t39;			*(Bim+__odx[0x1C])=__t19+__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t15+__t14)*ISRT2;			__it =(__t15-__t14)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 + __t25*s16;			__t25=__t25*c16 - __t24*s16;	__t24=__rt;
		__rt =__t34*s16 + __t35*c16;			__it =__t35*s16 - __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x02])=__t04+__t24;			*(Bim+__odx[0x02])=__t05+__t25;
		*(__B+__odx[0x12])=__t04-__t24;			*(Bim+__odx[0x12])=__t05-__t25;
		*(__B+__odx[0x0A])=__t14+__t35;			*(Bim+__odx[0x0A])=__t15-__t34;
		*(__B+__odx[0x1A])=__t14-__t35;			*(Bim+__odx[0x1A])=__t15+__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1C-__t1D)*ISRT2;			__it =(__t1C+__t1D)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 + __t2D*c16;			__t2D=__t2D*s16 - __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 + __t3D*s16;			__it =__t3D*c16 - __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		*(__B+__odx[0x06])=__t0C+__t2C;			*(Bim+__odx[0x06])=__t0D+__t2D;
		*(__B+__odx[0x16])=__t0C-__t2C;			*(Bim+__odx[0x16])=__t0D-__t2D;
		*(__B+__odx[0x0E])=__t1C+__t3D;			*(Bim+__odx[0x0E])=__t1D-__t3C;
		*(__B+__odx[0x1E])=__t1C-__t3D;			*(Bim+__odx[0x1E])=__t1D+__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 + __t13*s16;			__it =__t13*c16 - __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 + __t23*s32_1;	__t23=__t23*c32_1 - __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 + __t33*s32_3;	__it =__t33*c32_3 - __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x01])=__t02+__t22;			*(Bim+__odx[0x01])=__t03+__t23;
		*(__B+__odx[0x11])=__t02-__t22;			*(Bim+__odx[0x11])=__t03-__t23;
		*(__B+__odx[0x09])=__t12+__t33;			*(Bim+__odx[0x09])=__t13-__t32;
		*(__B+__odx[0x19])=__t12-__t33;			*(Bim+__odx[0x19])=__t13+__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 - __t1B*c16;			__it =__t1B*s16 + __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 + __t2B*c32_3;	__t2B=__t2B*s32_3 - __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 - __t3B*s32_1;	__it =__t3B*c32_1 + __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x05])=__t0A+__t2A;			*(Bim+__odx[0x05])=__t0B+__t2B;
		*(__B+__odx[0x15])=__t0A-__t2A;			*(Bim+__odx[0x15])=__t0B-__t2B;
		*(__B+__odx[0x0D])=__t1A+__t3B;			*(Bim+__odx[0x0D])=__t1B-__t3A;
		*(__B+__odx[0x1D])=__t1A-__t3B;			*(Bim+__odx[0x1D])=__t1B+__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 + __t17*c16;			__it =__t17*s16 - __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 + __t27*s32_3;	__t27=__t27*c32_3 - __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 - __t37*c32_1;	__it =__t37*s32_1 + __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x03])=__t06+__t26;			*(Bim+__odx[0x03])=__t07+__t27;
		*(__B+__odx[0x13])=__t06-__t26;			*(Bim+__odx[0x13])=__t07-__t27;
		*(__B+__odx[0x0B])=__t16+__t37;			*(Bim+__odx[0x0B])=__t17-__t36;
		*(__B+__odx[0x1B])=__t16-__t37;			*(Bim+__odx[0x1B])=__t17+__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 - __t1F*s16;			__it =__t1F*c16 + __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 + __t2F*c32_1;	__t2F=__t2F*s32_1 - __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 + __t3F*c32_3;	__it =__t3F*s32_3 - __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		*(__B+__odx[0x07])=__t0E+__t2E;			*(Bim+__odx[0x07])=__t0F+__t2F;
		*(__B+__odx[0x17])=__t0E-__t2E;			*(Bim+__odx[0x17])=__t0F-__t2F;
		*(__B+__odx[0x0F])=__t1E+__t3F;			*(Bim+__odx[0x0F])=__t1F-__t3E;
		*(__B+__odx[0x1F])=__t1E-__t3F;			*(Bim+__odx[0x1F])=__t1F+__t3E;
}


// Twiddleless Radix-32 DIF subtransform macro for use by larger radix-32*k twiddleless-DFT macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't fundamentally matter whether the in/outputs are local
// scalars or array locations, but our implementation of this macro is driven by our twiddleless-DFT scheme
// in which the DIF step has the power-of-2 component macros following a set of odd-radix ones, thus with:
//
//	o __A-inputs read from a block of contiguous local-allocated storage with unit index stride;
//	o __B-outputs written to an array with arbitrary index strides encoded in the __odx auxiliary array.
//
void RADIX_32_DIF_OOP(
	double *__A,			/*  Inputs: Base address plus 32 (index) offsets */
	double *__B, const int *__odx	/* Outputs: Base address plus 32 (index) offsets */
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Are = (double *)__A, *Aim = Are + 1;
	double *Bre = (double *)__B, *Bim = Bre + RE_IM_STRIDE;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(Are+0x00);			__t01=*(Aim+0x00);
		__rt =*(Are+0x20);			__it =*(Aim+0x20);
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(Are+0x10);			__t05=*(Aim+0x10);
		__rt =*(Are+0x30);			__it =*(Aim+0x30);
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02+__it;					__t07=__t03-__rt;
		__t02=__t02-__it;					__t03=__t03+__rt;

		__t08=*(Are+0x08);			__t09=*(Aim+0x08);
		__rt =*(Are+0x28);			__it =*(Aim+0x28);
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(Are+0x18);			__t0D=*(Aim+0x18);
		__rt =*(Are+0x38);			__it =*(Aim+0x38);
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A+__it;					__t0F=__t0B-__rt;
		__t0A=__t0A-__it;					__t0B=__t0B+__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04+__it;					__t0D=__t05-__rt;
		__t04=__t04-__it;					__t05=__t05+__rt;

		__rt =(__t0A-__t0B)*ISRT2;			__it =(__t0A+__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03-__it;
		__t02=__t02+__rt;					__t03=__t03+__it;

		__rt =(__t0E+__t0F)*ISRT2;			__it =(__t0F-__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(Are+0x04);			__t11=*(Aim+0x04);
		__rt =*(Are+0x24);			__it =*(Aim+0x24);
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(Are+0x14);			__t15=*(Aim+0x14);
		__rt =*(Are+0x34);			__it =*(Aim+0x34);
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12+__it;					__t17=__t13-__rt;
		__t12=__t12-__it;					__t13=__t13+__rt;

		__t18=*(Are+0x0c);			__t19=*(Aim+0x0c);
		__rt =*(Are+0x2c);			__it =*(Aim+0x2c);
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(Are+0x1c);			__t1D=*(Aim+0x1c);
		__rt =*(Are+0x3c);			__it =*(Aim+0x3c);
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A+__it;					__t1F=__t1B-__rt;
		__t1A=__t1A-__it;					__t1B=__t1B+__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14+__it;					__t1D=__t15-__rt;
		__t14=__t14-__it;					__t15=__t15+__rt;

		__rt =(__t1A-__t1B)*ISRT2;			__it =(__t1A+__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13-__it;
		__t12=__t12+__rt;					__t13=__t13+__it;

		__rt =(__t1E+__t1F)*ISRT2;			__it =(__t1F-__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(Are+0x02);			__t21=*(Aim+0x02);
		__rt =*(Are+0x22);			__it =*(Aim+0x22);
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(Are+0x12);			__t25=*(Aim+0x12);
		__rt =*(Are+0x32);			__it =*(Aim+0x32);
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22+__it;					__t27=__t23-__rt;
		__t22=__t22-__it;					__t23=__t23+__rt;

		__t28=*(Are+0x0a);			__t29=*(Aim+0x0a);
		__rt =*(Are+0x2a);			__it =*(Aim+0x2a);
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(Are+0x1a);			__t2D=*(Aim+0x1a);
		__rt =*(Are+0x3a);			__it =*(Aim+0x3a);
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A+__it;					__t2F=__t2B-__rt;
		__t2A=__t2A-__it;					__t2B=__t2B+__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24+__it;					__t2D=__t25-__rt;
		__t24=__t24-__it;					__t25=__t25+__rt;

		__rt =(__t2A-__t2B)*ISRT2;			__it =(__t2A+__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23-__it;
		__t22=__t22+__rt;					__t23=__t23+__it;

		__rt =(__t2E+__t2F)*ISRT2;			__it =(__t2F-__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(Are+0x06);			__t31=*(Aim+0x06);
		__rt =*(Are+0x26);			__it =*(Aim+0x26);
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(Are+0x16);			__t35=*(Aim+0x16);
		__rt =*(Are+0x36);			__it =*(Aim+0x36);
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32+__it;					__t37=__t33-__rt;
		__t32=__t32-__it;					__t33=__t33+__rt;

		__t38=*(Are+0x0e);			__t39=*(Aim+0x0e);
		__rt =*(Are+0x2e);			__it =*(Aim+0x2e);
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(Are+0x1e);			__t3D=*(Aim+0x1e);
		__rt =*(Are+0x3e);			__it =*(Aim+0x3e);
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A+__it;					__t3F=__t3B-__rt;
		__t3A=__t3A-__it;					__t3B=__t3B+__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34+__it;					__t3D=__t35-__rt;
		__t34=__t34-__it;					__t35=__t35+__rt;

		__rt =(__t3A-__t3B)*ISRT2;			__it =(__t3A+__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33-__it;
		__t32=__t32+__rt;					__t33=__t33+__it;

		__rt =(__t3E+__t3F)*ISRT2;			__it =(__t3F-__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(Bre+__odx[0x00])=__t00+__t20;		*(Bim+__odx[0x00])=__t01+__t21;
		*(Bre+__odx[0x01])=__t00-__t20;		*(Bim+__odx[0x01])=__t01-__t21;
		*(Bre+__odx[0x02])=__t10-__t31;		*(Bim+__odx[0x02])=__t11+__t30;
		*(Bre+__odx[0x03])=__t10+__t31;		*(Bim+__odx[0x03])=__t11-__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08+__t19;					__t08=__t08-__t19;
		__t19=__t09-__rt;					__t09=__t09+__rt;

		__rt =(__t28-__t29)*ISRT2;			__t29=(__t28+__t29)*ISRT2;	__t28=__rt;
		__rt =(__t39+__t38)*ISRT2;			__it =(__t39-__t38)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		*(Bre+__odx[0x04])=__t08+__t28;		*(Bim+__odx[0x04])=__t09+__t29;
		*(Bre+__odx[0x05])=__t08-__t28;		*(Bim+__odx[0x05])=__t09-__t29;
		*(Bre+__odx[0x06])=__t18-__t39;		*(Bim+__odx[0x06])=__t19+__t38;
		*(Bre+__odx[0x07])=__t18+__t39;		*(Bim+__odx[0x07])=__t19-__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t14-__t15)*ISRT2;			__it =(__t14+__t15)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 - __t25*s16;	__t25=__t25*c16 + __t24*s16;	__t24=__rt;
		__rt =__t34*s16 - __t35*c16;	__it =__t35*s16 + __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		*(Bre+__odx[0x08])=__t04+__t24;		*(Bim+__odx[0x08])=__t05+__t25;
		*(Bre+__odx[0x09])=__t04-__t24;		*(Bim+__odx[0x09])=__t05-__t25;
		*(Bre+__odx[0x0A])=__t14-__t35;		*(Bim+__odx[0x0A])=__t15+__t34;
		*(Bre+__odx[0x0B])=__t14+__t35;		*(Bim+__odx[0x0B])=__t15-__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1D+__t1C)*ISRT2;			__it =(__t1D-__t1C)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 - __t2D*c16;	__t2D=__t2D*s16 + __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 - __t3D*s16;	__it =__t3D*c16 + __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		*(Bre+__odx[0x0C])=__t0C+__t2C;		*(Bim+__odx[0x0C])=__t0D+__t2D;
		*(Bre+__odx[0x0D])=__t0C-__t2C;		*(Bim+__odx[0x0D])=__t0D-__t2D;
		*(Bre+__odx[0x0E])=__t1C-__t3D;		*(Bim+__odx[0x0E])=__t1D+__t3C;
		*(Bre+__odx[0x0F])=__t1C+__t3D;		*(Bim+__odx[0x0F])=__t1D-__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 - __t13*s16;	__it =__t13*c16 + __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 - __t23*s32_1;	__t23=__t23*c32_1 + __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 - __t33*s32_3;	__it =__t33*c32_3 + __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		*(Bre+__odx[0x10])=__t02+__t22;		*(Bim+__odx[0x10])=__t03+__t23;
		*(Bre+__odx[0x11])=__t02-__t22;		*(Bim+__odx[0x11])=__t03-__t23;
		*(Bre+__odx[0x12])=__t12-__t33;		*(Bim+__odx[0x12])=__t13+__t32;
		*(Bre+__odx[0x13])=__t12+__t33;		*(Bim+__odx[0x13])=__t13-__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 + __t1B*c16;	__it =__t1B*s16 - __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 - __t2B*c32_3;	__t2B=__t2B*s32_3 + __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 + __t3B*s32_1;	__it =__t3B*c32_1 - __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		*(Bre+__odx[0x14])=__t0A+__t2A;		*(Bim+__odx[0x14])=__t0B+__t2B;
		*(Bre+__odx[0x15])=__t0A-__t2A;		*(Bim+__odx[0x15])=__t0B-__t2B;
		*(Bre+__odx[0x16])=__t1A-__t3B;		*(Bim+__odx[0x16])=__t1B+__t3A;
		*(Bre+__odx[0x17])=__t1A+__t3B;		*(Bim+__odx[0x17])=__t1B-__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 - __t17*c16;	__it =__t17*s16 + __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 - __t27*s32_3;	__t27=__t27*c32_3 + __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 + __t37*c32_1;	__it =__t37*s32_1 - __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		*(Bre+__odx[0x18])=__t06+__t26;		*(Bim+__odx[0x18])=__t07+__t27;
		*(Bre+__odx[0x19])=__t06-__t26;		*(Bim+__odx[0x19])=__t07-__t27;
		*(Bre+__odx[0x1A])=__t16-__t37;		*(Bim+__odx[0x1A])=__t17+__t36;
		*(Bre+__odx[0x1B])=__t16+__t37;		*(Bim+__odx[0x1B])=__t17-__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 + __t1F*s16;	__it =__t1F*c16 - __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 - __t2F*c32_1;	__t2F=__t2F*s32_1 + __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 - __t3F*c32_3;	__it =__t3F*s32_3 + __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		*(Bre+__odx[0x1C])=__t0E+__t2E;		*(Bim+__odx[0x1C])=__t0F+__t2F;
		*(Bre+__odx[0x1D])=__t0E-__t2E;		*(Bim+__odx[0x1D])=__t0F-__t2F;
		*(Bre+__odx[0x1E])=__t1E-__t3F;		*(Bim+__odx[0x1E])=__t1F+__t3E;
		*(Bre+__odx[0x1F])=__t1E+__t3F;		*(Bim+__odx[0x1F])=__t1F-__t3E;
}


// Twiddleless Radix-32 DIT subtransform macro for use by larger radix-8*k macros.
// OOP = out of place, i.e. Assumes output locs != input locs.
// For the scalar-data versino of this macro it doesn't matter whether the in/outputs are local
// scalars or array locations, but the __A-naming indicates that this macro is intended to serve
// as a prototype for a SIMD/ASM macro in which the __A-inputs are read from an array with arbitrary
// index stride and the __t-outputs go into a block of contiguous local-allocated storage:
void RADIX_32_DIT_OOP(
	double *__A, const int *__idx,	/*  Inputs: Base address plus 32 (index) offsets */
	double *__B					/* Outputs: Base address plus 32 (index) offsets */
)
{
	double __rt,__it
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F;
	double *Are = (double *)__A, *Aim = Are + RE_IM_STRIDE;
	double *Bre = (double *)__B, *Bim = Bre + 1;
	/* Gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...	*/
	/* Each complex radix-8 subtransform needs 52 ADD, 4 MUL (not counting load/store/address-compute): */
	/*...Block 1: */
		__t00=*(__A+__idx[0x00]);			__t01=*(Aim+__idx[0x00]);
		__rt =*(__A+__idx[0x01]);			__it =*(Aim+__idx[0x01]);
		__t02=__t00-__rt;					__t03=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__t04=*(__A+__idx[0x02]);			__t05=*(Aim+__idx[0x02]);
		__rt =*(__A+__idx[0x03]);			__it =*(Aim+__idx[0x03]);
		__t06=__t04-__rt;					__t07=__t05-__it;
		__t04=__t04+__rt;					__t05=__t05+__it;

		__rt =__t04;						__it =__t05;
		__t04=__t00-__rt;					__t05=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t06;						__it =__t07;
		__t06=__t02-__it;					__t07=__t03+__rt;
		__t02=__t02+__it;					__t03=__t03-__rt;

		__t08=*(__A+__idx[0x04]);			__t09=*(Aim+__idx[0x04]);
		__rt =*(__A+__idx[0x05]);			__it =*(Aim+__idx[0x05]);
		__t0A=__t08-__rt;					__t0B=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__t0C=*(__A+__idx[0x06]);			__t0D=*(Aim+__idx[0x06]);
		__rt =*(__A+__idx[0x07]);			__it =*(Aim+__idx[0x07]);
		__t0E=__t0C-__rt;					__t0F=__t0D-__it;
		__t0C=__t0C+__rt;					__t0D=__t0D+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t08-__rt;					__t0D=__t09-__it;
		__t08=__t08+__rt;					__t09=__t09+__it;

		__rt =__t0E;						__it =__t0F;
		__t0E=__t0A-__it;					__t0F=__t0B+__rt;
		__t0A=__t0A+__it;					__t0B=__t0B-__rt;

		__rt =__t08;						__it =__t09;
		__t08=__t00-__rt;					__t09=__t01-__it;
		__t00=__t00+__rt;					__t01=__t01+__it;

		__rt =__t0C;						__it =__t0D;
		__t0C=__t04-__it;					__t0D=__t05+__rt;
		__t04=__t04+__it;					__t05=__t05-__rt;

		__rt =(__t0A+__t0B)*ISRT2;			__it =(__t0A-__t0B)*ISRT2;
		__t0A=__t02-__rt;					__t0B=__t03+__it;
		__t02=__t02+__rt;					__t03=__t03-__it;

		__rt =(__t0E-__t0F)*ISRT2;			__it =(__t0F+__t0E)*ISRT2;
		__t0E=__t06+__rt;					__t0F=__t07+__it;
		__t06=__t06-__rt;					__t07=__t07-__it;

	/*...Block 2:;*/
		__t10=*(__A+__idx[0x08]);			__t11=*(Aim+__idx[0x08]);
		__rt =*(__A+__idx[0x09]);			__it =*(Aim+__idx[0x09]);
		__t12=__t10-__rt;					__t13=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__t14=*(__A+__idx[0x0A]);			__t15=*(Aim+__idx[0x0A]);
		__rt =*(__A+__idx[0x0B]);			__it =*(Aim+__idx[0x0B]);
		__t16=__t14-__rt;					__t17=__t15-__it;
		__t14=__t14+__rt;					__t15=__t15+__it;

		__rt =__t14;						__it =__t15;
		__t14=__t10-__rt;					__t15=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t16;						__it =__t17;
		__t16=__t12-__it;					__t17=__t13+__rt;
		__t12=__t12+__it;					__t13=__t13-__rt;

		__t18=*(__A+__idx[0x0C]);			__t19=*(Aim+__idx[0x0C]);
		__rt =*(__A+__idx[0x0D]);			__it =*(Aim+__idx[0x0D]);
		__t1A=__t18-__rt;					__t1B=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__t1C=*(__A+__idx[0x0E]);			__t1D=*(Aim+__idx[0x0E]);
		__rt =*(__A+__idx[0x0F]);			__it =*(Aim+__idx[0x0F]);
		__t1E=__t1C-__rt;					__t1F=__t1D-__it;
		__t1C=__t1C+__rt;					__t1D=__t1D+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t18-__rt;					__t1D=__t19-__it;
		__t18=__t18+__rt;					__t19=__t19+__it;

		__rt =__t1E;						__it =__t1F;
		__t1E=__t1A-__it;					__t1F=__t1B+__rt;
		__t1A=__t1A+__it;					__t1B=__t1B-__rt;

		__rt =__t18;						__it =__t19;
		__t18=__t10-__rt;					__t19=__t11-__it;
		__t10=__t10+__rt;					__t11=__t11+__it;

		__rt =__t1C;						__it =__t1D;
		__t1C=__t14-__it;					__t1D=__t15+__rt;
		__t14=__t14+__it;					__t15=__t15-__rt;

		__rt =(__t1A+__t1B)*ISRT2;			__it =(__t1A-__t1B)*ISRT2;
		__t1A=__t12-__rt;					__t1B=__t13+__it;
		__t12=__t12+__rt;					__t13=__t13-__it;

		__rt =(__t1E-__t1F)*ISRT2;			__it =(__t1F+__t1E)*ISRT2;
		__t1E=__t16+__rt;					__t1F=__t17+__it;
		__t16=__t16-__rt;					__t17=__t17-__it;

	/*...Block 3: */
		__t20=*(__A+__idx[0x10]);			__t21=*(Aim+__idx[0x10]);
		__rt =*(__A+__idx[0x11]);			__it =*(Aim+__idx[0x11]);
		__t22=__t20-__rt;					__t23=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__t24=*(__A+__idx[0x12]);			__t25=*(Aim+__idx[0x12]);
		__rt =*(__A+__idx[0x13]);			__it =*(Aim+__idx[0x13]);
		__t26=__t24-__rt;					__t27=__t25-__it;
		__t24=__t24+__rt;					__t25=__t25+__it;

		__rt =__t24;						__it =__t25;
		__t24=__t20-__rt;					__t25=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t26;						__it =__t27;
		__t26=__t22-__it;					__t27=__t23+__rt;
		__t22=__t22+__it;					__t23=__t23-__rt;

		__t28=*(__A+__idx[0x14]);			__t29=*(Aim+__idx[0x14]);
		__rt =*(__A+__idx[0x15]);			__it =*(Aim+__idx[0x15]);
		__t2A=__t28-__rt;					__t2B=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__t2C=*(__A+__idx[0x16]);			__t2D=*(Aim+__idx[0x16]);
		__rt =*(__A+__idx[0x17]);			__it =*(Aim+__idx[0x17]);
		__t2E=__t2C-__rt;					__t2F=__t2D-__it;
		__t2C=__t2C+__rt;					__t2D=__t2D+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t28-__rt;					__t2D=__t29-__it;
		__t28=__t28+__rt;					__t29=__t29+__it;

		__rt =__t2E;						__it =__t2F;
		__t2E=__t2A-__it;					__t2F=__t2B+__rt;
		__t2A=__t2A+__it;					__t2B=__t2B-__rt;

		__rt =__t28;						__it =__t29;
		__t28=__t20-__rt;					__t29=__t21-__it;
		__t20=__t20+__rt;					__t21=__t21+__it;

		__rt =__t2C;						__it =__t2D;
		__t2C=__t24-__it;					__t2D=__t25+__rt;
		__t24=__t24+__it;					__t25=__t25-__rt;

		__rt =(__t2A+__t2B)*ISRT2;			__it =(__t2A-__t2B)*ISRT2;
		__t2A=__t22-__rt;					__t2B=__t23+__it;
		__t22=__t22+__rt;					__t23=__t23-__it;

		__rt =(__t2E-__t2F)*ISRT2;			__it =(__t2F+__t2E)*ISRT2;
		__t2E=__t26+__rt;					__t2F=__t27+__it;
		__t26=__t26-__rt;					__t27=__t27-__it;

	/*...Block 4: */
		__t30=*(__A+__idx[0x18]);			__t31=*(Aim+__idx[0x18]);
		__rt =*(__A+__idx[0x19]);			__it =*(Aim+__idx[0x19]);
		__t32=__t30-__rt;					__t33=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__t34=*(__A+__idx[0x1A]);			__t35=*(Aim+__idx[0x1A]);
		__rt =*(__A+__idx[0x1B]);			__it =*(Aim+__idx[0x1B]);
		__t36=__t34-__rt;					__t37=__t35-__it;
		__t34=__t34+__rt;					__t35=__t35+__it;

		__rt =__t34;						__it =__t35;
		__t34=__t30-__rt;					__t35=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t36;						__it =__t37;
		__t36=__t32-__it;					__t37=__t33+__rt;
		__t32=__t32+__it;					__t33=__t33-__rt;

		__t38=*(__A+__idx[0x1C]);			__t39=*(Aim+__idx[0x1C]);
		__rt =*(__A+__idx[0x1D]);			__it =*(Aim+__idx[0x1D]);
		__t3A=__t38-__rt;					__t3B=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__t3C=*(__A+__idx[0x1E]);			__t3D=*(Aim+__idx[0x1E]);
		__rt =*(__A+__idx[0x1F]);			__it =*(Aim+__idx[0x1F]);
		__t3E=__t3C-__rt;					__t3F=__t3D-__it;
		__t3C=__t3C+__rt;					__t3D=__t3D+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t38-__rt;					__t3D=__t39-__it;
		__t38=__t38+__rt;					__t39=__t39+__it;

		__rt =__t3E;						__it =__t3F;
		__t3E=__t3A-__it;					__t3F=__t3B+__rt;
		__t3A=__t3A+__it;					__t3B=__t3B-__rt;

		__rt =__t38;						__it =__t39;
		__t38=__t30-__rt;					__t39=__t31-__it;
		__t30=__t30+__rt;					__t31=__t31+__it;

		__rt =__t3C;						__it =__t3D;
		__t3C=__t34-__it;					__t3D=__t35+__rt;
		__t34=__t34+__it;					__t35=__t35-__rt;

		__rt =(__t3A+__t3B)*ISRT2;			__it =(__t3A-__t3B)*ISRT2;
		__t3A=__t32-__rt;					__t3B=__t33+__it;
		__t32=__t32+__rt;					__t33=__t33-__it;

		__rt =(__t3E-__t3F)*ISRT2;			__it =(__t3F+__t3E)*ISRT2;
		__t3E=__t36+__rt;					__t3F=__t37+__it;
		__t36=__t36-__rt;					__t37=__t37-__it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors: */
	/* Totals for the eight radix-4: 168 ADD, 72 MUL: */
	/*...Block 1: __t00,__t10,__t20,__t30	*/
		__rt =__t10;	__t10=__t00-__rt;	__t00=__t00+__rt;
		__it =__t11;	__t11=__t01-__it;	__t01=__t01+__it;

		__rt =__t30;	__t30=__t20-__rt;	__t20=__t20+__rt;
		__it =__t31;	__t31=__t21-__it;	__t21=__t21+__it;
	/* 16 ADD, 0 MUL: */
		*(Bre+0x00)=__t00+__t20;			*(Bim+0x00)=__t01+__t21;
		*(Bre+0x20)=__t00-__t20;			*(Bim+0x20)=__t01-__t21;
		*(Bre+0x10)=__t10+__t31;			*(Bim+0x10)=__t11-__t30;
		*(Bre+0x30)=__t10-__t31;			*(Bim+0x30)=__t11+__t30;

	/*...Block 5: __t08,__t18,__t28,__t38	*/
		__rt =__t18;
		__t18=__t08-__t19;					__t08=__t08+__t19;
		__t19=__t09+__rt;					__t09=__t09-__rt;

		__rt =(__t29+__t28)*ISRT2;			__t29=(__t29-__t28)*ISRT2;	__t28=__rt;
		__rt =(__t38-__t39)*ISRT2;			__it =(__t38+__t39)*ISRT2;
		__t38=__t28+__rt;					__t28=__t28-__rt;
		__t39=__t29+__it;					__t29=__t29-__it;
	/* 20 ADD, 4 MUL: */
		Bre += 0x08;						Bim += 0x08;
		*(Bre+0x00)=__t08+__t28;			*(Bim+0x00)=__t09+__t29;
		*(Bre+0x20)=__t08-__t28;			*(Bim+0x20)=__t09-__t29;
		*(Bre+0x10)=__t18+__t39;			*(Bim+0x10)=__t19-__t38;
		*(Bre+0x30)=__t18-__t39;			*(Bim+0x30)=__t19+__t38;

	/*...Block 3: __t04,__t14,__t24,__t34	*/
		__rt =(__t15+__t14)*ISRT2;			__it =(__t15-__t14)*ISRT2;
		__t14=__t04-__rt;					__t04=__t04+__rt;
		__t15=__t05-__it;					__t05=__t05+__it;

		__rt =__t24*c16 + __t25*s16;	__t25=__t25*c16 - __t24*s16;	__t24=__rt;
		__rt =__t34*s16 + __t35*c16;	__it =__t35*s16 - __t34*c16;
		__t34=__t24-__rt;					__t24=__t24+__rt;
		__t35=__t25-__it;					__t25=__t25+__it;
	/* 22 ADD, 10 MUL: */
		Bre -= 0x04;						Bim -= 0x04;
		*(Bre+0x00)=__t04+__t24;			*(Bim+0x00)=__t05+__t25;
		*(Bre+0x20)=__t04-__t24;			*(Bim+0x20)=__t05-__t25;
		*(Bre+0x10)=__t14+__t35;			*(Bim+0x10)=__t15-__t34;
		*(Bre+0x30)=__t14-__t35;			*(Bim+0x30)=__t15+__t34;

	/*...Block 7: __t0C,__t1C,__t2C,__t3C	*/
		__rt =(__t1C-__t1D)*ISRT2;			__it =(__t1C+__t1D)*ISRT2;
		__t1C=__t0C+__rt;					__t0C=__t0C-__rt;
		__t1D=__t0D+__it;					__t0D=__t0D-__it;

		__rt =__t2C*s16 + __t2D*c16;	__t2D=__t2D*s16 - __t2C*c16;	__t2C=__rt;
		__rt =__t3C*c16 + __t3D*s16;	__it =__t3D*c16 - __t3C*s16;
		__t3C=__t2C+__rt;					__t2C=__t2C-__rt;
		__t3D=__t2D+__it;					__t2D=__t2D-__it;
	/* 22 ADD, 10 MUL: */
		Bre += 0x08;						Bim += 0x08;
		*(Bre+0x00)=__t0C+__t2C;			*(Bim+0x00)=__t0D+__t2D;
		*(Bre+0x20)=__t0C-__t2C;			*(Bim+0x20)=__t0D-__t2D;
		*(Bre+0x10)=__t1C+__t3D;			*(Bim+0x10)=__t1D-__t3C;
		*(Bre+0x30)=__t1C-__t3D;			*(Bim+0x30)=__t1D+__t3C;

	/*...Block 2: __t02,__t12,__t22,__t32	*/
		__rt =__t12*c16 + __t13*s16;	__it =__t13*c16 - __t12*s16;
		__t12=__t02-__rt;					__t02=__t02+__rt;
		__t13=__t03-__it;					__t03=__t03+__it;

		__rt =__t22*c32_1 + __t23*s32_1;	__t23=__t23*c32_1 - __t22*s32_1;	__t22=__rt;
		__rt =__t32*c32_3 + __t33*s32_3;	__it =__t33*c32_3 - __t32*s32_3;
		__t32=__t22-__rt;					__t22=__t22+__rt;
		__t33=__t23-__it;					__t23=__t23+__it;
	/* 22 ADD, 12 MUL: */
		Bre -= 0x0a;						Bim -= 0x0a;
		*(Bre+0x00)=__t02+__t22;			*(Bim+0x00)=__t03+__t23;
		*(Bre+0x20)=__t02-__t22;			*(Bim+0x20)=__t03-__t23;
		*(Bre+0x10)=__t12+__t33;			*(Bim+0x10)=__t13-__t32;
		*(Bre+0x30)=__t12-__t33;			*(Bim+0x30)=__t13+__t32;

	/*...Block 6: __t0A,__t1A,__t2A,__t3A	*/
		__rt =__t1A*s16 - __t1B*c16;	__it =__t1B*s16 + __t1A*c16;
		__t1A=__t0A+__rt;					__t0A=__t0A-__rt;
		__t1B=__t0B+__it;					__t0B=__t0B-__it;

		__rt =__t2A*s32_3 + __t2B*c32_3;	__t2B=__t2B*s32_3 - __t2A*c32_3;	__t2A=__rt;
		__rt =__t3A*c32_1 - __t3B*s32_1;	__it =__t3B*c32_1 + __t3A*s32_1;
		__t3A=__t2A+__rt;					__t2A=__t2A-__rt;
		__t3B=__t2B+__it;					__t2B=__t2B-__it;
	/* 22 ADD, 12 MUL: */
		Bre += 0x08;						Bim += 0x08;
		*(Bre+0x00)=__t0A+__t2A;			*(Bim+0x00)=__t0B+__t2B;
		*(Bre+0x20)=__t0A-__t2A;			*(Bim+0x20)=__t0B-__t2B;
		*(Bre+0x10)=__t1A+__t3B;			*(Bim+0x10)=__t1B-__t3A;
		*(Bre+0x30)=__t1A-__t3B;			*(Bim+0x30)=__t1B+__t3A;

	/*...Block 4: __t06,__t16,__t26,__t36	*/
		__rt =__t16*s16 + __t17*c16;	__it =__t17*s16 - __t16*c16;
		__t16=__t06-__rt;					__t06=__t06+__rt;
		__t17=__t07-__it;					__t07=__t07+__it;

		__rt =__t26*c32_3 + __t27*s32_3;	__t27=__t27*c32_3 - __t26*s32_3;	__t26=__rt;
		__rt =__t36*s32_1 - __t37*c32_1;	__it =__t37*s32_1 + __t36*c32_1;
		__t36=__t26+__rt;					__t26=__t26-__rt;
		__t37=__t27+__it;					__t27=__t27-__it;
	/* 22 ADD, 12 MUL: */
		Bre -= 0x04;						Bim -= 0x04;
		*(Bre+0x00)=__t06+__t26;			*(Bim+0x00)=__t07+__t27;
		*(Bre+0x20)=__t06-__t26;			*(Bim+0x20)=__t07-__t27;
		*(Bre+0x10)=__t16+__t37;			*(Bim+0x10)=__t17-__t36;
		*(Bre+0x30)=__t16-__t37;			*(Bim+0x30)=__t17+__t36;

	/*...Block 8: __t0E,__t1E,__t2E,__t3E	*/
		__rt =__t1E*c16 - __t1F*s16;	__it =__t1F*c16 + __t1E*s16;
		__t1E=__t0E+__rt;					__t0E=__t0E-__rt;
		__t1F=__t0F+__it;					__t0F=__t0F-__it;

		__rt =__t2E*s32_1 + __t2F*c32_1;	__t2F=__t2F*s32_1 - __t2E*c32_1;	__t2E=__rt;
		__rt =__t3E*s32_3 + __t3F*c32_3;	__it =__t3F*s32_3 - __t3E*c32_3;
		__t3E=__t2E+__rt;					__t2E=__t2E-__rt;
		__t3F=__t2F+__it;					__t2F=__t2F-__it;
	/* 22 ADD, 12 MUL: */
		Bre += 0x08;						Bim += 0x08;
		*(Bre+0x00)=__t0E+__t2E;			*(Bim+0x00)=__t0F+__t2F;
		*(Bre+0x20)=__t0E-__t2E;			*(Bim+0x20)=__t0F-__t2F;
		*(Bre+0x10)=__t1E+__t3F;			*(Bim+0x10)=__t1F-__t3E;
		*(Bre+0x30)=__t1E-__t3F;			*(Bim+0x30)=__t1F+__t3E;
}

/************** RADIX-63 DIF/DIT: *****************************/

void RADIX_63_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 63 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 63 (index) offsets */
)
{
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;
	int l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
	const uint8 *iptr,
		dif_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x36,0x2d,0x24,0x1b,0x12,0x09,
			0x38,0x2f,0x26,0x1d,0x14,0x0b,0x02,
			0x31,0x28,0x1f,0x16,0x0d,0x04,0x3a,
			0x2a,0x21,0x18,0x0f,0x06,0x3c,0x33,
			0x23,0x1a,0x11,0x08,0x3e,0x35,0x2c,
			0x1c,0x13,0x0a,0x01,0x37,0x2e,0x25,
			0x15,0x0c,0x03,0x39,0x30,0x27,0x1e,
			0x0e,0x05,0x3b,0x32,0x29,0x20,0x17,
			0x07,0x3d,0x34,0x2b,0x22,0x19,0x10},
		dif_operm[64] = {	// ditto
			0x00,0x07,0x03,0x02,0x06,0x05,0x01,0x08,0x04,
			0x37,0x3e,0x3a,0x36,0x3d,0x39,0x38,0x3c,0x3b,
			0x32,0x2e,0x35,0x31,0x2d,0x34,0x30,0x2f,0x33,
			0x2a,0x29,0x25,0x2c,0x28,0x24,0x2b,0x27,0x26,
			0x1d,0x21,0x20,0x1c,0x23,0x1f,0x1b,0x22,0x1e,
			0x15,0x14,0x18,0x17,0x13,0x1a,0x16,0x12,0x19,
			0x10,0x0c,0x0b,0x0f,0x0e,0x0a,0x11,0x0d,0x09};
	const double	uc1 = .62348980185873353053,	/* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	/* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	/* cos(2u)	*/
					us2 = .97492791218182360702,	/* sin(2u)	*/
					uc3 =-.90096886790241912622,	/* cos(3u)	*/
					us3 = .43388373911755812050;	/* sin(3u)	*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[63], *tptr;
	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 9 radix-7 transforms:
	/*
	Twiddleless version arranges 9 sets of radix-7 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 63) 9 horizontally and 7 vertically. Display result of DIF/DIT input-scramble array in hex:

		00,36,2d,24,1b,12,09
		38,2f,26,1d,14,0b,02
		31,28,1f,16,0d,04,3a
		2a,21,18,0f,06,3c,33
		23,1a,11,08,3e,35,2c
		1c,13,0a,01,37,2e,25
		15,0c,03,39,30,27,1e
		0e,05,3b,32,29,20,17
		07,3d,34,2b,22,19,10
	*/
	tptr = t; iptr = dif_iperm;
	for(l = 0; l < 9; l++) {
		k0 = __idx[*iptr]; k1 = __idx[*(iptr+1)]; k2 = __idx[*(iptr+2)]; k3 = __idx[*(iptr+3)]; k4 = __idx[*(iptr+4)]; k5 = __idx[*(iptr+5)]; k6 = __idx[*(iptr+6)];
		RADIX_07_DFT(
			*(__A+k0),*(Aim+k0),*(__A+k1),*(Aim+k1),*(__A+k2),*(Aim+k2),*(__A+k3),*(Aim+k3),*(__A+k4),*(Aim+k4),*(__A+k5),*(Aim+k5),*(__A+k6),*(Aim+k6),
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
			tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
			uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
		);	tptr++; iptr += 7;
	}
	/*...and now do 7 radix-9 transforms. The required output permutation is

		00,07,03,02,06,05,01,08,04,
		37,3e,3a,36,3d,39,38,3c,3b,
		32,2e,35,31,2d,34,30,2f,33,
		2a,29,25,2c,28,24,2b,27,26,
		1d,21,20,1c,23,1f,1b,22,1e,
		15,14,18,17,13,1a,16,12,19,
		10,0c,0b,0f,0e,0a,11,0d,09.
	*/
	tptr = t; iptr = dif_operm;
	for(l = 0; l < 7; l++) {
		// When 63 is used to build a larger DFT radix (e.g. 1008), these indices will be permuted (nonmonotone), no simplification possible:
		k0 = __odx[*iptr]; k1 = __odx[*(iptr+1)]; k2 = __odx[*(iptr+2)]; k3 = __odx[*(iptr+3)]; k4 = __odx[*(iptr+4)]; k5 = __odx[*(iptr+5)]; k6 = __odx[*(iptr+6)]; k7 = __odx[*(iptr+7)]; k8 = __odx[*(iptr+8)];
		RADIX_09_DIF(
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
			*(__B+k0),*(Bim+k0),*(__B+k1),*(Bim+k1),*(__B+k2),*(Bim+k2),*(__B+k3),*(Bim+k3),*(__B+k4),*(Bim+k4),*(__B+k5),*(Bim+k5),*(__B+k6),*(Bim+k6),*(__B+k7),*(Bim+k7),*(__B+k8),*(Bim+k8),
			rt,it,re
		);	tptr += 9; iptr += 9;
	}
}

/***************/

void RADIX_63_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 63 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 63 (index) offsets */
)
{
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;
	int l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
	const uint8 *iptr,
		dit_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,
			0x32,0x31,0x30,0x2f,0x2e,0x2d,0x34,0x33,0x35,
			0x1d,0x1c,0x1b,0x22,0x21,0x23,0x1f,0x1e,0x20,
			0x10,0x0f,0x11,0x0d,0x0c,0x0e,0x0a,0x09,0x0b,
			0x37,0x36,0x38,0x3c,0x3e,0x3d,0x39,0x3b,0x3a,
			0x2a,0x2c,0x2b,0x27,0x29,0x28,0x24,0x26,0x25,
			0x15,0x17,0x16,0x12,0x14,0x13,0x1a,0x19,0x18,0},
		dit_operm[64] = {	// ditto
			0x00,0x24,0x09,0x2d,0x12,0x36,0x1b,
			0x0e,0x32,0x17,0x3b,0x20,0x05,0x29,
			0x1c,0x01,0x25,0x0a,0x2e,0x13,0x37,
			0x2a,0x0f,0x33,0x18,0x3c,0x21,0x06,
			0x38,0x1d,0x02,0x26,0x0b,0x2f,0x14,
			0x07,0x2b,0x10,0x34,0x19,0x3d,0x22,
			0x15,0x39,0x1e,0x03,0x27,0x0c,0x30,
			0x23,0x08,0x2c,0x11,0x35,0x1a,0x3e,
			0x31,0x16,0x3a,0x1f,0x04,0x28,0x0d,0};
	const double	uc1 = .62348980185873353053,	/* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	/* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	/* cos(2u)	*/
					us2 = .97492791218182360702,	/* sin(2u)	*/
					uc3 =-.90096886790241912622,	/* cos(3u)	*/
					us3 = .43388373911755812050;	/* sin(3u)	*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[63], *tptr;

	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 7 radix-9 transforms:
	/*
	Twiddleless version arranges 7 sets of radix-9 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 63) 9 horizontally and 7 vertically. Applying a further bit-reversal to that,
	Display result of Combined DIT input-scramble array in hex:

		00,02,01,08,07,06,05,04,03,
		32,31,30,2f,2e,2d,34,33,35,
		1d,1c,1b,22,21,23,1f,1e,20,
		10,0f,11,0d,0c,0e,0a,09,0b,
		37,36,38,3c,3e,3d,39,3b,3a,
		2a,2c,2b,27,29,28,24,26,25,
		15,17,16,12,14,13,1a,19,18.
	*/
	tptr = t; iptr = dit_iperm;
	for(l = 0; l < 7; l++) {
		// When 63 is used to build a larger DFT radix (e.g. 1008), these indices will be permuted (nonmonotone), no simplification possible:
		k0 = __idx[*iptr]; k1 = __idx[*(iptr+1)]; k2 = __idx[*(iptr+2)]; k3 = __idx[*(iptr+3)]; k4 = __idx[*(iptr+4)]; k5 = __idx[*(iptr+5)]; k6 = __idx[*(iptr+6)]; k7 = __idx[*(iptr+7)]; k8 = __idx[*(iptr+8)];
		RADIX_09_DIT(
			*(__A+k0),*(Aim+k0),*(__A+k1),*(Aim+k1),*(__A+k2),*(Aim+k2),*(__A+k3),*(Aim+k3),*(__A+k4),*(Aim+k4),*(__A+k5),*(Aim+k5),*(__A+k6),*(Aim+k6),*(__A+k7),*(Aim+k7),*(__A+k8),*(Aim+k8),
			tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
			rt,it,re
		);
		tptr += 9; iptr += 9;
	}
	/*...and now do 9 radix-7 transforms. The required output permutation is

		00,24,09,2d,12,36,1b,
		0e,32,17,3b,20,05,29,
		1c,01,25,0a,2e,13,37,
		2a,0f,33,18,3c,21,06,
		38,1d,02,26,0b,2f,14,
		07,2b,10,34,19,3d,22,
		15,39,1e,03,27,0c,30,
		23,08,2c,11,35,1a,3e,
		31,16,3a,1f,04,28,0d
	*/
	tptr = t; iptr = dit_operm;
	for(l = 0; l < 9; l++) {
		k0 = __odx[*iptr]; k1 = __odx[*(iptr+1)]; k2 = __odx[*(iptr+2)]; k3 = __odx[*(iptr+3)]; k4 = __odx[*(iptr+4)]; k5 = __odx[*(iptr+5)]; k6 = __odx[*(iptr+6)];
		RADIX_07_DFT(
			tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
			*(__B+k0),*(Bim+k0),*(__B+k1),*(Bim+k1),*(__B+k2),*(Bim+k2),*(__B+k3),*(Bim+k3),*(__B+k4),*(Bim+k4),*(__B+k5),*(Bim+k5),*(__B+k6),*(Bim+k6),
			uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
		);	tptr++; iptr += 7;
	}
}

// For power-of-2 radix > 32 we use a common macro for both the first-pass-of-2-pass-pow2
// and twiddleless-subtransform-to-be-combined-with-odd-radix cases:

/************** RADIX-64 DIF/DIT: *****************************/

void RADIX_64_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 64 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 64 (index) offsets */
)
{
	double __t[128];
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;

/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */

	//...Block 0: jt = j1;	jp = j2;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x00]),*(Aim+__idx[0x00]),*(__A+__idx[0x08]),*(Aim+__idx[0x08]),*(__A+__idx[0x10]),*(Aim+__idx[0x10]),*(__A+__idx[0x18]),*(Aim+__idx[0x18]),*(__A+__idx[0x20]),*(Aim+__idx[0x20]),*(__A+__idx[0x28]),*(Aim+__idx[0x28]),*(__A+__idx[0x30]),*(Aim+__idx[0x30]),*(__A+__idx[0x38]),*(Aim+__idx[0x38])
		,__t[0x00],__t[0x01],__t[0x02],__t[0x03],__t[0x04],__t[0x05],__t[0x06],__t[0x07],__t[0x08],__t[0x09],__t[0x0A],__t[0x0B],__t[0x0C],__t[0x0D],__t[0x0E],__t[0x0F]
	);
	//...Block 1: jt = j1 + p04;	jp = j2 + p04;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x04]),*(Aim+__idx[0x04]),*(__A+__idx[0x0c]),*(Aim+__idx[0x0c]),*(__A+__idx[0x14]),*(Aim+__idx[0x14]),*(__A+__idx[0x1c]),*(Aim+__idx[0x1c]),*(__A+__idx[0x24]),*(Aim+__idx[0x24]),*(__A+__idx[0x2c]),*(Aim+__idx[0x2c]),*(__A+__idx[0x34]),*(Aim+__idx[0x34]),*(__A+__idx[0x3c]),*(Aim+__idx[0x3c])
		,__t[0x10],__t[0x11],__t[0x12],__t[0x13],__t[0x14],__t[0x15],__t[0x16],__t[0x17],__t[0x18],__t[0x19],__t[0x1A],__t[0x1B],__t[0x1C],__t[0x1D],__t[0x1E],__t[0x1F]
	);
	//...Block 2: jt = j1 + p02;	jp = j2 + p02;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x02]),*(Aim+__idx[0x02]),*(__A+__idx[0x0a]),*(Aim+__idx[0x0a]),*(__A+__idx[0x12]),*(Aim+__idx[0x12]),*(__A+__idx[0x1a]),*(Aim+__idx[0x1a]),*(__A+__idx[0x22]),*(Aim+__idx[0x22]),*(__A+__idx[0x2a]),*(Aim+__idx[0x2a]),*(__A+__idx[0x32]),*(Aim+__idx[0x32]),*(__A+__idx[0x3a]),*(Aim+__idx[0x3a])
		,__t[0x20],__t[0x21],__t[0x22],__t[0x23],__t[0x24],__t[0x25],__t[0x26],__t[0x27],__t[0x28],__t[0x29],__t[0x2A],__t[0x2B],__t[0x2C],__t[0x2D],__t[0x2E],__t[0x2F]
	);
	//...Block 3: jt = j1 + p06;	jp = j2 + p06;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x06]),*(Aim+__idx[0x06]),*(__A+__idx[0x0e]),*(Aim+__idx[0x0e]),*(__A+__idx[0x16]),*(Aim+__idx[0x16]),*(__A+__idx[0x1e]),*(Aim+__idx[0x1e]),*(__A+__idx[0x26]),*(Aim+__idx[0x26]),*(__A+__idx[0x2e]),*(Aim+__idx[0x2e]),*(__A+__idx[0x36]),*(Aim+__idx[0x36]),*(__A+__idx[0x3e]),*(Aim+__idx[0x3e])
		,__t[0x30],__t[0x31],__t[0x32],__t[0x33],__t[0x34],__t[0x35],__t[0x36],__t[0x37],__t[0x38],__t[0x39],__t[0x3A],__t[0x3B],__t[0x3C],__t[0x3D],__t[0x3E],__t[0x3F]
	);
	//...Block 4: jt = j1 + p01;	jp = j2 + p01;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x01]),*(Aim+__idx[0x01]),*(__A+__idx[0x09]),*(Aim+__idx[0x09]),*(__A+__idx[0x11]),*(Aim+__idx[0x11]),*(__A+__idx[0x19]),*(Aim+__idx[0x19]),*(__A+__idx[0x21]),*(Aim+__idx[0x21]),*(__A+__idx[0x29]),*(Aim+__idx[0x29]),*(__A+__idx[0x31]),*(Aim+__idx[0x31]),*(__A+__idx[0x39]),*(Aim+__idx[0x39])
		,__t[0x40],__t[0x41],__t[0x42],__t[0x43],__t[0x44],__t[0x45],__t[0x46],__t[0x47],__t[0x48],__t[0x49],__t[0x4A],__t[0x4B],__t[0x4C],__t[0x4D],__t[0x4E],__t[0x4F]
	);
	//...Block 5: jt = j1 + p05;	jp = j2 + p05;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x05]),*(Aim+__idx[0x05]),*(__A+__idx[0x0d]),*(Aim+__idx[0x0d]),*(__A+__idx[0x15]),*(Aim+__idx[0x15]),*(__A+__idx[0x1d]),*(Aim+__idx[0x1d]),*(__A+__idx[0x25]),*(Aim+__idx[0x25]),*(__A+__idx[0x2d]),*(Aim+__idx[0x2d]),*(__A+__idx[0x35]),*(Aim+__idx[0x35]),*(__A+__idx[0x3d]),*(Aim+__idx[0x3d])
		,__t[0x50],__t[0x51],__t[0x52],__t[0x53],__t[0x54],__t[0x55],__t[0x56],__t[0x57],__t[0x58],__t[0x59],__t[0x5A],__t[0x5B],__t[0x5C],__t[0x5D],__t[0x5E],__t[0x5F]
	);
	//...Block 6: jt = j1 + p03;	jp = j2 + p03;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x03]),*(Aim+__idx[0x03]),*(__A+__idx[0x0b]),*(Aim+__idx[0x0b]),*(__A+__idx[0x13]),*(Aim+__idx[0x13]),*(__A+__idx[0x1b]),*(Aim+__idx[0x1b]),*(__A+__idx[0x23]),*(Aim+__idx[0x23]),*(__A+__idx[0x2b]),*(Aim+__idx[0x2b]),*(__A+__idx[0x33]),*(Aim+__idx[0x33]),*(__A+__idx[0x3b]),*(Aim+__idx[0x3b])
		,__t[0x60],__t[0x61],__t[0x62],__t[0x63],__t[0x64],__t[0x65],__t[0x66],__t[0x67],__t[0x68],__t[0x69],__t[0x6A],__t[0x6B],__t[0x6C],__t[0x6D],__t[0x6E],__t[0x6F]
	);
	//...Block 7: jt = j1 + p07;	jp = j2 + p07;
	RADIX_08_DIF_OOP(
		*(__A+__idx[0x07]),*(Aim+__idx[0x07]),*(__A+__idx[0x0f]),*(Aim+__idx[0x0f]),*(__A+__idx[0x17]),*(Aim+__idx[0x17]),*(__A+__idx[0x1f]),*(Aim+__idx[0x1f]),*(__A+__idx[0x27]),*(Aim+__idx[0x27]),*(__A+__idx[0x2f]),*(Aim+__idx[0x2f]),*(__A+__idx[0x37]),*(Aim+__idx[0x37]),*(__A+__idx[0x3f]),*(Aim+__idx[0x3f])
		,__t[0x70],__t[0x71],__t[0x72],__t[0x73],__t[0x74],__t[0x75],__t[0x76],__t[0x77],__t[0x78],__t[0x79],__t[0x7A],__t[0x7B],__t[0x7C],__t[0x7D],__t[0x7E],__t[0x7F]
	);

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

	/* Block 0: */
	// jt = j1;	jp = j2;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
	so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
	RADIX_08_DIF_OOP(
		__t[0x00],__t[0x01],__t[0x40],__t[0x41],__t[0x20],__t[0x21],__t[0x60],__t[0x61],__t[0x10],__t[0x11],__t[0x50],__t[0x51],__t[0x30],__t[0x31],__t[0x70],__t[0x71],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x07]),*(Bim+__odx[0x07])
	);
	/* Block 4: */
	__odx += 8;	// jt = j1 + p08;	jp = j2 + p08;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x08],__t[0x09],__t[0x18],__t[0x19],__t[0x28],__t[0x29],__t[0x38],__t[0x39],__t[0x48],__t[0x49],__t[0x58],__t[0x59],__t[0x68],__t[0x69],__t[0x78],__t[0x79],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16
	);
	/* Block 2: */
	__odx += 8;	// jt = j1 + p10;	jp = j2 + p10;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x04],__t[0x05],__t[0x14],__t[0x15],__t[0x24],__t[0x25],__t[0x34],__t[0x35],__t[0x44],__t[0x45],__t[0x54],__t[0x55],__t[0x64],__t[0x65],__t[0x74],__t[0x75],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
	);
	/* Block 6: */
	__odx += 8;	// jt = j1 + p18;	jp = j2 + p18;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x0C],__t[0x0D],__t[0x1C],__t[0x1D],__t[0x2C],__t[0x2D],__t[0x3C],__t[0x3D],__t[0x4C],__t[0x4D],__t[0x5C],__t[0x5D],__t[0x6C],__t[0x6D],__t[0x7C],__t[0x7D],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
	);
	/* Block 1: */
	__odx += 8;	// jt = j1 + p20;	jp = j2 + p20;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x02],__t[0x03],__t[0x12],__t[0x13],__t[0x22],__t[0x23],__t[0x32],__t[0x33],__t[0x42],__t[0x43],__t[0x52],__t[0x53],__t[0x62],__t[0x63],__t[0x72],__t[0x73],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
	);
	/* Block 5: */
	__odx += 8;	// jt = j1 + p28;	jp = j2 + p28;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x0A],__t[0x0B],__t[0x1A],__t[0x1B],__t[0x2A],__t[0x2B],__t[0x3A],__t[0x3B],__t[0x4A],__t[0x4B],__t[0x5A],__t[0x5B],__t[0x6A],__t[0x6B],__t[0x7A],__t[0x7B],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
	);
	/* Block 3: */
	__odx += 8;	// jt = j1 + p30;	jp = j2 + p30;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x06],__t[0x07],__t[0x16],__t[0x17],__t[0x26],__t[0x27],__t[0x36],__t[0x37],__t[0x46],__t[0x47],__t[0x56],__t[0x57],__t[0x66],__t[0x67],__t[0x76],__t[0x77],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
	);
	/* Block 7: */
	__odx += 8;	// jt = j1 + p38;	jp = j2 + p38;
	RADIX_08_DIF_TWIDDLE_OOP(
		__t[0x0E],__t[0x0F],__t[0x1E],__t[0x1F],__t[0x2E],__t[0x2F],__t[0x3E],__t[0x3F],__t[0x4E],__t[0x4F],__t[0x5E],__t[0x5F],__t[0x6E],__t[0x6F],__t[0x7E],__t[0x7F],
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x07]),*(Bim+__odx[0x07]),
		-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
	);
}

void RADIX_64_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 64 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 64 (index) offsets */
)
{
	double
		 __t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F
		,__t40,__t41,__t42,__t43,__t44,__t45,__t46,__t47,__t48,__t49,__t4A,__t4B,__t4C,__t4D,__t4E,__t4F
		,__t50,__t51,__t52,__t53,__t54,__t55,__t56,__t57,__t58,__t59,__t5A,__t5B,__t5C,__t5D,__t5E,__t5F
		,__t60,__t61,__t62,__t63,__t64,__t65,__t66,__t67,__t68,__t69,__t6A,__t6B,__t6C,__t6D,__t6E,__t6F
		,__t70,__t71,__t72,__t73,__t74,__t75,__t76,__t77,__t78,__t79,__t7A,__t7B,__t7C,__t7D,__t7E,__t7F;
	double *Aim = __A + __re_im_stride_in, *Bim = __B + __re_im_stride_out;
/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
	/*...Block 0: */
	// jt = j1;	jp = j2;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t00,__t01,__t02,__t03,__t04,__t05,__t06,__t07,__t08,__t09,__t0A,__t0B,__t0C,__t0D,__t0E,__t0F
	);
	/*...Block 1: */
	__idx += 8;	// jt = j1 + p08;	jp = j2 + p08;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t10,__t11,__t12,__t13,__t14,__t15,__t16,__t17,__t18,__t19,__t1A,__t1B,__t1C,__t1D,__t1E,__t1F
	);
	/*...Block 2: */
	__idx += 8;	// jt = j1 + p10;	jp = j2 + p10;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t20,__t21,__t22,__t23,__t24,__t25,__t26,__t27,__t28,__t29,__t2A,__t2B,__t2C,__t2D,__t2E,__t2F
	);
	/*...Block 3: */
	__idx += 8;	// jt = j1 + p18;	jp = j2 + p18;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t30,__t31,__t32,__t33,__t34,__t35,__t36,__t37,__t38,__t39,__t3A,__t3B,__t3C,__t3D,__t3E,__t3F
	);
	/*...Block 4: */
	__idx += 8;	// jt = j1 + p20;	jp = j2 + p20;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t40,__t41,__t42,__t43,__t44,__t45,__t46,__t47,__t48,__t49,__t4A,__t4B,__t4C,__t4D,__t4E,__t4F
	);
	/*...Block 5: */
	__idx += 8;	// jt = j1 + p28;	jp = j2 + p28;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t50,__t51,__t52,__t53,__t54,__t55,__t56,__t57,__t58,__t59,__t5A,__t5B,__t5C,__t5D,__t5E,__t5F
	);
	/*...Block 6: */
	__idx += 8;	// jt = j1 + p30;	jp = j2 + p30;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t60,__t61,__t62,__t63,__t64,__t65,__t66,__t67,__t68,__t69,__t6A,__t6B,__t6C,__t6D,__t6E,__t6F
	);
	/*...Block 7: */
	__idx += 8;	// jt = j1 + p38;	jp = j2 + p38;
	RADIX_08_DIT_OOP(
		*(__A+__idx[0x0]),*(Aim+__idx[0x0]),*(__A+__idx[0x1]),*(Aim+__idx[0x1]),*(__A+__idx[0x2]),*(Aim+__idx[0x2]),*(__A+__idx[0x3]),*(Aim+__idx[0x3]),*(__A+__idx[0x4]),*(Aim+__idx[0x4]),*(__A+__idx[0x5]),*(Aim+__idx[0x5]),*(__A+__idx[0x6]),*(Aim+__idx[0x6]),*(__A+__idx[0x7]),*(Aim+__idx[0x7])
		,__t70,__t71,__t72,__t73,__t74,__t75,__t76,__t77,__t78,__t79,__t7A,__t7B,__t7C,__t7D,__t7E,__t7F
	);

/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

	// Block 0: jt = j1;	jp = j2;
	/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
	RADIX_08_DIT_OOP(
		__t00,__t01,__t10,__t11,__t20,__t21,__t30,__t31,__t40,__t41,__t50,__t51,__t60,__t61,__t70,__t71,
		*(__B+__odx[0x00]),*(Bim+__odx[0x00]),*(__B+__odx[0x08]),*(Bim+__odx[0x08]),*(__B+__odx[0x10]),*(Bim+__odx[0x10]),*(__B+__odx[0x18]),*(Bim+__odx[0x18]),*(__B+__odx[0x20]),*(Bim+__odx[0x20]),*(__B+__odx[0x28]),*(Bim+__odx[0x28]),*(__B+__odx[0x30]),*(Bim+__odx[0x30]),*(__B+__odx[0x38]),*(Bim+__odx[0x38])
	);
	// Block 4: jt = j1 + p04;	jp = j2 + p04;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t08,__t09,__t18,__t19,__t28,__t29,__t38,__t39,__t48,__t49,__t58,__t59,__t68,__t69,__t78,__t79,
		*(__B+__odx[0x04]),*(Bim+__odx[0x04]),*(__B+__odx[0x24]),*(Bim+__odx[0x24]),*(__B+__odx[0x14]),*(Bim+__odx[0x14]),*(__B+__odx[0x34]),*(Bim+__odx[0x34]),*(__B+__odx[0x0c]),*(Bim+__odx[0x0c]),*(__B+__odx[0x2c]),*(Bim+__odx[0x2c]),*(__B+__odx[0x1c]),*(Bim+__odx[0x1c]),*(__B+__odx[0x3c]),*(Bim+__odx[0x3c]),
		0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16
	);
	// Block 2: jt = j1 + p02;	jp = j2 + p02;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t04,__t05,__t14,__t15,__t24,__t25,__t34,__t35,__t44,__t45,__t54,__t55,__t64,__t65,__t74,__t75,
		*(__B+__odx[0x02]),*(Bim+__odx[0x02]),*(__B+__odx[0x22]),*(Bim+__odx[0x22]),*(__B+__odx[0x12]),*(Bim+__odx[0x12]),*(__B+__odx[0x32]),*(Bim+__odx[0x32]),*(__B+__odx[0x0a]),*(Bim+__odx[0x0a]),*(__B+__odx[0x2a]),*(Bim+__odx[0x2a]),*(__B+__odx[0x1a]),*(Bim+__odx[0x1a]),*(__B+__odx[0x3a]),*(Bim+__odx[0x3a]),
		ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
	);
	// Block 6: jt = j1 + p06;	jp = j2 + p06;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t0C,__t0D,__t1C,__t1D,__t2C,__t2D,__t3C,__t3D,__t4C,__t4D,__t5C,__t5D,__t6C,__t6D,__t7C,__t7D,
		*(__B+__odx[0x06]),*(Bim+__odx[0x06]),*(__B+__odx[0x26]),*(Bim+__odx[0x26]),*(__B+__odx[0x16]),*(Bim+__odx[0x16]),*(__B+__odx[0x36]),*(Bim+__odx[0x36]),*(__B+__odx[0x0e]),*(Bim+__odx[0x0e]),*(__B+__odx[0x2e]),*(Bim+__odx[0x2e]),*(__B+__odx[0x1e]),*(Bim+__odx[0x1e]),*(__B+__odx[0x3e]),*(Bim+__odx[0x3e]),
		-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
	);
	// Block 1: jt = j1 + p01;	jp = j2 + p01;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t02,__t03,__t12,__t13,__t22,__t23,__t32,__t33,__t42,__t43,__t52,__t53,__t62,__t63,__t72,__t73,
		*(__B+__odx[0x01]),*(Bim+__odx[0x01]),*(__B+__odx[0x21]),*(Bim+__odx[0x21]),*(__B+__odx[0x11]),*(Bim+__odx[0x11]),*(__B+__odx[0x31]),*(Bim+__odx[0x31]),*(__B+__odx[0x09]),*(Bim+__odx[0x09]),*(__B+__odx[0x29]),*(Bim+__odx[0x29]),*(__B+__odx[0x19]),*(Bim+__odx[0x19]),*(__B+__odx[0x39]),*(Bim+__odx[0x39]),
		c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
	);
	// Block 5: jt = j1 + p05;	jp = j2 + p05;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t0A,__t0B,__t1A,__t1B,__t2A,__t2B,__t3A,__t3B,__t4A,__t4B,__t5A,__t5B,__t6A,__t6B,__t7A,__t7B,
		*(__B+__odx[0x05]),*(Bim+__odx[0x05]),*(__B+__odx[0x25]),*(Bim+__odx[0x25]),*(__B+__odx[0x15]),*(Bim+__odx[0x15]),*(__B+__odx[0x35]),*(Bim+__odx[0x35]),*(__B+__odx[0x0d]),*(Bim+__odx[0x0d]),*(__B+__odx[0x2d]),*(Bim+__odx[0x2d]),*(__B+__odx[0x1d]),*(Bim+__odx[0x1d]),*(__B+__odx[0x3d]),*(Bim+__odx[0x3d]),
		-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
	);
	// Block 3: jt = j1 + p03;	jp = j2 + p03;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t06,__t07,__t16,__t17,__t26,__t27,__t36,__t37,__t46,__t47,__t56,__t57,__t66,__t67,__t76,__t77,
		*(__B+__odx[0x03]),*(Bim+__odx[0x03]),*(__B+__odx[0x23]),*(Bim+__odx[0x23]),*(__B+__odx[0x13]),*(Bim+__odx[0x13]),*(__B+__odx[0x33]),*(Bim+__odx[0x33]),*(__B+__odx[0x0b]),*(Bim+__odx[0x0b]),*(__B+__odx[0x2b]),*(Bim+__odx[0x2b]),*(__B+__odx[0x1b]),*(Bim+__odx[0x1b]),*(__B+__odx[0x3b]),*(Bim+__odx[0x3b]),
		s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
	);
	// Block 7: jt = j1 + p07;	jp = j2 + p07;
	RADIX_08_DIT_TWIDDLE_OOP(
		__t0E,__t0F,__t1E,__t1F,__t2E,__t2F,__t3E,__t3F,__t4E,__t4F,__t5E,__t5F,__t6E,__t6F,__t7E,__t7F,
		*(__B+__odx[0x07]),*(Bim+__odx[0x07]),*(__B+__odx[0x27]),*(Bim+__odx[0x27]),*(__B+__odx[0x17]),*(Bim+__odx[0x17]),*(__B+__odx[0x37]),*(Bim+__odx[0x37]),*(__B+__odx[0x0f]),*(Bim+__odx[0x0f]),*(__B+__odx[0x2f]),*(Bim+__odx[0x2f]),*(__B+__odx[0x1f]),*(Bim+__odx[0x1f]),*(__B+__odx[0x3f]),*(Bim+__odx[0x3f]),
		-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
	);
}


/************** RADIX-128 DIF/DIT: *****************************/

void RADIX_128_DIF(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 128 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 128 (index) offsets */
)
{
	int i,j;
	struct complex t[128], *tptr;
	const struct complex *cd_ptr;
	double *Are = __A,*Aim = __A + __re_im_stride_in, *Bre = __B,*Bim = __B + __re_im_stride_out;
	int p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf;

// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:
	tptr = t;
	for(i = 0; i < 8; i++) {
		j = reverse8[i];	// A-array offsets processed in BR8 order = p[04261537]
		p0 = __idx[j+0x00]; p1 = __idx[j+0x08]; p2 = __idx[j+0x10]; p3 = __idx[j+0x18]; p4 = __idx[j+0x20]; p5 = __idx[j+0x28]; p6 = __idx[j+0x30]; p7 = __idx[j+0x38]; p8 = __idx[j+0x40]; p9 = __idx[j+0x48]; pa = __idx[j+0x50]; pb = __idx[j+0x58]; pc = __idx[j+0x60]; pd = __idx[j+0x68]; pe = __idx[j+0x70]; pf = __idx[j+0x78];
		RADIX_16_DIF(
			*(Are+p0),*(Aim+p0),*(Are+p1),*(Aim+p1),*(Are+p2),*(Aim+p2),*(Are+p3),*(Aim+p3),*(Are+p4),*(Aim+p4),*(Are+p5),*(Aim+p5),*(Are+p6),*(Aim+p6),*(Are+p7),*(Aim+p7),*(Are+p8),*(Aim+p8),*(Are+p9),*(Aim+p9),*(Are+pa),*(Aim+pa),*(Are+pb),*(Aim+pb),*(Are+pc),*(Aim+pc),*(Are+pd),*(Aim+pd),*(Are+pe),*(Aim+pe),*(Are+pf),*(Aim+pf),
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16
		);	tptr += 16;
	}

/*...and now do 16 radix-8 subtransforms w/internal twiddles - cf. radi1024_dif_pass1 for details. Process blocks in BR16 order: */
	tptr = t;	j = 0;
	/* Block 0 gets pulled out of parametrized-loop to take advantage of its unity-twiddles-ness: */
	// Remember, the twiddleless DIF also bit-reverses its outputs, so a_p* terms appear in BR-order here [swap index pairs 1/4 and 3/6]:
	RADIX_08_DIF_OOP(
		 tptr      ->re, tptr      ->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x60)->re,(tptr+0x60)->im,
		(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x70)->re,(tptr+0x70)->im,
		*(Bre+__odx[0x00]),*(Bim+__odx[0x00]),*(Bre+__odx[0x04]),*(Bim+__odx[0x04]),*(Bre+__odx[0x02]),*(Bim+__odx[0x02]),*(Bre+__odx[0x06]),*(Bim+__odx[0x06]),*(Bre+__odx[0x01]),*(Bim+__odx[0x01]),*(Bre+__odx[0x05]),*(Bim+__odx[0x05]),*(Bre+__odx[0x03]),*(Bim+__odx[0x03]),*(Bre+__odx[0x07]),*(Bim+__odx[0x07])
	); tptr++;
	/* Blocks 84c2a6e195d3b7f processed via loop: */
	cd_ptr = (struct complex *)DFT128_TWIDDLES[1];
	for(i = 1; i < 16; i++) {
		j += 8;	// jt = j1 + p08;	jp = j2 + p08;
		p0 = __odx[j]; p1 = __odx[j+1]; p2 = __odx[j+2]; p3 = __odx[j+3]; p4 = __odx[j+4]; p5 = __odx[j+5]; p6 = __odx[j+6]; p7 = __odx[j+7];
		RADIX_08_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
			*(Bre+p0),*(Bim+p0),*(Bre+p1),*(Bim+p1),*(Bre+p2),*(Bim+p2),*(Bre+p3),*(Bim+p3),*(Bre+p4),*(Bim+p4),*(Bre+p5),*(Bim+p5),*(Bre+p6),*(Bim+p6),*(Bre+p7),*(Bim+p7),
			cd_ptr->re,cd_ptr->im,(cd_ptr+1)->re,(cd_ptr+1)->im,(cd_ptr+2)->re,(cd_ptr+2)->im,(cd_ptr+3)->re,(cd_ptr+3)->im,(cd_ptr+4)->re,(cd_ptr+4)->im,(cd_ptr+5)->re,(cd_ptr+5)->im,(cd_ptr+6)->re,(cd_ptr+6)->im
		); tptr++; cd_ptr += 7;
	}
}

void RADIX_128_DIT(
	double *__A, const int *__idx, const int __re_im_stride_in,	/*  Inputs: Base address plus 128 (index) offsets */
	double *__B, const int *__odx, const int __re_im_stride_out	/* Outputs: Base address plus 128 (index) offsets */
)
{
	int i,j;
	struct complex t[128], *tptr;
	const struct complex *cd_ptr;
	double *Are = __A,*Aim = __A + __re_im_stride_in, *Bre = __B,*Bim = __B + __re_im_stride_out;
	int p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf;

// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in ascendingorder:
	tptr = t;
	for(i = 0; i < 8; i++) {
		j = i<<4;
		p0 = __idx[j]; p1 = __idx[j+0x1]; p2 = __idx[j+0x2]; p3 = __idx[j+0x3]; p4 = __idx[j+0x4]; p5 = __idx[j+0x5]; p6 = __idx[j+0x6]; p7 = __idx[j+0x7]; p8 = __idx[j+0x8]; p9 = __idx[j+0x9]; pa = __idx[j+0xa]; pb = __idx[j+0xb]; pc = __idx[j+0xc]; pd = __idx[j+0xd]; pe = __idx[j+0xe]; pf = __idx[j+0xf];
		RADIX_16_DIT(
			*(Are+p0),*(Aim+p0),*(Are+p1),*(Aim+p1),*(Are+p2),*(Aim+p2),*(Are+p3),*(Aim+p3),*(Are+p4),*(Aim+p4),*(Are+p5),*(Aim+p5),*(Are+p6),*(Aim+p6),*(Are+p7),*(Aim+p7),*(Are+p8),*(Aim+p8),*(Are+p9),*(Aim+p9),*(Are+pa),*(Aim+pa),*(Are+pb),*(Aim+pb),*(Are+pc),*(Aim+pc),*(Are+pd),*(Aim+pd),*(Are+pe),*(Aim+pe),*(Are+pf),*(Aim+pf),
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16
		);	tptr += 16;
	}

/*...and now do 16 radix-8 subtransforms w/internal twiddles - cf. radix128_dit_pass1 for details: */
	tptr = t;
	/* Block 0 gets pulled out of parametrized-loop to take advantage of its unity-twiddles-ness: */
	// Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here:
	RADIX_08_DIT_OOP(
		 tptr      ->re, tptr      ->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,
		(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
		*(Bre+__odx[0x00]),*(Bim+__odx[0x00]),*(Bre+__odx[0x10]),*(Bim+__odx[0x10]),*(Bre+__odx[0x20]),*(Bim+__odx[0x20]),*(Bre+__odx[0x30]),*(Bim+__odx[0x30]),*(Bre+__odx[0x40]),*(Bim+__odx[0x40]),*(Bre+__odx[0x50]),*(Bim+__odx[0x50]),*(Bre+__odx[0x60]),*(Bim+__odx[0x60]),*(Bre+__odx[0x70]),*(Bim+__odx[0x70])
	);
	/* Blocks 84c2a6e195d3b7f processed via loop: */
	cd_ptr = (struct complex *)DFT128_TWIDDLES[1];
	for(i = 1; i < 16; i++) {
		j = reverse16[i];
		tptr = t+j;
		p0 = __odx[j]; p1 = __odx[j+0x40]; p2 = __odx[j+0x20]; p3 = __odx[j+0x60]; p4 = __odx[j+0x10]; p5 = __odx[j+0x50]; p6 = __odx[j+0x30]; p7 = __odx[j+0x70];
		RADIX_08_DIT_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,
			*(Bre+p0),*(Bim+p0),*(Bre+p1),*(Bim+p1),*(Bre+p2),*(Bim+p2),*(Bre+p3),*(Bim+p3),*(Bre+p4),*(Bim+p4),*(Bre+p5),*(Bim+p5),*(Bre+p6),*(Bim+p6),*(Bre+p7),*(Bim+p7),
			cd_ptr->re,cd_ptr->im,(cd_ptr+1)->re,(cd_ptr+1)->im,(cd_ptr+2)->re,(cd_ptr+2)->im,(cd_ptr+3)->re,(cd_ptr+3)->im,(cd_ptr+4)->re,(cd_ptr+4)->im,(cd_ptr+5)->re,(cd_ptr+5)->im,(cd_ptr+6)->re,(cd_ptr+6)->im
		); cd_ptr += 7;
	}
}

/************** RADIX-256 DIF/DIT: *****************************/

void RADIX_256_DIF(
	// inputs: Base address plus 32 index offsets:
	double *__A,
	// Since i_array may either coincide with oarray [if using it for a complete radix-256 pass]
	// or a local-complex-data scratch array [if using for phase of a larger-radix DFT algo] and we need
	// this routine to work correctly with both scalar-complex and SIMD main-array data layouts, need a param
	// telling routine whether to use scalar-complex re/im stride of 1 or SIMD [2,4,8... depending on type]:
	const int __re_im_stride_in,
	int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
	int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
	// outputs: Base address plus 32 index offsets:
	double *__B,
	const int __re_im_stride_out,	// Similar as for in-array
	int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
	uint32 o_idx,	// Bitfield encoding the sequence of the o_offsets_lo sub-vectors to use for the radix-256 DFT's outputs
	int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
)
{
	struct complex t[256], *tptr;
	const double *addr, *addi;
	double *Are,*Aim, *Bre,*Bim;
	// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
	int i,j,nshift, *off_ptr;
	int p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf;
	int p00,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0;

// Gather the needed data and do 16 twiddleless length-16 subtransforms, with p-offsets in br8 order: 084c2a6e195d3b7f:
// NOTE that RADIX_16_DIF outputs are IN-ORDER rather than BR:
	p00 = i_offsets_hi[0x0];p10 = i_offsets_hi[0x1];p20 = i_offsets_hi[0x2];p30 = i_offsets_hi[0x3];p40 = i_offsets_hi[0x4];p50 = i_offsets_hi[0x5];p60 = i_offsets_hi[0x6];p70 = i_offsets_hi[0x7];p80 = i_offsets_hi[0x8];p90 = i_offsets_hi[0x9];pa0 = i_offsets_hi[0xa];pb0 = i_offsets_hi[0xb];pc0 = i_offsets_hi[0xc];pd0 = i_offsets_hi[0xd];pe0 = i_offsets_hi[0xe];pf0 = i_offsets_hi[0xf];

	tptr = t;
	for(i = 0; i < 16; i++) {
		j = reverse16[i];	// A-array offsets processed in BR16 order = p[084c2a6e195d3b7f]
		Are = __A + i_offsets_lo[j]; Aim = Are + __re_im_stride_in;
		RADIX_16_DIF(
			*(Are+p00),*(Aim+p00),*(Are+p10),*(Aim+p10),*(Are+p20),*(Aim+p20),*(Are+p30),*(Aim+p30),*(Are+p40),*(Aim+p40),*(Are+p50),*(Aim+p50),*(Are+p60),*(Aim+p60),*(Are+p70),*(Aim+p70),*(Are+p80),*(Aim+p80),*(Are+p90),*(Aim+p90),*(Are+pa0),*(Aim+pa0),*(Are+pb0),*(Aim+pb0),*(Are+pc0),*(Aim+pc0),*(Are+pd0),*(Aim+pd0),*(Are+pe0),*(Aim+pe0),*(Are+pf0),*(Aim+pf0),
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);	tptr += 16;
	}

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors: */

	// Extract index of the 16-element o_offsets_lo sub-vector to use for the current set of outputs:
	off_ptr = o_offsets_lo + ( (o_idx&0x3) << 4 );	// Low 2 bits of o_idx; loop below will use remaining 30 bits in ascending pairs
	p0  = off_ptr[0x0];p1  = off_ptr[0x1];p2  = off_ptr[0x2];p3  = off_ptr[0x3];p4  = off_ptr[0x4];p5  = off_ptr[0x5];p6  = off_ptr[0x6];p7  = off_ptr[0x7];p8  = off_ptr[0x8];p9  = off_ptr[0x9];pa  = off_ptr[0xa];pb  = off_ptr[0xb];pc  = off_ptr[0xc];pd  = off_ptr[0xd];pe  = off_ptr[0xe];pf  = off_ptr[0xf];

	// Block 0: Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
	tptr = t;
	Bre = __B + o_offsets_hi[0]; Bim = Bre + __re_im_stride_out;
	RADIX_16_DIF(
		tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		*(Bre+p0),*(Bim+p0),*(Bre+p1),*(Bim+p1),*(Bre+p2),*(Bim+p2),*(Bre+p3),*(Bim+p3),*(Bre+p4),*(Bim+p4),*(Bre+p5),*(Bim+p5),*(Bre+p6),*(Bim+p6),*(Bre+p7),*(Bim+p7),*(Bre+p8),*(Bim+p8),*(Bre+p9),*(Bim+p9),*(Bre+pa),*(Bim+pa),*(Bre+pb),*(Bim+pb),*(Bre+pc),*(Bim+pc),*(Bre+pd),*(Bim+pd),*(Bre+pe),*(Bim+pe),*(Bre+pf),*(Bim+pf),
		c16,s16
	);	tptr++;

	// Remaining 15 sets of macro calls done in loop:
	for(i = 1; i < 16; i++) {
		if(o_idx) {
			nshift = i+i;	// o_idx shift counts here run as >>2,4,...,30
			off_ptr = o_offsets_lo + ( ((o_idx>>nshift)&0x3) << 4 );
			p0 = off_ptr[0x0];p1 = off_ptr[0x1];p2 = off_ptr[0x2];p3 = off_ptr[0x3];p4 = off_ptr[0x4];p5 = off_ptr[0x5];p6 = off_ptr[0x6];p7 = off_ptr[0x7];p8 = off_ptr[0x8];p9 = off_ptr[0x9];pa = off_ptr[0xa];pb = off_ptr[0xb];pc = off_ptr[0xc];pd = off_ptr[0xd];pe = off_ptr[0xe];pf = off_ptr[0xf];
		}
		addr = DFT256_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array [whose row order already BRed above]
		Bre = __B + o_offsets_hi[i]; Bim = Bre + __re_im_stride_out;	// o_offsets_hi[] = p10,p20,...,pf0
		RADIX_16_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
			*(Bre+p0),*(Bim+p0),*(Bre+p1),*(Bim+p1),*(Bre+p2),*(Bim+p2),*(Bre+p3),*(Bim+p3),*(Bre+p4),*(Bim+p4),*(Bre+p5),*(Bim+p5),*(Bre+p6),*(Bim+p6),*(Bre+p7),*(Bim+p7),*(Bre+p8),*(Bim+p8),*(Bre+p9),*(Bim+p9),*(Bre+pa),*(Bim+pa),*(Bre+pb),*(Bim+pb),*(Bre+pc),*(Bim+pc),*(Bre+pd),*(Bim+pd),*(Bre+pe),*(Bim+pe),*(Bre+pf),*(Bim+pf),
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);	tptr++;
	}
}

void RADIX_256_DIT(
	// inputs: Base address plus 32 index offsets:
	double *__A,
	const int __re_im_stride_in,
	int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
	uint32 i_idx,	// Bitfield encoding the sequence of the i_offsets_lo sub-vectors to use for the radix-256 DFT's inputs
	int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
	// outputs: Base address plus 32 index offsets:
	double *__B,
	const int __re_im_stride_out,
	int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
	int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
)
{
	struct complex t[256], *tptr;
	const double *addr, *addi;
	double *Are,*Aim, *Bre,*Bim;
	// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
	int i,j,nshift, *off_ptr;
	// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
	int p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf;
	int p00,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0;

	// Extract index of the 16-element i_offsets_lo sub-vector to use for the current set of outputs:
	off_ptr = i_offsets_lo + ( (i_idx&0x3) << 4 );	// 16*(Low 2 bits of i_idx); loop below will use remaining 30 bits in ascending pairs
	p0  = off_ptr[0x0];p1  = off_ptr[0x1];p2  = off_ptr[0x2];p3  = off_ptr[0x3];p4  = off_ptr[0x4];p5  = off_ptr[0x5];p6  = off_ptr[0x6];p7  = off_ptr[0x7];p8  = off_ptr[0x8];p9  = off_ptr[0x9];pa  = off_ptr[0xa];pb  = off_ptr[0xb];pc  = off_ptr[0xc];pd  = off_ptr[0xd];pe  = off_ptr[0xe];pf  = off_ptr[0xf];

// Gather the needed data and do 16 twiddleless length-16 subtransforms, with p-offsets in-order:

	tptr = t;
	for(i = 0; i < 16; i++) {
		Are = __A + i_offsets_hi[i]; Aim = Are + __re_im_stride_in;	// i_offsets_hi[] = p10,p20,...,pf0
		RADIX_16_DIT(
			*(Are+p0),*(Aim+p0),*(Are+p1),*(Aim+p1),*(Are+p2),*(Aim+p2),*(Are+p3),*(Aim+p3),*(Are+p4),*(Aim+p4),*(Are+p5),*(Aim+p5),*(Are+p6),*(Aim+p6),*(Are+p7),*(Aim+p7),*(Are+p8),*(Aim+p8),*(Are+p9),*(Aim+p9),*(Are+pa),*(Aim+pa),*(Are+pb),*(Aim+pb),*(Are+pc),*(Aim+pc),*(Are+pd),*(Aim+pd),*(Are+pe),*(Aim+pe),*(Are+pf),*(Aim+pf),
			tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);	tptr += 16;
		if(i_idx) {	// vvv +2 here because this is setup for next loop pass
			nshift = i+i+2;	// i_idx shift counts here run as >>2,4,...,30
			off_ptr = i_offsets_lo + ( ((i_idx>>nshift)&0x3) << 4 );
			p0  = off_ptr[0x0];p1  = off_ptr[0x1];p2  = off_ptr[0x2];p3  = off_ptr[0x3];p4  = off_ptr[0x4];p5  = off_ptr[0x5];p6  = off_ptr[0x6];p7  = off_ptr[0x7];p8  = off_ptr[0x8];p9  = off_ptr[0x9];pa  = off_ptr[0xa];pb  = off_ptr[0xb];pc  = off_ptr[0xc];pd  = off_ptr[0xd];pe  = off_ptr[0xe];pf  = off_ptr[0xf];
		}
	}

/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
*/
	p00 = o_offsets_hi[0x0];p10 = o_offsets_hi[0x1];p20 = o_offsets_hi[0x2];p30 = o_offsets_hi[0x3];p40 = o_offsets_hi[0x4];p50 = o_offsets_hi[0x5];p60 = o_offsets_hi[0x6];p70 = o_offsets_hi[0x7];p80 = o_offsets_hi[0x8];p90 = o_offsets_hi[0x9];pa0 = o_offsets_hi[0xa];pb0 = o_offsets_hi[0xb];pc0 = o_offsets_hi[0xc];pd0 = o_offsets_hi[0xd];pe0 = o_offsets_hi[0xe];pf0 = o_offsets_hi[0xf];

	tptr = t;
	// Block 0: jt = j1;	jp = j2;
	Bre = __B + o_offsets_lo[0]; Bim = Bre + __re_im_stride_out;
	RADIX_16_DIT(
		tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
		*(Bre+p00),*(Bim+p00),*(Bre+p10),*(Bim+p10),*(Bre+p20),*(Bim+p20),*(Bre+p30),*(Bim+p30),*(Bre+p40),*(Bim+p40),*(Bre+p50),*(Bim+p50),*(Bre+p60),*(Bim+p60),*(Bre+p70),*(Bim+p70),*(Bre+p80),*(Bim+p80),*(Bre+p90),*(Bim+p90),*(Bre+pa0),*(Bim+pa0),*(Bre+pb0),*(Bim+pb0),*(Bre+pc0),*(Bim+pc0),*(Bre+pd0),*(Bim+pd0),*(Bre+pe0),*(Bim+pe0),*(Bre+pf0),*(Bim+pf0),
		c16,s16
	);

	// Remaining 15 sets of macro calls done in loop:
	for(i = 1; i < 16; i++) {
		j = reverse16[i];
		tptr = t + j;
		Bre = __B + o_offsets_lo[j]; Bim = Bre + __re_im_stride_out;	// o_offsets_lo[] = p[084c2a6e195d3b7f]
		addr = DFT256_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_16_DIT_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xd0)->re,(tptr+0xd0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0xf0)->re,(tptr+0xf0)->im,
			*(Bre+p00),*(Bim+p00),*(Bre+p10),*(Bim+p10),*(Bre+p20),*(Bim+p20),*(Bre+p30),*(Bim+p30),*(Bre+p40),*(Bim+p40),*(Bre+p50),*(Bim+p50),*(Bre+p60),*(Bim+p60),*(Bre+p70),*(Bim+p70),*(Bre+p80),*(Bim+p80),*(Bre+p90),*(Bim+p90),*(Bre+pa0),*(Bim+pa0),*(Bre+pb0),*(Bim+pb0),*(Bre+pc0),*(Bim+pc0),*(Bre+pd0),*(Bim+pd0),*(Bre+pe0),*(Bim+pe0),*(Bre+pf0),*(Bim+pf0),
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);
	}
}

// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC)

	/************** RADIX-63 DIF/DIT: *****************************/

	void SSE2_RADIX_63_DIF(
		const int init,	// Init consts (in 1-thread mode only!) and exit
		const int thr_id,
		vec_dbl *__A, const int *__idx,	/* Inputs : Base address plus 63 (index) offsets */
		vec_dbl *__B, const int *__odx	/* Outputs: Base address plus 63 (index) offsets */
	)
	{
		static int max_threads = 0;
		static vec_dbl *sc_arr = 0x0, *sc_ptr;	// Pad with 4 extra slots for scratch storage needed by SSE2_RADIX_07_DFT macro!!!
	  #ifdef MULTITHREAD
		static vec_dbl *__r0;					// Base address for discrete per-thread local stores
		// In || mode, only above base-pointer (shared by all threads) is static:
		vec_dbl *two,*one,	// Needed for FMA support
				*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,	// radix-7 DFT roots; Use dc0-ptr as "need to init local statics?" sentinel
				*cc1,*ss1,*cc2,*ss2,*cc3m1,*ss3,*cc4,*ss4;	// radix-9 DFT roots
	  #else
		static vec_dbl *two,*one,
				*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,
				*cc1,*ss1,*cc2,*ss2,*cc3m1,*ss3,*cc4,*ss4;	// 18 vec_dbl consts plus 4 padding slots = 22 slots alloc per thread
	  #endif
		// Unline e.g. the 64-DFTs in this file, the 63-DFT defines its own local storage for the
		// 2x63 vec_dbl intermediates needed by each thread. That storage is pointed to by tdat:
		vec_dbl *tdat, *tmp,
			*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,
			*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8;
		int l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
		const uint8 *iptr,
			dif_iperm[64] = {	// Only need 63, but pad to 64-bytes
				0x00,0x36,0x2d,0x24,0x1b,0x12,0x09,
				0x38,0x2f,0x26,0x1d,0x14,0x0b,0x02,
				0x31,0x28,0x1f,0x16,0x0d,0x04,0x3a,
				0x2a,0x21,0x18,0x0f,0x06,0x3c,0x33,
				0x23,0x1a,0x11,0x08,0x3e,0x35,0x2c,
				0x1c,0x13,0x0a,0x01,0x37,0x2e,0x25,
				0x15,0x0c,0x03,0x39,0x30,0x27,0x1e,
				0x0e,0x05,0x3b,0x32,0x29,0x20,0x17,
				0x07,0x3d,0x34,0x2b,0x22,0x19,0x10},
			dif_operm[64] = {	// ditto
				0x00,0x07,0x03,0x02,0x06,0x05,0x01,0x08,0x04,
				0x37,0x3e,0x3a,0x36,0x3d,0x39,0x38,0x3c,0x3b,
				0x32,0x2e,0x35,0x31,0x2d,0x34,0x30,0x2f,0x33,
				0x2a,0x29,0x25,0x2c,0x28,0x24,0x2b,0x27,0x26,
				0x1d,0x21,0x20,0x1c,0x23,0x1f,0x1b,0x22,0x1e,
				0x15,0x14,0x18,0x17,0x13,0x1a,0x16,0x12,0x19,
				0x10,0x0c,0x0b,0x0f,0x0e,0x0a,0x11,0x0d,0x09};
	  // AVX2 (i.e. FMA) means non-Nussbaumer radix-7, uses these sincos constants:
	  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
		// AVX2 (i.e. FMA) means non-Nussbaumer radix-7, uses these sincos constants:
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/
	  #else
		// Roots for radix-7 DFTs: SSE2/AVX versions assume LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation:
		const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
						cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
						cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
						cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
					/* Switch the sign of ss3 in these: */
						sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
						sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
						sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
						sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
	  #endif	// AVX2?
		// Roots for radix-9 DFTs:
		const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
						s   =  0.64278760968653932631,	/* sin(2*pi/9) */
						c2  =  0.17364817766693034887,	/* cos(2*u) */
						s2  =  0.98480775301220805936,	/* sin(2*u) */
						c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
						s3  =  0.86602540378443864677,	/* sin(3*u) */
						c4  = -0.93969262078590838404,	/* cos(4*u) */
						s4  =  0.34202014332566873307;	/* sin(4*u) */

		// If this is first time here, init pointers and associated data:
		if(thr_id == -1)	// Value of init stores #threads
		{
			if(init <= max_threads) {	// Previously inited with sufficient #threads
				ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
				return;
			}
			max_threads = init;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
			if(sc_arr) { free((void *)sc_arr); }
			// 126 slots for DFT-63 data, 22 for DFT-7,9 consts and DFT-7 pads, 4 to allow for alignment = 152:
			sc_arr = ALLOC_VEC_DBL(sc_arr, 152*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = ALIGN_VEC_DBL(sc_arr);
			ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		#ifdef MULTITHREAD
			__r0 = tdat = sc_ptr;	tmp = tdat + 126;
			two    = tmp + 0x0;
			one    = tmp + 0x1;
			tmp += 2;
			// Roots for radix-9 DFTs:	// Roots for radix-7 DFTs: MUST be followed by 4 alloc'ed padding slots needed by SSE2_RADIX_07_DFT macro!!! **************/
			cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
			ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
			cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
			ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
			cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
			ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
			cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
			ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
			for(l = 0; l < max_threads; ++l) {
			/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
			  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
				// FMA means non-Nussbaumer radix-7, use these sincos constants:
				VEC_DBL_INIT(dc0, uc1  );	VEC_DBL_INIT(ds0, us1);
				VEC_DBL_INIT(dc1, uc2  );	VEC_DBL_INIT(ds1, us2);
				VEC_DBL_INIT(dc2, uc3  );	VEC_DBL_INIT(ds2, us3);
				VEC_DBL_INIT(dc3, 0.0  );	VEC_DBL_INIT(ds3, 0.0);	// Unused in non-Nussbaumer mode
			  #else
				// Roots for radix-7 DFTs: cc2 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation
				VEC_DBL_INIT(dc0, cx0-1);	VEC_DBL_INIT(ds0, sx0);
				VEC_DBL_INIT(dc1, cx1  );	VEC_DBL_INIT(ds1, sx1);
				VEC_DBL_INIT(dc2, cx2  );	VEC_DBL_INIT(ds2, sx2);
				VEC_DBL_INIT(dc3, cx3  );	VEC_DBL_INIT(ds3, sx3);
			  #endif
				// Roots for radix-9 DFTs:
				VEC_DBL_INIT(cc1  , c	);	VEC_DBL_INIT(ss1, s );
				VEC_DBL_INIT(cc2  , c2  );	VEC_DBL_INIT(ss2, s2);
				VEC_DBL_INIT(cc3m1, c3m1);	VEC_DBL_INIT(ss3, s3);
				VEC_DBL_INIT(cc4  , c4  );	VEC_DBL_INIT(ss4, s4);
			/* Move on to next thread's local store */
				cc1   += 152;			dc0 += 152;
				ss1   += 152;			ds0 += 152;
				cc2   += 152;			dc1 += 152;
				ss2   += 152;			ds1 += 152;
				cc3m1 += 152;			dc2 += 152;
				ss3   += 152;			ds2 += 152;
				cc4   += 152;			dc3 += 152;
				ss4   += 152;			ds3 += 152;
			}
		#else
			tdat = sc_ptr;	tmp = tdat + 126;
			two    = tmp + 0x0;
			one    = tmp + 0x1;
			tmp += 2;
			// Roots for radix-9 DFTs:	// Roots for radix-7 DFTs: MUST be followed by 4 alloc'ed padding slots needed by SSE2_RADIX_07_DFT macro!!! **************/
			cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
			ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
			cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
			ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
			cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
			ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
			cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
			ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
		  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
			// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
			VEC_DBL_INIT(dc0, uc1  );	VEC_DBL_INIT(ds0, us1);
			VEC_DBL_INIT(dc1, uc2  );	VEC_DBL_INIT(ds1, us2);
			VEC_DBL_INIT(dc2, uc3  );	VEC_DBL_INIT(ds2, us3);
			VEC_DBL_INIT(dc3, 0.0  );	VEC_DBL_INIT(ds3, 0.0);	// Unused in non-Nussbaumer mode
		  #else
			// Roots for radix-7 DFTs: cc2 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation
			VEC_DBL_INIT(dc0, cx0-1);	VEC_DBL_INIT(ds0, sx0);
			VEC_DBL_INIT(dc1, cx1  );	VEC_DBL_INIT(ds1, sx1);
			VEC_DBL_INIT(dc2, cx2  );	VEC_DBL_INIT(ds2, sx2);
			VEC_DBL_INIT(dc3, cx3  );	VEC_DBL_INIT(ds3, sx3);
		  #endif
			// Roots for radix-9 DFTs:
			VEC_DBL_INIT(cc1  , c	);	VEC_DBL_INIT(ss1, s );
			VEC_DBL_INIT(cc2  , c2  );	VEC_DBL_INIT(ss2, s2);
			VEC_DBL_INIT(cc3m1, c3m1);	VEC_DBL_INIT(ss3, s3);
			VEC_DBL_INIT(cc4  , c4  );	VEC_DBL_INIT(ss4, s4);
		#endif
			return;
		} else {
			ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		tdat = __r0 + thr_id*152;	tmp = tdat + 126;
		two    = tmp + 0x0;
		one    = tmp + 0x1;
		tmp += 2;
		cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
		ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
		cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
		ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
		cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
		ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
		cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
		ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
	#else
		tdat = sc_ptr;
	#endif

		//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 9 radix-7 transforms:
		/*
		Twiddleless version arranges 9 sets of radix-7 DFT inputs as follows: 0 in upper left corner,
		decrement (mod 63) 9 horizontally and 7 vertically. Display result of DIF/DIT input-scramble array in hex:

			00,36,2d,24,1b,12,09
			38,2f,26,1d,14,0b,02
			31,28,1f,16,0d,04,3a
			2a,21,18,0f,06,3c,33
			23,1a,11,08,3e,35,2c
			1c,13,0a,01,37,2e,25
			15,0c,03,39,30,27,1e
			0e,05,3b,32,29,20,17
			07,3d,34,2b,22,19,10
		*/
		tmp = tdat; iptr = dif_iperm;
		for(l = 0; l < 9; l++) {
			k0 = __idx[*iptr]; k1 = __idx[*(iptr+1)]; k2 = __idx[*(iptr+2)]; k3 = __idx[*(iptr+3)]; k4 = __idx[*(iptr+4)]; k5 = __idx[*(iptr+5)]; k6 = __idx[*(iptr+6)];
			va0 = __A+k0; va1 = __A+k1; va2 = __A+k2; va3 = __A+k3; va4 = __A+k4; va5 = __A+k5; va6 = __A+k6;
			// Since there is no vec_cmplx type, ptr-offs double those of analogous scalar code:
			vb0 = tmp; vb1 = tmp+18; vb2 = tmp+36; vb3 = tmp+54; vb4 = tmp+72; vb5 = tmp+90; vb6 = tmp+108;
		  #ifdef USE_AVX2
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, dc0,two, vb0,vb1,vb2,vb3,vb4,vb5,vb6)
		  #else
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, dc0,     vb0,vb1,vb2,vb3,vb4,vb5,vb6)
		  #endif
			tmp += 2; iptr += 7;
		}
		/*...and now do 7 radix-9 transforms. The required output permutation is

			00,07,03,02,06,05,01,08,04,
			37,3e,3a,36,3d,39,38,3c,3b,
			32,2e,35,31,2d,34,30,2f,33,
			2a,29,25,2c,28,24,2b,27,26,
			1d,21,20,1c,23,1f,1b,22,1e,
			15,14,18,17,13,1a,16,12,19,
			10,0c,0b,0f,0e,0a,11,0d,09.
		*/
		tmp = tdat; iptr = dif_operm;
		for(l = 0; l < 7; l++) {
			// When 63 is used to build a larger DFT radix (e.g. 1008), these indices will be permuted (nonmonotone), no simplification possible:
			k0 = __odx[*iptr]; k1 = __odx[*(iptr+1)]; k2 = __odx[*(iptr+2)]; k3 = __odx[*(iptr+3)]; k4 = __odx[*(iptr+4)]; k5 = __odx[*(iptr+5)]; k6 = __odx[*(iptr+6)]; k7 = __odx[*(iptr+7)]; k8 = __odx[*(iptr+8)];
			// Since there is no vec_cmplx type, ptr-offs double those of analogous scalar code:
			va0 = tmp; va1 = tmp+2; va2 = tmp+4; va3 = tmp+6; va4 = tmp+8; va5 = tmp+10; va6 = tmp+12; va7 = tmp+14; va8 = tmp+16;
			vb0 = __B+k0; vb1 = __B+k1; vb2 = __B+k2; vb3 = __B+k3; vb4 = __B+k4; vb5 = __B+k5; vb6 = __B+k6; vb7 = __B+k7; vb8 = __B+k8;
			SSE2_RADIX_09_DIF(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,
				cc1,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8
			);	tmp += 18; iptr += 9;
		}
	}

	/***************/

	void SSE2_RADIX_63_DIT(
		const int init,	// Init consts (in 1-thread mode only!) and exit
		const int thr_id,
		vec_dbl *__A, const int *__idx,	/* Inputs : Base address plus 63 (index) offsets */
		vec_dbl *__B, const int *__odx	/* Outputs: Base address plus 63 (index) offsets */
	)
	{
		static int max_threads = 0;
		static vec_dbl *sc_arr = 0x0, *sc_ptr;	// Pad with 4 extra slots for scratch storage needed by SSE2_RADIX_07_DFT macro!!!
	  #ifdef MULTITHREAD
		static vec_dbl *__r0;					// Base address for discrete per-thread local stores
		// In || mode, only above base-pointer (shared by all threads) is static:
		vec_dbl *two,*one,	// Needed for FMA support
				*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,	// radix-7 DFT roots; Use dc0-ptr as "need to init local statics?" sentinel
				*cc1,*ss1,*cc2,*ss2,*cc3m1,*ss3,*cc4,*ss4;	// radix-9 DFT roots
	  #else
		static vec_dbl *two,*one,
				*dc0,*ds0,*dc1,*ds1,*dc2,*ds2,*dc3,*ds3,
				*cc1,*ss1,*cc2,*ss2,*cc3m1,*ss3,*cc4,*ss4;	// 18 vec_dbl consts plus 4 padding slots = 22 slots alloc per thread
	  #endif
		// Unline e.g. the 64-DFTs in this file, the 63-DFT defines its own local storage for the
		// 2x63 vec_dbl intermediates needed by each thread. That storage is pointed to by tdat:
		vec_dbl *tdat, *tmp,
			*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,
			*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8;
		int l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
		const uint8 *iptr,
			dit_iperm[64] = {	// Only need 63, but pad to 64-bytes
				0x00,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,
				0x32,0x31,0x30,0x2f,0x2e,0x2d,0x34,0x33,0x35,
				0x1d,0x1c,0x1b,0x22,0x21,0x23,0x1f,0x1e,0x20,
				0x10,0x0f,0x11,0x0d,0x0c,0x0e,0x0a,0x09,0x0b,
				0x37,0x36,0x38,0x3c,0x3e,0x3d,0x39,0x3b,0x3a,
				0x2a,0x2c,0x2b,0x27,0x29,0x28,0x24,0x26,0x25,
				0x15,0x17,0x16,0x12,0x14,0x13,0x1a,0x19,0x18,0},
			dit_operm[64] = {	// ditto
				0x00,0x24,0x09,0x2d,0x12,0x36,0x1b,
				0x0e,0x32,0x17,0x3b,0x20,0x05,0x29,
				0x1c,0x01,0x25,0x0a,0x2e,0x13,0x37,
				0x2a,0x0f,0x33,0x18,0x3c,0x21,0x06,
				0x38,0x1d,0x02,0x26,0x0b,0x2f,0x14,
				0x07,0x2b,0x10,0x34,0x19,0x3d,0x22,
				0x15,0x39,0x1e,0x03,0x27,0x0c,0x30,
				0x23,0x08,0x2c,0x11,0x35,0x1a,0x3e,
				0x31,0x16,0x3a,0x1f,0x04,0x28,0x0d,0};
	  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
		// FMA means non-Nussbaumer radix-7, use these sincos constants:
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/
	  #else
		// Roots for radix-7 DFTs: SSE2/AVX versions assume LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation:
		const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
						cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
						cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
						cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
					/* Switch the sign of ss3 in these: */
						sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
						sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
						sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
						sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
	  #endif	// AVX2?
		// Roots for radix-9 DFTs:
		const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
						s   =  0.64278760968653932631,	/* sin(2*pi/9) */
						c2  =  0.17364817766693034887,	/* cos(2*u) */
						s2  =  0.98480775301220805936,	/* sin(2*u) */
						c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
						s3  =  0.86602540378443864677,	/* sin(3*u) */
						c4  = -0.93969262078590838404,	/* cos(4*u) */
						s4  =  0.34202014332566873307;	/* sin(4*u) */

		// If this is first time here, init pointers and associated data:
		if(thr_id == -1)	// Value of init stores #threads
		{
			if(init <= max_threads) {	// Previously inited with sufficient #threads
				ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
				return;
			}
			max_threads = init;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
			if(sc_arr) { free((void *)sc_arr); }
			// 126 slots for DFT-63 data, 22 for DFT-7,9 consts and DFT-7 pads, 4 to allow for alignment = 152:
			sc_arr = ALLOC_VEC_DBL(sc_arr, 152*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = ALIGN_VEC_DBL(sc_arr);
			ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		#ifdef MULTITHREAD
			__r0 = tdat = sc_ptr;	tmp = tdat + 126;
			two    = tmp + 0x0;
			one    = tmp + 0x1;
			tmp += 2;
			// Roots for radix-9 DFTs:	// Roots for radix-7 DFTs: MUST be followed by 4 alloc'ed padding slots needed by SSE2_RADIX_07_DFT macro!!! **************/
			cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
			ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
			cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
			ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
			cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
			ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
			cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
			ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
			for(l = 0; l < max_threads; ++l) {
			/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
			  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
				// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
				VEC_DBL_INIT(dc0, uc1  );	VEC_DBL_INIT(ds0, us1);
				VEC_DBL_INIT(dc1, uc2  );	VEC_DBL_INIT(ds1, us2);
				VEC_DBL_INIT(dc2, uc3  );	VEC_DBL_INIT(ds2, us3);
				VEC_DBL_INIT(dc3, 0.0  );	VEC_DBL_INIT(ds3, 0.0);	// Unused in non-Nussbaumer mode
			  #else
				// Roots for radix-7 DFTs: cc2 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation
				VEC_DBL_INIT(dc0, cx0-1);	VEC_DBL_INIT(ds0, sx0);
				VEC_DBL_INIT(dc1, cx1  );	VEC_DBL_INIT(ds1, sx1);
				VEC_DBL_INIT(dc2, cx2  );	VEC_DBL_INIT(ds2, sx2);
				VEC_DBL_INIT(dc3, cx3  );	VEC_DBL_INIT(ds3, sx3);
			  #endif
				// Roots for radix-9 DFTs:
				VEC_DBL_INIT(cc1  , c	);	VEC_DBL_INIT(ss1, s );
				VEC_DBL_INIT(cc2  , c2  );	VEC_DBL_INIT(ss2, s2);
				VEC_DBL_INIT(cc3m1, c3m1);	VEC_DBL_INIT(ss3, s3);
				VEC_DBL_INIT(cc4  , c4  );	VEC_DBL_INIT(ss4, s4);
			/* Move on to next thread's local store */
				cc1   += 152;			dc0 += 152;
				ss1   += 152;			ds0 += 152;
				cc2   += 152;			dc1 += 152;
				ss2   += 152;			ds1 += 152;
				cc3m1 += 152;			dc2 += 152;
				ss3   += 152;			ds2 += 152;
				cc4   += 152;			dc3 += 152;
				ss4   += 152;			ds3 += 152;
			}
		#else
			tdat = sc_ptr;	tmp = tdat + 126;
			two    = tmp + 0x0;
			one    = tmp + 0x1;
			tmp += 2;
			// Roots for radix-9 DFTs:	// Roots for radix-7 DFTs: MUST be followed by 4 alloc'ed padding slots needed by SSE2_RADIX_07_DFT macro!!! **************/
			cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
			ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
			cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
			ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
			cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
			ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
			cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
			ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
			/* These remain fixed within each per-thread local store: */
			VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
		  #if defined(USE_AVX2) || defined( USE_ARM_V8_SIMD)
			// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
			VEC_DBL_INIT(dc0, uc1  );	VEC_DBL_INIT(ds0, us1);
			VEC_DBL_INIT(dc1, uc2  );	VEC_DBL_INIT(ds1, us2);
			VEC_DBL_INIT(dc2, uc3  );	VEC_DBL_INIT(ds2, us3);
			VEC_DBL_INIT(dc3, 0.0  );	VEC_DBL_INIT(ds3, 0.0);	// Unused in non-Nussbaumer mode
		  #else
			// Roots for radix-7 DFTs: cc2 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation
			VEC_DBL_INIT(dc0, cx0-1);	VEC_DBL_INIT(ds0, sx0);
			VEC_DBL_INIT(dc1, cx1  );	VEC_DBL_INIT(ds1, sx1);
			VEC_DBL_INIT(dc2, cx2  );	VEC_DBL_INIT(ds2, sx2);
			VEC_DBL_INIT(dc3, cx3  );	VEC_DBL_INIT(ds3, sx3);
		  #endif
			// Roots for radix-9 DFTs:
			VEC_DBL_INIT(cc1  , c	);	VEC_DBL_INIT(ss1, s );
			VEC_DBL_INIT(cc2  , c2  );	VEC_DBL_INIT(ss2, s2);
			VEC_DBL_INIT(cc3m1, c3m1);	VEC_DBL_INIT(ss3, s3);
			VEC_DBL_INIT(cc4  , c4  );	VEC_DBL_INIT(ss4, s4);
		#endif
			return;
		} else {
			ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		tdat = __r0 + thr_id*152;	tmp = tdat + 126;
		two    = tmp + 0x0;
		one    = tmp + 0x1;
		tmp += 2;
		cc1    = tmp + 0x0;			dc0    = tmp + 0x8;
		ss1    = tmp + 0x1;			ds0    = tmp + 0x9;
		cc2    = tmp + 0x2;			dc1    = tmp + 0xa;
		ss2    = tmp + 0x3;			ds1    = tmp + 0xb;
		cc3m1  = tmp + 0x4;			dc2    = tmp + 0xc;
		ss3    = tmp + 0x5;			ds2    = tmp + 0xd;
		cc4    = tmp + 0x6;			dc3    = tmp + 0xe;
		ss4    = tmp + 0x7;			ds3    = tmp + 0xf;
	#else
		tdat = sc_ptr;
	#endif

		//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 7 radix-9 transforms:
		/*
		Twiddleless version arranges 7 sets of radix-9 DFT inputs as follows: 0 in upper left corner,
		decrement (mod 63) 9 horizontally and 7 vertically. Applying a further bit-reversal to that,
		Display result of Combined DIT input-scramble array in hex:

			00,02,01,08,07,06,05,04,03,
			32,31,30,2f,2e,2d,34,33,35,
			1d,1c,1b,22,21,23,1f,1e,20,
			10,0f,11,0d,0c,0e,0a,09,0b,
			37,36,38,3c,3e,3d,39,3b,3a,
			2a,2c,2b,27,29,28,24,26,25,
			15,17,16,12,14,13,1a,19,18.
		*/
		tmp = tdat; iptr = dit_iperm;
		for(l = 0; l < 7; l++) {
			// When 63 is used to build a larger DFT radix (e.g. 1008), these indices will be permuted (nonmonotone), no simplification possible:
			k0 = __idx[*iptr]; k1 = __idx[*(iptr+1)]; k2 = __idx[*(iptr+2)]; k3 = __idx[*(iptr+3)]; k4 = __idx[*(iptr+4)]; k5 = __idx[*(iptr+5)]; k6 = __idx[*(iptr+6)]; k7 = __idx[*(iptr+7)]; k8 = __idx[*(iptr+8)];
			va0 = __A+k0; va1 = __A+k1; va2 = __A+k2; va3 = __A+k3; va4 = __A+k4; va5 = __A+k5; va6 = __A+k6; va7 = __A+k7; va8 = __A+k8;
			// Since there is no vec_cmplx type, ptr-offs double those of analogous scalar code:
			vb0 = tmp; vb1 = tmp+2; vb2 = tmp+4; vb3 = tmp+6; vb4 = tmp+8; vb5 = tmp+10; vb6 = tmp+12; vb7 = tmp+14; vb8 = tmp+16;
			SSE2_RADIX_09_DIT(
				va0,va1,va2,va3,va4,va5,va6,va7,va8,
				cc1,
				vb0,vb1,vb2,vb3,vb4,vb5,vb6,vb7,vb8
			);	tmp += 18; iptr += 9;
		}
		/*...and now do 9 radix-7 transforms. The required output permutation is

			00,24,09,2d,12,36,1b,
			0e,32,17,3b,20,05,29,
			1c,01,25,0a,2e,13,37,
			2a,0f,33,18,3c,21,06,
			38,1d,02,26,0b,2f,14,
			07,2b,10,34,19,3d,22,
			15,39,1e,03,27,0c,30,
			23,08,2c,11,35,1a,3e,
			31,16,3a,1f,04,28,0d
		*/
		tmp = tdat; iptr = dit_operm;
		for(l = 0; l < 9; l++) {
			k0 = __odx[*iptr]; k1 = __odx[*(iptr+1)]; k2 = __odx[*(iptr+2)]; k3 = __odx[*(iptr+3)]; k4 = __odx[*(iptr+4)]; k5 = __odx[*(iptr+5)]; k6 = __odx[*(iptr+6)];
			// Since there is no vec_cmplx type, ptr-offs double those of analogous scalar code:
			va0 = tmp; va1 = tmp+18; va2 = tmp+36; va3 = tmp+54; va4 = tmp+72; va5 = tmp+90; va6 = tmp+108;
			vb0 = __B+k0; vb1 = __B+k1; vb2 = __B+k2; vb3 = __B+k3; vb4 = __B+k4; vb5 = __B+k5; vb6 = __B+k6;
		  #ifdef USE_AVX2
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, dc0,two, vb0,vb1,vb2,vb3,vb4,vb5,vb6)
		  #else
			SSE2_RADIX_07_DFT(va0,va1,va2,va3,va4,va5,va6, dc0,     vb0,vb1,vb2,vb3,vb4,vb5,vb6)
		  #endif
			tmp += 2; iptr += 7;
		}
	}

	/************** RADIX-64 DIF/DIT: *****************************/

	void SSE2_RADIX_64_DIF(
		const int init,	// Init consts (in 1-thread mode only!) and exit
		const int thr_id,
		// Need to know if the DIF-64 is standalone within a contig block of 64 vec_cmplx data
		// (as occurs if it's a standalone or part of a larger transform of length N = odd*64), or
		// part of a larger power-of-2 transform of length N = 2^k > 64, in which we need to adjust data strides:
		const int pow2_stride_shift,	// set = trailz(N) - trailz(64)
		// Inputs: Base address plus index offsets:
		double *__A, const int *i_offsets,
		// Intermediates-storage pointer:
		vec_dbl*r00,
		// Outputs: Base address plus index offsets:
		double *__B, const int *o_offsets
	)
	{
		// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
		// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
		const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
		double dtmp;
		const int k = pow2_stride_shift;
		int scale, k1,k2,k3,k4,k5,k6,k7;
		static int max_threads = 0;
		static vec_dbl *sc_arr = 0x0, *sc_ptr;	// Pad with 4 extra slots for scratch storage needed by SSE2_RADIX_07_DFT macro!!!
	  #ifdef MULTITHREAD
		static vec_dbl *__r0;	// Base address for discrete per-thread local stores - alloc 9x2e vec_dbl slots per thread
		// In || mode, only above base-pointer (shared by all threads) is static:
		vec_dbl *two,*one,*sqrt2,*isrt2, *cc0,*ss0,
			 *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
			*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7;
	  #else
		static vec_dbl *two,*one,*sqrt2,*isrt2, *cc0,*ss0,
			 *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
			*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7;
	  #endif
		// Intermediates pointers:
		vec_dbl *tmp, *v0,*v1,*v2,*v3,*v4,*v5,*v6,*v7;
		/* Addresses into array sections */
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j;
		const int *off_ptr;
	  #if COMPACT_OBJ	// This seems slower than unlooped, so mke latter the default
		vec_dbl *w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIF_TWIDDLE_OOP
		// ptr-offsets w.r.to base-ptr nisrt2:
		const uint8 twid_offsets[112] = {
			0x03,0x02, 0x01,0x01, 0x00,0x01, 0x17,0x18, 0x1b,0x17, 0x18,0x17, 0x1a,0x18,// ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
			0x01,0x01, 0x17,0x18, 0x18,0x17, 0x0b,0x0c, 0x24,0x23, 0x23,0x24, 0x0c,0x0b,// isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
			0x00,0x01, 0x18,0x17, 0x1a,0x1b, 0x23,0x24, 0x0e,0x0c ,0x0f,0x0b, 0x27,0x26,// nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
			0x17,0x18, 0x0b,0x0c, 0x23,0x24, 0x05,0x06, 0x1d,0x1e, 0x11,0x12, 0x29,0x2a,// cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
			0x1b,0x17, 0x24,0x23, 0x0e,0x0c, 0x1d,0x1e, 0x2c,0x2a, 0x06,0x05, 0x14,0x15,// nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
			0x18,0x17, 0x23,0x24, 0x0f,0x0b, 0x11,0x12, 0x06,0x05, 0x2a,0x29, 0x21,0x1d,// ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
			0x1a,0x18, 0x0c,0x0b, 0x27,0x26, 0x29,0x2a, 0x14,0x15, 0x21,0x1d, 0x06,0x08 // ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
		}, *bptr;
	  #endif

		// If this is first time here, init pointers and associated data:
		if(thr_id == -1)	// Value of init stores #threads
		{
			if(init <= max_threads) {	// Previously inited with sufficient #threads
				ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
				return;
			}
			max_threads = init;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
			if(sc_arr) { free((void *)sc_arr); }
			sc_arr = ALLOC_VEC_DBL(sc_arr, 0x32*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = ALIGN_VEC_DBL(sc_arr);
			ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		#ifdef MULTITHREAD
			__r0 = tmp = sc_ptr;
			// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
			two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
			one     = tmp + 1;
			sqrt2	= tmp + 2;
		//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
			tmp += 4;
			nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies the
			 isrt2	= tmp + 0x01;	// same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
			 cc0	= tmp + 0x02;
			 ss0	= tmp + 0x03;
			 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
			 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
			 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
			 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
			 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
			 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
			 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
			 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
			 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
			 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
			 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
			 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
			 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
			 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
			for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
				tmp = sqrt2+1;	// Unnamed copy of isrt2
				// Alternate "rounded the other way" copies of sqrt2,isrt2:
				dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
				dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
			//	VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(tmp, ISRT2);
				VEC_DBL_INIT(nisrt2,-dtmp);
				VEC_DBL_INIT( isrt2, dtmp);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
				VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);	//	tmp =  cc0-1; ASSERT(tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");	Disable to allow "round down" variant
				VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc2, c32_1);		VEC_DBL_INIT( ss2, s32_1);		tmp =  cc2-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc4, c16  );		VEC_DBL_INIT( ss4, s16  );		tmp =  cc4-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc6, c32_3);		VEC_DBL_INIT( ss6, s32_3);		tmp =  cc6-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc2,-c32_1);		VEC_DBL_INIT(nss2,-s32_1);		tmp = ncc2-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc4,-c16  );		VEC_DBL_INIT(nss4,-s16  );		tmp = ncc4-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc6,-c32_3);		VEC_DBL_INIT(nss6,-s32_3);		tmp = ncc6-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, dtmp);
			/* Move on to next thread's local store */
				 two += 0x32;
				 one += 0x32;
				 sqrt2 += 0x32;
				nisrt2 += 0x32;
				 isrt2 += 0x32;
				 cc0 += 0x32;
				 ss0 += 0x32;
				 cc1 += 0x32;		ncc1 += 0x32;
				 ss1 += 0x32;		nss1 += 0x32;
				 cc2 += 0x32;		ncc2 += 0x32;
				 ss2 += 0x32;		nss2 += 0x32;
				 cc3 += 0x32;		ncc3 += 0x32;
				 ss3 += 0x32;		nss3 += 0x32;
				 cc4 += 0x32;		ncc4 += 0x32;
				 ss4 += 0x32;		nss4 += 0x32;
				 cc5 += 0x32;		ncc5 += 0x32;
				 ss5 += 0x32;		nss5 += 0x32;
				 cc6 += 0x32;		ncc6 += 0x32;
				 ss6 += 0x32;		nss6 += 0x32;
				 cc7 += 0x32;		ncc7 += 0x32;
				 ss7 += 0x32;		nss7 += 0x32;
			}
		#else
			tmp = sc_ptr;
			// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
			two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
			one     = tmp + 1;
			sqrt2	= tmp + 2;
		//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
			tmp += 4;
			/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
			radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
			and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
			SSE2_RADIX8_DIF_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
			need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
			in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
			we get it at (vec_dbl*)cc-1.
			/Stupidity */
			nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies the
			 isrt2	= tmp + 0x01;	// same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
			 cc0	= tmp + 0x02;
			 ss0	= tmp + 0x03;
		// [copy isrt2]	= tmp + 0x04;	// [copy isrt2]	= tmp + 0x07;
			 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
			 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
		// [copy isrt2]	= tmp + 0x0a;	// [copy isrt2]	= tmp + 0x0d;
			 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
			 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
		// [copy isrt2]	= tmp + 0x10;	// [copy isrt2]	= tmp + 0x13;
			 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
			 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
		// [copy isrt2]	= tmp + 0x16;	// [copy isrt2]	= tmp + 0x19;
			 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
			 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
		// [copy isrt2]	= tmp + 0x1c;	// [copy isrt2]	= tmp + 0x1f;
			 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
			 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
		// [copy isrt2]	= tmp + 0x22;	// [copy isrt2]	= tmp + 0x25;
			 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
			 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
		// [copy isrt2]	= tmp + 0x28;	// [copy isrt2]	= tmp + 0x2b;
			 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
			 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
			/* These remain fixed: */
			VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
			tmp = sqrt2+1;	// Unnamed copy of isrt2
			// Alternate "rounded the other way" copies of sqrt2,isrt2:
			dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
			dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
		//	VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(tmp, ISRT2);
			VEC_DBL_INIT(nisrt2,-dtmp);
			VEC_DBL_INIT( isrt2, dtmp);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
			VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);	//	tmp =  cc0-1; ASSERT(tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");	Disable to allow "round down" variant
			VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc2, c32_1);		VEC_DBL_INIT( ss2, s32_1);		tmp =  cc2-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc4, c16  );		VEC_DBL_INIT( ss4, s16  );		tmp =  cc4-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc6, c32_3);		VEC_DBL_INIT( ss6, s32_3);		tmp =  cc6-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc2,-c32_1);		VEC_DBL_INIT(nss2,-s32_1);		tmp = ncc2-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc4,-c16  );		VEC_DBL_INIT(nss4,-s16  );		tmp = ncc4-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc6,-c32_3);		VEC_DBL_INIT(nss6,-s32_3);		tmp = ncc6-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, dtmp);
		#endif
		//	fprintf(stderr, "Init SSE2_RADIX_64_DIF with max_threads = %d\n",max_threads);
			return;
		} else {
			ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		tmp = __r0 + thr_id*0x32;
		// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
		two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one     = tmp + 1;
		sqrt2	= tmp + 2;
	//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
		tmp += 4;
		nisrt2	= tmp + 0x00;
		 isrt2	= tmp + 0x01;
		 cc0	= tmp + 0x02;
		 ss0	= tmp + 0x03;
		 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
		 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
		 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
		 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
		 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
		 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
		 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
		 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
		 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
		 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
		 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
		 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
		 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
		 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
	#endif

	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7;
		OFF1 = 0x100;
		OFF2 = 0x200;
		OFF3 = 0x300;
		OFF4 = 0x400;
		OFF5 = 0x500;
		OFF6 = 0x600;
		OFF7 = 0x700;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x100
		#define OFF2	4*0x200
		#define OFF3	4*0x300
		#define OFF4	4*0x400
		#define OFF5	4*0x500
		#define OFF6	4*0x600
		#define OFF7	4*0x700
	  #elif defined(USE_AVX)
		#define OFF1	2*0x100
		#define OFF2	2*0x200
		#define OFF3	2*0x300
		#define OFF4	2*0x400
		#define OFF5	2*0x500
		#define OFF6	2*0x600
		#define OFF7	2*0x700
	   #else
		#define OFF1	0x100
		#define OFF2	0x200
		#define OFF3	0x300
		#define OFF4	0x400
		#define OFF5	0x500
		#define OFF6	0x600
		#define OFF7	0x700
	   #endif
		if(i_offsets) {	/* This is geared toward radix = odd*64 DIF DFT: */
			scale = i_offsets[1];
		//	scale = i_offsets[1]<<1; *** should handle this in the i_offs *****	// Need 2x offset because of vec_dbl -> vec_cmplx implied cast
			// Apr 2014: Generalized-index scheme replaces original fixed __A-offsets OFF[1-7]*[pow2 scaling] with strides
			// taken from i_offsets array. Due to the way the radix-64 DFTs are used to build larger pow2 and non-pow2 DFTs
			// there is never an issue of irregular (i.e. not simple multiple of the basc stride in i_offsets[1]) strides,
			// just one of whether that basic 'unit' stride will amount to one vec_dbl pair or a larger stride. That means
			// we only need a very small sampling of the i_offsets data - in fact just i_offsets[1-8] - to infer the rest:
			k1 = i_offsets[0x8]<<L2_SZ_VD; k2 = k1+k1; k3 = (k1<<1)+k1;
			j = k1<<2;	k4 = j; k5 = k1+j; k6 = k2+j; k7 = k3+j;
		} else {		/* ...and this is geared toward radix = 2^k * 64: */
			scale = 1<<(1 + k);
			k1 = OFF1<<k; k2 = OFF2<<k; k3 = OFF3<<k; k4 = OFF4<<k; k5 = OFF5<<k; k6 = OFF6<<k; k7 = OFF7<<k;
		}
		// Outs are BRed:
		v0 = r00; v1 = v0+8; v2 = v0+4; v3 = v0+12; v4 = v0+2; v5 = v0+10; v6 = v0+6; v7 = v0+14;
		for(i = 0; i < 8; i++) {
			j = reverse8[i];	// i=1/j=4 should land us in the middle of the first 1/8-chunk of the contiguous in-array __A
							// For k = 0 __A has (64 vec_cmplex) = (128 vec_dbl) elements, thus i=1/j=4 => j*scale = 128/8 = 16.
			tmp = (vec_dbl*)__A+j*scale;	// __A-offsets are processed in BR8 order
		  #ifdef USE_AVX2
			SSE2_RADIX8_DIF_0TWIDDLE_B(	// Use B-version of macro, which takes the i-strides as intvars rather than literal bytes:
				tmp,k1,k2,k3,k4,k5,k6,k7,
				v0,v1,v2,v3,v4,v5,v6,v7, isrt2,two
			);
		  #else
			SSE2_RADIX8_DIF_0TWIDDLE_B(	// Use B-version of macro, which takes the i-strides as intvars rather than literal bytes:
				tmp,k1,k2,k3,k4,k5,k6,k7,
				v0,v1,v2,v3,v4,v5,v6,v7, isrt2
			);
		  #endif
			v0 += 16; v1 += 16; v2 += 16; v3 += 16; v4 += 16; v5 += 16; v6 += 16; v7 += 16;
		}

	/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

	/* Block 0: */
		tmp = r00;
		off_ptr = o_offsets; add0 = __B+off_ptr[0];add1 = __B+off_ptr[1];add2 = __B+off_ptr[2];add3 = __B+off_ptr[3];add4 = __B+off_ptr[4];add5 = __B+off_ptr[5];add6 = __B+off_ptr[6];add7 = __B+off_ptr[7];
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
		so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIF_0TWIDDLE(
			r00, OFF4,OFF2,OFF6,OFF1,OFF5,OFF3,OFF7,
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two
		);
	  #else
		SSE2_RADIX8_DIF_0TWIDDLE(
			r00, OFF4,OFF2,OFF6,OFF1,OFF5,OFF3,OFF7,
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2
		);
	  #endif
	  #if COMPACT_OBJ

		bptr = twid_offsets;	// Byte-pointer loops over slots in twid_offsets[] byte-array
		for(i = 1, j = 0; i < 8; i++, j += 14) {
			tmp = r00 + (reverse8[i] << 1);
			off_ptr += 8;
			w[0x0] = nisrt2 + *bptr++; w[0x1] = nisrt2 + *bptr++; w[0x2] = nisrt2 + *bptr++; w[0x3] = nisrt2 + *bptr++; w[0x4] = nisrt2 + *bptr++; w[0x5] = nisrt2 + *bptr++; w[0x6] = nisrt2 + *bptr++; w[0x7] = nisrt2 + *bptr++; w[0x8] = nisrt2 + *bptr++; w[0x9] = nisrt2 + *bptr++; w[0xa] = nisrt2 + *bptr++; w[0xb] = nisrt2 + *bptr++; w[0xc] = nisrt2 + *bptr++; w[0xd] = nisrt2 + *bptr++;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				tmp,OFF1,
				__B,off_ptr,	// Do the output-address pointer arithmetic inside the macro
				twid_ptrs, two
			);
		}

	  #else	// COMPACT_OBJ = False:
			// Each of these rhi-ptrs points to next 8 vec_cmplx= 16 vec_dbl data sets:
		vec_dbl *r10 = r00+0x10,*r20 = r00+0x20,*r30 = r00+0x30,*r40 = r00+0x40,*r50 = r00+0x50,*r60 = r00+0x60,*r70 = r00+0x70;
	// Remaining 7 macro calls need explicitly named BRed iptrs with stride 16, use tmp(in place of iarg r00),r40,r20,r60,r10,r50,r30,r70:
		#error SSE2_RADIX8_DIF_TWIDDLE_OOP in this section need twiddles collected into w[]-array as above!
	/* Block 4: */
		tmp += 8;	// r08
		off_ptr += 8;
	//...and another (needed by FMA version of macro)) kludge for the 30-arg limit: put copies of (vec_dbl)2.0,SQRT2 into the first 2 of each set of outputs.
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4, two
		);
	/* Block 2: */
		tmp -= 4;	// r04
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2, two
		);
	/* Block 6: */
		tmp += 8;	// r0C
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6, two
		);
	/* Block 1: */
		tmp -= 10;	// r02
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7, two
		);
	/* Block 5: */
		tmp += 8;	// r0A
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3, two
		);
	/* Block 3: */
		tmp -= 4;	// r06
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5, two
		);
	/* Block 7: */
		tmp += 8;	// r0E
		off_ptr += 8;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
			tmp,OFF1,
			__B,off_ptr,
			ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1, two
		);

	  #endif

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
		#undef OFF5
		#undef OFF6
		#undef OFF7
	  #endif
	}

	void SSE2_RADIX_64_DIT(
		const int init,	// Init consts (in 1-thread mode only!) and exit
		const int thr_id,
		// Inputs: Base address plus index offsets:
		double *__A, const int *i_offsets,
		// Intermediates-storage pointer:
		vec_dbl*r00,
		// Outputs: Base address plus index offsets:
		vec_dbl*__B, const int *o_offsets
	)
	{
		const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
		double dtmp;
		static int max_threads = 0;
		static vec_dbl *sc_arr = 0x0, *sc_ptr;	// Pad with 4 extra slots for scratch storage needed by SSE2_RADIX_07_DFT macro!!!
	  #ifdef MULTITHREAD
		static vec_dbl *__r0;	// Base address for discrete per-thread local stores - alloc 9x2e vec_dbl slots per thread
		// In || mode, only above base-pointer (shared by all threads) is static:
		vec_dbl *two,*one,*sqrt2,*isrt2, *cc0,*ss0,
			 *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
			*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7;
	  #else
		static vec_dbl *two,*one,*sqrt2,*isrt2, *cc0,*ss0,
			 *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
			*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7;
	  #endif
		// Intermediates pointers:
		vec_dbl *tmp, *v0,*v1,*v2,*v3,*v4,*v5,*v6,*v7,
			// Each of these rhi-ptrs points to next 8 vec_cmplx= 16 vec_dbl data sets:
			*r10 = r00+0x10,*r20 = r00+0x20,*r30 = r00+0x30,*r40 = r00+0x40,*r50 = r00+0x50,*r60 = r00+0x60,*r70 = r00+0x70,
			*w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIF_TWIDDLE_OOP
		/* Addresses into array sections */
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j;
		const int *off_ptr;

		// If this is first time here, init pointers and associated data:
		if(thr_id == -1)	// Value of init stores #threads
		{
			if(init <= max_threads) {	// Previously inited with sufficient #threads
				ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
				return;
			}
			max_threads = init;
		#ifndef COMPILER_TYPE_GCC
			ASSERT(NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
		#endif
			if(sc_arr) { free((void *)sc_arr); }
			sc_arr = ALLOC_VEC_DBL(sc_arr, 0x32*max_threads);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			sc_ptr = ALIGN_VEC_DBL(sc_arr);
			ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		#ifdef MULTITHREAD
			__r0 = tmp = sc_ptr;
			// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
			two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
			one     = tmp + 1;
			sqrt2	= tmp + 2;
		//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
			tmp += 4;
			nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies the
			 isrt2	= tmp + 0x01;	// same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
			 cc0	= tmp + 0x02;
			 ss0	= tmp + 0x03;
			 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
			 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
			 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
			 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
			 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
			 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
			 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
			 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
			 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
			 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
			 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
			 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
			 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
			 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
			for(i = 0; i < max_threads; ++i) {
			/* These remain fixed within each per-thread local store: */
				VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
				tmp = sqrt2+1;	// Unnamed copy of isrt2
				// Alternate "rounded the other way" copies of sqrt2,isrt2:
				dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
				dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
			//	VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(tmp, ISRT2);
				VEC_DBL_INIT(nisrt2,-dtmp);
				VEC_DBL_INIT( isrt2, dtmp);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
				VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);	//	tmp =  cc0-1; ASSERT(tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");	Disable to allow "round down" variant
				VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc2, c32_1);		VEC_DBL_INIT( ss2, s32_1);		tmp =  cc2-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc4, c16  );		VEC_DBL_INIT( ss4, s16  );		tmp =  cc4-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc6, c32_3);		VEC_DBL_INIT( ss6, s32_3);		tmp =  cc6-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc2,-c32_1);		VEC_DBL_INIT(nss2,-s32_1);		tmp = ncc2-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc4,-c16  );		VEC_DBL_INIT(nss4,-s16  );		tmp = ncc4-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc6,-c32_3);		VEC_DBL_INIT(nss6,-s32_3);		tmp = ncc6-1; VEC_DBL_INIT(tmp, dtmp);
				VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, dtmp);
			/* Move on to next thread's local store */
				 two += 0x32;
				 one += 0x32;
				 sqrt2 += 0x32;
				nisrt2 += 0x32;
				 isrt2 += 0x32;
				 cc0 += 0x32;
				 ss0 += 0x32;
				 cc1 += 0x32;		ncc1 += 0x32;
				 ss1 += 0x32;		nss1 += 0x32;
				 cc2 += 0x32;		ncc2 += 0x32;
				 ss2 += 0x32;		nss2 += 0x32;
				 cc3 += 0x32;		ncc3 += 0x32;
				 ss3 += 0x32;		nss3 += 0x32;
				 cc4 += 0x32;		ncc4 += 0x32;
				 ss4 += 0x32;		nss4 += 0x32;
				 cc5 += 0x32;		ncc5 += 0x32;
				 ss5 += 0x32;		nss5 += 0x32;
				 cc6 += 0x32;		ncc6 += 0x32;
				 ss6 += 0x32;		nss6 += 0x32;
				 cc7 += 0x32;		ncc7 += 0x32;
				 ss7 += 0x32;		nss7 += 0x32;
			}
		#else
			tmp = sc_ptr;
			// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
			two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
			one     = tmp + 1;
			sqrt2	= tmp + 2;
		//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
			tmp += 4;
			/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
			radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
			and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
			SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
			need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
			in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
			we get it at (vec_dbl*)cc-1.
			/Stupidity */
			nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies the
			 isrt2	= tmp + 0x01;	// same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
			 cc0	= tmp + 0x02;
			 ss0	= tmp + 0x03;
		// [copy isrt2]	= tmp + 0x04;	// [copy isrt2]	= tmp + 0x07;
			 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
			 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
		// [copy isrt2]	= tmp + 0x0a;	// [copy isrt2]	= tmp + 0x0d;
			 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
			 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
		// [copy isrt2]	= tmp + 0x10;	// [copy isrt2]	= tmp + 0x13;
			 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
			 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
		// [copy isrt2]	= tmp + 0x16;	// [copy isrt2]	= tmp + 0x19;
			 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
			 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
		// [copy isrt2]	= tmp + 0x1c;	// [copy isrt2]	= tmp + 0x1f;
			 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
			 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
		// [copy isrt2]	= tmp + 0x22;	// [copy isrt2]	= tmp + 0x25;
			 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
			 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
		// [copy isrt2]	= tmp + 0x28;	// [copy isrt2]	= tmp + 0x2b;
			 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
			 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
			/* These remain fixed: */
			VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	tmp = sqrt2+1;
			tmp = sqrt2+1;	// Unnamed copy of isrt2
			// Alternate "rounded the other way" copies of sqrt2,isrt2:
			dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
			dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
		//	VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(tmp, ISRT2);
			VEC_DBL_INIT(nisrt2,-dtmp);
			VEC_DBL_INIT( isrt2, dtmp);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
			VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);	//	tmp =  cc0-1; ASSERT(tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");	Disable to allow "round down" variant
			VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc2, c32_1);		VEC_DBL_INIT( ss2, s32_1);		tmp =  cc2-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc4, c16  );		VEC_DBL_INIT( ss4, s16  );		tmp =  cc4-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc6, c32_3);		VEC_DBL_INIT( ss6, s32_3);		tmp =  cc6-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc2,-c32_1);		VEC_DBL_INIT(nss2,-s32_1);		tmp = ncc2-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc4,-c16  );		VEC_DBL_INIT(nss4,-s16  );		tmp = ncc4-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc6,-c32_3);		VEC_DBL_INIT(nss6,-s32_3);		tmp = ncc6-1; VEC_DBL_INIT(tmp, dtmp);
			VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, dtmp);
		#endif
		//	fprintf(stderr, "Init SSE2_RADIX_64_DIT with max_threads = %d\n",max_threads);
			return;
		} else {
			ASSERT(sc_arr != 0, "This function requires an initial Init-consts-mode call (in 1-thread mode only) before use!");
		}	/* end of inits */

		/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT((uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		tmp = __r0 + thr_id*0x32;
		// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
		two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one     = tmp + 1;
		sqrt2	= tmp + 2;
	//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
		tmp += 4;
		nisrt2	= tmp + 0x00;
		 isrt2	= tmp + 0x01;
		 cc0	= tmp + 0x02;
		 ss0	= tmp + 0x03;
		 cc1	= tmp + 0x05;		ncc1	= tmp + 0x08;
		 ss1	= tmp + 0x06;		nss1	= tmp + 0x09;
		 cc2	= tmp + 0x0b;		ncc2	= tmp + 0x0e;
		 ss2	= tmp + 0x0c;		nss2	= tmp + 0x0f;
		 cc3	= tmp + 0x11;		ncc3	= tmp + 0x14;
		 ss3	= tmp + 0x12;		nss3	= tmp + 0x15;
		 cc4	= tmp + 0x17;		ncc4	= tmp + 0x1a;
		 ss4	= tmp + 0x18;		nss4	= tmp + 0x1b;
		 cc5	= tmp + 0x1d;		ncc5	= tmp + 0x20;
		 ss5	= tmp + 0x1e;		nss5	= tmp + 0x21;
		 cc6	= tmp + 0x23;		ncc6	= tmp + 0x26;
		 ss6	= tmp + 0x24;		nss6	= tmp + 0x27;
		 cc7	= tmp + 0x29;		ncc7	= tmp + 0x2c;
		 ss7	= tmp + 0x2a;		nss7	= tmp + 0x2d;
	#endif

	/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
	// Because of roots-sign diffs between SSE2 and scalar macros, 1/7 2/6 3/5 swapped in DIT_0TWIDDLE outputs!

		off_ptr = i_offsets;
		tmp = r00;
		for(i = 0; i < 8; i++) {
			add0 = __A+off_ptr[0];add1 = __A+off_ptr[1];add2 = __A+off_ptr[2];add3 = __A+off_ptr[3];add4 = __A+off_ptr[4];add5 = __A+off_ptr[5];add6 = __A+off_ptr[6];add7 = __A+off_ptr[7];
		#ifdef USE_AVX2
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2,two)
		#else
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, tmp, isrt2)
		#endif
			tmp += 16;
			off_ptr += 8;
		}

	/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details.
	Use the same positive-power roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies.
	Twiddles occur in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,4,2,6,1,5,3,7].
	*/
		// s1p[00,08,10,18,20,28,30,38]:
		off_ptr = o_offsets;
	/*
	printf("off_ptr[0x8] = %2x\n",off_ptr[0x8]);	<*** 2x larger than the 0x10-mults in the 'else' below *>**
	exit(0);
	*/
		// Apr 2014: Generalized-index scheme replaces original fixed __B-offsets 0x[0-7]0 with strides taken from
		// o_offsets array. Due to the way the radix-64 DFTs are used to build larger pow2 and non-pow2 DFTs there
		// is never an issue of irregular (i.e. not simple multiple of the basic stride in o_offsets[1]) strides,
		// just one of whether that basic 'unit' stride will amount to one vec_dbl pair or a larger stride. That means
		// we only need a very small sampling of the o_offsets data - in fact just o_offsets[1-8] - to infer the rest:
		j = off_ptr[8];	v0 = __B; v1 = __B+j; v2 = __B+(j<<1); v3 = __B+j+(j<<1);
		uint64 dit_o_ptr_stride = (uint64)v1-(uint64)v0;
		j <<= 2;		v4 = v0+j; v5 = v1+j; v6 = v2+j; v7 = v3+j;
	// Block 0: All unity twiddles:
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			r00,r10,r20,r30,r40,r50,r60,r70,
			v0,v7,v6,v5,v4,v3,v2,v1, isrt2,two
		);
	  #else
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			r00,r10,r20,r30,r40,r50,r60,r70,
			v0,v7,v6,v5,v4,v3,v2,v1, isrt2
		);
	  #endif
		tmp = r00;

	// Note: Assumed 1-before 1st of the 14 sincos args in each call to SSE2_RADIX8_DIT_TWIDDLE_OOP is the basic isrt2 arg
	// used for radix-8. This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
	//...and another kludge for the 30-arg limit: put a copy of (vec_dbl)2.0 into the first of each set of outputs:
	   #ifdef USE_ARM_V8_SIMD
		uint32 OFF1 = 0x100;
	  #elif defined(USE_AVX512)
		#define OFF1	0x400
	  #elif defined(USE_AVX)
		#define OFF1	0x200
	  #else
		#define OFF1	0x100
	   #endif
	// Block 4: jt = j1 + p04;	jp = j2 + p04;
		// Note our r** pointers have indices that run 2x faster than the s1p** ones, e.g. r08-r00 == s1p04-s1p00:
		j = +8;	tmp += j;	// tmp = r08
		j = off_ptr[4]; v0 += j;	// s1p04
		i = 0; w[i++]=ss0; w[i++]=cc0; w[i++]= isrt2; w[i++]=isrt2; w[i++]= nisrt2; w[i++]=isrt2; w[i++]= cc4; w[i++]=ss4; w[i++]= nss4; w[i++]=cc4; w[i++]= ss4; w[i++]=cc4; w[i++]= ncc4; w[i++]=ss4;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 2: jt = j1 + p02;	jp = j2 + p02;
		j = +4;	tmp += j;	// tmp = r0C
		j = off_ptr[2]; v0 -= j;	// s1p02
		i = 0; w[i++]=isrt2; w[i++]=isrt2; w[i++]=cc4; w[i++]=ss4; w[i++]=ss4; w[i++]=cc4; w[i++]=cc2; w[i++]=ss2; w[i++]=ss6; w[i++]=cc6; w[i++]=cc6; w[i++]=ss6; w[i++]=ss2; w[i++]=cc2;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 6: jt = j1 + p06;	jp = j2 + p06;
		j = -8;	tmp += j;	// tmp = r04
		j = off_ptr[4]; v0 += j;	// s1p06
		i = 0; w[i++]=nisrt2; w[i++]=isrt2; w[i++]=ss4; w[i++]=cc4; w[i++]=ncc4; w[i++]=nss4; w[i++]=cc6; w[i++]=ss6; w[i++]=ncc2; w[i++]=ss2; w[i++]=nss2; w[i++]=cc2; w[i++]=nss6; w[i++]=ncc6;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 1: jt = j1 + p01;	jp = j2 + p01;
		j =+10;	tmp += j;	// tmp = r0E
		j = off_ptr[5]; v0 -= j;	// s1p01
		i = 0; w[i++]=cc4; w[i++]=ss4; w[i++]=cc2; w[i++]=ss2; w[i++]=cc6; w[i++]=ss6; w[i++]=cc1; w[i++]=ss1; w[i++]=cc5; w[i++]=ss5; w[i++]=cc3; w[i++]=ss3; w[i++]=cc7; w[i++]=ss7;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 5: jt = j1 + p05;	jp = j2 + p05;
		j = -8;	tmp += j;	// tmp = r06
		j = off_ptr[4]; v0 += j;	// s1p05
		i = 0; w[i++]=nss4; w[i++]=cc4; w[i++]=ss6; w[i++]=cc6; w[i++]=ncc2; w[i++]=ss2; w[i++]=cc5; w[i++]=ss5; w[i++]=ncc7; w[i++]=ss7; w[i++]=ss1; w[i++]=cc1; w[i++]=ncc3; w[i++]=nss3;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 3: jt = j1 + p03;	jp = j2 + p03;
		j = +4;	tmp += j;	// tmp = r0A
		j = off_ptr[2]; v0 -= j;	// s1p03
		i = 0; w[i++]=ss4; w[i++]=cc4; w[i++]=cc6; w[i++]=ss6; w[i++]=nss2; w[i++]=cc2; w[i++]=cc3; w[i++]=ss3; w[i++]=ss1; w[i++]=cc1; w[i++]=ss7; w[i++]=cc7; w[i++]=nss5; w[i++]=cc5;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	// Block 7: jt = j1 + p07;	jp = j2 + p07;
		j = -8;	tmp += j;	// tmp = r02
		j = off_ptr[4]; v0 += j;	// s1p07
		i = 0; w[i++]=ncc4; w[i++]=ss4; w[i++]=ss2; w[i++]=cc2; w[i++]=nss6; w[i++]=ncc6; w[i++]=cc7; w[i++]=ss7; w[i++]=ncc3; w[i++]=nss3; w[i++]=nss5; w[i++]=cc5; w[i++]=ss1; w[i++]=ncc1;
		SSE2_RADIX8_DIT_TWIDDLE_OOP(
			tmp,OFF1,
			v0,dit_o_ptr_stride,
			twid_ptrs, two
		);
	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
	  #endif
	}

	/************** RADIX-128 DIF/DIT: *****************************/

	void SSE2_RADIX128_DIF(
		// Input pointer: Base ptr of 16 local-mem:
		vec_dbl*__A,
		// Intermediates-storage pointer:
		vec_dbl*tmp,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
		// Outputs: Base address plus 30 index offsets:
		double *__B,
		int *o_offsets	// Array storing output index offsets
	)
	{
		// Intermediates pointers:
		vec_dbl *tm0,*tm1,
			*w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIF_TWIDDLE_OOP
		/* Addresses into array sections */
		double *addr,*add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7;
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j, *off_ptr;
		int p0,p1,p2,p3,p4,p5,p6,p7;

	// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:
	// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:

	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4,OFF5,OFF6,OFF7;
		OFF1 = 0x100;
		OFF2 = 0x200;
		OFF3 = 0x300;
		OFF4 = 0x400;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x100
		#define OFF2	4*0x200
		#define OFF3	4*0x300
		#define OFF4	4*0x400
	  #elif defined(USE_AVX)
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF3	0x600
		#define OFF4	0x800
	  #else
		#define OFF1	0x100
		#define OFF2	0x200
		#define OFF3	0x300
		#define OFF4	0x400
	  #endif
		// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:
		// Use that [isrt2] ptr heads up the needed isrt2/cc16/ss16 pointer-triplet
		// (actually a quartet in AVX2 mode, since there need sqrt2 preceding isrt2)
		tm1 = tmp;	// Outputs of 16-DFTs into local scratch storage
	//	printf("B: Pass 1 inputs:\n");
		for(i = 0; i < 8; i++) {
			j = reverse8[i]<<1;	// __A-offsets are processed in BR8 order
			tm0 = __A+j;
	//	for(j = 0; j < 16; j += 2) { printf("%2x: %15.5f,%15.5f, %15.5f,%15.5f\n",i*16+j,(tm0+j)->d0,(tm0+j)->d1,(tm0+j+1)->d0,(tm0+j+1)->d1); };
			SSE2_RADIX16_DIF_0TWIDDLE_B(tm0,OFF1,OFF2,OFF3,OFF4, isrt2,two, tm1);
			tm1 += 32;
		}
	//	exit(0);
	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors: */
	//	printf("B: Pass 2 inputs:\n");
	//	for(i = 0; i < 256; i += 2) { printf("%2x: %15.5f,%15.5f, %15.5f,%15.5f\n",i,(tmp+i)->d0,(tmp+i)->d1,(tmp+i+1)->d0,(tmp+i+1)->d1); };
	//	exit(0);
	   #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
		OFF2 = 0x400;
		OFF3 = 0x600;
		OFF4 = 0x800;
		OFF5 = 0xa00;
		OFF6 = 0xc00;
		OFF7 = 0xe00;
	   #endif
		/* Block 0: */
		off_ptr = o_offsets;
		p0 = off_ptr[0x0];p1 = off_ptr[0x1];p2 = off_ptr[0x2];p3 = off_ptr[0x3];p4 = off_ptr[0x4];p5 = off_ptr[0x5];p6 = off_ptr[0x6];p7 = off_ptr[0x7];
		addr = __B; add0 = addr+p0; add1 = addr+p1; add2 = addr+p2; add3 = addr+p3; add4 = addr+p4; add5 = addr+p5; add6 = addr+p6; add7 = addr+p7;
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
		so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
		SSE2_RADIX8_DIF_0TWIDDLE(
		#ifdef USE_ARM_V8_SIMD
			tmp, OFF4,OFF2,OFF6,OFF1,OFF5,OFF3,OFF7,	//	tmp,r40,r20,r60,r10,r50,r30,r70
		#elif defined(USE_AVX512)
			tmp,0x2000,0x1000,0x3000,0x0800,0x2800,0x1800,0x3800,
		#elif defined(USE_AVX)
			tmp,0x1000,0x0800,0x1800,0x0400,0x1400,0x0c00,0x1c00,
		#else	// USE_SSE2
			tmp, 0x800, 0x400, 0xc00, 0x200, 0xa00, 0x600, 0xe00,
		#endif
		#ifdef USE_AVX2
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2,two
		#else
			add0,add1,add2,add3,add4,add5,add6,add7, isrt2
		#endif
		);

	   #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
	  #elif defined(USE_AVX512)
		#define OFF1	0x800
	  #elif defined(USE_AVX)
		#define OFF1	0x400
	  #else
		#define OFF1	0x200
	   #endif
		// Remaining 15 sets of macro calls done in loop:
		for(i = 1; i < 16; i++) {
			off_ptr += 8;
			j = reverse16[i];
			// I-ptrs and sincos-ptrs; Twid-offsets are multiples of 21 vec_dbl; +1 to point to cos [middle] term of each [isrt2,c,s] triplet
			tm1 = tmp + (i<<1); tm0 = twid0 + (j<<4)+(j<<2)+j;
			w[0x0] = tm0+ 1; w[0x1] = tm0+ 2;
			w[0x2] = tm0+ 4; w[0x3] = tm0+ 5;
			w[0x4] = tm0+ 7; w[0x5] = tm0+ 8;
			w[0x6] = tm0+10; w[0x7] = tm0+11;
			w[0x8] = tm0+13; w[0x9] = tm0+14;
			w[0xa] = tm0+16; w[0xb] = tm0+17;
			w[0xc] = tm0+19; w[0xd] = tm0+20;
			// O-ptrs: a[j1] offset here = p08,p10,p18,...
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				tm1,OFF1,
				__B,off_ptr,	// Do the output-address pointer arithmetic inside the macro
				twid_ptrs, two
			);
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
	  #endif
	}

	void SSE2_RADIX128_DIT(
		// Inputs: Base address plus index offsets:
		double *__A,
		int *i_offsets,	// Array storing input index offsets
		// Intermediates-storage pointer:
		vec_dbl*tmp,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
		// Output pointer: Base ptr of 16 local-mem:
		vec_dbl*__B
	)
	{
		// Intermediates pointers:
		vec_dbl *tm0,*tm1,*tm2,
			*u0,*u1,*u2,*u3,*u4,*u5,*u6,*u7, *v0,*v1,*v2,*v3,*v4,*v5,*v6,*v7,
			*w[14], **twid_ptrs = &(w[0]);	// Array of 14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DIT_TWIDDLE_OOP
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j, *off_ptr;

	/* Gather the needed data (128 64-bit complex, i.e. 256 64-bit reals) and do 8 twiddleless length-16 subtransforms: */

	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x20;
		OFF2 = 0x40;
		OFF3 = 0x60;
		OFF4 = 0x80;
	  #elif defined(USE_AVX512)
		#define OFF1	4*0x20
		#define OFF2	4*0x40
		#define OFF3	4*0x60
		#define OFF4	4*0x80
	  #elif defined(USE_AVX)
		#define OFF1	2*0x20
		#define OFF2	2*0x40
		#define OFF3	2*0x60
		#define OFF4	2*0x80
	  #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
	  #endif

		// Start of isrt2/cc16/ss16 triplet needed for radix-16 SSE2 DFT macros [use instead of isrt2 here]
		// Use that [isrt2] ptr heads up the needed isrt2/cc16/ss16 pointer-triplet
		// (actually a quartet in AVX2 mode, since there need sqrt2 preceding isrt2)
		tm1 = tmp;	// Outputs of 16-DFTs into local scratch storage
		off_ptr = i_offsets;
		for(i = 0; i < 8; i++) {
			// I-ptrs: a[j1] offset here = p08,p10,p18,...
			SSE2_RADIX16_DIT_0TWIDDLE(
				__A,off_ptr,
				isrt2,two,
				tm1,OFF1,OFF2,OFF3,OFF4
			);
			tm1 += 32;
			off_ptr += 16;
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	/*...and now do 16 radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

	/* Block 0: */
		tm0 = tmp; tm1 = __B;
		u0 = tm0; u1 = tm0 + 0x20; u2 = tm0 + 0x40; u3 = tm0 + 0x60; u4 = tm0 + 0x80; u5 = tm0 + 0xa0; u6 = tm0 + 0xc0; u7 = tm0 + 0xe0;
		v0 = tm1; v1 = tm1 + 0xe0; v2 = tm1 + 0xc0; v3 = tm1 + 0xa0; v4 = tm1 + 0x80; v5 = tm1 + 0x60; v6 = tm1 + 0x40; v7 = tm1 + 0x20;
	  #ifdef USE_AVX2
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			u0,u1,u2,u3,u4,u5,u6,u7,
			v0,v1,v2,v3,v4,v5,v6,v7,
			isrt2,two
		);
	  #else
		SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
			u0,u1,u2,u3,u4,u5,u6,u7,
			v0,v1,v2,v3,v4,v5,v6,v7,
			isrt2
		);
	  #endif

	// Note: 1st of the 14 sincos args in each call to SSE2_RADIX8_DIT_TWIDDLE_OOP is the basic isrt2 needed for
	// radix-8. This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.

	   #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
	  #elif defined(USE_AVX512)
		#define OFF1	0x800
	  #elif defined(USE_AVX)
		#define OFF1	0x400
	  #else
		#define OFF1	0x200
	   #endif
		uint64 dit_o_ptr_stride = (uint64)OFF1;	// Dec 2020: Do we need to generalize this as in Radix-64 DIT, as called e.g. by radix-960 carry routine?
		// Remaining 15 sets of macro calls done in loop:
		for(i = 1; i < 16; i++) {
			j = reverse16[i];
			// Twid-offsets are multiples of 21 vec_dbl; +1 to point to cos [middle] term of each [isrt2,c,s] triplet
			/*** Unlike 128-DIF, need added ptr tm2 here to avid overwriting tmp-intermdiates ptr: ***/
			tm0 = tmp + (j<<1); tm1 = __B + (j<<1); tm2 = twid0 + (j<<4)+(j<<2)+j;
			// 7 [c,s] pointer-pairs, +=3 between pairs; c0 = tm2+1 to point to cos [middle] term of each [isrt2,c,s] triplet:
			w[0x0] = tm2+ 1; w[0x1] = tm2+ 2;
			w[0x2] = tm2+ 4; w[0x3] = tm2+ 5;
			w[0x4] = tm2+ 7; w[0x5] = tm2+ 8;
			w[0x6] = tm2+10; w[0x7] = tm2+11;
			w[0x8] = tm2+13; w[0x9] = tm2+14;
			w[0xa] = tm2+16; w[0xb] = tm2+17;
			w[0xc] = tm2+19; w[0xd] = tm2+20;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				tm0,OFF1,
				tm1,dit_o_ptr_stride,
				twid_ptrs, two
			);
		}
	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
	  #endif
	}

	/************** RADIX-256 DIF/DIT: *****************************/

	void SSE2_RADIX256_DIF(
		// Input pointer: Base ptr of 16 local-mem:
		vec_dbl*__A,
		// Intermediates-storage pointer:
		vec_dbl*r00,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
		// Outputs: Base address plus 30 index offsets:
		double *__B,
		int *o_offsets_lo,	// Array storing  low parts of output index offsets in 16 slots
		uint32 o_idx,	// Bitfield encoding the sequence of the o_offsets_lo sub-vectors to use for the radix-256 DFT's outputs
		int *o_offsets_hi	// Array storing high parts of output index offsets in 16 slots
	)
	{
		// Intermediates pointers:
		vec_dbl *tm0,*tm1,*tm2;
		/* Addresses into array sections */
		double *addr;
	  #ifdef USE_AVX2
		double *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7,*add8,*add9,*adda,*addb,*addc,*addd,*adde,*addf;
	  #endif
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j,nshift, *off_ptr;
	  #ifdef USE_AVX2
		int p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf;
	  #endif
		ASSERT(o_idx != 0x0, "Null o_idx pointer in SSE2_RADIX256_DIF!");
	// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x200;
		OFF2 = 0x400;
		OFF3 = 0x600;
		OFF4 = 0x800;
	  #elif defined(USE_AVX512)
		#define OFF1	0x200*4
		#define OFF2	0x400*4
		#define OFF3	0x600*4
		#define OFF4	0x800*4
	  #elif defined(USE_AVX)
		#define OFF1	0x200*2
		#define OFF2	0x400*2
		#define OFF3	0x600*2
		#define OFF4	0x800*2
	  #else
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF3	0x600
		#define OFF4	0x800
	  #endif
		tm1 = r00;
		for(i = 0; i < 16; i++) {
			j = reverse16[i]<<1;	// __A-offsets are processed in BR16 order
			tm0 = __A+j;
			SSE2_RADIX16_DIF_0TWIDDLE_B(tm0,OFF1,OFF2,OFF3,OFF4, isrt2,two, tm1);
			tm1 += 32;
		}

	/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors: */

	// Block 0: has all-unity twiddles
		// Extract index of the 16-element o_offsets_lo sub-vector to use for the current set of outputs -
		// the 16 p-offsets p0-pf will be in some permuted order starting at off_ptr:
		off_ptr = o_offsets_lo + ( (o_idx&0x3) << 4 );	// Low 2 bits of o_idx; loop below will use remaining 30 bits in ascending pairs
		tm1 = r00;
		addr = __B + o_offsets_hi[0];
		SSE2_RADIX16_DIF_TWIDDLE_OOP(
			tm1,OFF1,OFF4,
			addr,off_ptr,
			isrt2,twid0
		);	tm1 += 2;

	  #ifdef USE_AVX2

		// Due to tangent-twiddles scheme and resulting singularity of tangent(arg(I)) = 1/0,
		// only last 14 of the 15 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3:
	// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
		off_ptr = o_offsets_lo + ((o_idx<<2)&0x30);	// Shorthand for 16*((o_idx>>2)&0x3)
		j = 16;	// = 8<<1; Mimics elided i = 1 pass of length-14 loop below
		tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl - this one points to twid8 [of twid0-f]
		addr = __B + o_offsets_hi[1];
		SSE2_RADIX16_DIF_TWIDDLE_OOP(
			tm1,OFF1,OFF4,
			addr,off_ptr,
			isrt2,tm2
		);	tm1 += 2;

		// Remaining 14 sets of macro calls done in loop:
		for(i = 2; i < 16; i++) {
			nshift = i+i;	// o_idx shift counts here run as >>2,4,...,30
			off_ptr = o_offsets_lo + ( ((o_idx>>nshift)&0x3) << 4 );
			p0  = off_ptr[0x0];p1  = off_ptr[0x1];p2  = off_ptr[0x2];p3  = off_ptr[0x3];p4  = off_ptr[0x4];p5  = off_ptr[0x5];p6  = off_ptr[0x6];p7  = off_ptr[0x7];p8  = off_ptr[0x8];p9  = off_ptr[0x9];pa  = off_ptr[0xa];pb  = off_ptr[0xb];pc  = off_ptr[0xc];pd  = off_ptr[0xd];pe  = off_ptr[0xe];pf  = off_ptr[0xf];
			j = reverse16[i]<<1;
			tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			addr = __B + o_offsets_hi[i];	// o_offsets_hi[] = p10,p20,...,pf0
			add0 = addr+p0; add1 = addr+p1; add2 = addr+p2; add3 = addr+p3; add4 = addr+p4; add5 = addr+p5; add6 = addr+p6; add7 = addr+p7;
				add8 = addr+p8; add9 = addr+p9; adda = addr+pa; addb = addr+pb; addc = addr+pc; addd = addr+pd; adde = addr+pe; addf = addr+pf;
			SSE2_RADIX16_DIF_FMA_OOP(
				tm1,OFF1,OFF2,OFF3,OFF4, add0,add1,add2,add3,add4,add5,add6,add7,add8,add9,adda,addb,addc,addd,adde,addf, tm2
			);	tm1 += 2;
		}

	  #else	// Non-FMA version:

		// Remaining 15 sets of macro calls done in loop:
		for(i = 1; i < 16; i++) {
			nshift = i+i;	// o_idx shift counts here run as >>2,4,...,30
			off_ptr = o_offsets_lo + ( ((o_idx>>nshift)&0x3) << 4 );
			j = reverse16[i]<<1;
			tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			addr = __B + o_offsets_hi[i];
			SSE2_RADIX16_DIF_TWIDDLE_OOP(
				tm1,OFF1,OFF4,
				addr,off_ptr,
				isrt2,tm2
			);	tm1 += 2;
		}

	  #endif	// FMA/AVX2 ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif
	}

	void SSE2_RADIX256_DIT(
		// Inputs: Base address plus index offsets:
		double *__A,
		int *i_offsets_lo,	// Array storing  low parts of input index offsets in 16 slots
		uint32 i_idx,	// Bitfield encoding the sequence of the i_offsets_lo sub-vectors to use for the radix-256 DFT's inputs
		int *i_offsets_hi,	// Array storing high parts of input index offsets in 16 slots
		// Intermediates-storage pointer:
		vec_dbl*r00,
		// Pointers to base-roots data and first of 16 twiddle vectors:
		vec_dbl*isrt2, vec_dbl*two, vec_dbl*twid0,
		// Output pointer: Base ptr of 16 local-mem:
		vec_dbl*__B
	)
	{
		// Intermediates pointers:
		vec_dbl *tm0,*tm1,*tm2;
		/* Addresses into array sections */
		double *addr;
		// Index-offset names here reflect original unpermuted inputs, but the math also works for permuted ones:
		int i,j,nshift;
		const int *off_ptr;

	/* Gather the needed data (256 64-bit complex, i.e. 512 64-bit reals) and do 8 twiddleless length-16 subtransforms: */
	  #ifdef USE_ARM_V8_SIMD
		uint32 OFF1,OFF2,OFF3,OFF4;
		OFF1 = 0x20;
		OFF2 = 0x40;
		OFF3 = 0x60;
		OFF4 = 0x80;
	  #elif defined(USE_AVX512)
		#define OFF1	0x20*4
		#define OFF2	0x40*4
		#define OFF3	0x60*4
		#define OFF4	0x80*4
	  #elif defined(USE_AVX)
		#define OFF1	0x20*2
		#define OFF2	0x40*2
		#define OFF3	0x60*2
		#define OFF4	0x80*2
	  #else
		#define OFF1	0x20
		#define OFF2	0x40
		#define OFF3	0x60
		#define OFF4	0x80
	  #endif

		// Extract index of the 16-element i_offsets_lo sub-vector to use for the current set of outputs:
		off_ptr = i_offsets_lo + ( (i_idx&0x3) << 4 );	// 16*(Low 2 bits of i_idx); loop below will use remaining 30 bits in ascending pairs

	// Gather the needed data and do 16 twiddleless length-16 subtransforms, with p-offsets in-order:

		tm0 = r00;
		for(i = 0; i < 16; i++) {
			addr = __A + i_offsets_hi[i];	// i_offsets_hi[] = p10,p20,...,pf0
			SSE2_RADIX16_DIT_0TWIDDLE(
				addr,off_ptr,
				isrt2,two,
				tm0,OFF1,OFF2,OFF3,OFF4
			); tm0 += 32;
			if(i_idx) {	// vvv +2 here because this is setup for next loop pass
				nshift = i+i+2;	// i_idx shift counts here run as >>2,4,...,30
				off_ptr = i_offsets_lo + ( ((i_idx>>nshift)&0x3) << 4 );
			}
		}

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif

	/*...and now do 16 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f].
	*/
	  #ifdef USE_ARM_V8_SIMD
		OFF1 = 0x200;
		OFF2 = 0x400;
		OFF3 = 0x600;
		OFF4 = 0x800;
	  #elif defined(USE_AVX512)
		#define OFF1	0x200*4
		#define OFF2	0x400*4
		#define OFF3	0x600*4
		#define OFF4	0x800*4
	  #elif defined(USE_AVX)
		#define OFF1	0x200*2
		#define OFF2	0x400*2
		#define OFF3	0x600*2
		#define OFF4	0x800*2
	  #else
		#define OFF1	0x200
		#define OFF2	0x400
		#define OFF3	0x600
		#define OFF4	0x800
	  #endif

	// Block 0: All unity twiddles:
		tm1 = __B;
		// Jan 2020: For reduced-#args version of macro, need base-pointer r00 and array of 16
		// double-array-equivalent index offsets corr. to the address strides between r00 and r10,r20...,
		// which are 32 vec_dbl = 32*RE_IM_STRIDE doubles apart:
		const int dd = RE_IM_STRIDE<<5, offs[16] = {0,dd,2*dd,3*dd,4*dd,5*dd,6*dd,7*dd,8*dd,9*dd,10*dd,11*dd,12*dd,13*dd,14*dd,15*dd};
		off_ptr = &(offs[0]);
		SSE2_RADIX16_DIT_0TWIDDLE(
			r00,off_ptr,
			isrt2,two,
			tm1,OFF1,OFF2,OFF3,OFF4
		);

	  #ifdef USE_AVX2

		// Due to tangent-twiddles scheme and resulting singularity of tangent(arg(I)) = 1/0,
		// only last 14 of the 15 with-twiddles DFTs allow use of FMA-based macros under Intel AVX2/FMA3:
	// Block 8: BR twiddles = {  I.{},  C^ 8,-~C^ 8,  C^ 4,*~C^ 4, *C^ 4,-~C^ 4,  C^ 2,*~C^ 2, *C^ 6,-~C^ 6,  C^ 6,*~C^ 6, *C^ 2,-~C^ 2}
		j = 8<<1;	// Mimics elided i = 1 pass of length-14 loop below
		tm0 = r00 + j; tm1 = __B+j; tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
		SSE2_RADIX16_DIT_TWIDDLE_OOP(
			tm0,OFF1,OFF2,OFF3,OFF4, tm1,OFF1,OFF2,OFF3,OFF4, isrt2, tm2
		);
		// Remaining 14 sets of macro calls done in loop:
		for(i = 2; i < 16; i++) {
			j = reverse16[i]<<1;	// __B-offsets are processed in BR16 order
			tm0 = r00 + j; tm1 = __B+j; tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			SSE2_RADIX16_DIT_FMA_OOP(
				tm0,OFF1,OFF2,OFF3,OFF4, tm1,OFF1,OFF2,OFF3,OFF4, tm2
			);
		}

	  #else	// Non-FMA version:

		// Remaining 15 sets of macro calls done in loop:
		for(i = 1; i < 16; i++) {
			j = reverse16[i]<<1;	// __B-offsets are processed in BR16 order
			tm0 = r00 + j; tm1 = __B+j; tm2 = twid0 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			SSE2_RADIX16_DIT_TWIDDLE_OOP(
				tm0,OFF1,OFF2,OFF3,OFF4, tm1,OFF1,OFF2,OFF3,OFF4, isrt2, tm2
			);
		}

	  #endif	// FMA/AVX2 ?

	  #ifndef USE_ARM_V8_SIMD
		#undef OFF1
		#undef OFF2
		#undef OFF3
		#undef OFF4
	  #endif
	}

#endif	// USE_SSE2 = False, or non-64-bit-GCC:

