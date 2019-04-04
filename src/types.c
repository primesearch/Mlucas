/*******************************************************************************
*                                                                              *
*   (C) 1997-2016 by Ernst W. Mayer.                                           *
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

#include "types.h"

/* Useful extern constants to export: */

/* Multiword ints have word significance increasing from left to right: */

/* 5/04/2005: uint96/160s are really uint128/192s with upper 32 bits zero: */
const uint96  NIL96  = {(uint64)0, (uint32)0};
const uint96  ONE96  = {(uint64)1, (uint32)0};
const uint96  TWO96  = {(uint64)2, (uint32)0};

const uint128 NIL128 = {(uint64)0, (uint64)0};
const uint128 ONE128 = {(uint64)1, (uint64)0};
const uint128 TWO128 = {(uint64)2, (uint64)0};

const uint160 NIL160 = {(uint64)0, (uint64)0, (uint32)0};
const uint160 ONE160 = {(uint64)1, (uint64)0, (uint32)0};
const uint160 TWO160 = {(uint64)2, (uint64)0, (uint32)0};

const uint192 NIL192 = {(uint64)0, (uint64)0, (uint64)0};
const uint192 ONE192 = {(uint64)1, (uint64)0, (uint64)0};
const uint192 TWO192 = {(uint64)2, (uint64)0, (uint64)0};

const uint256 NIL256 = {(uint64)0, (uint64)0, (uint64)0, (uint64)0};
const uint256 ONE256 = {(uint64)1, (uint64)0, (uint64)0, (uint64)0};
const uint256 TWO256 = {(uint64)2, (uint64)0, (uint64)0, (uint64)0};

const uint512 NIL512 = {(uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0};
const uint512 ONE512 = {(uint64)1, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0};
const uint512 TWO512 = {(uint64)2, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0, (uint64)0};

/* Binary predicates for use of stdlib qsort(): */
int ncmp_int(const void * a, const void * b)	// Default-int compare predicate
{
	return ( *(int*)a - *(int*)b );
}

int ncmp_uint32(const void * a, const void * b)	// Mnemonic: "Numeric CoMPare of UINT32 data"
{
	uint32 diff = *(uint32*)a - *(uint32*)b;
	uint32 borrow = 1 - ((diff > *(uint32*)a) << 1);	// -1 if (a < b), +1 otherwise
	// If (diff > a) == 1, had a borrow, i.e. a < b, return -1.
	// Otherwise return 0 if diff == 0, +1 if diff != 0. Can roll all 3 possibilities into one expression:
	return ( borrow & -(diff != 0) );
	/*
	a < b: bw = -1, (diff != 0) = 1, -() = -1 ===> -1 & -1 = -1
	a = b: bw = +1, (diff != 0) = 0, -() =  0 ===> +1 &  0 =  0
	a > b: bw = +1, (diff != 0) = 1, -() = -1 ===> +1 & -1 = +1
	*/
}

int ncmp_sint32(const void * a, const void * b)
{
	return ( *(sint32*)a - *(sint32*)b );
}

int ncmp_uint64(const void * a, const void * b)
{
	uint64 diff = *(uint64*)a - *(uint64*)b;
	uint64 borrow = 1 - ((diff > *(uint64*)a) << 1);	// -1 if (a < b), +1 otherwise
	return ( borrow & -(diff != 0) );
}

int ncmp_sint64(const void * a, const void * b)
{
	return ( *(sint64*)a - *(sint64*)b );
}

