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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef radix4032_included
#define radix4032_included

	#include "types.h"

	// Low parts [p0-3f] of output-index perms - need 64 bytes (one init row below) for each radix-64 DFT:
	const uint8 dif64_oidx_lo[4032] = {
		0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,
		0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,
		0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,
		0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,
		0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,
		0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,
		0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,
		0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,
		0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,
		0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,
		0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,
		0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,
		0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,
		0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,
		0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,
		0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,
		0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,
		0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,
		0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
		0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,
		0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,
		0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,
		0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,
		0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,
		0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,
		0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,
		0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,
		0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,
		0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,
		0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,
		0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,
		0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,
		0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,
		0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,
		0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,
		0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,
		0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
		0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,
		0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,
		0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,
		0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,
		0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,
		0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,
		0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,
		0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,
		0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,
		0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,
		0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,
		0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,
		0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,
		0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,
		0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,
		0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,
		0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,
		0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,0x2f,0x2e,0x2c,0x2d,0x28,0x29,0x2a,0x2b,0x37,0x36,0x34,0x35,0x30,0x31,0x32,0x33,0x3f,0x3e,0x3c,0x3d,0x38,0x39,0x3a,0x3b,0x17,0x16,0x14,0x15,0x10,0x11,0x12,0x13,0x1f,0x1e,0x1c,0x1d,0x18,0x19,0x1a,0x1b,0x0f,0x0e,0x0c,0x0d,0x08,0x09,0x0a,0x0b,0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
		0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,0x05,0x04,0x07,0x06,0x03,0x02,0x00,0x01,0x19,0x18,0x1b,0x1a,0x1d,0x1c,0x1f,0x1e,0x15,0x14,0x17,0x16,0x13,0x12,0x10,0x11,0x29,0x28,0x2b,0x2a,0x2d,0x2c,0x2f,0x2e,0x25,0x24,0x27,0x26,0x23,0x22,0x20,0x21,0x39,0x38,0x3b,0x3a,0x3d,0x3c,0x3f,0x3e,0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,
		0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,0x32,0x33,0x31,0x30,0x36,0x37,0x35,0x34,0x22,0x23,0x21,0x20,0x26,0x27,0x25,0x24,0x2a,0x2b,0x29,0x28,0x2e,0x2f,0x2d,0x2c,0x02,0x03,0x01,0x00,0x06,0x07,0x05,0x04,0x0a,0x0b,0x09,0x08,0x0e,0x0f,0x0d,0x0c,0x12,0x13,0x11,0x10,0x16,0x17,0x15,0x14,0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,
		0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,0x1b,0x1a,0x18,0x19,0x1f,0x1e,0x1c,0x1d,0x0b,0x0a,0x08,0x09,0x0f,0x0e,0x0c,0x0d,0x07,0x06,0x04,0x05,0x00,0x01,0x02,0x03,0x33,0x32,0x30,0x31,0x37,0x36,0x34,0x35,0x3b,0x3a,0x38,0x39,0x3f,0x3e,0x3c,0x3d,0x2b,0x2a,0x28,0x29,0x2f,0x2e,0x2c,0x2d,0x27,0x26,0x24,0x25,0x20,0x21,0x22,0x23,
		0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,0x21,0x20,0x23,0x22,0x25,0x24,0x27,0x26,0x3e,0x3f,0x3d,0x3c,0x39,0x38,0x3b,0x3a,0x31,0x30,0x33,0x32,0x35,0x34,0x37,0x36,0x1e,0x1f,0x1d,0x1c,0x19,0x18,0x1b,0x1a,0x11,0x10,0x13,0x12,0x15,0x14,0x17,0x16,0x01,0x00,0x03,0x02,0x05,0x04,0x07,0x06,0x09,0x08,0x0b,0x0a,0x0d,0x0c,0x0f,0x0e,
		0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00,0x0c,0x0d,0x0e,0x0f,0x0a,0x0b,0x09,0x08,0x14,0x15,0x16,0x17,0x12,0x13,0x11,0x10,0x1c,0x1d,0x1e,0x1f,0x1a,0x1b,0x19,0x18,0x24,0x25,0x26,0x27,0x22,0x23,0x21,0x20,0x2c,0x2d,0x2e,0x2f,0x2a,0x2b,0x29,0x28,0x34,0x35,0x36,0x37,0x32,0x33,0x31,0x30,0x3c,0x3d,0x3e,0x3f,0x3a,0x3b,0x39,0x38,
		0x35,0x34,0x37,0x36,0x33,0x32,0x30,0x31,0x3d,0x3c,0x3f,0x3e,0x3b,0x3a,0x38,0x39,0x2d,0x2c,0x2f,0x2e,0x2b,0x2a,0x28,0x29,0x23,0x22,0x20,0x21,0x27,0x26,0x24,0x25,0x0d,0x0c,0x0f,0x0e,0x0b,0x0a,0x08,0x09,0x03,0x02,0x00,0x01,0x07,0x06,0x04,0x05,0x1d,0x1c,0x1f,0x1e,0x1b,0x1a,0x18,0x19,0x13,0x12,0x10,0x11,0x17,0x16,0x14,0x15,
		0x1a,0x1b,0x19,0x18,0x1e,0x1f,0x1d,0x1c,0x16,0x17,0x15,0x14,0x11,0x10,0x13,0x12,0x06,0x07,0x05,0x04,0x01,0x00,0x03,0x02,0x0e,0x0f,0x0d,0x0c,0x09,0x08,0x0b,0x0a,0x3a,0x3b,0x39,0x38,0x3e,0x3f,0x3d,0x3c,0x36,0x37,0x35,0x34,0x31,0x30,0x33,0x32,0x26,0x27,0x25,0x24,0x21,0x20,0x23,0x22,0x2e,0x2f,0x2d,0x2c,0x29,0x28,0x2b,0x2a,
		0x20,0x21,0x22,0x23,0x24,0x25,0x26,0x27,0x28,0x29,0x2a,0x2b,0x2c,0x2d,0x2e,0x2f,0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x3b,0x3c,0x3d,0x3e,0x3f,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,0x04,0x05,0x06,0x07,0x02,0x03,0x01,0x00
	};

#endif	/* #ifndef radix4032_included */
