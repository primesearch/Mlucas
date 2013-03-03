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

#include "util.h"

#if(defined(CPU_IS_X86) || defined(CPU_IS_IA64) || defined(CPU_IS_X86_64))

	/* cpuid(func,ax,bx,cx,dx) macro: put whatever function number you want in for func,
	and put in 4 32-bit ints to store the output values of eax, ebx, ecx, and edx, respectively.

	See the wonderfully brief note copied from the Intel ISR, below for a complete description
	of the various types of information returned, for different input values of the eax register (func).
	The most common one (and the one we use in our SSEx detection) is func = 1.
	*/

  #if(defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))

	#define CPUID(func,ax,bx,cx,dx)\
		__asm__ __volatile__ ("cpuid":\
	"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func));

/*
	uint32	get_cpuid()
	{
		uint32 retval = 0;
		asm (
		"movl     $1,%%eax	\n\t"
		"cpuid            	\n\t"
		"mov      %%edx,%0	\n\t"
	: : "r" (retval) : "eax","ebx","ecx","edx","cc" );
		return retval;
	}
*/

  #elif(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_ICC))

	/* For MSVC and ICC, Replace with fancier call based on __cpuid intrisic: */
#if 0
	#define CPUID(func,a,b,c,d)\
	{\
		__asm	mov	eax, func\
		__asm	cpuid\
		__asm	mov	a, eax\
		__asm	mov	b, ebx\
		__asm	mov	c, ecx\
		__asm	mov	d, edx\
	}
#endif

	#include <intrin.h>

	int CPUInfo[4] = {-1};

	#define CPUID(func,a,b,c,d)\
	{\
		__cpuid(CPUInfo, func);	\
		a = CPUInfo[0];	\
		b = CPUInfo[1];	\
		c = CPUInfo[2];	\
		d = CPUInfo[3];	\
	}

	const char* szFeatures[] =
	{
		"x87 FPU On Chip",
		"Virtual-8086 Mode Enhancement",
		"Debugging Extensions",
		"Page Size Extensions",
		"Time Stamp Counter",
		"RDMSR and WRMSR Support",
		"Physical Address Extensions",
		"Machine Check Exception",
		"CMPXCHG8B Instruction",
		"APIC On Chip",
		"Unknown1",
		"SYSENTER and SYSEXIT",
		"Memory Type Range Registers",
		"PTE Global Bit",
		"Machine Check Architecture",
		"Conditional Move/Compare Instruction",
		"Page Attribute Table",
		"Page Size Extension",
		"Processor Serial Number",
		"CFLUSH Extension",
		"Unknown2",
		"Debug Store",
		"Thermal Monitor and Clock Ctrl",
		"MMX Technology",
		"FXSAVE/FXRSTOR",
		"SSE Extensions",
		"SSE2 Extensions",
		"Self Snoop",
		"Hyper-threading Technology",
		"Thermal Monitor",
		"Unknown4",
		"Pend. Brk. EN."
	};

	void	cpu_details(void)
	{
		char CPUString[0x20];
		char CPUBrandString[0x40];
		int nSteppingID = 0;
		int nModel = 0;
		int nFamily = 0;
		int nProcessorType = 0;
		int nExtendedmodel = 0;
		int nExtendedfamily = 0;
		int nBrandIndex = 0;
		int nCLFLUSHcachelinesize = 0;
		int nAPICPhysicalID = 0;
		int nFeatureInfo = 0;
		int nCacheLineSize = 0;
		int nL2Associativity = 0;
		int nCacheSizeK = 0;
		unsigned	nIds, nExIds, i;
		int	bSSE3NewInstructions = FALSE;
		int	bMONITOR_MWAIT = FALSE;
		int	bCPLQualifiedDebugStore = FALSE;
		int	bThermalMonitor2 = FALSE;

		/*
		// __cpuid with an InfoType argument of 0 returns the number of
		// valid Ids in CPUInfo[0] and the CPU identification string in
		// the other three array elements. The CPU identification string is
		// not in linear order. The code below arranges the information
		// in a human readable form.
		*/
		__cpuid(CPUInfo, 0);
		nIds = CPUInfo[0];
		memset(CPUString, 0, sizeof(CPUString));
		*((int*)CPUString) = CPUInfo[1];
		*((int*)(CPUString+4)) = CPUInfo[3];
		*((int*)(CPUString+8)) = CPUInfo[2];

		/* Get the information associated with each valid Id */
		for (i=0; i<=nIds; ++i)
		{
			__cpuid(CPUInfo, i);
		#if 0
			printf("\nFor InfoType %d\n", i);
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
		#endif
			/* Interpret CPU feature information. */
			if(i == 1)
			{
				nSteppingID = CPUInfo[0] & 0xf;
				nModel = (CPUInfo[0] >> 4) & 0xf;
				nFamily = (CPUInfo[0] >> 8) & 0xf;
				nProcessorType = (CPUInfo[0] >> 12) & 0x3;
				nExtendedmodel = (CPUInfo[0] >> 16) & 0xf;
				nExtendedfamily = (CPUInfo[0] >> 20) & 0xff;
				nBrandIndex = CPUInfo[1] & 0xff;
				nCLFLUSHcachelinesize = ((CPUInfo[1] >> 8) & 0xff) * 8;
				nAPICPhysicalID = (CPUInfo[1] >> 24) & 0xff;
				bSSE3NewInstructions = (CPUInfo[2] & 0x1);
				bMONITOR_MWAIT = (CPUInfo[2] & 0x8);
				bCPLQualifiedDebugStore = (CPUInfo[2] & 0x10);
				bThermalMonitor2 = (CPUInfo[2] & 0x100);
				nFeatureInfo = CPUInfo[3];
			}
		}

		/* Calling __cpuid with 0x80000000 as the InfoType argument
		gets the number of valid extended IDs.
		 */
		__cpuid(CPUInfo, 0x80000000);
		nExIds = CPUInfo[0];
		memset(CPUBrandString, 0, sizeof(CPUBrandString));

		/* Get the information associated with each extended ID. */
		for (i=0x80000000; i<=nExIds; ++i)
		{
			__cpuid(CPUInfo, i);
			printf("\nFor InfoType %x\n", i);
			printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
			printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
			printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
			printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);

			/* Interpret CPU brand string and cache information. */
			if(i == 0x80000002)
				memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
			else if(i == 0x80000003)
				memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
			else if(i == 0x80000004)
				memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
			else if(i == 0x80000006)
			{
				nCacheLineSize = CPUInfo[2] & 0xff;
				nL2Associativity = (CPUInfo[2] >> 12) & 0xf;
				nCacheSizeK = (CPUInfo[2] >> 16) & 0xffff;
			}
		}

		/* Display all the information in user-friendly format. */

		printf("\n\nCPU String: %s\n", CPUString);

		if(nIds >= 1)
		{
			if(nSteppingID)
				printf("Stepping ID = %d\n", nSteppingID);
			if(nModel)
				printf("Model = %d\n", nModel);
			if(nFamily)
				printf("Family = %d\n", nFamily);
			if(nProcessorType)
				printf("Processor Type = %d\n", nProcessorType);
			if(nExtendedmodel)
				printf("Extended model = %d\n", nExtendedmodel);
			if(nExtendedfamily)
				printf("Extended family = %d\n", nExtendedfamily);
			if(nBrandIndex)
				printf("Brand Index = %d\n", nBrandIndex);
			if(nCLFLUSHcachelinesize)
				printf("CLFLUSH cache line size = %d\n", nCLFLUSHcachelinesize);
			if(nAPICPhysicalID)
				printf("APIC Physical ID = %d\n", nAPICPhysicalID);

			if(nFeatureInfo || bSSE3NewInstructions || bMONITOR_MWAIT || bCPLQualifiedDebugStore || bThermalMonitor2)
			{
				printf("\nThe following features are supported:\n");

				if(bSSE3NewInstructions)
					printf("\tSSE3 New Instructions\n");
				if(bMONITOR_MWAIT)
					printf("\tMONITOR/MWAIT\n");
				if(bCPLQualifiedDebugStore)
					printf("\tCPL Qualified Debug Store\n");
				if(bThermalMonitor2)
					printf("\tThermal Monitor 2\n");

				i = 0;
				nIds = 1;
				while (i < (sizeof(szFeatures)/sizeof(const char*)))
				{
					if(nFeatureInfo & nIds)
					{
						printf("\t");
						printf(szFeatures[i]);
						printf("\n");
					}

					nIds <<= 1;
					++i;
				}
			}
		}

		if(nExIds >= 0x80000004)
			printf("\nCPU Brand String: %s\n", CPUBrandString);

		if(nExIds >= 0x80000006)
		{
			printf("Cache Line Size = %d\n", nCacheLineSize);
			printf("L2 Associativity = %d\n", nL2Associativity);
			printf("Cache Size = %dK\n\n", nCacheSizeK);
		}
	}

	/*
	Obsolete version showing all the junky hacks needed for pre-CPUID-compliant CPUs:
	uint32	get_cpuid()
	{
		uint32 retval;
		{
			// Detection of pre-80286/80286/386+ processors //
			__asm		mov     ax, 7202h	// set bits 12-14 and clear bit 15 //
			__asm		push    ax
			__asm		popf
			__asm		pushf
			__asm		pop     ax

			__asm		test    ah, 0f0h
			__asm		js      end			// bit 15 of FLAGS is set on pre-286 //
			__asm		jz      end			// bits 12..15 of FLAGS are clear on 286 //

			// it's a 80386 or higher processor //
			// Detection of 80386/80486(w/out CPUID)/80486+(CPUID compliant) //
			__asm		pushfd
			__asm		pop     eax
			__asm		mov     edx, eax
			__asm		xor     eax, 00240000h	// flip bits 18 (AC) and 21 (ID) //
			__asm		push    eax
			__asm		popfd
			__asm		pushfd
			__asm		pop     eax

			__asm		xor     eax, edx	// check if both bits didn't toggle //
			__asm		jz      end			// 386 //
			__asm		shr     eax, 19		// check if only bit 18 toggled //
			__asm		jz      end			// 80486_without_CPUID //
			__asm		mov     eax, 1
			__asm		cpuid
			__asm		jmp     en
			__asm	end:xor     edx,edx
			__asm	en:	mov     retval,edx
		};
		return retval;
	}
	*/
  #else

/*	#error get_cpuid() only supported for ia32 under GCC, MSVC, and Intel C!	*/
	#define CPUID(x,a,b,c,d)	/* */

  #endif

	/* Calling CPUID with EAX = 0 returns maximum input value of EAX register for CPUID call in EAX, vendor string <"Your compnay Logo Here!!!> in EBX:EDX:ECX: */
	void	get_cpu(void)
	{
		uint32 i,j,shift;
		uint32 a,b,c,d;
		/* Need 12+3*16 char for extended cpuid, add 4 more for terminator and " : " padding */
		char cpu_str[65] = "0123456789ab : f0123456789abcdef0123456789abcdef0123456789abcdef";

		CPUID(0,a,b,c,d);
		for(i = 0; i < 4; i++)
		{
			shift = i*8;
			cpu_str[i   ] = (char)((b >> shift) & 0xff);
			cpu_str[i+ 4] = (char)((d >> shift) & 0xff);
			cpu_str[i+ 8] = (char)((c >> shift) & 0xff);
		}

		CPUID(0x80000002,a,b,c,d);
		for(i = 0; i < 4; i++)
		{
			j = i+15;
			shift = i*8;
			cpu_str[j   ] = (char)((a >> shift) & 0xff);
			cpu_str[j+ 4] = (char)((b >> shift) & 0xff);
			cpu_str[j+ 8] = (char)((c >> shift) & 0xff);
			cpu_str[j+12] = (char)((d >> shift) & 0xff);
		}

		CPUID(0x80000003,a,b,c,d);
		for(i = 0; i < 4; i++)
		{
			j = i+31;
			shift = i*8;
			cpu_str[j   ] = (char)((a >> shift) & 0xff);
			cpu_str[j+ 4] = (char)((b >> shift) & 0xff);
			cpu_str[j+ 8] = (char)((c >> shift) & 0xff);
			cpu_str[j+12] = (char)((d >> shift) & 0xff);
		}

		CPUID(0x80000004,a,b,c,d);
		for(i = 0; i < 4; i++)
		{
			j = i+47;
			shift = i*8;
			cpu_str[j   ] = (char)((a >> shift) & 0xff);
			cpu_str[j+ 4] = (char)((b >> shift) & 0xff);
			cpu_str[j+ 8] = (char)((c >> shift) & 0xff);
			cpu_str[j+12] = (char)((d >> shift) & 0xff);
		}

		fprintf(stderr, "CPU Type = %s\n",cpu_str);
	}

	/* SSE  is bit 25 of EDX returned by calling CPUID with input EAX = 1: */
	uint32	has_sse()
	{
		uint32 a,b,c,d;
		CPUID(1,a,b,c,d);

		if(d & 0x02000000)
			return 1;
		else
			return 0;
	}

	/* SSE2 is bit 26 of EDX returned by calling CPUID with input EAX = 1: */
	uint32	has_sse2()
	{
		uint32 a,b,c,d;
		CPUID(1,a,b,c,d);

		if(d & 0x04000000)
			return 1;
		else
			return 0;
	}

	/* SSE3 is bit  0 of ECX returned by calling CPUID with input EAX = 1: */
	uint32	has_sse3()
	{
		uint32 a,b,c,d;
		CPUID(1,a,b,c,d);

		if(c & 0x00000001)
			return 1;
		else
			return 0;
	}

#endif

