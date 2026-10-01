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
*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
*  more details.                                                               *
*                                                                              *
*  You should have received a copy of the GNU General Public License along     *
*  with this program; see the file GPL.txt.  If not, you may view one at       *
*  http://www.fsf.org/licenses/licenses.html, or obtain one by writing to the  *
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,        *
*  MA 02111-1307, USA.                                                         *
*                                                                              *
*******************************************************************************/

/* Run-time access to the Windows 7 / Server 2008 R2 processor-group API.

makemake.sh targets Vista (-D_WIN32_WINNT=0x0600) for msys/cygwin builds, and below 0x0601 the
mingw-w64 headers do not declare SetThreadGroupAffinity, GetActiveProcessorCount or
GetActiveProcessorGroupCount. Raising the target instead is not an option: calling them directly
makes them static imports, and a binary carrying those cannot be loaded at all on Vista - a worse
outcome than losing the feature. Resolving them through GetProcAddress gives one binary that uses
processor groups wherever they exist and still starts on Vista.

GROUP_AFFINITY itself is declared even at 0x0600, so only the entry points need typedefs. Include
this after <windows.h>.
*/
#ifndef win_procgroup_h_included
#define win_procgroup_h_included

#if defined(OS_TYPE_WINDOWS) || defined(__MINGW32__)

  #ifndef ALL_PROCESSOR_GROUPS
	#define ALL_PROCESSOR_GROUPS 0xffff
  #endif

typedef BOOL  (WINAPI *PFN_SetThreadGroupAffinity     )(HANDLE, const GROUP_AFFINITY*, PGROUP_AFFINITY);
typedef DWORD (WINAPI *PFN_GetActiveProcessorCount    )(WORD);
typedef WORD  (WINAPI *PFN_GetActiveProcessorGroupCount)(void);

extern PFN_SetThreadGroupAffinity		pSetThreadGroupAffinity;
extern PFN_GetActiveProcessorCount		pGetActiveProcessorCount;
extern PFN_GetActiveProcessorGroupCount	pGetActiveProcessorGroupCount;

/* Resolve the above. Idempotent, but call it only from the main thread, before any worker is
spawned - the workers read the pointers and must never race the writes. */
void win7_procgroup_init(void);

#endif	// OS_TYPE_WINDOWS || __MINGW32__

#endif	// win_procgroup_h_included
