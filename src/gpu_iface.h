/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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

// Thanks to Jason Papadopoulos for the original version of the GPU interface ... this is now so
// heavily modified by me that any resemblance to the original in the nontrivial details should be
// considered coincidental, and any faults strictly mine.

#ifndef gpu_iface_h_included
#define gpu_iface_h_included

#ifndef GPU_IFACE_STANDALONE
	// Non-standalone build assumes the non-main functions in this file will serve as GPU diagnostics
	// for an Mlucas or Mfactor build, so require same compile flag as for the other sources in such a build:
	#ifndef USE_GPU
		#error Compilation of any source file using a gpu-specific header requires the user-defined preprocessor flag USE_GPU
	#endif

	#include "masterdefs.h"
	#include "types.h"
#else
	#include <stdio.h>
	typedef int int32;
#endif

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_types.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_GPU 16

typedef struct cudaDeviceProp gpu_info_t;
/*
cudaDeviceProp struct members:

int 	canMapHostMemory 	Device can map host memory with cudaHostAlloc/cudaHostGetDevicePointer.
int 	clockRate			Clock frequency in kilohertz.
int 	computeMode			Compute mode (See cudaComputeMode).
int 	deviceOverlap		Device can concurrently copy memory and execute a kernel.
int 	integrated			Device is integrated as opposed to discrete.
int 	kernelExecTimeoutEnabled 	Specified whether there is a run time limit on kernels.
int 	major				Major compute capability.
int 	minor				Minor compute capability.
int 	maxGridSize [3] 	Maximum size of each dimension of a grid.
int 	maxThreadsDim [3] 	Maximum size of each dimension of a block.
int 	maxThreadsPerBlock 	Maximum number of threads per block.
size_t 	memPitch			Maximum pitch in bytes allowed by memory copies.
int 	multiProcessorCount Number of multiprocessors on device.
char 	name [256]			ASCII string identifying device.
int 	regsPerBlock		32-bit registers available per block
size_t 	sharedMemPerBlock 	Shared memory available per block in bytes.
size_t 	textureAlignment 	Alignment requirement for textures.
size_t 	totalConstMem		Constant memory available on device in bytes.
size_t 	totalGlobalMem		Global memory available on device in bytes.
int 	warpSize			Warp size in threads.
*/

typedef struct {
	int32 num_gpu;
	gpu_info_t gpu_info[MAX_GPU];
} gpu_config_t;

char * cuGetErrorMessage(CUresult result);

void gpu_init(gpu_config_t *config);

#define CUDA_TRY(func) \
	{ 			 				\
		CUresult status = func;				\
		if (status != CUDA_SUCCESS) {			\
			printf("error (line %d): %s\n", __LINE__,\
				cuGetErrorMessage(status));	\
			exit(-1);				\
		}						\
	}

#define CUDA_ALIGN_PARAM(offset, pow2align) \
	(offset) = ((offset) + (pow2align) - 1) & ~((pow2align) - 1)

#ifdef __cplusplus
}
#endif

#endif /* !gpu_iface_h_included_ */

