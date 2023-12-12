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

#include "gpu_iface.h"

#ifdef __CUDACC__
	#warning using nvcc
	#ifdef __CUDA_ARCH__
		#warning device code trajectory
		#if __CUDA_ARCH__ > 120
			#warning compiling with double precision
		#else
			#warning compiling with single precision
		#endif
	#else
		#warning nvcc host code trajectory
	#endif
#else
	#warning non-nvcc code trajectory
#endif

#ifndef OS_BITS
	#error Bitness not defined!
#elif OS_BITS == 32
	#warning compiling in 32-bit mode
#elif OS_BITS == 64
	#warning compiling in 64-bit mode
#else
	#error Bitness defined but not supported!
#endif

// 50 Ways to say "Houston, we have a problem":
char *
cuGetErrorMessage(CUresult result)
{
	switch (result) {
	case CUDA_SUCCESS: return "CUDA_SUCCESS";
	case CUDA_ERROR_INVALID_VALUE: return "CUDA_ERROR_INVALID_VALUE";
	case CUDA_ERROR_OUT_OF_MEMORY: return "CUDA_ERROR_OUT_OF_MEMORY";
	case CUDA_ERROR_NOT_INITIALIZED: return "CUDA_ERROR_NOT_INITIALIZED";
	case CUDA_ERROR_DEINITIALIZED: return "CUDA_ERROR_DEINITIALIZED";
	case CUDA_ERROR_NO_DEVICE: return "CUDA_ERROR_NO_DEVICE";
	case CUDA_ERROR_INVALID_DEVICE: return "CUDA_ERROR_INVALID_DEVICE";
	case CUDA_ERROR_INVALID_IMAGE: return "CUDA_ERROR_INVALID_IMAGE";
	case CUDA_ERROR_INVALID_CONTEXT: return "CUDA_ERROR_INVALID_CONTEXT";
	case CUDA_ERROR_CONTEXT_ALREADY_CURRENT: return "CUDA_ERROR_CONTEXT_ALREADY_CURRENT";
	case CUDA_ERROR_MAP_FAILED: return "CUDA_ERROR_MAP_FAILED";
	case CUDA_ERROR_UNMAP_FAILED: return "CUDA_ERROR_UNMAP_FAILED";
	case CUDA_ERROR_ARRAY_IS_MAPPED: return "CUDA_ERROR_ARRAY_IS_MAPPED";
	case CUDA_ERROR_ALREADY_MAPPED: return "CUDA_ERROR_ALREADY_MAPPED";
	case CUDA_ERROR_NO_BINARY_FOR_GPU: return "CUDA_ERROR_NO_BINARY_FOR_GPU";
	case CUDA_ERROR_ALREADY_ACQUIRED: return "CUDA_ERROR_ALREADY_ACQUIRED";
	case CUDA_ERROR_NOT_MAPPED: return "CUDA_ERROR_NOT_MAPPED";
	case CUDA_ERROR_INVALID_SOURCE: return "CUDA_ERROR_INVALID_SOURCE";
	case CUDA_ERROR_FILE_NOT_FOUND: return "CUDA_ERROR_FILE_NOT_FOUND";
	case CUDA_ERROR_INVALID_HANDLE: return "CUDA_ERROR_INVALID_HANDLE";
	case CUDA_ERROR_NOT_FOUND: return "CUDA_ERROR_NOT_FOUND";
	case CUDA_ERROR_NOT_READY: return "CUDA_ERROR_NOT_READY";
	case CUDA_ERROR_LAUNCH_FAILED: return "CUDA_ERROR_LAUNCH_FAILED";
	case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES: return "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES";
	case CUDA_ERROR_LAUNCH_TIMEOUT: return "CUDA_ERROR_LAUNCH_TIMEOUT";
	case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING: return "CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING";
	case CUDA_ERROR_UNKNOWN: return "CUDA_ERROR_UNKNOWN";
	default: return "CUDA: unexpected error";
	}
}

// Read information on all available GPUs into input arg:
void
gpu_init(gpu_config_t *gpu_config)
{
	int32 device, nskip = 0;
	memset(gpu_config, 0, sizeof(gpu_config_t));

//	CUDA_TRY(cudaGetDeviceCount(&gpu_config->num_gpu))	*** error: a value of type "cudaError_t" cannot be used to initialize an entity of type "CUresult"
	cudaGetDeviceCount(&gpu_config->num_gpu);
	for (device = 0; device < (int32)gpu_config->num_gpu; device++)
	{
		// Get pointer to info for [device]th GPU having the minimum required capability:
		gpu_info_t *info = gpu_config->gpu_info + device - nskip;
//		CUDA_TRY(cudaGetDeviceProperties(info, device))	*** error: a value of type "cudaError_t" cannot be used to initialize an entity of type "CUresult"
		cudaGetDeviceProperties(info, device);
		if(info->major < 2) {
			printf("GPU #%d compute capability %d.%d is less than min-supported 2.x ... ignoring this device.\n",device,info->major,info->minor);
			++nskip;
		}
		// Note: Devices with cc = 2.x have (32 + 16*x) shader cores per multiprocessor (At least for x = 0 and 1 ... may need table for this
	}
	gpu_config->num_gpu -= nskip;
	return;
}

#ifdef GPU_IFACE_STANDALONE
	int main(int argc, char *argv[])
	{
		gpu_config_t gpu_config;
		gpu_info_t ginfo;
		int32 igpu;

		gpu_init(&gpu_config);
		if (gpu_config.num_gpu > 0) {
			printf("Detected %u CUDA-enabled GPU devices.\n", gpu_config.num_gpu);
			for(igpu = 0; igpu < gpu_config.num_gpu; ++igpu) {
				ginfo = gpu_config.gpu_info[igpu];
				printf("GPU #%u: %s v%u.%u\n", igpu, ginfo.name, ginfo.major, ginfo.minor);
				printf("clock_speed = %u MHz\n", ginfo.clockRate/1000);
				printf("num_compute_units = %u\n", ginfo.multiProcessorCount);
				printf("constant_mem_size = %u\n", ginfo.totalConstMem);
				printf("shared_mem_size = %u\n", ginfo.sharedMemPerBlock);
				printf("global_mem_size = %u\n", ginfo.totalGlobalMem);
				printf("registers_per_block = %u\n", ginfo.regsPerBlock);
				printf("max_threads_per_block = %u\n", ginfo.maxThreadsPerBlock);
				printf("can_overlap = %u\n", ginfo.deviceOverlap);
				printf("concurrent_kernels = %u\n", ginfo.concurrentKernels);
				printf("warp_size = %u\n", ginfo.warpSize);
				printf("max_thread_dim[3] = [%u,%u,%u]\n", ginfo.maxThreadsDim[0], ginfo.maxThreadsDim[1], ginfo.maxThreadsDim[2]);
				printf("max_grid_size[3] = [%u,%u,%u]\n", ginfo.maxGridSize[0], ginfo.maxGridSize[1], ginfo.maxGridSize[2]);
			}
			exit(0);
		} else {
			printf("ERROR: No CUDA-enabled GPUs found\n");
			exit(-1);
		}
	}
#endif

