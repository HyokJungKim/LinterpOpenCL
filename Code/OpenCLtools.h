#pragma once

#ifndef OPENCLTOOLS_H_
#define OPENCLTOOLS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.h>

void PrintErrorCL(cl_int inError)
{
  std::string tempString = std::to_string(inError);
  std::cout << "OpenCL Error: ";

  switch (inError) {
  case CL_SUCCESS: tempString = "SUCCESS"; break;
  case CL_DEVICE_NOT_FOUND: tempString = "CL_DEVICE_NOT_FOUND"; break;
  case CL_DEVICE_NOT_AVAILABLE: tempString = "CL_DEVICE_NOT_AVAILABLE"; break;
  case CL_COMPILER_NOT_AVAILABLE: tempString = "CL_COMPILER_NOT_AVAILABLE"; break;
  case CL_MEM_OBJECT_ALLOCATION_FAILURE: tempString = "CL_MEM_OBJECT_ALLOCATION_FAILURE"; break;
  case CL_OUT_OF_RESOURCES: tempString =  "CL_OUT_OF_RESOURCES"; break;
  case CL_OUT_OF_HOST_MEMORY: tempString =  "CL_OUT_OF_HOST_MEMORY"; break;
  case CL_PROFILING_INFO_NOT_AVAILABLE: tempString =  "CL_PROFILING_INFO_NOT_AVAILABLE"; break;
  case CL_MEM_COPY_OVERLAP: tempString =  "CL_MEM_COPY_OVERLAP"; break;
  case CL_IMAGE_FORMAT_MISMATCH: tempString =  "CL_IMAGE_FORMAT_MISMATCH"; break;
  case CL_IMAGE_FORMAT_NOT_SUPPORTED: tempString =  "CL_IMAGE_FORMAT_NOT_SUPPORTED"; break;
  case CL_BUILD_PROGRAM_FAILURE: tempString =  "CL_BUILD_PROGRAM_FAILURE"; break;
  case CL_MAP_FAILURE: tempString =  "CL_MAP_FAILURE"; break;
  case CL_MISALIGNED_SUB_BUFFER_OFFSET: tempString =  "CL_MISALIGNED_SUB_BUFFER_OFFSET"; break;
  case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST: tempString =  "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST"; break;
  case CL_COMPILE_PROGRAM_FAILURE: tempString =  "CL_COMPILE_PROGRAM_FAILURE"; break;
  case CL_LINKER_NOT_AVAILABLE: tempString =  "CL_LINKER_NOT_AVAILABLE"; break;
  case CL_LINK_PROGRAM_FAILURE: tempString =  "CL_LINK_PROGRAM_FAILURE"; break;
  case CL_DEVICE_PARTITION_FAILED: tempString =  "CL_DEVICE_PARTITION_FAILED"; break;
  case CL_KERNEL_ARG_INFO_NOT_AVAILABLE: tempString =  "CL_KERNEL_ARG_INFO_NOT_AVAILABLE"; break;

      // compile-time errors
  case CL_INVALID_VALUE: tempString =  "CL_INVALID_VALUE"; break;
  case CL_INVALID_DEVICE_TYPE: tempString =  "CL_INVALID_DEVICE_TYPE"; break;
  case CL_INVALID_PLATFORM: tempString =  "CL_INVALID_PLATFORM"; break;
  case CL_INVALID_DEVICE: tempString =  "CL_INVALID_DEVICE"; break;
  case CL_INVALID_CONTEXT: tempString =  "CL_INVALID_CONTEXT"; break;
  case CL_INVALID_QUEUE_PROPERTIES: tempString =  "CL_INVALID_QUEUE_PROPERTIES"; break;
  case CL_INVALID_COMMAND_QUEUE: tempString =  "CL_INVALID_COMMAND_QUEUE"; break;
  case CL_INVALID_HOST_PTR: tempString =  "CL_INVALID_HOST_PTR"; break;
  case CL_INVALID_MEM_OBJECT: tempString =  "CL_INVALID_MEM_OBJECT";
  case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: tempString =  "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"; break;
  case CL_INVALID_IMAGE_SIZE: tempString =  "CL_INVALID_IMAGE_SIZE"; break;
  case CL_INVALID_SAMPLER: tempString =  "CL_INVALID_SAMPLER"; break;
  case CL_INVALID_BINARY: tempString =  "CL_INVALID_BINARY"; break;
  case CL_INVALID_BUILD_OPTIONS: tempString =  "CL_INVALID_BUILD_OPTIONS"; break;
  case CL_INVALID_PROGRAM: tempString =  "CL_INVALID_PROGRAM"; break;
  case CL_INVALID_PROGRAM_EXECUTABLE: tempString =  "CL_INVALID_PROGRAM_EXECUTABLE"; break;
  case CL_INVALID_KERNEL_NAME: tempString =  "CL_INVALID_KERNEL_NAME"; break;
  case CL_INVALID_KERNEL_DEFINITION: tempString =  "CL_INVALID_KERNEL_DEFINITION"; break;
  case CL_INVALID_KERNEL: tempString =  "CL_INVALID_KERNEL"; break;
  case CL_INVALID_ARG_INDEX: tempString =  "CL_INVALID_ARG_INDEX"; break;
  case CL_INVALID_ARG_VALUE: tempString =  "CL_INVALID_ARG_VALUE"; break;
  case CL_INVALID_ARG_SIZE: tempString =  "CL_INVALID_ARG_SIZE"; break;
  case CL_INVALID_KERNEL_ARGS: tempString =  "CL_INVALID_KERNEL_ARGS"; break;
  case CL_INVALID_WORK_DIMENSION: tempString =  "CL_INVALID_WORK_DIMENSION"; break;
  case CL_INVALID_WORK_GROUP_SIZE: tempString =  "CL_INVALID_WORK_GROUP_SIZE"; break;
  case CL_INVALID_WORK_ITEM_SIZE: tempString =  "CL_INVALID_WORK_ITEM_SIZE"; break;
  case CL_INVALID_GLOBAL_OFFSET: tempString =  "CL_INVALID_GLOBAL_OFFSET"; break;
  case CL_INVALID_EVENT_WAIT_LIST: tempString =  "CL_INVALID_EVENT_WAIT_LIST"; break;
  case CL_INVALID_EVENT: tempString =  "CL_INVALID_EVENT"; break;
  case CL_INVALID_OPERATION: tempString =  "CL_INVALID_OPERATION"; break;
  case CL_INVALID_GL_OBJECT: tempString =  "CL_INVALID_GL_OBJECT"; break;
  case CL_INVALID_BUFFER_SIZE: tempString =  "CL_INVALID_BUFFER_SIZE"; break;
  case CL_INVALID_MIP_LEVEL: tempString =  "CL_INVALID_MIP_LEVEL"; break;
  case CL_INVALID_GLOBAL_WORK_SIZE: tempString =  "CL_INVALID_GLOBAL_WORK_SIZE"; break;
  case CL_INVALID_PROPERTY: tempString =  "CL_INVALID_PROPERTY"; break;
  case CL_INVALID_IMAGE_DESCRIPTOR: tempString =  "CL_INVALID_IMAGE_DESCRIPTOR"; break;
  case CL_INVALID_COMPILER_OPTIONS: tempString =  "CL_INVALID_COMPILER_OPTIONS"; break;
  case CL_INVALID_LINKER_OPTIONS: tempString =  "CL_INVALID_LINKER_OPTIONS"; break;
  case CL_INVALID_DEVICE_PARTITION_COUNT: tempString =  "CL_INVALID_DEVICE_PARTITION_COUNT"; break;

  default: tempString =  "Unknown OpenCL error number: " + tempString; break;
  }

  std::cout << tempString << std::endl;
}

cl_context CreateContext()
{
  cl_int errNum;
  cl_uint numPlatforms;
  cl_platform_id firstPlatformId;
  cl_context context = NULL;

  errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
  if (errNum != CL_SUCCESS || numPlatforms <= 0) {
      std::cerr << "Failed to find any OpenCL platforms." << std::endl;
      return NULL;
  }

  cl_context_properties contextProperties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)firstPlatformId, 0 };
  context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU, NULL, NULL, &errNum);

  if (errNum != CL_SUCCESS) { // No GPU Found
      std::cout << "Could not create GPU context, trying CPU..." << std::endl;
      context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU, NULL, NULL, &errNum);

      if (errNum != CL_SUCCESS) {
          std::cout << "Failed to create an OpenCL GPU or CPU context." << std::endl;
          return NULL;
      }
  }

  return context;
}

void showDeviceInfo(cl_device_id inDevice, cl_uint &ComputeUnits, bool PrintBool) {
  cl_int errNum;
  size_t InfoBufferSize = -1;

  if (PrintBool) {
    std::cout << std::endl;
    std::cout << "-------------------------- Selected Device Info --------------------------" << std::endl;
    // Vendor Name
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_VENDOR, 0, NULL, &InfoBufferSize);
    char* DeviceVendor = (char*)malloc(InfoBufferSize);
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_VENDOR, InfoBufferSize, DeviceVendor, &InfoBufferSize);
    std::cout << "Vendor Name: " << DeviceVendor << std::endl;
    free(DeviceVendor);

    // Device Name
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_NAME, 0, NULL, &InfoBufferSize);
    char* DeviceName = (char*)malloc(InfoBufferSize);
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_NAME, InfoBufferSize, DeviceName, &InfoBufferSize);
    std::cout << "Device Name: " << DeviceName << std::endl;
    free(DeviceName);

    // Device OpenCL Version
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_VERSION, 0, NULL, &InfoBufferSize);
    char* DeviceProfile = (char*)malloc(InfoBufferSize);
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_VERSION, InfoBufferSize, DeviceProfile, &InfoBufferSize);
    std::cout << "OpenCL Version: " << DeviceProfile << std::endl;
    free(DeviceProfile);

    // Global memory size
    cl_ulong memsize;

    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &memsize, NULL);
    std::cout << "Global Memory Size in Megabytes: " << memsize / 1000000 << std::endl;

    // Local memory size
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &memsize, NULL);
    std::cout << "Local Memory Size in Kilobytes: " << memsize / 1024 << std::endl;

    // Max Compute Units
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &ComputeUnits, NULL);
    std::cout << "Max Compute Units: " << ComputeUnits << std::endl;

    // Kernel Work Group Size - 1
		size_t MaxLocalSizeTotal;
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &MaxLocalSizeTotal, NULL);
    std::cout << "Kernel Work Group (Local) Total Size across Dimensions: " << MaxLocalSizeTotal << std::endl;

    cl_uint tempUint;
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &tempUint, NULL);
    std::cout << "Max Dimensions: " << tempUint << std::endl;

    // Kernel Work Group Size - 2
    size_t *MaxLocalSizeEach = new size_t[tempUint];
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*tempUint, MaxLocalSizeEach, NULL);
    for (int ii = 0; ii < (int)tempUint; ii++) {
        std::cout << "  -- Dimension " << ii << " Max Size: " << MaxLocalSizeEach[ii] << std::endl;
    }
    delete[] MaxLocalSizeEach;

    // Double Precision Configuration
    cl_device_fp_config fc;
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(cl_device_fp_config), &fc, NULL);
    std::cout << "Double Precision Configuration: " << std::endl;
		if (fc & CL_FP_DENORM) { std::cout << "  - (1) denorms are supported.\n"; }
		else { std::cout << "  - (1) denorms are NOT supported.\n"; }

		if (fc & CL_FP_INF_NAN) { std::cout << "  - (2) INF and NaNs are supported.\n"; }
		else { std::cout << "  - (2) INF and NaNs are NOT supported.\n"; }

		if (fc & CL_FP_ROUND_TO_NEAREST) { std::cout << "  - (3) round to nearest even rounding mode supported.\n"; }
		else { std::cout << "  - (3) round to nearest even rounding mode NOT supported.\n"; }

		if (fc & CL_FP_ROUND_TO_ZERO) { std::cout << "  - (4) round to zero rounding mode supported.\n"; }
		else { std::cout << "  - (4) round to zero rounding mode NOT supported.\n"; }

		if (fc & CL_FP_ROUND_TO_INF) { std::cout << "  - (5) round to +ve and -ve infinity rounding modes supported.\n"; }
		else { std::cout << "  - (5) round to +ve and -ve infinity rounding modes NOT supported.\n"; }

		if (fc & CL_FP_FMA) { std::cout << "  - (6) IEEE754 - 2008 fused multiply - add is supported.\n"; }
		else { std::cout << "  - (6) IEEE754-2008 fused multiply-add is NOT supported.\n"; }
		
		// Preferred vector types
		cl_uint vector_size;
		errNum = clGetDeviceInfo(inDevice, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint), &vector_size, NULL);
		std::cout << "Preferred FP64 Vector Size: " << vector_size << std::endl;

		errNum = clGetDeviceInfo(inDevice, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(cl_uint), &vector_size, NULL);
		std::cout << "Preferred FP32 Vector Size: " << vector_size << std::endl;

		// Native Instruction Set
		errNum = clGetDeviceInfo(inDevice, CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE, sizeof(cl_uint), &vector_size, NULL);
		std::cout << "Native Instruction Set FP64 Vector Size: " << vector_size << std::endl;

		errNum = clGetDeviceInfo(inDevice, CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT, sizeof(cl_uint), &vector_size, NULL);
		std::cout << "Native Instruction Set FP32 Vector Size: " << vector_size << std::endl;
		std::cout << "--------------------------------------------------------------------------\n\n";
  }
  else {
    // Max Compute Units
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &ComputeUnits, NULL);


    // Kernel Work Group Size - 1
		size_t MaxLocalSizeTotal;
    errNum = clGetDeviceInfo(inDevice, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &MaxLocalSizeTotal, NULL);
  }
}

cl_command_queue CreateCommandQueue(cl_context context, cl_device_id &device,
	int device_num) {
  cl_int errNum;
  cl_device_id *devices;
  cl_command_queue commandQueue = NULL;
  size_t deviceBufferSize = -1;

  errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);

  if (errNum != CL_SUCCESS) {
    PrintErrorCL(errNum);
    return NULL;
  }

  if (deviceBufferSize <= 0)
  {
    std::cerr << "No devices available." << std::endl;
    return NULL;
  }

  devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];

  errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);

  if (errNum != CL_SUCCESS)
  {
    std::cerr << "Failed to get device IDs" << std::endl;
    return NULL;
  }

  // WARNING
  // Nvidia GPUs only support OpenCL of version 1.2 or lower / This commandQueue is obsolete with OpenCL 2.0 or higher
  // USE clCreateCommandQueueWithProperties INSTEAD FOR NEWER OPENCL
  commandQueue = clCreateCommandQueueWithProperties(context, devices[device_num], 0, NULL);
  if (commandQueue == NULL) {
    delete[] devices;
    std::cerr << "Failed to create commandQueue for device 0" << std::endl;
    return NULL;
  }

  device = devices[device_num];

  delete[] devices;

  // false: do not print output;
	cl_uint ComputeUnits;
  showDeviceInfo(device, ComputeUnits, true);

  return commandQueue;
}

cl_program CreateProgram(cl_context context, cl_device_id device, 
	const char* fileName, const char* compileoption) {
  cl_int errNum;
  cl_program program;

  std::ifstream kernelFile(fileName, std::ios::in);
  if (!kernelFile.is_open()) {
    std::cerr << "`file for reading: " << fileName << std::endl;
    return NULL;
  }

  std::ostringstream oss;
  oss << kernelFile.rdbuf();

  std::string srcStdStr = oss.str();

  const char *srcStr = srcStdStr.c_str();

  program = clCreateProgramWithSource(context, 1, (const char**)&srcStr, NULL, NULL);

  if (program == NULL) {
    std::cerr << "Failed to create CL program from source." << std::endl;
    return NULL;
  }

  errNum = clBuildProgram(program, 0, NULL, compileoption, NULL, NULL);

  if (errNum != CL_SUCCESS) {
    char buildLog[16384];
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);

    std::cerr << "Error in kernel: " << std::endl;
    std::cerr << buildLog;
    clReleaseProgram(program);
    return NULL;
  }

  return program;
}

/*
	Fetch CPU
*/

cl_context CreateContext_CPU()
{
	cl_int errNum;
	cl_uint numPlatforms;
	cl_platform_id firstPlatformId;
	cl_context context = NULL;

	errNum = clGetPlatformIDs(2, &firstPlatformId, &numPlatforms);
	if (errNum != CL_SUCCESS || numPlatforms <= 0) {
		std::cerr << "Failed to find any OpenCL platforms." << std::endl;
		return NULL;
	}

	cl_context_properties contextProperties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)firstPlatformId, 0 };
	context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU, NULL, NULL, &errNum);

	if (errNum != CL_SUCCESS) {
		std::cout << "Could not create CPU context..." << std::endl;
	}

	return context;
}

#endif