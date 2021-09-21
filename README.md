# LinterpOpenCL

This is a header file only package doing linear interpolation *designed for scientific applications* programmed with [OpenCL](https://www.khronos.org/opencl/) to utilize GPUs. Both CPU (without OpenCL) and GPU versions are available. Let us start with a simple example!



**Step 1.** Initialize OpenCL environment
```cpp
    cl_context clContext = nullptr;
    cl_device_id clDevice = nullptr;
    cl_command_queue clQueue = nullptr;
    cl_program clProgram = nullptr;
    
    initializeCL(clContext, clDevice, clQueue, clProgram, LINTERP_FILE); // LINTERP_FILE: path to `linterp_ext.cl'
```

**Step 2.** Setup the problem
