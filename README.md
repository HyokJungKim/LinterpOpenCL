# LinterpOpenCL

This is a header file only package doing linear interpolation *designed for scientific applications* programmed with [OpenCL](https://www.khronos.org/opencl/). Let us start with a simple example!

**Step 1. ** Initialize OpenCL environment
'''cpp
    cl_context clContext = nullptr;
    cl_device_id clDevice = nullptr;
    cl_command_queue clQueue = nullptr;
    cl_program clProgram = nullptr;
'''
