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
```cpp
    int Nd = 4; // Size of the dimension

    vd ubb {1.0, 1.0, 1.0, 1.0}; // Upper bounds
    vd lbb {-1.0, 0.5, 0.0, -0.5}; // Lower bounds
    vd scaler{0.4, 0.3, 0.8, 1.0}; // Scaler in the grid
    vi Ngrids {16, 16, 16, 16}; // Number of grids in each dimension as interpolation nodes
    vi Ninterp {16, 32, 16, 32}; // Number of points in each dimention to interpolate
    
    v2d o_grids(Nd); // Interpolation nodes
    v2d p_grids(Nd); // Number of points to evaluate

    for (int nn = 0; nn < Nd; nn++) {
        o_grids[nn] = linspacex(ubb[nn], lbb[nn], Ngrids[nn], scaler[nn]); // linspacex 
        p_grids[nn] = linspacex(ubb[nn], lbb[nn], Ninterp[nn], scaler[nn]);
    }
```
* The function *linspace* produces equally spaced grids between upper and lower bounds with desired number of points.
* The function *linspacex* uses scaled grids between upper and lower bounds where *i*-th grid point is computed by <img src="https://render.githubusercontent.com/render/math?math=l_i %2B \left( \frac{i-1}{N-1} \right)"> is converted as <img src="https://render.githubusercontent.com/render/math?math=m_i">
