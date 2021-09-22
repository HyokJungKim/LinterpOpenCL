# LinterpOpenCL

This is a header file only package doing linear interpolation *designed for scientific applications* programmed with [OpenCL](https://www.khronos.org/opencl/) to utilize GPUs. Both CPU (without OpenCL) and GPU versions are available. Let us start with a simple example!

## Simple Example

**Step 0.** Include the necessary header
```cpp
#include "linterpGPU.h"

int main() {
    
    // All codes mentioned in the steps below goes here.
    
    return 0;
}
```

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
* The function *linspacex* uses scaled grids between upper (*ub*) and lower (*lb*) bounds where *i*-th grid point based on total of *N* number of points and scaler *s* is computed by <img src="https://render.githubusercontent.com/render/math?math=lb %2B (ub - lb) \left( \frac{i-1}{N-1} \right)^{1.0/s}">. Choosing *s* between zero and one makes the grids concentrated near zero, which is useful in some applications.
* *p_grids* and *Ninterp* are related to points to interpolate in this example. In actual application they do not have to be supplied in this way. (See **Step 4.** for details)

**Step 3.** Evaluate the function and create an *linterpGPU* object
```cpp
    // Cartesian product of grids at each dimension
    // Output is NN by Nd where NN is the product of all members in Ngrids
    v2d test_grids = cartesian<double>(o_grids);

    // Calculate grids based on user supplied function testfun
    // Output is a vector with NN elements
    vd test_f = testfun(test_grids);
    
    // Initialize linear interpolation object using GPU 
    linterpGPU testobj_GPU(test_f, o_grids, scaler, lbb, ubb, clContext, clQueue, clProgram);
```

**Step 4.** Interpolate the desired grids
```cpp
    v2d interp_grids = cartesian<double>(p_grids); // Make cartesian product of points to interpolate
    vd interp_grids_flat = flattenmat(interp_grids); // Flatten the matrix (Fortran ordered)
    vd resultvec_GPU = testobj_GPU.interpN(interp_grids_flat); // The resulting vector
```
Of course, the grids to interpolate does not have to be supplied by *catersian* or *flattenmat* function. Any flattened vector which has dimension of *NN* (total number of grids to interpolate) times *Nd* (size of the dimension) inside the lower and upper bounds specified before should be okay.

## To-do List
* Coming soon!
