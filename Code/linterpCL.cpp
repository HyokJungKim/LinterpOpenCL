//#include "linterpCPU.h"
#include "linterpGPU.h"
//#include "OpenCLtools.h"

vd testfun(const v2d& grids) {
    int NN = static_cast<int>(grids.size());
    vd outvec(NN, 0.0);

    for (int nn = 0; nn < NN; nn++) {
        outvec[nn] = std::exp(grids[nn][0]) + std::sqrt(grids[nn][1] + grids[nn][2]) * grids[nn][3];
    }

    return outvec;
}

vd testfun2(const v2d& grids) {
    int NN = static_cast<int>(grids.size());
    vd outvec(NN, 0.0);

    for (int nn = 0; nn < NN; nn++) {
        outvec[nn] = 0.5* std::exp(grids[nn][0]) + std::sqrt(grids[nn][1] + grids[nn][2]) * grids[nn][3];
    }

    return outvec;
}

int main() {
    /*
     * OpenCL settings
    */
    cl_context clContext = nullptr;
    cl_device_id clDevice = nullptr;
    cl_command_queue clQueue = nullptr;
    cl_program clProgram = nullptr;

    initializeCL(clContext, clDevice, clQueue, clProgram, LINTERP_FILE);

    int Nd = 4; // Size of the dimension

    vd ubb {1.0, 1.0, 1.0, 1.0}; // Upper bounds
    vd lbb {-0.5, 0.0, 0.0, 0.0}; // Lower bounds
    vd scaler{0.4, 0.3, 0.8, 1.0}; // Scaler in the grid
    vi Ngrids {16, 16, 16, 16}; // Number of grids in each dimension as interpolation nodes
    vi Ninterp {32, 32, 32, 32}; // Number of points in each dimention to evaluate

    v2d o_grids(Nd); // Interpolation nodes
    v2d p_grids(Nd); // Number of points to evaluate

    for (int nn = 0; nn < Nd; nn++) {
        o_grids[nn] = linspacex(lbb[nn], ubb[nn], Ngrids[nn], scaler[nn]);
        p_grids[nn] = linspacex(lbb[nn], ubb[nn], Ninterp[nn], scaler[nn]);
    }

    // Cartesian product of grids at each dimension
    // Output is NN by Nd where NN is the product of all members in Ngrids
    v2d test_grids = cartesian<double>(o_grids);

    // Calculate grids based on user supplied function testfun
    // Output is a vector with NN elements
    vd test_f = testfun(test_grids);

    linterpCPU testobj(test_f, o_grids, scaler, lbb, ubb);

    v2d interp_grids = cartesian<double>(p_grids);
    vd interp_grids_flat = flattenmat(interp_grids); // Fortran ordered matrix

    clock_t t1, t2;
    t1 = clock();
    vd resultvec = testobj.interpN(interp_grids);
    t2 = clock();

    printf("CPU: multilinear interp: %d clocks, %f sec\n",
           static_cast<int>(t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

    printf("\n\n\n");

    // Initialize linear interpolation object using GPU
    linterpGPU testobj_GPU(test_f, o_grids, scaler, lbb, ubb, clContext, clQueue, clProgram);

    t1 = clock();
    vd resultvec_GPU;

    for (int ii = 0; ii < 1; ii++) {
        resultvec_GPU = testobj_GPU.interpN(interp_grids_flat);
    }
    t2 = clock();

    printf("GPU: multilinear interp: %d clocks, %f sec\n",
           static_cast<int>(t2-t1), ((double)(t2 - t1)) / CLOCKS_PER_SEC);

    double error_size = 0.0;
    for (int ii = 0; ii < static_cast<int>(resultvec.size()); ii++) {
        error_size += std::fabs(resultvec[ii] - resultvec_GPU[ii]);

        if (std::isnan(resultvec_GPU[ii])) {
          printf("%d: %f %f %f %f\n", ii, interp_grids[ii][0], interp_grids[ii][1],
            interp_grids[ii][2], interp_grids[ii][3]);
        }
    }
    std::cout << "Error Size (CPU vs. GPU): " << error_size << "\n";

    vd test_f2 = testfun(interp_grids);
    auto divis = static_cast<double>(resultvec.size());
    error_size = 0.0;
    for (int ii = 0; ii < static_cast<int>(resultvec.size()); ii++) {
        error_size += std::fabs((resultvec[ii] - test_f2[ii]) / test_f2[ii]) / divis;
    }
    std::cout << "Error Size (CPU vs. True): " << error_size << "\n";

    vd test_f3 = testfun2(test_grids);

    testobj_GPU.update_fvals(test_f3);

    vd resultvec_GPU2 = testobj_GPU.interpN(interp_grids_flat);

    testobj_GPU.clearGPU();

    clReleaseProgram(clProgram);
    clReleaseCommandQueue(clQueue);
    clReleaseDevice(clDevice);
    clReleaseContext(clContext);


    return 0;
}
