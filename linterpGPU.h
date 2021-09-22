//
// Created by hyokzzang on 10/27/20.
//

#include "linterpCPU.h"
#include "OpenCLtools.h"

#ifndef LINTERPCL_LINTERPGPU_H
#define LINTERPCL_LINTERPGPU_H

void initializeCL(cl_context &clContext, cl_device_id& clDevice,
                  cl_command_queue& clQueue, cl_program& clProgram,
                  const std::string& in_string) {

    clContext = CreateContext();
    clQueue = CreateCommandQueue(clContext, clDevice, 0);

    //std::string compileoption = "-cl-fast-relaxed-math";
    //std::string compileoption = "-w -g";
    std::string compileoption = "";
    clProgram = CreateProgram(clContext, clDevice, in_string.c_str(), compileoption.c_str());
};

void makeidxlist(const int totcol, const int Nd, const vi& tempi, vi& idxlist) {
    for (int ii = 0; ii < totcol; ii++) {
        for (int jj = 0; jj < Nd; jj++) {
            idxlist[ii * Nd + jj] = (ii / tempi[jj]) % 2;
        }
    }
};

template <class T> T multiples(const T in_number, const T in_mult) {
    return in_mult * ((in_number - 1) / in_mult + 1);
}

template <class T> std::vector<T> flattenmat(const std::vector<std::vector<T>>& in_vectors) {
    int Nrow = static_cast<int>(in_vectors.size());
    int Ncol = static_cast<int>(in_vectors[0].size());

    std::vector<T> outvec(Nrow*Ncol);

    for (int ii = 0; ii < Nrow; ii++) {
        for (int jj = 0; jj < Ncol; jj++) {
            outvec[ii * Ncol + jj] = in_vectors[ii][jj];
        }
    }
    return outvec;
}

template <class T> std::vector<T> flattenmat_fortran(const std::vector<std::vector<T>>& in_vectors) {
    int Nrow = static_cast<int>(in_vectors.size());
    int Ncol = static_cast<int>(in_vectors[0].size());

    std::vector<T> outvec(Nrow*Ncol);

    for (int ii = 0; ii < Nrow; ii++) {
        for (int jj = 0; jj < Ncol; jj++) {
            outvec[jj * Nrow + ii] = in_vectors[ii][jj];
        }
    }
    return outvec;
}

class linterpGPU {
private:
    int Nd = 1, totcol = 1; // default

    cl_context clContext = nullptr;
    cl_command_queue clQueue = nullptr;
    cl_program clProgram = nullptr;
    cl_kernel clKernel = nullptr;

    cl_mem mem_grids = nullptr, mem_o_grids = nullptr, mem_grid_sizes = nullptr, mem_idx_info_i = nullptr;
    cl_mem mem_idx_info_d = nullptr, mem_fvals = nullptr, mem_outvec = nullptr, mem_idxlist = nullptr;
    vi idx_info_i, idxlist;
    vd o_grids, grid_sizes, idx_info_d, fvals;

    size_t global_t[1]{0};
    size_t local_t[1] = {256};

    cl_int errNum = 0; // default

    // How many grid points each worker processes?
    // Small: many cores are working at the same time but many repetitive initializations of local memory
    // Large: pros and cons are opposite of small
    int workerN = 16;

public:
    cl_mem mem_outvec_keep = nullptr, mem_grids_keep = nullptr;
    vd lbb, ubb;

    void setCL(cl_context& in_clContext, cl_command_queue& in_clQueue, cl_program& in_clProgram) {
        clContext = in_clContext;
        clQueue = in_clQueue;
        clProgram = in_clProgram;
        clKernel = clCreateKernel(clProgram, "linterpN", nullptr);
    }
    
    linterpGPU() = default;

    linterpGPU(const vd& in_fvals, const v2d& in_o_grids, const vd& in_scaler,
               const vd& in_lbb, const vd& in_ubb,
               cl_context& in_clContext, cl_command_queue& in_clQueue, cl_program& in_clProgram) {
        // Make kernel and set OpenCL variables
        setCL(in_clContext, in_clQueue, in_clProgram);

        o_grids = vd(256, 0.0);
        grid_sizes = vd(256, 0.0);
        idx_info_i = vi(32, 0); // only up to 64 allowed?
        idx_info_d = vd(32, 0); // only up to 64 allowed?

        fvals = in_fvals;

        lbb = in_lbb;
        ubb = in_ubb;

        Nd = static_cast<int>(in_o_grids.size());

        int Ntot = 1;
        int Ncount = 0;
        vi tempi(Nd);
        totcol = 1;

        for (int ii = 0; ii < Nd; ii++) {
            idx_info_i[ii] = Ntot;
            idx_info_i[Nd + ii] = Ncount;
            idx_info_i[Nd * 2 + ii] = static_cast<int>(in_o_grids[ii].size());
            tempi[ii] = totcol;

            Ntot *= idx_info_i[Nd * 2 + ii];
            Ncount += idx_info_i[Nd * 2 + ii];
            totcol *= 2;

            idx_info_d[ii] = lbb[ii];
            idx_info_d[Nd + ii] = ubb[ii];
            idx_info_d[2*Nd + ii] = in_scaler[ii];

            for (int jj = 0; jj < idx_info_i[Nd * 2 + ii]; jj++) {
                o_grids[idx_info_i[Nd + ii] + jj] = in_o_grids[ii][jj];
            }

            for (int jj = 0; jj < idx_info_i[Nd * 2 + ii] - 1; jj++) {
                grid_sizes[idx_info_i[Nd + ii] + jj] = std::abs(in_o_grids[ii][jj+1] - in_o_grids[ii][jj]);
            }
        }

        int mult = multiples<int>(totcol * Nd, 256);
        idxlist = vi(mult, 0);
        makeidxlist(totcol, Nd, tempi, idxlist);

        // Assign OpenCL memory objects
        mem_o_grids = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(double) * o_grids.size(),
                                        o_grids.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_grid_sizes = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     sizeof(double) * grid_sizes.size(),
                                        grid_sizes.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idxlist = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     sizeof(int) * idxlist.size(),
                                     idxlist.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idx_info_i = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(int) * idx_info_i.size(),
                                        idx_info_i.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idx_info_d = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(double) * idx_info_d.size(),
                                        idx_info_d.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_fvals = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   sizeof(double) * fvals.size(),
                                   fvals.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        // Set Kernels
        errNum = clSetKernelArg(clKernel, 1, sizeof(cl_mem), &mem_o_grids);
        errNum = clSetKernelArg(clKernel, 2, sizeof(cl_mem), &mem_grid_sizes);

        errNum &= clSetKernelArg(clKernel, 3, sizeof(cl_mem), &mem_idxlist);
        errNum &= clSetKernelArg(clKernel, 4, sizeof(cl_mem), &mem_idx_info_i);
        errNum &= clSetKernelArg(clKernel, 5, sizeof(cl_mem), &mem_idx_info_d);
        errNum &= clSetKernelArg(clKernel, 6, sizeof(int), (void*)&Nd);
        errNum &= clSetKernelArg(clKernel, 7, sizeof(int), (void*)&totcol);
        errNum &= clSetKernelArg(clKernel, 8, sizeof(int), (void*)&workerN);
        errNum &= clSetKernelArg(clKernel, 9, sizeof(cl_mem), &mem_fvals);

        // Local memory
        errNum &= clSetKernelArg(clKernel, 11, sizeof(cl_double) * 256, nullptr);
        errNum &= clSetKernelArg(clKernel, 12, sizeof(cl_double) * 256, nullptr);
        errNum &= clSetKernelArg(clKernel, 13, sizeof(cl_int) * 32, nullptr);
        errNum &= clSetKernelArg(clKernel, 14, sizeof(cl_double) * 32, nullptr);
        errNum &= clSetKernelArg(clKernel, 15, sizeof(cl_int) * mult, nullptr);
        if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }
    }
    
    vd interpN(vd& in_grids);
    void interpN(vd& out_rslt, vd& in_grids);
    void interpN_keep(vd& in_grids);
    void update_fvals(vd& in_fvals);
    void clearGPU();
};

vd linterpGPU::interpN(vd &in_grids) {
    global_t[0] = multiples<size_t>(in_grids.size() / static_cast<size_t>(Nd),
                                    256 * static_cast<size_t>(workerN));
    in_grids.resize(global_t[0] * static_cast<size_t>(Nd), 0.0); // Make it a multiple of 256
    vd outvec(global_t[0]);
    global_t[0] /= static_cast<size_t>(workerN);

    // Assign OpenCL memory objects
    mem_grids = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_grids.size(),
                               in_grids.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    mem_outvec = clCreateBuffer(clContext, CL_MEM_WRITE_ONLY,
                                sizeof(double) * outvec.size(), nullptr, &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    // Set Kernels
    errNum = clSetKernelArg(clKernel, 0, sizeof(cl_mem), &mem_grids);
    errNum |= clSetKernelArg(clKernel, 10, sizeof(cl_mem), &mem_outvec);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }

    // Run kernel
    errNum = clEnqueueNDRangeKernel(clQueue, clKernel, 1, nullptr,
                                    global_t, local_t, 0,
                                    nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error running clKernel. \n"; PrintErrorCL(errNum); }

    // Read data from CPU
    errNum = clEnqueueReadBuffer(clQueue, mem_outvec, CL_TRUE, 0, sizeof(double) * outvec.size(),
                                 outvec.data(), 0, nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error reading buffer. \n"; PrintErrorCL(errNum); }

    clReleaseMemObject(mem_grids);
    clReleaseMemObject(mem_outvec);

    return outvec;
};

void linterpGPU::interpN(vd& outvec, vd &in_grids) {
    global_t[0] = multiples<size_t>(in_grids.size() / static_cast<size_t>(Nd),
                                    256 * static_cast<size_t>(workerN));
    in_grids.resize(global_t[0] * static_cast<size_t>(Nd), 0.0); // Make it a multiple of 256
    global_t[0] /= static_cast<size_t>(workerN);

    // Assign OpenCL memory objects
    mem_grids = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_grids.size(),
                               in_grids.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    mem_outvec = clCreateBuffer(clContext, CL_MEM_WRITE_ONLY,
                                sizeof(double) * outvec.size(), nullptr, &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    // Set Kernels
    errNum = clSetKernelArg(clKernel, 0, sizeof(cl_mem), &mem_grids);
    errNum |= clSetKernelArg(clKernel, 10, sizeof(cl_mem), &mem_outvec);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }

    // Run kernel
    errNum = clEnqueueNDRangeKernel(clQueue, clKernel, 1, nullptr,
                                    global_t, local_t, 0,
                                    nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error running clKernel. \n"; PrintErrorCL(errNum); }

    // Read data from CPU
    errNum = clEnqueueReadBuffer(clQueue, mem_outvec, CL_TRUE, 0, sizeof(double) * outvec.size(),
                                 outvec.data(), 0, nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error reading buffer. \n"; PrintErrorCL(errNum); }

    clReleaseMemObject(mem_grids);
    clReleaseMemObject(mem_outvec);
};

void linterpGPU::interpN_keep(vd &in_grids) {
    global_t[0] = multiples<size_t>(in_grids.size() / static_cast<size_t>(Nd),
                                    256 * static_cast<size_t>(workerN));
    in_grids.resize(global_t[0] * static_cast<size_t>(Nd), 0.0); // Make it a multiple of 256
    vd outvec(global_t[0]);
    global_t[0] /= static_cast<size_t>(workerN);

    // Assign OpenCL memory objects
    mem_grids_keep = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_grids.size(),
                               in_grids.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    mem_outvec_keep = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS,
                                sizeof(double) * outvec.size(), nullptr, &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    // Set Kernels
    errNum = clSetKernelArg(clKernel, 0, sizeof(cl_mem), &mem_grids_keep);
    errNum |= clSetKernelArg(clKernel, 10, sizeof(cl_mem), &mem_outvec_keep);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }

    // Run kernel
    errNum = clEnqueueNDRangeKernel(clQueue, clKernel, 1, nullptr,
                                    global_t, local_t, 0,
                                    nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error running clKernel. \n"; PrintErrorCL(errNum); }
};

void linterpGPU::update_fvals(vd& in_fvals) {
    clReleaseMemObject(mem_fvals);

    mem_fvals = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_fvals.size(), in_fvals.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting buffer. \n"; PrintErrorCL(errNum); }

    errNum = clSetKernelArg(clKernel, 9, sizeof(cl_mem), &mem_fvals);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting Kernel arguments. \n"; PrintErrorCL(errNum); }
}

void linterpGPU::clearGPU() {
    errNum = clReleaseMemObject(mem_o_grids);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_o_grids\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_grid_sizes);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_grid_sizes\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idx_info_i);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idx_info_i\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idx_info_d);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idx_info_d\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_fvals);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_fvals\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idxlist);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idxlist\n"; PrintErrorCL(errNum); }

    errNum = clReleaseKernel(clKernel);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. Kernel\n"; PrintErrorCL(errNum); }
}

/*
    Linterp Pre-set the Size of the Output matrix (GPU)
*/
class linterpGPU_preset {
private:
    int Nd = 1, totcol = 1; // default

    cl_context clContext = nullptr;
    cl_command_queue clQueue = nullptr;
    cl_program clProgram = nullptr;
    cl_kernel clKernel = nullptr;

    cl_mem mem_grids = nullptr, mem_o_grids = nullptr, mem_grid_sizes = nullptr, mem_idx_info_i = nullptr;
    cl_mem mem_idx_info_d = nullptr, mem_fvals = nullptr, mem_idxlist = nullptr;
    vi idx_info_i, idxlist;
    vd o_grids, grid_sizes, idx_info_d;

    size_t global_t[1]{0};
    size_t local_t[1] = {256};

    cl_int errNum = 0; // default

    int outvec_size;

    // How many grid points each worker processes?
    // Small: many cores are working at the same time but many repetitive initializations of local memory
    // Large: pros and cons are opposite of small
    int workerN = 16;

public:
    cl_mem mem_outvec = nullptr, mem_outvec_keep = nullptr, mem_grids_keep = nullptr;
    vd outvec, fvals, lbb, ubb;

    void setCL(cl_context& in_clContext, cl_command_queue& in_clQueue, cl_program& in_clProgram) {
        clContext = in_clContext;
        clQueue = in_clQueue;
        clProgram = in_clProgram;
        clKernel = clCreateKernel(clProgram, "linterpN", nullptr);
    }

    linterpGPU_preset() = default;

    linterpGPU_preset& operator = (const linterpGPU_preset& inObj);

    linterpGPU_preset(const vd& in_fvals, const v2d& in_o_grids, const vd& in_scaler,
                      const vd& in_lbb, const vd& in_ubb,
                      cl_context& in_clContext, cl_command_queue& in_clQueue, cl_program& in_clProgram,
                      const int in_outvec_size) {
        // Make kernel and set OpenCL variables
        setCL(in_clContext, in_clQueue, in_clProgram);

        o_grids = vd(256, 0.0);
        grid_sizes = vd(256, 0.0);
        idx_info_i = vi(32, 0);
        idx_info_d = vd(32, 0);

        fvals = in_fvals;

        outvec_size = in_outvec_size;
        outvec = vd(outvec_size);

        lbb = in_lbb;
        ubb = in_ubb;

        Nd = static_cast<int>(in_o_grids.size());

        int Ntot = 1;
        int Ncount = 0;
        vi tempi(Nd);
        totcol = 1;

        for (int ii = 0; ii < Nd; ii++) {
            idx_info_i[ii] = Ntot;
            idx_info_i[Nd + ii] = Ncount;
            idx_info_i[Nd * 2 + ii] = static_cast<int>(in_o_grids[ii].size());
            tempi[ii] = totcol;

            Ntot *= idx_info_i[Nd * 2 + ii];
            Ncount += idx_info_i[Nd * 2 + ii];
            totcol *= 2;

            idx_info_d[ii] = lbb[ii];
            idx_info_d[Nd + ii] = ubb[ii];
            idx_info_d[Nd*2 + ii] = in_scaler[ii];

            for (int jj = 0; jj < idx_info_i[Nd * 2 + ii]; jj++) {
                o_grids[idx_info_i[Nd + ii] + jj] = in_o_grids[ii][jj];
            }

            for (int jj = 0; jj < idx_info_i[Nd * 2 + ii] - 1; jj++) {
                grid_sizes[idx_info_i[Nd + ii] + jj] = std::abs(in_o_grids[ii][jj+1] - in_o_grids[ii][jj]);
            }
        }

        int mult = multiples<int>(totcol * Nd, 256);
        idxlist = vi(mult, 0);
        makeidxlist(totcol, Nd, tempi, idxlist);

        // Assign OpenCL memory objects
        mem_o_grids = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(double) * o_grids.size(),
                                     o_grids.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_grid_sizes = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     sizeof(double) * grid_sizes.size(),
                                        grid_sizes.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idxlist = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                     sizeof(int) * idxlist.size(),
                                     idxlist.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idx_info_i = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(int) * idx_info_i.size(),
                                        idx_info_i.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_idx_info_d = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                        sizeof(double) * idx_info_d.size(),
                                        idx_info_d.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_fvals = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   sizeof(double) * fvals.size(),
                                   fvals.data(), &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        mem_outvec = clCreateBuffer(clContext, CL_MEM_WRITE_ONLY,
                                sizeof(double) * outvec_size, nullptr, &errNum);
        if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

        // Set Kernels
        errNum = clSetKernelArg(clKernel, 1, sizeof(cl_mem), &mem_o_grids);
        errNum = clSetKernelArg(clKernel, 2, sizeof(cl_mem), &mem_grid_sizes);

        errNum &= clSetKernelArg(clKernel, 3, sizeof(cl_mem), &mem_idxlist);
        errNum &= clSetKernelArg(clKernel, 4, sizeof(cl_mem), &mem_idx_info_i);
        errNum &= clSetKernelArg(clKernel, 5, sizeof(cl_mem), &mem_idx_info_d);
        errNum &= clSetKernelArg(clKernel, 6, sizeof(int), (void*)&Nd);
        errNum &= clSetKernelArg(clKernel, 7, sizeof(int), (void*)&totcol);
        errNum &= clSetKernelArg(clKernel, 8, sizeof(int), (void*)&workerN);
        errNum &= clSetKernelArg(clKernel, 9, sizeof(cl_mem), &mem_fvals);
        errNum &= clSetKernelArg(clKernel, 10, sizeof(cl_mem), &mem_outvec);
        // Local memory
        errNum &= clSetKernelArg(clKernel, 11, sizeof(cl_double) * 256, nullptr);
        errNum &= clSetKernelArg(clKernel, 12, sizeof(cl_double) * 256, nullptr);
        errNum &= clSetKernelArg(clKernel, 13, sizeof(cl_int) * 32, nullptr);
        errNum &= clSetKernelArg(clKernel, 14, sizeof(cl_double) * 32, nullptr);
        errNum &= clSetKernelArg(clKernel, 15, sizeof(cl_int) * mult, nullptr);
        if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }
    }
    
    vd interpN(vd& in_grids);
    //void interpN_keep(vd& in_grids);
    void clearGPU();
    void update_fvals(vd&);
};

linterpGPU_preset& linterpGPU_preset::operator = (const linterpGPU_preset& inObj) {
    Nd = inObj.Nd;
    totcol = inObj.totcol;

    clContext = inObj.clContext;
    clQueue = inObj.clQueue;
    clProgram = inObj.clProgram;
    clKernel = inObj.clKernel;

    mem_grids = inObj.mem_grids;
    mem_o_grids = inObj.mem_o_grids;
    mem_grid_sizes = inObj.mem_grid_sizes;
    mem_idx_info_i = inObj.mem_idx_info_i;
    mem_idx_info_d = inObj.mem_idx_info_d;
    mem_fvals = inObj.mem_fvals;
    mem_idxlist = inObj.mem_idxlist;

    idx_info_i = inObj.idx_info_i;
    idxlist = inObj.idxlist;

    o_grids = inObj.o_grids;
    grid_sizes = inObj.grid_sizes;
    idx_info_d = inObj.idx_info_d;
    fvals = inObj.fvals;

    global_t[0] = inObj.global_t[0];
    local_t[0] = inObj.local_t[0];

    lbb = inObj.lbb;
    ubb = inObj.ubb;

    errNum = 0; // default

    outvec_size = inObj.outvec_size;

    // How many grid points each worker processes?
    // Small: many cores are working at the same time but many repetitive initializations of local memory
    // Large: pros and cons are opposite of small
    workerN = inObj.workerN;

    mem_outvec = inObj.mem_outvec;
    mem_outvec_keep = inObj.mem_outvec_keep;
    mem_grids_keep = inObj.mem_grids_keep;
    outvec = inObj.outvec;

    return *this;
}

void linterpGPU_preset::update_fvals(vd& in_fvals) {
    fvals = in_fvals;

    clReleaseMemObject(mem_fvals);

    mem_fvals = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_fvals.size(),
                               in_fvals.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }

    errNum = clSetKernelArg(clKernel, 9, sizeof(cl_mem), &mem_fvals);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }
}

vd linterpGPU_preset::interpN(vd &in_grids) {
    global_t[0] = multiples<size_t>(in_grids.size() / static_cast<size_t>(Nd),
                                    256 * static_cast<size_t>(workerN));
    in_grids.resize(global_t[0] * static_cast<size_t>(Nd), 0.0); // Make it a multiple of 256

    global_t[0] /= static_cast<size_t>(workerN);

    // Assign OpenCL memory objects
    mem_grids = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_grids.size(),
                               in_grids.data(), &errNum);
    if (errNum != CL_SUCCESS) { std::cerr << "Error creating buffer. \n"; PrintErrorCL(errNum); }
                               
    // Set Kernels
    errNum = clSetKernelArg(clKernel, 0, sizeof(cl_mem), &mem_grids);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }

    // Run kernel
    errNum = clEnqueueNDRangeKernel(clQueue, clKernel, 1, nullptr,
                                    global_t, local_t, 0,
                                    nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error running clKernel. \n"; PrintErrorCL(errNum); }

    // Read data from CPU
    errNum = clEnqueueReadBuffer(clQueue, mem_outvec, CL_TRUE, 0, sizeof(double) * outvec.size(),
                                 outvec.data(), 0, nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error reading buffer. \n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_grids);
    if (errNum != CL_SUCCESS) { std::cerr << "Error releasing memory. \n"; PrintErrorCL(errNum); }


    return outvec;
};

/*
void linterpGPU_preset::interpN_keep(vd &in_grids) {
    global_t[0] = multiples<size_t>(in_grids.size() / static_cast<size_t>(Nd),
                                    256 * static_cast<size_t>(workerN));
    in_grids.resize(global_t[0] * static_cast<size_t>(Nd), 0.0); // Make it a multiple of 256
    vd outvec_pre(global_t[0]);
    global_t[0] /= static_cast<size_t>(workerN);

    // Assign OpenCL memory objects
    mem_grids_keep = clCreateBuffer(clContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                               sizeof(double) * in_grids.size(),
                               in_grids.data(), nullptr);

    mem_outvec_keep = clCreateBuffer(clContext, CL_MEM_READ_WRITE | CL_MEM_HOST_NO_ACCESS,
                                sizeof(double) * outvec.size(), nullptr, nullptr);

    // Set Kernels
    errNum = clSetKernelArg(clKernel, 0, sizeof(cl_mem), &mem_grids_keep);
    errNum |= clSetKernelArg(clKernel, 9, sizeof(cl_mem), &mem_outvec_keep);
    if (errNum != CL_SUCCESS) { std::cerr << "Error setting kernel arguments. \n"; PrintErrorCL(errNum); }

    // Run kernel
    errNum = clEnqueueNDRangeKernel(clQueue, clKernel, 1, nullptr,
                                    global_t, local_t, 0,
                                    nullptr, nullptr);
    if (errNum != CL_SUCCESS) { std::cerr << "Error running clKernel. \n"; PrintErrorCL(errNum); }
};
 */

void linterpGPU_preset::clearGPU() {
    errNum = clReleaseMemObject(mem_o_grids);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_o_grids\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_grid_sizes);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_grid_sizes\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idx_info_i);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idx_info_i\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idx_info_d);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idx_info_d\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_fvals);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_fvals\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_idxlist);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_idxlist\n"; PrintErrorCL(errNum); }

    errNum = clReleaseMemObject(mem_outvec);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. mem_outvec\n"; PrintErrorCL(errNum); }

    errNum = clReleaseKernel(clKernel);
    if (errNum != CL_SUCCESS) { std::cerr << "Release failure. clKernel\n"; PrintErrorCL(errNum); }
}

#endif //LINTERPCL_LINTERPGPU_H
