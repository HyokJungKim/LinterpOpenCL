// Minimum size of grid points: 2
// Maximum number of points added over all dimensions: 256
// Maximum dimension: 6
// Local work size is fixed to 256
// Memory coalesced verion has no significant performance improvement
__kernel void linterpN(__global const double* in_grids,
    __global const double* in_o_grids, __global const double* in_grid_sizes,
     __global const int* in_idxlist, 
     __global const int* in_idx_info_i, __global const double* in_idx_info_d,
    const int Nd, const int totcol, const int workerN, 
    __global const double* fvals, __global double* outvec,
    __local double* o_grids, __local double* grid_sizes,
    __local int* idx_info_i, __local double* idx_info_d, __local int* idxlist) {
    
    // caution: definition of this variable changes after the local memory fence
    int tempid = get_local_id(0);

    // Allocate local memories (3 operations too much?)
    o_grids[tempid] = in_o_grids[tempid];
    grid_sizes[tempid] = in_grid_sizes[tempid];

    if (tempid < 64) {
        if (tempid < 32) {
            idx_info_i[tempid] = in_idx_info_i[tempid];
        }
        else {
            idx_info_d[tempid-32] = in_idx_info_d[tempid-32];
        }
    }
    
    // combination of grids to interpolate
    for (int ii = 0; ii < (totcol * Nd - 1) / 256 + 1; ii++) {
        idxlist[ii*256 + tempid] = in_idxlist[ii*256 + tempid];
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    /*
        Structure of idx_info_i (compare with linterpCPU.h)
        idx_info_i = [ N_indexer_f, N_indexer_s, Nd_sizes ]
        idx_info_d = [ lbb, ubb, scaler ]
    */
    int tempint[4]; // Maybe 6 is too much (lower this if you want)
    double gridvec[4]; // Maybe 6 is too much (lower this if you want)

    int tempidx_w, tempidx_f;
    double tempval, wght, returnval;
    
    for (int tt = 0; tt < workerN; tt++) {
        returnval = 0.0;
        tempid = get_global_id(0) * workerN + tt ; // global id!

        for (int jj = 0; jj < Nd; jj++) {
            gridvec[jj] = in_grids[tempid * Nd + jj];
            tempval = (gridvec[jj]-idx_info_d[jj]) / (idx_info_d[Nd+jj]-idx_info_d[jj]);
            tempval = pow(tempval, idx_info_d[jj+Nd*2]);
            tempval *= (double)(idx_info_i[Nd*2+jj]-1);
            tempint[jj] = min((int)(tempval + 1e-10), idx_info_i[Nd*2+jj] - 2);
            tempint[jj] = max(tempint[jj], 0);
        }

        for (int nn = 0; nn < totcol; nn++) {
            wght = 1.0;
            tempidx_f = 0;

            for (int kk = 0; kk < Nd; kk++) {
                tempidx_w = tempint[kk] + idxlist[nn*Nd + kk];
                tempval = gridvec[kk] - o_grids[idx_info_i[Nd + kk] + tempidx_w];
                tempval = fabs(tempval) / grid_sizes[idx_info_i[Nd + kk] + tempint[kk]];
                tempval = 1.0 - tempval; // Maybe clamp here

                wght *= tempval;
                tempidx_f += idx_info_i[kk] * tempidx_w;
            }
            
            returnval += wght * fvals[tempidx_f];
        }

        outvec[tempid] = returnval;
    }
}