//
// CPU version of the linear interpolation
//

#ifndef LINTERPCL_LINTERPCL_H
#define LINTERPCL_LINTERPCL_H
#endif //LINTERPCL_LINTERPCL_H

#include "AllSettings.h"

void makeidxlist(const int totcol, const int Nd, const vi& tempi, v2i& idxlist) {
    for (int ii = 0; ii < totcol; ii++) {
        for (int jj = 0; jj < Nd; jj++) {
            idxlist[ii][jj] = (ii / tempi[jj]) % 2;
        }
    }
}

/*
 * CPU version for reference
 */
// Consider using the outer product function I made before...
class linterpCPU {
private:
    int Nd, Ntot, totcol;
    vi N_indexer, Nd_sizes;
    v2i idxlist;
    vd fvals, scaler, ubb, lbb; // function values
    v2d o_grids, grid_sizes; // original grid points

public:
    linterpCPU(const vd& in_fvals, const v2d& in_o_grids, const vd& in_scaler,
               const vd& in_lbb, const vd& in_ubb) {
        fvals = in_fvals;
        o_grids = in_o_grids;
        Nd = static_cast<int>(in_o_grids.size());
        Nd_sizes = vi(Nd);
        Ntot = 1;
        totcol = 1;
        grid_sizes = v2d(Nd);
        vi tempi(Nd);
        N_indexer = vi(Nd);

        for (int ii = 0; ii < Nd; ii++) {
            Nd_sizes[ii] = static_cast<int>(o_grids[ii].size());
            tempi[ii] = totcol;
            N_indexer[ii] = Ntot;
            Ntot *= Nd_sizes[ii];
            totcol *= 2;
            grid_sizes[ii] = vd(Nd_sizes[ii] - 1);

            for(int jj = 0; jj < Nd_sizes[ii] - 1; jj++) {
                grid_sizes[ii][jj] = std::fabs(in_o_grids[ii][jj+1] - in_o_grids[ii][jj]);
            }
        }
        idxlist = v2i(totcol, vi(Nd));
        makeidxlist(totcol, Nd, tempi, idxlist);
        scaler = in_scaler;
        ubb = in_ubb;
        lbb = in_lbb;
    };

    v2i grid2idx(const v2d& in_grids);
    vd interpN(const v2d& in_grids);
};

v2i linterpCPU::grid2idx(const v2d& in_grids) {
    int TT = static_cast<int>(in_grids.size());
    v2i outMat(TT, vi(Nd));
    double tempval;

    for (int ii = 0; ii < TT; ii++) {
        for (int nn = 0; nn < Nd; nn++) {
            tempval = std::pow((in_grids[ii][nn] - lbb[nn]) / (ubb[nn] - lbb[nn]), scaler[nn]);
            tempval *= static_cast<double>(Nd_sizes[nn] - 1);
            outMat[ii][nn] = std::min(static_cast<int>(tempval + 1e-10), Nd_sizes[nn] - 2);
        }
    }
    return outMat;
};

vd linterpCPU::interpN(const v2d &in_grids) {
    int tempidx_f, tempidx_w, tempi, TT = static_cast<int>(in_grids.size());
    v2i search_idx = grid2idx(in_grids);
    double wght, tempval;
    vd outMat(TT, 0.0);

    for (int tt = 0; tt < TT; tt++) {
        for (int nn = 0; nn < totcol; nn++) {
            tempidx_f = 0;
            wght = 1.0;

            for (int kk = 0; kk < Nd; kk++) {
                tempi = search_idx[tt][kk];
                tempidx_w = tempi + idxlist[nn][kk];
                tempval = 1.0 - std::abs(in_grids[tt][kk] - o_grids[kk][tempidx_w]) / grid_sizes[kk][tempi];
                wght *= tempval;
                tempidx_f += N_indexer[kk] * tempidx_w;
            }
            outMat[tt] += wght * fvals[tempidx_f];
        }
    }
    return outMat;
};