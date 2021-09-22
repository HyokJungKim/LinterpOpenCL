//
// Created by hyokzzang on 10/27/20.
//

#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <ctime>
#include <string>

#ifndef LINTERPCL_ALLSETTINGS_H
#define LINTERPCL_ALLSETTINGS_H

typedef std::vector<int> vi;
typedef std::vector<vi> v2i;
typedef std::vector<v2i> v3i;
typedef std::vector<double> vd;
typedef std::vector<vd> v2d;
typedef std::vector<v2d> v3d;

//#define CL_TARGET_OPENCL_VERSION 210

std::string LINTERP_FILE = "/media/hyokzzang/linux_hdd/CloudStation/Codes/LinterpCL/linterp_ext.cl";
//std::string LINTERP_FILE = "C:/Users/hyokz/CloudStation/Codes/LinterpCL/";
//std::string LINTERP_FILE = "D:/CloudStation/Codes/LinterpCL/linterp.cl";

vd linspace(const double lb, const double ub, const int NN) {
    vd outvec(NN, 0.0);
    double divid = (ub - lb) / (static_cast<double>(NN) - 1.0);

    for (int nn = 0; nn < NN; nn++) {
        outvec[nn] = lb + static_cast<double>(nn) * divid;
    }

    return outvec;
}

vd linspacex(const double lb, const double ub, const int NN, const double scalex) {
    vd outvec(NN, 0.0);
    double divid = 1.0 / (static_cast<double>(NN) - 1.0);

    for (int nn = 0; nn < NN; nn++) {
        outvec[nn] = std::pow(static_cast<double>(nn) * divid, 1.0 / scalex) * (ub - lb) + lb;
    }

    return outvec;
}

template<class T> std::vector<std::vector<T>> cartesian(const std::vector<std::vector<T>>& in_vectors) {
    int total_cols = static_cast<int>(in_vectors.size());
    int total_rows = 1;
    vi Divisors(total_cols, 0);
    vi Indexors(total_cols, 1);

    for (int cc = 0; cc < total_cols; cc++) {
        Divisors[cc] = static_cast<int>(in_vectors[cc].size());
        total_rows *= Divisors[cc];

        if (cc > 0) {
            Indexors[cc] = Indexors[cc - 1] * Divisors[cc - 1];
        }
    }

    std::vector<std::vector<T>> out_vec(total_rows, std::vector<T>(total_cols));
    int tempidx;

    for (int rr = 0; rr < total_rows; rr++) {
        for (int cc = 0; cc < total_cols; cc++) {
            tempidx = (rr / Indexors[cc]) % Divisors[cc];
            out_vec[rr][cc] = in_vectors[cc][tempidx];
        }
    }

    return out_vec;
}

// Version without initializing the out_vec several times
template<class T> void cartesian(const std::vector<std::vector<T>>& in_vectors, 
  std::vector<std::vector<T>>& out_vec) {
  int total_cols = static_cast<int>(in_vectors.size());
  int total_rows = 1;
  vi Divisors(total_cols, 0);
  vi Indexors(total_cols, 1);

  for (int cc = 0; cc < total_cols; cc++) {
    Divisors[cc] = static_cast<int>(in_vectors[cc].size());
    total_rows *= Divisors[cc];

    if (cc > 0) {
      Indexors[cc] = Indexors[cc - 1] * Divisors[cc - 1];
    }
  }

  int tempidx;

  for (int rr = 0; rr < total_rows; rr++) {
    for (int cc = 0; cc < total_cols; cc++) {
      tempidx = (rr / Indexors[cc]) % Divisors[cc];
      out_vec[rr][cc] = in_vectors[cc][tempidx];
    }
  }
}

#endif //LINTERPCL_ALLSETTINGS_H