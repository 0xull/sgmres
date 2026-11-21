#pragma once

#include "sgmres/types.hpp"

namespace sgmres {
    // One step of the Arnoldi process. On entry the first j + 1 columns of q are orthonormal. On
    // return column j of h holds the coefficients of A * q_j along them, h(j + 1, j) holds the norm
    // of what is left, and column j + 1 of q is that leftover normalized.
    // 
    // q needs at least j + 2 columns, h at least j + 2 rows and j + 1 columns.
    // 
    // Returns false when A * q_j already lies in span(q_0, ..., q_j) to working precision. The Krylov
    // space is then invariant, h(j + 1, j) is set to zero and column j + 1 of q is left untouched.
    bool arnoldi_step(const SparseMatrix& a, Matrix& q, Matrix& h, Index j);
} // namespace sgmres
