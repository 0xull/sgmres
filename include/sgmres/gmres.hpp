#pragma once

#include "sgmres/solve_result.hpp"
#include "sgmres/types.hpp"

namespace sgmres {

struct GmresOptions {
  Index restart = 50;  // basis vectors per cycle
  Index max_cycles = 200;
  double tol = 1e-8;  // stop once ||f - Ax|| <= tol * ||f||
};

// Restarted GMRES with full Arnoldi. This is the baseline sGMRES gets measured against.
SolveResult gmres(const SparseMatrix& a, const Vector& f, const Vector& x0,
                  const GmresOptions& options = {});

}  // namespace sgmres