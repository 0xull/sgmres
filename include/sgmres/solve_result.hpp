#pragma once

#include <vector>

#include "sgmres/types.hpp"

namespace sgmres {

struct SolveResult {
  Vector x;
  double relative_residual = 0.0;  // ||f - Ax|| / ||f||, computed from x rather than estimated
  Index iterations = 0;            // basis vectors built, one product with A each
  Index cycles = 0;
  bool converged = false;
  std::vector<double> history;  // relative residual estimate after each iteration
};

}  // namespace sgmres