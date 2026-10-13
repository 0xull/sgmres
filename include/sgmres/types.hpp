#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sgmres {

using Index = Eigen::Index;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

// Row-major so y = A * x streams each row once. Every matvec in the solver goes through this.
using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

}  // namespace sgmres
