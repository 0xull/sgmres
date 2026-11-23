#include "sgmres/gmres.hpp"

#include <cassert>
#include <cmath>
#include <vector>

#include "sgmres/arnoldi.hpp"

namespace sgmres {

namespace {

// Rotation in the (i, i + 1) plane that zeroes the lower of two entries.
struct Givens {
  double c = 1.0;
  double s = 0.0;

  static Givens zeroing(double a, double b) {
    const double r = std::hypot(a, b);
    if (r == 0.0) return {};
    return {a / r, b / r};
  }

  void apply(double& x, double& y) const {
    const double top = c * x + s * y;
    y = -s * x + c * y;
    x = top;
  }
};

}  // namespace

SolveResult gmres(const SparseMatrix& a, const Vector& f, const Vector& x0,
                  const GmresOptions& options) {
  assert(options.restart > 0);
  const Index n = a.rows();
  const Index m = options.restart;
  const double f_norm = f.norm();

  SolveResult result;
  result.x = x0;
  if (f_norm == 0.0) {
    result.x.setZero();
    result.converged = true;
    return result;
  }
  const double target = options.tol * f_norm;

  Matrix q(n, m + 1);
  Matrix h(m + 1, m);
  Vector g(m + 1);
  std::vector<Givens> rotations(static_cast<std::size_t>(m));

  Vector r = f - a * result.x;
  double r_norm = r.norm();

  while (r_norm > target && result.cycles < options.max_cycles) {
    q.col(0) = r / r_norm;
    g.setZero();
    g(0) = r_norm;

    // Each new column of H is rotated into upper triangular form as it arrives. After step j
    // the least-squares residual is |g(j + 1)|, so convergence is known without forming x.
    Index steps = 0;
    for (Index j = 0; j < m; ++j) {
      const bool extended = arnoldi_step(a, q, h, j);
      for (Index i = 0; i < j; ++i) rotations[i].apply(h(i, j), h(i + 1, j));
      rotations[j] = Givens::zeroing(h(j, j), h(j + 1, j));
      rotations[j].apply(h(j, j), h(j + 1, j));
      rotations[j].apply(g(j), g(j + 1));

      steps = j + 1;
      ++result.iterations;
      result.history.push_back(std::abs(g(j + 1)) / f_norm);
      if (!extended || std::abs(g(j + 1)) <= target) break;
    }

    const Vector y =
        h.topLeftCorner(steps, steps).triangularView<Eigen::Upper>().solve(g.head(steps));
    result.x.noalias() += q.leftCols(steps) * y;

    // The estimate drifts from the true residual in floating point, so restart from the truth.
    r = f - a * result.x;
    r_norm = r.norm();
    ++result.cycles;
  }

  result.relative_residual = r_norm / f_norm;
  result.converged = r_norm <= target;
  return result;
}

}  // namespace sgmres