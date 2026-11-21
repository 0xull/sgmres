#include "sgmres/arnoldi.hpp"

#include <cassert>
#include <limits>

namespace sgmres {
bool arnoldi_step(const SparseMatrix& a, Matrix& q, Matrix& h, Index j) {
  assert(q.cols() >= j + 2 && h.rows() >= j + 2 && h.cols() >= j + 1);

  Vector w = a * q.col(j);
  const double scale = w.norm();
  const auto basis = q.leftCols(j + 1);

  // Classical Gram-Schmidt, run twice. One pass loses orthogonality once w is close to the
  // span of the basis, a second pass restores it for the price of one more sweep.
  Vector coeffs = basis.transpose() * w;
  w.noalias() -= basis * coeffs;
  const Vector correction = basis.transpose() * w;
  w.noalias() -= basis * correction;
  coeffs += correction;

  h.col(j).head(j + 1) = coeffs;

  // Anything below a few ulps of the original vector is rounding noise, not a new direction.
  const double norm = w.norm();
  if (norm <= 10 * std::numeric_limits<double>::epsilon() * scale) {
    h(j + 1, j) = 0.0;
    return false;
  }

  h(j + 1, j) = norm;
  q.col(j + 1) = w / norm;
  return true;
}
}  // namespace sgmres
