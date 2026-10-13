#include "sgmres/problems.hpp"

#include <random>
#include <vector>

namespace sgmres {

namespace {

using Triplet = Eigen::Triplet<double>;

SparseMatrix from_triplets(Index n, const std::vector<Triplet>& entries) {
  SparseMatrix a(n, n);
  a.setFromTriplets(entries.begin(), entries.end());
  a.makeCompressed();
  return a;
}

// Shared by both grid operators. east is the coupling to (i, j + 1), every other neighbour
// gets -1.
SparseMatrix five_point(Index m, double diag, double east) {
  const Index n = m * m;
  std::vector<Triplet> entries;
  entries.reserve(static_cast<std::size_t>(5 * n));

  for (Index i = 0; i < m; ++i) {
    for (Index j = 0; j < m; ++j) {
      const Index row = i * m + j;
      entries.emplace_back(row, row, diag);
      if (i > 0) entries.emplace_back(row, row - m, -1.0);
      if (i + 1 < m) entries.emplace_back(row, row + m, -1.0);
      if (j > 0) entries.emplace_back(row, row - 1, -1.0);
      if (j + 1 < m) entries.emplace_back(row, row + 1, east);
    }
  }
  return from_triplets(n, entries);
}

}  // namespace

SparseMatrix tridiagonal(Index n, double lower, double diag, double upper) {
  std::vector<Triplet> entries;
  entries.reserve(static_cast<std::size_t>(3 * n));

  for (Index i = 0; i < n; ++i) {
    if (i > 0) entries.emplace_back(i, i - 1, lower);
    entries.emplace_back(i, i, diag);
    if (i + 1 < n) entries.emplace_back(i, i + 1, upper);
  }
  return from_triplets(n, entries);
}

SparseMatrix laplacian_2d(Index m) { return five_point(m, 4.0, -1.0); }

SparseMatrix convection_diffusion_2d(Index m, double convection, double shift) {
  return five_point(m, 4.0 + shift, -1.0 - convection);
}

Vector random_normal(Index n, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal(0.0, 1.0);

  Vector v(n);
  for (Index i = 0; i < n; ++i) v[i] = normal(rng);
  return v;
}

}  // namespace sgmres
