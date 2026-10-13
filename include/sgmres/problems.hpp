#pragma once

#include <cstdint>

#include "sgmres/types.hpp"

namespace sgmres {

// n x n with constant diagonals. tridiagonal(8, -1, 3, 1) is the worked example from the talk.
SparseMatrix tridiagonal(Index n, double lower, double diag, double upper);

// 5-point Laplacian on an m x m grid with zero Dirichlet boundary, n = m * m. Symmetric
// positive definite.
SparseMatrix laplacian_2d(Index m);

// laplacian_2d plus an extra -convection on each east neighbour, which breaks symmetry, plus
// shift * I. A shift of 1 moves the spectrum away from zero and makes GMRES converge quickly.
SparseMatrix convection_diffusion_2d(Index m, double convection, double shift = 0.0);

// Entries drawn from N(0, 1). Fixed seeds keep runs comparable across solvers.
Vector random_normal(Index n, std::uint64_t seed);

}  // namespace sgmres
