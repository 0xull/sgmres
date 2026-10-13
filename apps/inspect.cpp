// Builds a test matrix, reports its shape and symmetry, and times y = A * x. The matvec time is
// the unit every solver cost gets compared against later.

#include <cstdio>
#include <cstdlib>
#include <string_view>

#include "sgmres/problems.hpp"
#include "sgmres/stopwatch.hpp"

namespace {

void usage(const char* argv0) {
  std::fprintf(stderr,
               "usage: %s <problem> [size] [reps]\n"
               "  problem  tridiag | laplacian | convdiff\n"
               "  size     n for tridiag, grid side m otherwise (default 300)\n"
               "  reps     number of matvecs to time (default 100)\n",
               argv0);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    usage(argv[0]);
    return 2;
  }

  const std::string_view problem = argv[1];
  const sgmres::Index size = argc > 2 ? std::atol(argv[2]) : 300;
  const int reps = argc > 3 ? std::atoi(argv[3]) : 100;
  if (size <= 0 || reps <= 0) {
    usage(argv[0]);
    return 2;
  }

  sgmres::SparseMatrix a;
  if (problem == "tridiag") {
    a = sgmres::tridiagonal(size, -1.0, 3.0, 1.0);
  } else if (problem == "laplacian") {
    a = sgmres::laplacian_2d(size);
  } else if (problem == "convdiff") {
    a = sgmres::convection_diffusion_2d(size, 0.5);
  } else {
    usage(argv[0]);
    return 2;
  }

  const sgmres::SparseMatrix at = a.transpose();
  const bool symmetric = (a - at).norm() == 0.0;

  const sgmres::Vector x = sgmres::random_normal(a.cols(), 1);
  sgmres::Vector y(a.rows());

  y.noalias() = a * x;  // first touch, keep it out of the timing
  sgmres::Stopwatch clock;
  for (int r = 0; r < reps; ++r) y.noalias() = a * x;
  const double matvec_ms = 1e3 * clock.seconds() / reps;

  std::printf("problem    %.*s\n", static_cast<int>(problem.size()), problem.data());
  std::printf("n          %ld\n", static_cast<long>(a.rows()));
  std::printf("nnz        %ld (%.2f per row)\n", static_cast<long>(a.nonZeros()),
              static_cast<double>(a.nonZeros()) / static_cast<double>(a.rows()));
  std::printf("symmetric  %s\n", symmetric ? "yes" : "no");
  std::printf("matvec     %.3f ms (mean of %d)\n", matvec_ms, reps);
  std::printf("|A x|      %.6e\n", y.norm());
  return 0;
}
