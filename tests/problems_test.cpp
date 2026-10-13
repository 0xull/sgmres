#include "sgmres/problems.hpp"

#include <cmath>
#include <cstdio>

namespace {

int failures = 0;

#define CHECK(cond)                                                                 \
  do {                                                                              \
    if (!(cond)) {                                                                  \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond); \
      ++failures;                                                                   \
    }                                                                               \
  } while (0)

#define CHECK_NEAR(a, b, tol) CHECK(std::abs((a) - (b)) <= (tol))

bool is_symmetric(const sgmres::SparseMatrix& a) {
  const sgmres::SparseMatrix at = a.transpose();
  return (a - at).norm() == 0.0;
}

// Numbers from the 8 x 8 example worked by hand in the talk.
void hand_example() {
  const sgmres::SparseMatrix a = sgmres::tridiagonal(8, -1.0, 3.0, 1.0);
  const sgmres::Vector b1 = sgmres::Vector::Ones(8) / std::sqrt(8.0);
  const sgmres::Vector m1 = a * b1;

  CHECK_NEAR(m1[0], 4.0 / std::sqrt(8.0), 1e-14);
  for (int i = 1; i < 7; ++i) CHECK_NEAR(m1[i], 3.0 / std::sqrt(8.0), 1e-14);
  CHECK_NEAR(m1[7], 2.0 / std::sqrt(8.0), 1e-14);
  CHECK_NEAR(b1.dot(m1), 3.0, 1e-14);
}

void laplacian() {
  const sgmres::SparseMatrix a = sgmres::laplacian_2d(4);
  CHECK(a.rows() == 16);
  CHECK(a.nonZeros() == 5 * 16 - 4 * 4);
  CHECK(is_symmetric(a));
  CHECK(a.coeff(5, 5) == 4.0);
  CHECK(a.coeff(5, 6) == -1.0);
  CHECK(a.coeff(3, 4) == 0.0);  // end of one grid row does not touch the start of the next

  CHECK(sgmres::laplacian_2d(100).nonZeros() == 5 * 10000 - 4 * 100);
}

void convection_diffusion() {
  const sgmres::SparseMatrix a = sgmres::convection_diffusion_2d(4, 0.5);
  CHECK(!is_symmetric(a));
  CHECK(a.nonZeros() == sgmres::laplacian_2d(4).nonZeros());
  CHECK(a.coeff(5, 6) == -1.5);
  CHECK(a.coeff(6, 5) == -1.0);

  const sgmres::SparseMatrix shifted = sgmres::convection_diffusion_2d(4, 0.5, 1.0);
  CHECK(shifted.coeff(5, 5) == 5.0);
}

void random_vectors() {
  const sgmres::Vector a = sgmres::random_normal(100000, 7);
  const sgmres::Vector b = sgmres::random_normal(100000, 7);
  const sgmres::Vector c = sgmres::random_normal(100000, 8);
  CHECK(a == b);
  CHECK(a != c);

  const double mean = a.mean();
  const double var = (a.array() - mean).square().mean();
  CHECK(std::abs(mean) < 0.02);
  CHECK(std::abs(var - 1.0) < 0.02);
}

}  // namespace

int main() {
  hand_example();
  laplacian();
  convection_diffusion();
  random_vectors();

  if (failures != 0) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::puts("all checks passed");
  return 0;
}
