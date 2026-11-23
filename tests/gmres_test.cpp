#include "sgmres/gmres.hpp"

#include <cmath>

#include "check.hpp"
#include "sgmres/problems.hpp"

namespace {

// One cycle of two steps on the 8 x 8 example from the talk. The exact answer there is
// x = 12/37 + (-0.1, 0, ..., 0, 0.1) with residual sqrt(22.2) / 37 = 0.12734.
void hand_example() {
  const sgmres::SparseMatrix a = sgmres::tridiagonal(8, -1.0, 3.0, 1.0);
  const sgmres::Vector f = sgmres::Vector::Ones(8);
  const sgmres::Vector x0 = sgmres::Vector::Zero(8);

  const auto result = sgmres::gmres(a, f, x0, {.restart = 2, .max_cycles = 1, .tol = 1e-14});
  CHECK(result.iterations == 2);
  CHECK(result.cycles == 1);
  CHECK(!result.converged);

  CHECK_NEAR(result.x[0], 12.0 / 37.0 - 0.1, 1e-14);
  for (int i = 1; i < 7; ++i) CHECK_NEAR(result.x[i], 12.0 / 37.0, 1e-14);
  CHECK_NEAR(result.x[7], 12.0 / 37.0 + 0.1, 1e-14);

  const double residual = std::sqrt(22.2) / 37.0;
  CHECK_NEAR(result.relative_residual * f.norm(), residual, 1e-14);

  // Each rotation scales the residual by its sine, so the history is a running product.
  const double s0 = 0.5 / std::sqrt(9.25);
  const double s1 = (std::sqrt(3.0) / 2.0) / std::sqrt(10.0);
  CHECK(result.history.size() == 2);
  CHECK_NEAR(result.history[0], s0, 1e-14);
  CHECK_NEAR(result.history[1], s0 * s1, 1e-14);
}

// With as many steps as unknowns the Krylov space is the whole space, so GMRES is exact.
void full_space() {
  const sgmres::SparseMatrix a = sgmres::tridiagonal(8, -1.0, 3.0, 1.0);
  const sgmres::Vector f = sgmres::Vector::Ones(8);

  const auto result =
      sgmres::gmres(a, f, sgmres::Vector::Zero(8), {.restart = 8, .max_cycles = 1, .tol = 1e-12});
  CHECK(result.converged);
  CHECK(result.relative_residual < 1e-12);
}

void restarted() {
  const sgmres::SparseMatrix a = sgmres::convection_diffusion_2d(16, 0.5, 1.0);
  const sgmres::Vector f = sgmres::random_normal(a.rows(), 11);

  const auto result = sgmres::gmres(a, f, sgmres::Vector::Zero(a.rows()),
                                    {.restart = 10, .max_cycles = 100, .tol = 1e-10});
  CHECK(result.converged);
  CHECK(result.cycles > 1);
  CHECK(result.relative_residual <= 1e-10);
  CHECK((f - a * result.x).norm() <= 1e-10 * f.norm());

  // GMRES never lets the residual grow, within a cycle or across a restart.
  for (std::size_t i = 1; i < result.history.size(); ++i) {
    CHECK(result.history[i] <= result.history[i - 1] * (1 + 1e-8));
  }
}

void starts_converged() {
  const sgmres::SparseMatrix a = sgmres::convection_diffusion_2d(8, 0.5, 1.0);
  const sgmres::Vector x = sgmres::random_normal(a.rows(), 2);
  const sgmres::Vector f = a * x;

  const auto result = sgmres::gmres(a, f, x);
  CHECK(result.converged);
  CHECK(result.iterations == 0);
}

}  // namespace

int main() {
  hand_example();
  full_space();
  restarted();
  starts_converged();
  return sgmres::test::finish();
}
