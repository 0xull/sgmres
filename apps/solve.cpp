// Solves Ax = f for one of the test problems with a random f and reports how it went.

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

#include "sgmres/gmres.hpp"
#include "sgmres/problems.hpp"
#include "sgmres/stopwatch.hpp"

namespace {

void usage(const char* argv0) {
  std::fprintf(stderr,
               "usage: %s <problem> <size> [options]\n"
               "  --restart m      basis vectors per cycle (default 50)\n"
               "  --tol t          relative residual to stop at (default 1e-8)\n"
               "  --max-cycles c    give up after c cycles (default 200)\n"
               "  --shift s        solve with A + s * I\n",
               argv0);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    usage(argv[0]);
    return 2;
  }

  const std::string_view problem = argv[1];
  const sgmres::Index size = std::atol(argv[2]);
  sgmres::GmresOptions options;
  double shift = 0.0;

  for (int i = 3; i + 1 < argc; i += 2) {
    const std::string_view flag = argv[i];
    const char* value = argv[i + 1];
    if (flag == "--restart")
      options.restart = std::atol(value);
    else if (flag == "--tol")
      options.tol = std::atof(value);
    else if (flag == "--max-cycles")
      options.max_cycles = std::atol(value);
    else if (flag == "--shift")
      shift = std::atof(value);
    else {
      usage(argv[0]);
      return 2;
    }
  }

  if (argc % 2 == 0 || size <= 0 || options.restart <= 0 || options.max_cycles <= 0) {
    usage(argv[0]);
    return 2;
  }

  sgmres::SparseMatrix a;
  try {
    a = sgmres::make_problem(problem, size, shift);
  } catch (const std::invalid_argument& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 2;
  }

  const sgmres::Vector f = sgmres::random_normal(a.rows(), 1);
  const sgmres::Vector x0 = sgmres::Vector::Zero(a.cols());

  sgmres::Stopwatch clock;
  const sgmres::SolveResult result = sgmres::gmres(a, f, x0, options);
  const double seconds = clock.seconds();

  std::printf("problem     %.*s, shift %g\n", static_cast<int>(problem.size()), problem.data(),
              shift);
  std::printf("n           %ld\n", static_cast<long>(a.rows()));
  std::printf("method      gmres(%ld)\n", static_cast<long>(options.restart));
  std::printf("iterations  %ld in %ld cycles\n", static_cast<long>(result.iterations),
              static_cast<long>(result.cycles));
  std::printf("residual    %.3e (relative)\n", result.relative_residual);
  std::printf("converged   %s\n", result.converged ? "yes" : "no");
  std::printf("time        %.3f s\n", seconds);

  return result.converged ? 0 : 1;
}
