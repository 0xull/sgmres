#pragma once

#include <cmath>
#include <cstdio>

namespace sgmres::test {
inline int failure = 0;

inline int finish() {
  if (failure != 0) {
    std::fprintf(stderr, "%d check(s) failed\n", failure);
    return 1;
  }
  std::puts("all checks passed");
  return 0;
}

}  // namespace sgmres::test

#define CHECK(cond)                                                                 \
  do {                                                                              \
    if (!(cond)) {                                                                  \
      std::fprintf(stderr, "%s:%d: CHECK(%s) failed\n", __FILE__, __LINE__, #cond); \
      ++::sgmres::test::failure;                                                    \
    }                                                                               \
  } while (0)

#define CHECK_NEAR(a, b, tol) CHECK(std::abs((a) - (b)) <= (tol))
