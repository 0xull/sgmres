# sgmres

A C++ implementation of sketched GMRES from Nakatsukasa and Tropp, *Fast and Accurate Randomized
Algorithms for Linear Systems and Eigenvalue Problems*, SIAM J. Matrix Anal. Appl. 45(2), 2024
([arXiv:2111.00113](https://arxiv.org/abs/2111.00113)).

I presented the linear-systems half of the paper to friends. The [slides](https://docs.google.com/presentation/d/1kwKwuFp0f44CgoTUpiEOIHIjKhg2pvIi2IuPbbTLwPc/edit?usp=sharing) and the [recording](https://drive.google.com/file/d/1xj4FvTcFOlwTL1ZIN6z47QcOCVCvMbPj/view?usp=sharing) are both public.

So far this has the build, the test matrices and a matvec timer. The solvers come next.

## Build

Needs CMake 3.20+, a C++20 compiler and Eigen 3.4.

```sh
# macOS
brew install cmake eigen

# Debian / Ubuntu
sudo apt install cmake g++ libeigen3-dev
```

```sh
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Run

```sh
./build/inspect convdiff 300
```

`inspect <problem> [size] [reps]` builds a matrix and times `y = A * x`. Problems are `tridiag`
(the 8 x 8 example from the talk is `tridiag 8`), `laplacian` and `convdiff`. For the grid
problems `size` is the grid side, so `300` gives n = 90,000.

## Layout

```
include/sgmres/   public headers
src/              library sources
apps/             command line tools
tests/            ctest executables
```
