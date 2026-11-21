#include "sgmres/arnoldi.hpp"

#include <cmath>

#include "check.hpp"
#include "sgmres/problems.hpp"
#include "sgmres/types.hpp"

namespace {
    // The first two Arnoldi steps of the 8 x 8 example from the talk.
    void hand_example() {
        const sgmres::SparseMatrix a = sgmres::tridiagonal(8, -1.0, 3.0, 1.0);
        sgmres::Matrix q = sgmres::Matrix::Zero(8, 3);
        sgmres::Matrix h = sgmres::Matrix::Zero(3, 2);
        q.col(0).setConstant(1.0 / std::sqrt(8.0));

        CHECK(sgmres::arnoldi_step(a, q, h, 0));
        CHECK_NEAR(h(0, 0), 3.0, 1e-14);
        CHECK_NEAR(h(1, 0), 0.5, 1e-14);

        sgmres::Vector b2 = sgmres::Vector::Zero(8);
        b2[0] = 1.0 / std::sqrt(2.0);
        b2[7] = -1.0 / std::sqrt(2.0);
        CHECK((q.col(1) - b2).norm() < 1e-14);

        CHECK(sgmres::arnoldi_step(a, q, h, 1));
        CHECK_NEAR(h(0, 1), -0.5, 1e-14);
        CHECK_NEAR(h(1, 1), 3.0, 1e-14);
        CHECK_NEAR(h(2, 1), std::sqrt(3.0) / 2.0, 1e-14);
    }

    // A * Q_m = Q_{m+1} * H with orthonormal Q, on a nonsymmetric matrix.
    void arnoldi_relation() {
        const sgmres::SparseMatrix a = sgmres::convection_diffusion_2d(10, 0.5);
        const sgmres::Index n = a.rows();
        const sgmres::Index m = 30;

        sgmres::Matrix q(n, m + 1);
        sgmres::Matrix h = sgmres::Matrix::Zero(m + 1, m);
        const sgmres::Vector r = sgmres::random_normal(n, 3);
        q.col(0) = r / r.norm();
        for (sgmres::Index j = 0; j < m; ++j) CHECK(sgmres::arnoldi_step(a, q, h, j));

        const sgmres::Matrix gram = q.transpose() * q;
        CHECK((gram - sgmres::Matrix::Identity(m + 1, m + 1)).norm() < 1e-12);

        const sgmres::Matrix aq = a * q.leftCols(m);
        CHECK((aq - q * h).norm() < 1e-12 * aq.norm());
    }

    // A = 2I maps every vector onto itself, so the Krylov space is one-dimensional and the first
    // step has nothing new to add.
    void invariant_space() {
        const sgmres::SparseMatrix a = sgmres::tridiagonal(6, 0.0, 2.0, 0.0);
        sgmres::Matrix q = sgmres::Matrix::Zero(6, 2);
        sgmres::Matrix h = sgmres::Matrix::Zero(2, 1);
        const sgmres::Vector r = sgmres::random_normal(6, 5);
        q.col(0) = r / r.norm();

        CHECK(!sgmres::arnoldi_step(a, q, h, 0));
        CHECK_NEAR(h(0, 0), 2.0, 1e-14);
        CHECK(h(1, 0) == 0.0);
    }
} // namespace

int main() {
    hand_example();
    arnoldi_relation();
    invariant_space();
    return sgmres::test::finish();
}
