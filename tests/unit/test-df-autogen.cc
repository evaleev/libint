#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/dfbs_generator.h>
#include <libint2/engine.h>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "catch.hpp"
#include "fixture.h"

Eigen::Tensor<double, 3> compute_eri3(const BasisSet &obs,
                                      const BasisSet &dfbs) {
  const auto nshells = obs.size();
  const auto nshells_df = dfbs.size();
  const auto &unitshell = libint2::Shell::unit();

  // construct the 2-electron 3-center repulsion integrals engine
  // since the code assumes (xx|xs) braket, and Engine/libint only produces
  // (xs|xx), use 4-center engine
  libint2::Engine engine(libint2::Operator::coulomb,
                         std::max(obs.max_nprim(), dfbs.max_nprim()),
                         std::max(obs.max_l(), dfbs.max_l()), 0);
  engine.set(BraKet::xs_xx);
  const auto shell2bf = obs.shell2bf();
  const auto shell2bf_df = dfbs.shell2bf();
  const auto nbf = obs.nbf();
  const auto ndf = dfbs.nbf();

  Eigen::Tensor<double, 3> result(ndf, nbf, nbf);
  const auto &buf = engine.results();
  for (auto s1 = 0; s1 != nshells_df; ++s1) {
    const auto bf1 = shell2bf_df[s1];
    const auto n1 = dfbs[s1].size();
    for (auto s2 = 0; s2 != nshells; ++s2) {
      const auto bf2 = shell2bf[s2];
      const auto n2 = obs[s2].size();
      for (auto s3 = 0; s3 != nshells;
           ++s3) {  // only loop over unique ao shells
        const auto bf3 = shell2bf[s3];
        const auto n3 = obs[s3].size();
        engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(
            dfbs[s1], unitshell, obs[s2], obs[s3]);
        size_t fij = 0;
        for (auto f = 0; f != n1; ++f) {
          for (auto i = 0; i != n2; ++i) {
            for (auto j = 0; j != n3; ++j) {
              result(bf1 + f, bf2 + i, bf3 + j) = buf[0][fij];
              ++fij;
            }
          }
        }
      }
    }
  }

  return result;
}

Eigen::MatrixXd compute_eri2(const BasisSet &dfbs) {
  const auto ndf = dfbs.nbf();
  Eigen::MatrixXd result(ndf, ndf);
  result.fill(0.0);
  libint2::Engine engine(libint2::Operator::coulomb, dfbs.max_nprim(),
                         dfbs.max_l(), 0);
  engine.set(BraKet::xs_xs);
  engine.set(libint2::ScreeningMethod::Conservative);
  const auto shell2bf_df = dfbs.shell2bf();
  const auto &buf = engine.results();
  for (size_t s1 = 0; s1 != dfbs.size(); ++s1) {
    const size_t bf1 = shell2bf_df[s1];
    const size_t n1 = dfbs[s1].size();
    for (size_t s2 = 0; s2 != dfbs.size(); ++s2) {
      const size_t bf2 = shell2bf_df[s2];
      const size_t n2 = dfbs[s2].size();
      engine.compute(dfbs[s1], dfbs[s2]);
      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
    }
  }
  return result;
}

Eigen::Tensor<double, 4> compute_eri4(const BasisSet &obs) {
  const auto nshells = obs.size();
  libint2::Engine engine(libint2::Operator::coulomb,
                         std::max(obs.max_nprim(), obs.max_nprim()),
                         std::max(obs.max_l(), obs.max_l()), 0);
  engine.set(BraKet::xx_xx);
  const auto shell2bf = obs.shell2bf();
  const auto nbf = obs.nbf();

  Eigen::Tensor<double, 4> result(nbf, nbf, nbf, nbf);
  const auto &buf = engine.results();
  for (auto s1 = 0; s1 != nshells; ++s1) {
    const auto bf1 = shell2bf[s1];
    const auto n1 = obs[s1].size();
    // TODO: loop over all unique quartets: s1 >= s2 >= s3 >= s4
    for (auto s2 = 0; s2 != nshells; ++s2) {
      const auto bf2 = shell2bf[s2];
      const auto n2 = obs[s2].size();
      for (auto s3 = 0; s3 != nshells; ++s3) {
        const auto bf3 = shell2bf[s3];
        const auto n3 = obs[s3].size();
        for (auto s4 = 0; s4 != nshells; ++s4) {
          const auto bf4 = shell2bf[s4];
          const auto n4 = obs[s4].size();
          engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
              obs[s1], obs[s2], obs[s3], obs[s4]);
          size_t ijkl = 0;
          for (auto i = 0; i != n1; ++i) {
            for (auto j = 0; j != n2; ++j) {
              for (auto k = 0; k != n3; ++k) {
                for (auto l = 0; l != n4; ++l) {
                  result(bf1 + i, bf2 + j, bf3 + k, bf4 + l) = buf[0][ijkl];
                  ++ijkl;
                }
              }
            }
          }
        }
      }
    }
  }

  return result;
}

TEST_CASE("DFBS-Generator", "[dfbs-generator]") {
  // Will use Neon as a test case
  libint2::Atom atom;
  atom.atomic_number = 2;
  atom.x = 0.0;
  atom.y = 0.0;
  atom.z = 0.0;

  // Using STO-3G basis set for Neon
  libint2::BasisSet bs("STO-3G", {atom});
  const auto nao = bs.nbf();

  // Use automatic DF basis set generator with tight cholesky threshold
  libint2::DFBasisSetGenerator dfbs_generator(bs.shells(), 1e-4);
  const auto dfbs = dfbs_generator.reduced_basis();
  const auto ndf = dfbs.nbf();

  // Compute 2-center integrals
  const auto eri2 = compute_eri2(dfbs);

  // Compute 3-center integrals
  const auto eri3 = compute_eri3(bs, dfbs);

  // Compute 4-center integrals
  const auto eri4 = compute_eri4(bs);

  // Compute L_inv
  Eigen::LLT<Eigen::MatrixXd> V_LLt(eri2);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ndf, ndf);
  auto L = V_LLt.matrixL();
  Eigen::MatrixXd Linv_t = L.solve(I).transpose();
  auto eri2_inv = Linv_t * Linv_t.transpose();

  // reconstruct 4-center integrals and compare with original
  double norm = 0.0;
  Eigen::Tensor<double, 4> eri4_df(nao, nao, nao, nao);
  for (size_t i = 0; i < nao; i++) {
    for (size_t j = 0; j < nao; j++) {
      for (size_t k = 0; k < nao; k++) {
        for (size_t l = 0; l < nao; l++) {
          for (size_t X = 0; X < ndf; X++) {
            for (size_t Y = 0; Y < ndf; Y++) {
              eri4_df(i, j, k, l) +=
                  eri3(X, i, j) * eri2_inv(X, Y) * eri3(Y, k, l);
            }
          }
          norm += std::pow(eri4_df(i, j, k, l) - eri4(i, j, k, l), 2);
        }
      }
    }
  }
  REQUIRE_NOTHROW(norm < 1e-10);
}
