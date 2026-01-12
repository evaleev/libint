/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_soad_fock_h_
#define _libint2_src_lib_libint_soad_fock_h_

#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/chemistry/sto3g_atomic_density.h>

#include <Eigen/Dense>

namespace libint2 {

/// @brief computes Superposition-Of-Atomic-Densities guess for the molecular
/// 1-RDM

/// The atomic densities are represented by the STO-3G AOs;
/// occupies subshells by smearing electrons (of neutral ground-state atoms)
/// evenly over the subshell components
/// @param[in] atoms list of atoms
/// @return the 1-RDM in the STO-3G AO basis
/// @return the 1-RDM, normalized to the # of electrons/2
Eigen::MatrixXd compute_soad(const std::vector<Atom> &atoms) {
  // compute number of atomic orbitals
  size_t nao = 0;
  for (const auto &atom : atoms) {
    const auto Z = atom.atomic_number;
    nao += libint2::sto3g_num_ao(Z);
  }

  // compute the minimal basis density
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nao, nao);
  size_t ao_offset = 0;  // first AO of this atom
  for (const auto &atom : atoms) {
    const auto Z = atom.atomic_number;
    const auto &occvec = libint2::sto3g_ao_occupation_vector(Z);
    for (const auto &occ : occvec) {
      D(ao_offset, ao_offset) = occ;
      ++ao_offset;
    }
  }

  return D * 0.5;  // we use densities normalized to # of electrons/2
}

Eigen::MatrixXd compute_2body_fock_general(const BasisSet &obs,
                                           const Eigen::MatrixXd &D,
                                           const BasisSet &D_bs,
                                           bool D_is_shelldiagonal,
                                           double precision) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto n_D = D_bs.nbf();
  assert(D.cols() == D.rows() && D.cols() == n_D);

  Eigen::MatrixXd G = Eigen::MatrixXd ::Zero(n, n);

  // construct the 2-electron repulsion integrals engine
  libint2::Engine engine(libint2::Operator::coulomb,
                         std::max(obs.max_nprim(), D_bs.max_nprim()),
                         std::max(obs.max_l(), D_bs.max_l()), 0);
  engine.set_precision(precision);
  auto shell2bf = obs.shell2bf();
  auto shell2bf_D = D_bs.shell2bf();

  const auto &buf = engine.results();

  // loop over permutationally-unique set of shells
  for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
    auto bf1_first = shell2bf[s1];  // first basis function in this shell
    auto n1 = obs[s1].size();       // number of basis functions in this shell

    for (auto s2 = 0; s2 <= s1; ++s2) {
      auto bf2_first = shell2bf[s2];
      auto n2 = obs[s2].size();

      for (auto s3 = 0; s3 < D_bs.size(); ++s3) {
        auto bf3_first = shell2bf_D[s3];
        auto n3 = D_bs[s3].size();

        auto s4_begin = D_is_shelldiagonal ? s3 : 0;
        auto s4_fence = D_is_shelldiagonal ? s3 + 1 : D_bs.size();

        for (auto s4 = s4_begin; s4 != s4_fence; ++s4, ++s1234) {
          auto bf4_first = shell2bf_D[s4];
          auto n4 = D_bs[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of
          // the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

          if (s3 >= s4) {
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s1234_deg = s12_deg * s34_deg;
            // auto s1234_deg = s12_deg;
            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], obs[s2], D_bs[s3], D_bs[s4]);
            const auto *buf_1234 = buf[0];
            if (buf_1234 != nullptr) {
              for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = buf_1234[f1234];
                      const auto value_scal_by_deg = value * s1234_deg;
                      G(bf1, bf2) += 2.0 * D(bf3, bf4) * value_scal_by_deg;
                    }
                  }
                }
              }
            }
          }

          engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
              obs[s1], D_bs[s3], obs[s2], D_bs[s4]);
          const auto *buf_1324 = buf[0];
          if (buf_1324 == nullptr)
            continue;  // if all integrals screened out, skip to next quartet

          for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for (auto f3 = 0; f3 != n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_1324[f1324];
                  const auto value_scal_by_deg = value * s12_deg;
                  G(bf1, bf2) -= D(bf3, bf4) * value_scal_by_deg;
                }
              }
            }
          }
        }
      }
    }
  }

  // symmetrize the result and return
  return 0.5 * (G + G.transpose());
}

Eigen::MatrixXd compute_soad_fock(const BasisSet &obs, const BasisSet &minbs,
                                  const std::vector<Atom> &atoms) {
  // use SOAD as guess for the density matrix
  auto D = compute_soad(atoms);
  // compute fock with guess density
  auto F = compute_2body_fock_general(
      obs, D, minbs, true /* SOAD_D_is_shelldiagonal */,
      std::numeric_limits<double>::epsilon()  // this is cheap, no reason
                                              // to be cheaper
  );
  return F;
}

}  // namespace libint2

#endif  //_libint2_src_lib_libint_soad_fock_h_
