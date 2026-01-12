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

#ifndef _libint2_src_lib_libint_dfbs_generator_h_
#define _libint2_src_lib_libint_dfbs_generator_h_

#include <libint2.h>
#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/boys.h>
#include <libint2/pivoted_cholesky.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

namespace libint2 {

namespace detail {

/// @brief Uncontracts a set of shells
/// @param[in] shells a vector of shells to be uncontracted
/// @return a vector of uncontracted shells
inline std::vector<Shell> uncontract(const std::vector<Shell> &shells) {
  std::vector<Shell> primitive_shells;
  for (const auto &contracted_shell : shells) {
    for (size_t p = 0; p < contracted_shell.nprim(); p++) {
      const auto prim_shell = contracted_shell.extract_primitive(p, true);
      // if dealing with generally contracted basis (e.g., cc-pvxz) represented
      // as a segmented basis need to remove duplicates
      if (std::find(primitive_shells.begin(), primitive_shells.end(),
                    prim_shell) == primitive_shells.end())
        primitive_shells.emplace_back(std::move(prim_shell));
    }
  }
  return primitive_shells;
}

/// @brief returns \f$ \Gamma(x) \f$
inline double gamma_function(const double x) { return std::tgamma(x); }

/// @brief Computes an effective exponent for a product of two primitive
/// gaussians as as described in J. Chem. Theory Comput. 2021, 17, 6886−6900.
/// \f$ \alpha_{eff} = \left[ \frac{\Gamma(L+2)\Gamma(l_\mu +l_\nu +
/// \frac{3}{2})}{\Gamma(L+ \frac{3}{2})\Gamma(l_\mu +l_\nu + 2)} \right]^2
/// (\alpha_\mu + \alpha_\nu) \f$
/// @param shell1 first primitive shell
/// @param shell2 second primitive shell
/// @param L total angular momentum of product function
/// @return effective exponent of product function
inline double alpha_eff(const Shell &shell1, const Shell &shell2, const int L) {
  const auto alpha1 = shell1.alpha[0];
  const auto alpha2 = shell2.alpha[0];
  const auto l1 = shell1.contr[0].l;
  const auto l2 = shell2.contr[0].l;
  const auto prefactor =
      std::pow((gamma_function(L + 2.) * gamma_function(l1 + l2 + 1.5)) /
                   (gamma_function(l1 + l2 + 2.) * gamma_function(L + 1.5)),
               2.);
  return prefactor * (alpha1 + alpha2);
}

/// @brief Creates a set of product functions from a set of primitive shells.
/// Each pair of primitive shells produces a set of product functions with an
/// effective exponent (\f$ \alpha_{eff} \f$ as described above) and angular
/// momentum ranging from \f$ |l1-l2| \leq L \leq l1+l2 \f$ as described in J.
/// Chem. Theory Comput. 2021, 17, 6886−6900.
/// @param primitive_shells set of primitive shells
/// @return set of product functions
inline std::vector<Shell> product_functions(
    const std::vector<Shell> &primitive_shells) {
  std::vector<Shell> product_functions;
  for (size_t i = 0; i < primitive_shells.size(); ++i) {
    for (size_t j = 0; j <= i; ++j) {
      const auto li = primitive_shells[i].contr[0].l;
      const auto lj = primitive_shells[j].contr[0].l;
      for (auto L = std::abs(li - lj); L <= li + lj; L++) {
        const auto alpha = libint2::svector<double>(
            {alpha_eff(primitive_shells[i], primitive_shells[j], L)});
        libint2::svector<Shell::Contraction> contr_;
        Shell::Contraction contr1;
        contr1.l = L;
        contr1.pure = true;  // libint2 needs solid harmonics for 2c2b integrals
        contr1.coeff = {1.0};
        contr_.push_back(contr1);
        assert(primitive_shells[i].O == primitive_shells[j].O);
        const auto shell = Shell(alpha, contr_, primitive_shells[i].O);
        if (std::find(product_functions.begin(), product_functions.end(),
                      shell) == product_functions.end())
          product_functions.emplace_back(shell);
      }
    }
  }
  return product_functions;
}

/// @brief returns a hash map of shell indices to basis function indices
inline std::vector<size_t> map_shell_to_basis_function(
    const std::vector<libint2::Shell> &shells) {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto &&shell : shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

/// @brief Computes the Coulomb matrix \f$ (\mu|rij^{-1}|\nu) \f$  for a set of
/// shells
/// @param shells set of shells
/// @return Coulomb matrix
inline Eigen::MatrixXd compute_coulomb_matrix(
    const std::vector<Shell> &shells) {
  const auto n = nbf(shells);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
  using libint2::Engine;
  const auto oper = libint2::Operator::coulomb;
  const auto braket = libint2::BraKet::xs_xs;
  Engine engine(oper, max_nprim(shells), max_l(shells), 0,
                std::numeric_limits<scalar_type>::epsilon(),
                libint2::default_params(oper), braket);
  engine.set(ScreeningMethod::Conservative);
  const auto shell2bf = map_shell_to_basis_function(shells);
  const auto &buf = engine.results();
  for (size_t s1 = 0; s1 != shells.size(); ++s1) {
    auto bf1 = shell2bf[s1];
    auto n1 = shells[s1].size();
    for (size_t s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();
      engine.compute(shells[s1], shells[s2]);
      Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }
  return result;
}

/// @brief Splits a set of shells by angular momentum
/// @param shells set of shells
/// @return vector of vectors of shells split by angular momentum L. The i-th
///  vector contains all shells with angular momentum i
inline std::vector<std::vector<Shell>> split_by_L(
    const std::vector<Shell> &shells) {
  int lmax = max_l(shells);
  std::vector<std::vector<Shell>> sorted_shells;
  sorted_shells.resize(lmax + 1);
  for (auto &&shell : shells) {
    auto l = shell.contr[0].l;
    sorted_shells[l].push_back(shell);
  }
  return sorted_shells;
}

/// @brief Computes the reduced set of product functions via pivoted Cholesky
/// decomposition as described in J. Chem. Theory Comput. 2021, 17, 6886−6900.
/// Cholesky decomposition is performed on the Coulomb matrix of the product
/// functions and the pivot indices are sorted in ascending order of column sums
/// of the Coulomb matrix. The reduced set of product functions is then
/// constructed by selecting the product functions corresponding to the pivot
/// indices.
/// @param shells set of shells
/// @param cholesky_threshold threshold for choosing a product function via
/// pivoted Cholesky decomposition
/// @return reduced set of product functions
inline std::vector<Shell> pivoted_cholesky_in_L(
    const std::vector<Shell> &shells, const double cholesky_threshold) {
  const auto n = shells.size();  // number of shells
  std::vector<size_t>
      shell_indices;  // hash map of basis function indices to shell indices
  const auto L = shells[0].contr[0].l;  // all shells must have same L
  const auto nbf = libint2::nbf(
      shells);  // total number of basis functions in vector of shells
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < 2 * L + 1;
         ++j)  // 2L+1 since libint2 strictly uses solid harmonics for 2c2b
      // integrals
      shell_indices.push_back(i);
  }
  assert(shell_indices.size() == nbf);
  const auto C = compute_coulomb_matrix(shells);
  std::vector<size_t> pivot(nbf);
  for (auto i = 0; i < nbf; ++i) {
    pivot[i] = i;
  }
  // set pivot indices in ascending order of off diagonal elements of Coulomb
  // matrix see Phys. Rev. A 101, 032504 (Accurate reproduction of strongly
  // repulsive interatomic potentials)
  Eigen::MatrixXd C_off_diag = C;
  auto col_sum = C_off_diag.colwise().sum();
  // sort pivot indices in ascending order of column sums
  std::sort(pivot.begin(), pivot.end(), [&col_sum](size_t i1, size_t i2) {
    return col_sum[i1] < col_sum[i2];
  });
  // compute Cholesky decomposition
  const auto reduced_pivots = pivoted_cholesky(C, cholesky_threshold, pivot);

  std::vector<Shell> reduced_shells;
  for (size_t i = 0; i < reduced_pivots.size(); ++i) {
    // check if the reduced shell is already in reduced shells
    if (std::find(reduced_shells.begin(), reduced_shells.end(),
                  shells[shell_indices[reduced_pivots[i]]]) ==
        reduced_shells.end())
      reduced_shells.push_back(shells[shell_indices[reduced_pivots[i]]]);
  }
  return reduced_shells;
}
}  // namespace detail

/// @brief This class produces density fitting basis sets for an atom from
/// products of AO basis functions and eliminates linearly dependent functions
/// via pivoted Cholesky decomposition see: J. Chem. Theory Comput. 2021, 17,
/// 6886−6900 (Straightforward and Accurate Automatic Auxiliary Basis Set
/// Generation for Molecular Calculations with Atomic Orbital Basis Sets)
class DFBasisSetGenerator {
 public:
  /// @brief constructor for DFBasisSetGenerator class, generates density
  /// fitting basis set from products of AO basis functions of and atom see: J.
  /// Chem. Theory Comput. 2021, 17, 6886−6900 (Straightforward and Accurate
  /// Automatic Auxiliary Basis Set Generation for Molecular Calculations with
  /// Atomic Orbital Basis Sets)
  /// @param obs_name name of basis set for AO functions
  /// @param atoms vector of atoms
  /// @param cholesky_threshold threshold for threshold for  pivoted Cholesky
  /// decomposition to be performed on the Coulomb matrix of the product
  /// functions
  DFBasisSetGenerator(std::string obs_name, const Atom &atom,
                      const double cholesky_threshold = 1e-7) {
    // get AO basis shells for each atom
    const auto atom_bs = BasisSet(obs_name, {atom});
    const auto obs_shells = atom_bs.shells();
    // get primitive shells from AO functions
    const auto primitive_shells = detail::uncontract(obs_shells);
    // compute candidate shells
    candidate_shells_ = detail::product_functions(primitive_shells);
    cholesky_threshold_ = cholesky_threshold;
  }

  /// @brief constructor for DFBasisSetGenerator class, generates density
  /// fitting basis set from products of AO shells provided by user
  /// @param shells vector of vector of shells for each atom
  /// @param cholesky_threshold threshold for  pivoted Cholesky decomposition to
  /// be performed on the Coulomb matrix of the product functions
  DFBasisSetGenerator(const std::vector<Shell> &shells,
                      const double cholesky_threshold = 1e-7) {
    const auto primitive_shells = detail::uncontract(shells);
    candidate_shells_ = detail::product_functions(primitive_shells);
    cholesky_threshold_ = cholesky_threshold;
  }

  DFBasisSetGenerator() = default;

  ~DFBasisSetGenerator() = default;

  /// @brief returns the candidate shells (full set of product functions)
  std::vector<Shell> candidate_shells() { return candidate_shells_; }

  /// @brief returns a set of shells reduced via pivoted Cholesky
  /// decomposition of the Coulomb matrix of the product functions for each
  /// angular momentum L as described in J. Chem. Theory Comput. 2021, 17,
  /// 6886−6900.
  std::vector<Shell> reduced_shells() {
    if (reduced_shells_computed_)
      return reduced_shells_;
    else {
      const auto candidate_splitted_in_L =
          detail::split_by_L(candidate_shells_);
      for (size_t i = 0; i < candidate_splitted_in_L.size(); ++i) {
        std::vector<Shell> reduced_shells_L;
        if (candidate_splitted_in_L[i].size() > 1)
          reduced_shells_L = detail::pivoted_cholesky_in_L(
              candidate_splitted_in_L[i], cholesky_threshold_);
        else
          reduced_shells_L = candidate_splitted_in_L[i];
        reduced_shells_.insert(reduced_shells_.end(), reduced_shells_L.begin(),
                               reduced_shells_L.end());
      }
      reduced_shells_computed_ = true;
    }
    return reduced_shells_;
  }

  /// @brief returns a BasisSet of shells reduced via pivoted Cholesky
  /// decomposition of the Coulomb matrix of the product functions for each
  /// angular momentum L as described in J. Chem. Theory Comput. 2021, 17,
  /// 6886−6900.
  const BasisSet reduced_basis() { return BasisSet(reduced_shells()); }

 private:
  double cholesky_threshold_;
  std::vector<Shell> candidate_shells_;  // full set of product functions
  std::vector<Shell> reduced_shells_;    // reduced set of product functions
  bool reduced_shells_computed_ = false;
};

}  // namespace libint2

#endif /* _libint2_src_lib_libint_dfbs_generator_h_ */
