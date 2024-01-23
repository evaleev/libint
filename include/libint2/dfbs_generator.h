/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
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
#include <libint2/soad_fock.h>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <tuple>

namespace libint2 {

namespace detail {

inline double fd_occupation(double e) { return 1. / (1. + std::exp(e)); }

/// computes a map for contraction vector in a shell
/// @param[in] shell libint shell of which a map is to be computed
/// @return vector of lower bound indices to shell contraction vector
inline std::vector<size_t> shell_hashmap(const libint2::Shell &shell) {
  std::vector<size_t> shell_map;
  auto contr = shell.contr;
  size_t n = 0;
  for (const auto &a : contr) {
    shell_map.push_back(n);
    n += a.size();
  }
  return shell_map;
}

/// computes contraction matrix for an uncontracted shell to a contracted shell
/// @param[in] unc_shell a uncontracted shell
/// @param[in] contr_shell a contracted shell
/// @return a contraction matrix from uncontracted shell to contracted shell
inline Eigen::MatrixXd shell_contraction_matrix(
    const libint2::Shell &unc_shell, const libint2::Shell &contr_shell) {
  Eigen::MatrixXd result(unc_shell.size(), contr_shell.size());
  result.fill(0.0);
  // check to see if the shell belongs to same center
  if (unc_shell.O == contr_shell.O) {
    auto contr_exps = contr_shell.alpha;
    auto unc_exp = unc_shell.alpha[0];
    // check if the uncontracted primitive exponent present in the contracted
    // shell
    if (std::find(contr_exps.begin(), contr_exps.end(), unc_exp) !=
        contr_exps.end()) {
      // get index of the exponent in contracted shell
      std::int64_t exponent_index = 0;
      for (auto i = 0; i < contr_exps.size(); i++) {
        if (unc_exp == contr_exps[i]) exponent_index = i;
      }
      // create hashmaps for all contractions in contracted shell
      auto unc_hashmap = shell_hashmap(unc_shell);
      auto contr_hashmap = shell_hashmap(contr_shell);
      // iterate through contractions in contracted shell
      for (auto p1 = 0; p1 < unc_shell.contr.size(); p1++) {
        for (auto p2 = 0; p2 < contr_shell.contr.size(); p2++) {
          auto unc_contr = unc_shell.contr[p1];
          auto contr_contr = contr_shell.contr[p2];
          auto n1 = unc_contr.size();
          auto n2 = contr_contr.size();
          // check to see if the contraction belongs to same am (l)
          if (unc_contr.l == contr_contr.l) {
            Eigen::MatrixXd block(unc_contr.size(), contr_contr.size());
            block.fill(0.0);
            // fill diagonal elements with transformation coefficient
            for (auto diag = 0; diag < contr_contr.size(); diag++) {
              block(diag, diag) =
                  contr_contr.coeff[exponent_index] / unc_contr.coeff[0];
            }
            result.block(unc_hashmap[p1], contr_hashmap[p2], n1, n2) = block;
          }
        }
      }
    }
    return result;
  }

  else
    return result;
}

inline Eigen::MatrixXd contraction_matrix(const BasisSet &primitive_basis,
                                          const BasisSet &contracted_basis) {
  const auto &unc_shell2bf = primitive_basis.shell2bf();
  const auto &contr_shell2bf = contracted_basis.shell2bf();
  const auto &unc_shells = primitive_basis.shells();
  const auto &contr_shells = contracted_basis.shells();
  Eigen::MatrixXd result(primitive_basis.nbf(), contracted_basis.nbf());
  result.fill(0.0);
  for (auto s1 = 0; s1 < unc_shells.size(); s1++) {
    for (auto s2 = 0; s2 < contr_shells.size(); s2++) {
      auto n1 = unc_shells[s1].size();
      auto n2 = contr_shells[s2].size();
      result.block(unc_shell2bf[s1], contr_shell2bf[s2], n1, n2) =
          shell_contraction_matrix(unc_shells[s1], contr_shells[s2]);
    }
  }

  return result;
}

/// Computes uncontracted shells from a vector of shells
/// @param[in] cluster a vector of shells
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

/// @brief returns \Gamma(x)  of x
inline double gamma_function(const double x) { return std::tgamma(x); }

/// @brief return effective exponent of product of two primitive shells
/// @param shell1 first shell
/// @param shell2 second shell
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

/// @brief creates a set of product functions from a set of primitive shells
/// @param primitive_shells set of primitive shells
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

inline std::vector<std::tuple<Shell, double>> product_functions(
    const std::vector<std::tuple<Shell, double>> &shells_weights) {
  std::vector<std::tuple<Shell, double>> product_functions;
  for (size_t i = 0; i < shells_weights.size(); ++i) {
    for (size_t j = 0; j <= i; ++j) {
      const auto si = get<0>(shells_weights[i]);
      const auto sj = get<0>(shells_weights[j]);
      const auto li = si.contr[0].l;
      const auto lj = sj.contr[0].l;
      for (auto L = std::abs(li - lj); L <= li + lj; L++) {
        const auto alpha = libint2::svector<double>({alpha_eff(si, sj, L)});
        libint2::svector<Shell::Contraction> contr_;
        Shell::Contraction contr1;
        contr1.l = L;
        contr1.pure = true;  // libint2 needs solid harmonics for 2c2b integrals
        contr1.coeff = {1.0};
        contr_.push_back(contr1);
        assert(si.O == sj.O);
        const auto shell = Shell(alpha, contr_, si.O);
        const auto weight =
            get<1>(shells_weights[i]) * get<1>(shells_weights[j]);
        std::tuple<Shell, double> shell_weight{shell, weight};
        if (std::find(product_functions.begin(), product_functions.end(),
                      shell_weight) == product_functions.end())
          product_functions.emplace_back(shell_weight);
      }
    }
  }
  return product_functions;
}

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

/// @brief computes 2 indexed integrals for an operator from a set of shells
/// @param op operator
/// @param shells set of shells
/// @param atoms vector of atoms
/// @return Coulomb matrix
inline Eigen::MatrixXd compute_2indexed_ints(const Operator &op,
                                             const std::vector<Shell> &shells,
                                             const std::vector<Atom> atoms) {
  const auto n = nbf(shells);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
  using libint2::Engine;
  Engine engine(op, max_nprim(shells), max_l(shells));
  if (op == libint2::Operator::coulomb) {
    engine.set(BraKet::xs_xs);
    engine.set(ScreeningMethod::Conservative);
  }
  if (op == libint2::Operator::nuclear) {
    engine.

        set_params(libint2::make_point_charges(atoms));
  }
  const auto shell2bf = detail::map_shell_to_basis_function(shells);
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

/// @brief Sorts a vector of shells by angular momentum
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

inline std::vector<std::vector<std::tuple<Shell, double>>> split_by_L(
    const std::vector<std::tuple<Shell, double>> &shells_weights_vec,
    int lmax) {
  std::vector<std::vector<std::tuple<Shell, double>>> result;
  result.resize(lmax + 1);
  for (auto &&candidate : shells_weights_vec) {
    auto l = get<0>(candidate).contr[0].l;
    result[l].push_back(candidate);
  }
  return result;
}

/// @brief computes the reduced set of product functions via pivoted Cholesky
/// decomposition
/// @param shells_weights set of tuples of shells and their weights
/// @param cholesky_threshold threshold for choosing a product function via
/// pivoted Cholesky decomposition
/// @return reduced set of product functions

inline std::vector<Shell> shell_pivoted_cholesky(
    const std::vector<std::tuple<Shell, double>> &shells_weights,
    const double cholesky_threshold) {
  std::vector<Shell> shells;
  std::vector<double> weights;

  for (auto &&shell_weight : shells_weights) {
    const auto shell = get<0>(shell_weight);
    shells.push_back(shell);
    const auto L = shell.contr[0].l;
    for (size_t l = 0; l < 2 * L + 1; ++l) {
      weights.push_back(get<1>(shell_weight));
    }
  }

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
  Atom dummy_atom;
  auto C =
      compute_2indexed_ints(libint2::Operator::coulomb, shells, {dummy_atom});

  for (size_t i = 0; i < C.rows(); ++i) {
    for (size_t j = 0; j < C.cols(); ++j) {
      C(i, j) *= std::sqrt(weights[i] * weights[j]);
    }
  }

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

/// @brief class produces density fitting basis sets from products of AO basis
/// functions eliminates linearly dependent functions via pivoted Cholesky
/// decomposition see: J. Chem. Theory Comput. 2021, 17, 6886−6900
/// (Straightforward and Accurate Automatic Auxiliary Basis Set Generation for
/// Molecular Calculations with Atomic Orbital Basis Sets)
class DFBasisSetGenerator {
 public:
  /// @brief constructor for DFBS generator class, generates density fitting
  /// basis set from products of AO basis functions see: J. Chem. Theory Comput.
  /// 2021, 17, 6886−6900 (Straightforward and Accurate Automatic Auxiliary
  /// Basis Set Generation for Molecular Calculations with Atomic Orbital Basis
  /// Sets)
  /// @param obs_name name of basis set for AO functions
  /// @param atoms vector of atoms
  /// @param cholesky_threshold threshold for choosing a product functions via
  /// pivoted Cholesky decomposition
  DFBasisSetGenerator(std::string obs_name, const Atom &atom,
                      const double cholesky_threshold = 1e-7,
                      bool use_weights = false,
                      std::string minbs_name = "MINI") {
    // get AO basis shells for each atom
    atom_ = atom;
    obs_ = BasisSet(obs_name, {atom_});
    minbs_ = BasisSet(minbs_name, {atom_});
    cholesky_threshold_ = cholesky_threshold;
    use_weights_ = use_weights;
    if (use_weights_) compute_weights();
    candidates_ = candidates();
  }
  /// @brief constructor for DFBS generator class, generates density fitting
  /// basis set from products of AO shells provided by user
  /// @param cluster vector of vector of shells for each atom
  /// @param cholesky_threshold threshold for choosing a product functions via
  /// pivoted Cholesky decomposition
  DFBasisSetGenerator(const std::vector<Shell> &shells, const Atom &atom,
                      const double cholesky_threshold = 1e-7,
                      bool use_weights = false,
                      const std::vector<Shell> &mini_shells = {}) {
    atom_ = atom;
    obs_ = BasisSet(shells);
    minbs_ = BasisSet(mini_shells);
    cholesky_threshold_ = cholesky_threshold;
    use_weights_ = use_weights;
    if (use_weights_) {
      compute_weights();
    }
    candidates_ = candidates();
  }

  DFBasisSetGenerator() = default;

  ~DFBasisSetGenerator() = default;

  /// @brief returns the reduced shells (reduced set of product functions)
  /// computed via pivoted Cholesky decomposition
  std::vector<Shell> reduced_shells() {
    if (reduced_shells_computed_)
      return reduced_shells_;
    else {
      const auto candidates_in_L =
          detail::split_by_L(candidates_, 2 * max_l(obs_.shells()));
      for (size_t i = 0; i < candidates_in_L.size(); ++i) {
        std::vector<Shell> reduced_shells_L;
        if (candidates_in_L[i].size() > 1) {
          reduced_shells_L = detail::shell_pivoted_cholesky(
              candidates_in_L[i], cholesky_threshold_);
        } else
          reduced_shells_L = {get<0>(candidates_in_L[i][0])};

        std::cout << "Number of shells in L = " << i << " is "
                  << reduced_shells_L.size() << std::endl;
        reduced_shells_.insert(reduced_shells_.end(), reduced_shells_L.begin(),
                               reduced_shells_L.end());
      }
      reduced_shells_computed_ = true;
      return reduced_shells_;
    }
  }

  /// @brief returns a set of weights of primitive obs shells
  void compute_weights() {
    const auto T = detail::compute_2indexed_ints(libint2::Operator::kinetic,
                                                 obs_.shells(), {atom_});
    const auto V = detail::compute_2indexed_ints(libint2::Operator::nuclear,
                                                 obs_.shells(), {atom_});
    auto H = T + V;
    auto F_soad = compute_soad_fock(obs_, minbs_, {atom_});
    auto F = H + F_soad;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(F);
    auto C = es.eigenvectors();
    auto E = es.eigenvalues();
    auto primitive_bs = BasisSet(detail::uncontract(obs_.shells()));
    const auto contraction_matrix =
        detail::contraction_matrix(primitive_bs, obs_);
    auto C_unc = contraction_matrix * C;

    const auto nocc = atom_.atomic_number / 2;
    auto C_occ = C_unc.leftCols(nocc);
    auto D = C_occ * C_occ.transpose();

    // compute square of coefficients
    //      Eigen::MatrixXd C_unc_sqr(C_unc.rows(),C_unc.cols());
    //      for(size_t i=0;i<C_unc.rows();++i){
    //          for(size_t j=0;j<C_unc.cols();++j){
    //              C_unc_sqr(i,j) = C_unc(i,j)*C_unc(i,j);
    //          }
    //      }
    //
    //      std::vector<double> prim_func_weights;
    //
    //      for(size_t i=0;i<C_unc.rows();++i){
    //          double max_weight = 0.0;
    //          for(size_t j=0;j<C_unc.cols();++j){
    //              max_weight = std::max(max_weight,detail::fd_occupation(E(j))
    //              * C_unc_sqr(i,j));
    //          }
    //          prim_func_weights.push_back(max_weight);
    //      }

    auto prim_shells = primitive_bs.shells();
    size_t nfns = 0;
    for (size_t i = 0; i < prim_shells.size(); ++i) {
      double max_weight = 0.0;
      for (size_t j = 0; j < prim_shells[i].size(); ++j) {
        // max_weight = std::max(max_weight, prim_func_weights[nfns]);
        max_weight = std::max(max_weight, D.diagonal()(nfns));
        nfns++;
      }
      prim_obs_weights_.push_back(max_weight);
    }
  }

  const std::vector<std::tuple<Shell, double>> candidates() {
    const auto prim_shells = detail::uncontract(obs_.shells());

    // create tuple of shells and their weights
    std::vector<std::tuple<Shell, double>> prim_shell_weights;
    for (size_t i = 0; i < prim_shells.size(); ++i) {
      std::tuple<Shell, double> shell_weight;
      get<0>(shell_weight) = prim_shells[i];
      if (use_weights_)
        get<1>(shell_weight) = prim_obs_weights_[i];
      else
        get<1>(shell_weight) = 1.;
      prim_shell_weights.push_back(shell_weight);
    }
    const auto candidates = detail::product_functions(prim_shell_weights);
    return candidates;  // return candidates once done
  }

  /// @brief returns the reduced basis set (reduced set of product
  /// functions)
  /// computed via pivoted Cholesky decomposition
  const BasisSet reduced_basis() { return BasisSet(reduced_shells()); }

 private:
  BasisSet obs_;
  BasisSet minbs_;
  Atom atom_;
  double cholesky_threshold_;
  std::vector<std::tuple<Shell, double>>
      candidates_;  // full set of product functions and their weights
  std::vector<Shell> reduced_shells_;  // reduced set of product functions
  bool reduced_shells_computed_ = false;
  bool use_weights_;
  std::vector<double> prim_obs_weights_;
};

}  // namespace libint2

#endif /* _libint2_src_lib_libint_dfbs_generator_h_ */
