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

#include <algorithm>
#include <libint2.h>
#include <math.h>
#include <libint2/shell.h>
#include <libint2/basis.h>
#include <libint2/atom.h>
#include <libint2/boys.h>
#include <libint2/pivoted_cholesky.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace libint2 {

    namespace datail {

        /// @brief returns \Gamma(x)  of x
        double gamma_function(const double &x) {
            return std::tgamma(x);
        }

        /// @brief return effective exponent of product of two primitive shells
        /// @param shell1 first shell
        /// @param shell2 second shell
        /// @param L total angular momentum of product function
        /// @return effective exponent of product function
        double alpha_eff(const Shell &shell1, const Shell &shell2, const int &L) {
            auto alpha1 = shell1.alpha[0];
            auto alpha2 = shell2.alpha[0];
            auto l1 = shell1.contr[0].l;
            auto l2 = shell2.contr[0].l;
            auto prefactor = std::pow((gamma_function(L + 2.) * gamma_function(l1 + l2 + 1.5)) /
                                      (gamma_function(l1 + l2 + 2.) * gamma_function(L + 1.5)),
                                      2.);
            return prefactor * (alpha1 + alpha2);
        }

        /// @brief creates a set of product functions from a set of primitive shells
        /// @param primitive_shells set of primitive shells
        std::vector<Shell> product_functions(const std::vector<Shell> &primitive_shells) {
            std::vector<Shell> product_functions;
            for (auto i = 0; i < primitive_shells.size(); ++i) {
                for (auto j = 0; j <= i; ++j) {
                    auto li = primitive_shells[i].contr[0].l;
                    auto lj = primitive_shells[j].contr[0].l;
                    for (auto L = std::abs(li - lj); L <= li + lj; L++) {
                        auto alpha = libint2::svector<double>({alpha_eff(primitive_shells[i], primitive_shells[j], L)});
                        libint2::svector<Shell::Contraction> contr_;
                        Shell::Contraction contr1;
                        contr1.l = L;
                        contr1.pure = true; // libint2 needs solid harmonics for 2c2b integrals
                        contr1.coeff = {1.0};
                        contr_.push_back(contr1);
                        assert(primitive_shells[i].O == primitive_shells[j].O);
                        auto shell = Shell(alpha, contr_, primitive_shells[i].O);
                        if (std::find(product_functions.begin(), product_functions.end(),
                                      shell) == product_functions.end())
                            product_functions.emplace_back(shell);
                    }
                }
            }
            return product_functions;
        }

        /// @brief creates a set of candidate product shells from a set of primitive shells
        /// @param primitive_shells set of primitive shells
        /// @return set of candidate product shells
        std::vector<std::vector<Shell>> candidate_functions(const std::vector<std::vector<Shell>> &primitive_shells) {
            std::vector<std::vector<Shell>> candidate_functions;
            for (auto i = 0; i < primitive_shells.size(); ++i) {
                candidate_functions.push_back(product_functions(primitive_shells[i]));
            }
            return candidate_functions;
        }

        /// @brief returns a hash map of shell indices to basis function indices
        std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell> &shells) {
            std::vector<size_t> result;
            result.reserve(shells.size());

            size_t n = 0;
            for (auto shell: shells) {
                result.push_back(n);
                n += shell.size();
            }

            return result;
        }

        /// @brief computes the Coulomb matrix (\mu|rij^{-1}|\nu)  for a set of shells
        /// @param shells set of shells
        /// @return Coulomb matrix
        Eigen::MatrixXd compute_coulomb_matrix(const std::vector<Shell> &shells) {
            const auto n = nbf(shells);
            Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
            using libint2::Engine;
            Engine engine(libint2::Operator::coulomb, max_nprim(shells), max_l(shells), 0);
            engine.set(BraKet::xs_xs);
            engine.set(ScreeningMethod::Conservative);
            auto shell2bf = map_shell_to_basis_function(shells);
            const auto &buf = engine.results();
            for (auto s1 = 0; s1 != shells.size(); ++s1) {
                auto bf1 = shell2bf[s1];
                auto n1 = shells[s1].size();
                for (auto s2 = 0; s2 <= s1; ++s2) {
                    auto bf2 = shell2bf[s2];
                    auto n2 = shells[s2].size();
                    engine.compute(shells[s1], shells[s2]);
                    Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
                    result.block(bf1, bf2, n1, n2) = buf_mat;
                    if (s1 != s2)
                        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
                }
            }
            return result;
        }

        /// @brief Sorts a vector of shells by angular momentum
        std::vector<std::vector<Shell>> split_by_L(const std::vector<Shell> &shells) {
            int lmax = max_l(shells);
            std::vector<std::vector<Shell>> sorted_shells;
            sorted_shells.resize(lmax + 1);
            for (auto shell: shells) {
                auto l = shell.contr[0].l;
                sorted_shells[l].push_back(shell);
            }
            return sorted_shells;
        }

        /// @brief computes the reduced set of product functions via pivoted Cholesky decomposition
        /// @param shells set of shells
        /// @param cholesky_threshold threshold for choosing a product function via pivoted Cholesky decomposition
        /// @return reduced set of product functions
        std::vector<Shell> shell_pivoted_cholesky(const std::vector<Shell> shells, const double cholesky_threshold) {

            auto n = shells.size(); // number of shells
            std::vector<size_t> shell_indices; // hash map of basis function indices to shell indices
            auto L = shells[0].contr[0].l; // all shells must have same L
            auto nbf = libint2::nbf(shells);  // total number of basis functions in vector of shells
            for (auto i = 0; i < n; ++i) {
                for (auto j = 0; j < 2 * L + 1; ++j)  // 2L+1 since libint2 strictly uses solid harmonics for 2c2b integrals
                    shell_indices.push_back(i);
            }
            assert(shell_indices.size() == nbf);
            auto C = compute_coulomb_matrix(shells);
            std::vector<size_t> pivot(nbf);
            for(auto i=0;i<nbf;++i){
                pivot[i] = i;
            }
            // set pivot indices in ascending order of off diagonal elements of Coulomb matrix
            // see Phys. Rev. A 101, 032504 (Accurate reproduction of strongly repulsive interatomic potentials)
            Eigen::MatrixXd C_off_diag = C;
            auto col_sum = C_off_diag.colwise().sum();
            // sort pivot indices in ascending order of column sums
            std::sort(pivot.begin(), pivot.end(), [&col_sum](size_t i1, size_t i2) { return col_sum[i1] < col_sum[i2]; });
            // compute Cholesky decomposition
            auto reduced_pivots = pivoted_cholesky(C, cholesky_threshold, pivot);

            std::vector<Shell> reduced_shells;
            for (auto i = 0; i < reduced_pivots.size(); ++i) {
                // check if the reduced shell is already in reduced shells
                if (std::find(reduced_shells.begin(), reduced_shells.end(),
                              shells[shell_indices[reduced_pivots[i]]]) == reduced_shells.end())
                reduced_shells.push_back(shells[shell_indices[reduced_pivots[i]]]);
            }
            return reduced_shells;
        }

    }// namespace detail

    class DFBasisSetGenerator {
    public:
        /// @brief constructor for DFBS generator class, generates density fitting basis set from products of AO basis functions
        /// see: J. Chem. Theory Comput. 2021, 17, 6886−6900 (Straightforward and Accurate Automatic Auxiliary Basis Set Generation for Molecular Calculations with Atomic Orbital Basis Sets)
        /// @param obs_name name of basis set for AO functions
        /// @param atoms vector of atoms
        /// @param cholesky_threshold threshold for choosing a product functions via pivoted Cholesky decomposition
        DFBasisSetGenerator(std::string obs_name,
                            const std::vector<Atom> &atoms, const double cholesky_thershold = 1e-7) : obs_name_(
                std::move(obs_name)), atoms_(std::move(atoms)) {
            std::vector<std::vector<Shell>> obs_shell_vec;
            std::vector<std::vector<Shell>> primitive_cluster;
            // get AO basis shells for each atom
            for (auto atom: atoms) {
                auto atom_bs = BasisSet(obs_name_, {atom});
                obs_shell_vec.emplace_back(atom_bs.shells());
            }
            // get primitive shells from AO functions
            for (auto obs_shells: obs_shell_vec) {
                primitive_cluster.emplace_back(uncontract(obs_shells));
            }

            //compute candidate shells
            candidate_shells_ = datail::candidate_functions(primitive_cluster);
            cholesky_threshold_ = cholesky_thershold;
        }

        /// @brief constructor for DFBS generator class, generates density fitting basis set from products of AO shells provided by user
        /// @param cluster vector of vector of shells for each atom
        /// @param cholesky_threshold threshold for choosing a product functions via pivoted Cholesky decomposition
        DFBasisSetGenerator(std::vector<std::vector<Shell>> cluster, const double cholesky_thershold = 1e-7) {
            std::vector<std::vector<Shell>> primitive_cluster;
            for (auto i = 0; i < cluster.size(); ++i) {
                primitive_cluster.emplace_back(uncontract(cluster[i]));
            }
            candidate_shells_ = datail::candidate_functions(primitive_cluster);
            cholesky_threshold_ = cholesky_thershold;
        }

        DFBasisSetGenerator() = default;

        ~DFBasisSetGenerator() = default;

        /// @brief returns the candidate shells (full set of product functions)
        std::vector<std::vector<Shell>> candidate_shells() {
            return candidate_shells_;
        }

        /// @brief returns the candidate basis set (full set of product functions)
        /// @warning generates huge and heavily linearly dependent basis sets
        BasisSet product_basis(){
            std::vector<Shell> product_shells;
            for(auto shells: candidate_shells_){
                product_shells.insert(product_shells.end(), shells.begin(), shells.end());
            }
            return BasisSet(std::move(product_shells));
        }

        /// @brief returns the candidate shells sorted by angular momentum
        std::vector<std::vector<std::vector<Shell>>> candidates_splitted_in_L() {
            std::vector<std::vector<std::vector<Shell>>> sorted_shells;
            for (auto shells: candidate_shells_) {
                sorted_shells.push_back(datail::split_by_L(shells));
            }
            return sorted_shells;
        }

        /// @brief returns the reduced shells (reduced set of product functions) computed via pivoted Cholesky decomposition
        std::vector<std::vector<Shell>> reduced_shells() {
            if (reduced_shells_.size() != 0)
                return reduced_shells_;
            else {
                auto candidate_splitted_in_L = candidates_splitted_in_L();
                for (auto i = 0; i < candidate_splitted_in_L.size(); ++i) {
                    std::vector<Shell> atom_shells;
                    for (auto j = 0; j < candidate_splitted_in_L[i].size(); ++j) {
                        auto reduced_shells = datail::shell_pivoted_cholesky(candidate_splitted_in_L[i][j],
                                                                             cholesky_threshold_);
                        atom_shells.insert(atom_shells.end(), reduced_shells.begin(), reduced_shells.end());
                    }
                    reduced_shells_.push_back(atom_shells);
                }
            }
            return reduced_shells_;
        }

        /// @brief returns the reduced basis set (reduced set of product functions) computed via pivoted Cholesky decomposition
        BasisSet reduced_basis() {
            auto reduced_cluster = reduced_shells();
            std::vector<Shell> reduced_shells;
            for (auto shells: reduced_cluster) {
                reduced_shells.insert(reduced_shells.end(), shells.begin(), shells.end());
            }
            return BasisSet(std::move(reduced_shells));
        }


    private:
        std::string obs_name_;  //name of AO basis set
        std::vector<Atom> atoms_; //vector of atoms
        double cholesky_threshold_;
        std::vector<std::vector<Shell>> candidate_shells_;  //full set of product functions
        std::vector<std::vector<Shell>> reduced_shells_;    //reduced set of product functions

    };

} // namespace libint2

#endif /* _libint2_src_lib_libint_dfbs_generator_h_ */
