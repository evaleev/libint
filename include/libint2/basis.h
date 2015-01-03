/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_basis_h_
#define _libint2_src_lib_libint_basis_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <vector>

#include <libint2.h>
#include <libint2/shell.h>
#include <libint2/atom.h>

namespace libint2 {

  /// BasisSet is a slightly decorated \c std::vector of \c libint2::Shell objects.
  class BasisSet : public std::vector<libint2::Shell> {
    public:
      BasisSet() : name_(""), nbf_(-1), max_nprim_(0), max_l_(-1) {}
      BasisSet(const BasisSet&) = default;
      BasisSet(BasisSet&&) = default;
      ~BasisSet() = default;

      BasisSet(std::string name,
               const std::vector<Atom>& atoms) : name_(std::move(name)) {

        // read in the library file contents
        std::string file_dot_g94(SRCDATADIR);
        file_dot_g94 += "/" + canonicalize_name(name_) + ".g94";
        std::vector<std::vector<libint2::Shell>> ref_shells = read_g94_basis_library(file_dot_g94);

        // for each atom find the corresponding basis
        for(auto a=0; a<atoms.size(); ++a) {

          auto Z = atoms[a].atomic_number;
          if (ref_shells[Z].empty())
            throw std::string("did not find the basis for this Z in ") + file_dot_g94;

          for(auto s: ref_shells[Z]) {
            this->push_back(std::move(s));
            this->back().move({{atoms[a].x, atoms[a].y, atoms[a].z}});
          }

        }

        // technical step: rescale contraction coefficients to include primitive normalization coefficients
        for(auto& s: *this) {
          s.renorm();
        }

        init();
      }

      /// @return the number of basis functions in the basis; -1 if uninitialized
      long nbf() const {
        return nbf_;
      }
      /// @return the maximum number of primitives in a contracted Shell, i.e. maximum contraction length; 0 if uninitialized
      size_t max_nprim() const {
        return max_nprim_;
      }
      /// @return the maximum angular momentum of a contraction; -1 if uninitialized
      size_t max_l() const {
        return max_l_;
      }
      /// @return the map from shell index to index of the first basis function from this shell
      /// \note basis functions are ordered as shells, i.e. shell2bf[i] >= shell2bf[j] iff i >= j
      const std::vector<size_t>& shell2bf() const {
        return shell2bf_;
      }

    private:
      std::string name_;
      long nbf_;
      size_t max_nprim_;
      int max_l_;
      std::vector<size_t> shell2bf_;

      void init() {
        nbf_ = nbf(*this);
        max_nprim_ = max_nprim(*this);
        max_l_ = max_l(*this);
        shell2bf_ = compute_shell2bf(*this);
      }

      struct canonicalizer {
          char operator()(char c) {
            char cc = ::tolower(c);
            switch (cc) {
              case '/': cc = 'I'; break;
            }
            return cc;
          }
      };

      static std::string canonicalize_name(const std::string& name) {
        auto result = name;
        std::transform(name.begin(), name.end(),
                       result.begin(), BasisSet::canonicalizer());
        return result;
      }

    public:

      static std::vector<std::vector<libint2::Shell>> read_g94_basis_library(std::string file_dot_g94) {

        std::cout << "Will read basis set from " << file_dot_g94 << std::endl;
        std::ifstream is(file_dot_g94);
        assert(is.good());
        std::vector<std::vector<libint2::Shell>> ref_shells(118); // 118 = number of chemical elements

        std::string comment, rest;

        // skip till first basis
        while(std::getline(is, comment) && comment != "****") {
        }

        size_t Z;
        auto nextbasis = true, nextshell = false;
        while(std::getline(is, comment) && comment != "") {
          if (comment == "****") {
            nextbasis = true;
            nextshell = false;
            continue;
          }
          if (nextbasis) {
            nextbasis = false;
            std::istringstream iss(comment);
            std::string elemsymbol;
            iss >> elemsymbol >> rest;

            bool found = false;
            using libint2::chemistry::element_info;
            for(const auto& e: element_info) {
              if (strcaseequal(e.symbol, elemsymbol)) {
                Z = e.Z;
                found = true;
                break;
              }
            }
            if (not found) {
              std::ostringstream oss;
              oss << "in file " << file_dot_g94
                  << " found G94 basis set for element symbol \""
                  << elemsymbol << "\", not found in Periodic Table.";
              throw std::runtime_error(oss.str());
            }

            nextshell = true;
            continue;
          }
          if (nextshell) {
            std::istringstream iss(comment);
            std::string amlabel;
            unsigned nprim;
            iss >> amlabel >> nprim >> rest;
            if (amlabel != "SP" && amlabel != "sp") {
              assert(amlabel.size() == 1);
              auto l = Shell::am_symbol_to_l(amlabel[0]);
              std::vector<double> exps;
              std::vector<double> coeffs;
              for(auto p = 0; p!=nprim; ++p) {
                std::getline(is, comment);
                std::istringstream iss(comment);
                double e, c;
                iss >> e >> c;
                exps.emplace_back(e);
                coeffs.emplace_back(c);
              }
              auto pure = l>1;
              ref_shells[Z].push_back(
                  libint2::Shell{
                    exps,
                    {
                      {l, pure, coeffs}
                    },
                    {{0.0, 0.0, 0.0}}
                  }
              );
            }
            else { // split the SP shells
              std::vector<double> exps;
              std::vector<double> coeffs_s, coeffs_p;
              for(auto p = 0; p!=nprim; ++p) {
                std::getline(is, comment);
                std::istringstream iss(comment);
                double e, c1, c2;
                iss >> e >> c1 >> c2;
                exps.emplace_back(e);
                coeffs_s.emplace_back(c1);
                coeffs_p.emplace_back(c2);
              }
              ref_shells[Z].push_back(
                  libint2::Shell{exps,
                {
                 {0, false, coeffs_s}
                },
                {{0.0, 0.0, 0.0}}
              }
              );
              ref_shells[Z].push_back(
                  libint2::Shell{ exps,
                {
                 {1, false, coeffs_p}
                },
                {{0.0, 0.0, 0.0}}
              }
              );
            }
          }
        }

        return ref_shells;
      }

      static size_t nbf(const std::vector<libint2::Shell>& shells) {
        size_t n = 0;
        for (const auto& shell: shells)
          n += shell.size();
        return n;
      }

      static size_t max_nprim(const std::vector<libint2::Shell>& shells) {
        size_t n = 0;
        for (auto shell: shells)
          n = std::max(shell.nprim(), n);
        return n;
      }

      static int max_l(const std::vector<libint2::Shell>& shells) {
        int l = 0;
        for (auto shell: shells)
          for (auto c: shell.contr)
            l = std::max(c.l, l);
        return l;
      }

      static std::vector<size_t> compute_shell2bf(const std::vector<libint2::Shell>& shells) {
        std::vector<size_t> result;
        result.reserve(shells.size());

        size_t n = 0;
        for (auto shell: shells) {
          result.push_back(n);
          n += shell.size();
        }

        return result;
      }

  }; // BasisSet

} // namespace libint2

#endif /* _libint2_src_lib_libint_basis_h_ */
