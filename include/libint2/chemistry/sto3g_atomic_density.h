/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef INCLUDE_LIBINT2_CHEMISTRY_STO3GATOMICDENSITY_H_
#define INCLUDE_LIBINT2_CHEMISTRY_STO3GATOMICDENSITY_H_

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace libint2 {

namespace detail {

/// computes orbital occupation numbers for a subshell of size \c size created
/// by smearing
/// no more than \c ne electrons (corresponds to spherical averaging)
///
/// @param[in,out] occvec occupation vector, increments by \c size on return
/// @param[in] size the size of the subshell
/// @param[in,out] ne the number of electrons, on return contains the number of
/// "remaining" electrons
template <typename Real>
void subshell_occvec(Real*& occvec, size_t size, size_t& ne) {
  const auto ne_alloc = (ne > 2 * size) ? 2 * size : ne;
  ne -= ne_alloc;
  // # of electrons / orbital compute as precisely as possible
  const double ne_per_orb = (ne_alloc % size == 0)
                                ? static_cast<Real>(ne_alloc / size)
                                : (static_cast<Real>(ne_alloc)) / size;
  for (size_t f = 0; f != size; ++f) occvec[f] = ne_per_orb;
  occvec += size;
}

}  // namespace detail

/// @param[in] Z the atomic number of the element
/// @throw if Z > 53
/// @return the number of STO-3G AOs for the element with atomic number \c Z
inline size_t sto3g_num_ao(size_t Z) {
  size_t nao;
  if (Z == 1 || Z == 2)  // H, He
    nao = 1;
  else if (Z <= 10)  // Li - Ne
    nao = 5;         // 2p is included even for Li and Be
  else if (Z <= 18)  // Na - Ar
    nao = 9;         // 3p is included even for Na and Mg
  else if (Z <= 20)  // K, Ca
    nao = 13;        // 4p is included
  else if (Z <= 36)  // Sc - Kr
    nao = 18;
  else if (Z <= 38)  // Rb, Sr
    nao = 22;        // 5p is included
  else if (Z <= 53)  // Y - I
    nao = 27;
  else
    throw std::invalid_argument{
        "STO-3G basis is not defined for elements with Z > 53"};
  return nao;
}

/// @brief computes average orbital occupancies in the ground state of a neutral
///        atoms
/// @throw if Z > 53
/// @return occupation vector corresponding to the ground state electronic
///         configuration of a neutral atom with atomic number \c Z
///         corresponding to the orbital ordering in STO-3G basis
template <typename Real = double>
const std::vector<Real>& sto3g_ao_occupation_vector(size_t Z) {
  static std::vector<Real> occvec(27, 0.0);

  using detail::subshell_occvec;

  occvec.resize(sto3g_num_ao(Z));
  auto* occs_ptr = &occvec[0];
  auto& occs = occs_ptr;

  size_t num_of_electrons = Z;  // # of electrons to allocate

  // neutral atom electronic configurations from NIST:
  // http://www.nist.gov/pml/data/images/illo_for_2014_PT_1.PNG
  subshell_occvec(occs, 1, num_of_electrons);    // 1s
  if (Z > 2) {                                   // Li+
    subshell_occvec(occs, 1, num_of_electrons);  // 2s
    subshell_occvec(occs, 3, num_of_electrons);  // 2p
  }
  if (Z > 10) {                                  // Na+
    subshell_occvec(occs, 1, num_of_electrons);  // 3s
    subshell_occvec(occs, 3, num_of_electrons);  // 3p
  }
  if (18 < Z && Z <= 36) {  // K .. Kr
    // NB 4s is singly occupied in K, Cr, and Cu
    size_t num_of_4s_electrons = (Z == 19 || Z == 24 || Z == 29) ? 1 : 2;
    num_of_electrons -= num_of_4s_electrons;
    subshell_occvec(occs, 1, num_of_4s_electrons);  // 4s

    size_t num_of_4p_electrons =
        std::min(static_cast<decltype(Z)>(6), (Z > 30) ? Z - 30 : 0);
    num_of_electrons -= num_of_4p_electrons;
    subshell_occvec(occs, 3, num_of_4p_electrons);  // 4p

    subshell_occvec(occs, 5, num_of_electrons);  // 3d
  }
  if (36 < Z && Z <= 53) {  // Rb .. I
    // 3d4s4p are fully occupied ...
    subshell_occvec(occs, 1, num_of_electrons);  // 4s
    subshell_occvec(occs, 3, num_of_electrons);  // 4p

    // NB 5s is singly occupied in Rb, Nb, Mo, Ru, Rh, and Ag
    size_t num_of_5s_electrons =
        (Z == 37 || Z == 41 || Z == 42 || Z == 44 || Z == 45 || Z == 47) ? 1
                                                                         : 2;
    num_of_electrons -= num_of_5s_electrons;
    subshell_occvec(occs, 1, num_of_5s_electrons);  // 5s

    size_t num_of_5p_electrons =
        std::min(static_cast<decltype(Z)>(6), (Z > 48) ? Z - 48 : 0);
    num_of_electrons -= num_of_5p_electrons;
    subshell_occvec(occs, 3, num_of_5p_electrons);  // 5p

    subshell_occvec(occs, 5, num_of_electrons);  // 3d
    subshell_occvec(occs, 5, num_of_electrons);  // 4d
  }

  return occvec;
}

}  // namespace libint2

#endif  // INCLUDE_LIBINT2_CHEMISTRY_STO3GATOMICDENSITY_H_
