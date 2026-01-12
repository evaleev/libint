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

#ifndef LIBINT2_PYTHON_BASIS_H
#define LIBINT2_PYTHON_BASIS_H

#include <libint2/basis.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <boost/container/small_vector.hpp>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <typename T, std::size_t N>
struct type_caster<boost::container::small_vector<T, N> >
    : list_caster<boost::container::small_vector<T, N>, T> {};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

#include <tuple>
#include <vector>

namespace py = pybind11;

namespace libint2::python::basis {

template <class... Args>
std::unique_ptr<BasisSet> make_basis(const Args &...args) {
  // py::print("make_basis:", args...);
  auto basis = std::make_unique<BasisSet>(args...);
  assert(!basis->empty());
  assert(basis->nbf());
  return basis;
}

template <>
inline std::unique_ptr<BasisSet> make_basis(
    const std::map<int, std::vector<Shell> > &element_basis,
    const std::vector<Atom> &atoms) {
  std::vector<Shell> basis;
  for (auto a : atoms) {
    auto Z = a.atomic_number;
    auto it = element_basis.find(Z);
    if (it == element_basis.end() || it->second.empty()) {
      throw std::logic_error("Missing element Z=" + std::to_string(Z));
    }
    for (Shell s : it->second) {
      s.move({{a.x, a.y, a.z}});
      basis.push_back(s);
    }
  }
  return make_basis(basis);
}

inline size_t __len__(const BasisSet &b) { return b.size(); }

inline const Shell &__getitem__(const BasisSet &b, size_t idx) {
  if (idx >= b.size()) throw py::index_error();
  return b.at(idx);
}

inline auto __iter__(const BasisSet &b) {
  return py::make_iterator(b.begin(), b.end());
}

inline auto enumerate(const BasisSet &bs) {
  std::vector<std::tuple<const Shell &, int> > basis;
  const auto &bf = bs.shell2bf();
  for (size_t i = 0; i < bs.size(); ++i) {
    basis.emplace_back(bs[i], bf[i]);
  }
  return basis;
}

template <class S>
S slice(const std::tuple<const Shell &, int> &s) {
  int nbf = std::get<0>(s).size();
  int idx = std::get<1>(s);
  auto step = std::integral_constant<int, 1>{};
  return S{idx, idx + nbf, step};
}

}  // namespace libint2::python::basis

#endif /* LIBINT2_PYTHON_BASIS_H */
