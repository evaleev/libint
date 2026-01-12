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

#include <libint2/chemistry/sto3g_atomic_density.h>

#include <Eigen/Dense>
#include <libint2.hpp>
#include <tuple>
#include <vector>

#if defined(_MSC_VER)
#include <BaseTsd.h>
// handles ssize_t in pybind11/numpy.h
typedef SSIZE_T ssize_t;
#endif

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "basis.h"

namespace py = pybind11;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

namespace libint2::python::engine {

size_t num_threads();
void set_num_threads(size_t num_threads);

template <class... Params>
inline Engine make_engine(Operator op, const BraKet *braket = nullptr,
                          int L = LIBINT2_MAX_AM, int K = 10) {
  Engine engine(op, K, L, 0, 1e-15);
  if (braket) engine.set(*braket);
  // engine.set_params(params)...;
  return engine;
}

template <Operator Operator>
Engine make_engine() {
  return make_engine(Operator);
}

template <Operator Operator, class Params>
Engine make_engine(const Params *params) {
  auto engine = make_engine<Operator>();
  if (params) engine.set_params(*params);
  return engine;
}

template <class Shells>
py::object compute2(Engine &engine, const Shells &A, const Shells &B);

template <class Shells>
py::object compute3(Engine &engine, const Shells &A, const Shells &B,
                    const Shells &C);

template <class Shells>
py::object compute4(Engine &engine, const Shells &A, const Shells &B,
                    const Shells &C, const Shells &D);

Matrix compute_1body_ints(Engine &engine, const BasisSet &basis);

Matrix compute_2body_fock(Engine &engine, const Eigen::Ref<const Matrix> &D,
                          const BasisSet &basis);

}  // namespace libint2::python::engine

namespace libint2 {
namespace python {

using Double3 = std::array<double, 3>;
using Primitive = std::tuple<double, double>;
using PointCharge = std::pair<double, Double3>;

template <class... Args>
Atom make_atom(Args... args) = delete;

template <>
inline Atom make_atom(int Z, Double3 R) {
  return Atom{Z, R[0], R[1], R[2]};
}

template <>
inline Atom make_atom(py::object atom) {
  auto [Z, R] = py::cast<std::tuple<int, Double3> >(atom);
  return make_atom(Z, R);
}

Shell make_shell(int L, std::vector<Primitive> &primitives,
                 const Double3 &center, bool pure) {
  libint2::svector<Shell::real_t> exponents, coeff;
  for (auto [a, C] : primitives) {
    exponents.push_back(a);
    coeff.push_back(C);
  }
  assert(!exponents.empty() && exponents.size() == coeff.size());
  return Shell(exponents, {{L, pure, coeff}}, center);
}

auto alpha(const Shell &s) {
  auto a = s.alpha;
  return std::vector<double>(a.begin(), a.end());
}

auto coeffs(const Shell &s) {
  auto c = s.contr[0].coeff;
  return std::vector<double>(c.begin(), c.end());
}

std::vector<double> coeffs_normalized(const Shell &s) {
  std::vector<double> c;
  for (size_t i = 0; i < s.nprim(); ++i) {
    c.push_back(s.coeff_normalized(0, i));
  }
  return c;
}

PYBIND11_MODULE(libint2, m) {
  py::enum_<SHGShellOrdering>(m, "SHGShellOrdering")
      .value("Standard", libint2::SHGShellOrdering_Standard)
      .value("Gaussian", libint2::SHGShellOrdering_Gaussian)
      .value("MOLDEN", libint2::SHGShellOrdering_Gaussian);

  libint2::initialize();

  m.attr("MAX_AM") = LIBINT2_MAX_AM;

  py::class_<Atom>(m, "Atom")
      .def(py::init(&make_atom<int, Double3>))
      .def(py::init(&make_atom<py::object>));
  py::implicitly_convertible<py::object, Atom>();

  py::class_<Shell>(m, "Shell")
      .def(py::init(&make_shell), py::arg("L"), py::arg("primitives"),
           py::arg("center") = Double3{0, 0, 0}, py::arg("pure") = true)
      .def(py::init([](const std::tuple<int, std::vector<Primitive> > &s) {
        auto [L, prims] = s;
        return make_shell(L, prims, Double3{0, 0, 0}, true);
      }))
      .def("size", &Shell::size)
      .def_property_static(
          "unit_normalization",
          [](py::object) { return Shell::do_enforce_unit_normalization(); },
          [](py::object, bool f) { Shell::do_enforce_unit_normalization(f); })
      .def_property(
          "pure",
          [](const Shell &s) {
            for (auto &c : s.contr) {
              if (!c.pure) return false;
            }
            return true;
          },
          [](Shell &s, bool pure) {
            for (auto &c : s.contr) {
              c.pure = pure;
            }
          })
      .def_property_readonly("coeffs", &coeffs)
      .def_property_readonly("alpha", &alpha)
      .def_property_readonly("coeffs_normalized", &coeffs_normalized)
      .def_property_readonly_static("unit",
                                    [](py::object) { return Shell::unit(); });
  py::implicitly_convertible<py::tuple, Shell>();

  using basis::make_basis;
  py::class_<BasisSet>(m, "BasisSet")
      .def(py::init(&make_basis<std::vector<Shell> >))
      .def(py::init(
          &make_basis<std::map<int, std::vector<Shell> >, std::vector<Atom> >))
      .def(py::init(&make_basis<std::string, std::vector<Atom>, bool>),
           py::arg("name"), py::arg("atoms"),
           py::arg("throw_if_no_match") = true)
      .def("__len__", &basis::__len__)
      // Essential: keep object alive while reference/iterator exists
      .def("__getitem__", basis::__getitem__,
           py::return_value_policy::reference, py::keep_alive<0, 1>())
      .def("__iter__", &basis::__iter__, py::keep_alive<0, 1>())
      .def_property("pure", nullptr, &BasisSet::set_pure)
      .def_property_readonly("nbf", [](const BasisSet &bs) { return bs.nbf(); })
      .def_property_readonly("functions",
                             [](const BasisSet &bs) { return bs.shell2bf(); });
  py::implicitly_convertible<py::list, BasisSet>();

  py::enum_<Operator>(m, "Operator")
      .value("overlap", Operator::overlap)
      .value("nuclear", Operator::nuclear)
      .value("coulomb", Operator::coulomb)
      .value("kinetic", Operator::kinetic)
      .def_static("rank", [](Operator op) { return libint2::rank(op); });

  py::enum_<BraKet>(m, "BraKet")
      .value("XX", BraKet::x_x)
      .value("XXXX", BraKet::xx_xx)
      .value("XXXS", BraKet::xx_xs)
      .value("XSXX", BraKet::xs_xx)
      .value("XSXS", BraKet::xs_xs);

  py::class_<Engine>(m, "Engine")
      .def(py::init())
      .def(py::init(&engine::make_engine<void>), py::arg("oper"),
           py::arg("braket") = nullptr, py::arg("L") = LIBINT2_MAX_AM,
           py::arg("K") = 10)
      .def_property("oper", &Engine::oper,
                    [](Engine &e, Operator op) { e.set(op); })
      .def_property("braket", &Engine::braket,
                    [](Engine &e, BraKet braket) { e.set(braket); })
      .def("set_params",
           &Engine::set_params<std::vector<std::pair<double, Double3> > >)
      .def("compute", &engine::compute2<Shell>)
      .def("compute", &engine::compute3<Shell>)
      .def("compute", &engine::compute4<Shell>)
      .def("compute", &engine::compute2<BasisSet>)
      .def("compute", &engine::compute3<BasisSet>)
      .def("compute", &engine::compute4<BasisSet>)
      .def("compute_1body_ints", &engine::compute_1body_ints, py::arg("basis"))
      .def("compute_2body_fock", &engine::compute_2body_fock, py::arg("basis"),
           py::arg("density"))
      .def_property_static(
          "num_threads", [](py::object) { return engine::num_threads(); },
          [](py::object, int num_threads) {
            engine::set_num_threads(num_threads);
          });

  m.def("kinetic", &engine::make_engine<Operator::kinetic>);
  m.def("overlap", &engine::make_engine<Operator::overlap>);
  m.def("coulomb", &engine::make_engine<Operator::coulomb>);
  m.def("nuclear",
        &engine::make_engine<Operator::nuclear, std::vector<PointCharge> >);

  m.def("sto3g_num_ao", &libint2::sto3g_num_ao);
  m.def("sto3g_ao_occupation_vector",
        &libint2::sto3g_ao_occupation_vector<double>);

  m.def("solid_harmonics_ordering", &libint2::solid_harmonics_ordering);
  m.def("set_solid_harmonics_ordering", &libint2::set_solid_harmonics_ordering);

  using SolidHarmonicsCoefficients =
      libint2::solidharmonics::SolidHarmonicsCoefficients<double>;

  py::class_<SolidHarmonicsCoefficients>(m, "SolidHarmonicsCoefficients")
      .def_static("coefficient", SolidHarmonicsCoefficients::coeff,
                  py::arg("l"), py::arg("m"), py::arg("lx"), py::arg("ly"),
                  py::arg("lz"));
}

}  // namespace python
}  // namespace libint2
