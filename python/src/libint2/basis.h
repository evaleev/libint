#ifndef LIBINT2_PYTHON_BASIS_H
#define LIBINT2_PYTHON_BASIS_H

#include <libint2/basis.h>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace libint2::python::basis {

template<class ... Args>
std::unique_ptr<BasisSet> make_basis(const Args& ... args) {
  //py::print("make_basis:", args...);
  auto basis = std::make_unique<BasisSet>(args...);
  assert(!basis->empty());
  assert(basis->nbf());
  return basis;
}

template<>
inline std::unique_ptr<BasisSet> make_basis(
  const std::map<int, std::vector<Shell> > &element_basis,
  const std::vector<Atom> &atoms)
{
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

inline size_t __len__(const BasisSet &b) {
  return b.size();
}

inline const Shell& __getitem__(const BasisSet &b, size_t idx) {
  if (idx >= b.size()) throw py::index_error();
  return b.at(idx);
}

inline auto __iter__(const BasisSet &b) {
  return py::make_iterator(b.begin(), b.end());
}

inline auto enumerate(const BasisSet& bs) {
  std::vector< std::tuple<const Shell&,int> > basis;
  const auto &bf = bs.shell2bf();
  for (size_t i = 0; i < bs.size(); ++i) {
    basis.emplace_back(bs[i], bf[i]);
  }
  return basis;
}

template<class S>
S slice(const std::tuple<const Shell&,int> &s) {
  int nbf = std::get<0>(s).size();
  int idx = std::get<1>(s);
  auto step = std::integral_constant<int,1>{};
  return S{idx, idx+nbf, step};
}

} // libint2::python::basis


#endif /* LIBINT2_PYTHON_BASIS_H */
