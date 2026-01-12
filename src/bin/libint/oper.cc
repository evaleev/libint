/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <exception.h>
#include <oper.h>
#include <rr.h>

#include <cmath>
#include <stdexcept>

using namespace std;
using namespace libint2;

////////////

std::string R12_k_G12_Descr::label_(int K, bool contracted) {
  ostringstream oss;
  oss << "r12^" << K << " * " << (contracted ? "G12" : "g12");
  return oss.str();
}

std::string R12_k_G12_Descr::symbol_(int K, bool contracted) {
  ostringstream oss;
  oss << "r12_" << (K < 0 ? "minus_" : "") << std::abs(K) << "_"
      << (contracted ? "G12" : "g12");
  return oss.str();
}

////////////

unsigned int R12k_R12l_G12_Descr::key() const {
  unsigned int k =
      (((((K_[0] * kmax + K_[1]) * kmax + K_[2]) * kmax + L_[0]) * kmax +
        L_[1]) *
           kmax +
       L_[2]);
  return k;
}

std::string R12k_R12l_G12_Descr::label_(const IntVec3& K, const IntVec3& L,
                                        bool contracted) {
  ostringstream oss;
  oss << "(r12x^" << K[0] << "*r12y^" << K[1] << "*r12z^" << K[2] << ")*"
      << "(r12x^" << L[0] << "*r12y^" << L[1] << "*r12z^" << L[2] << ") * "
      << (contracted ? "G12" : "g12");
  return oss.str();
}

std::string R12k_R12l_G12_Descr::symbol_(const IntVec3& K, const IntVec3& L,
                                         bool contracted) {
  ostringstream oss;
  oss << "r12x" << K[0] << "_r12y" << K[1] << "_r12z" << K[2] << "__"
      << "r12x" << L[0] << "_r12y" << L[1] << "_r12z" << L[2] << "__"
      << (contracted ? "G12" : "g12");
  return oss.str();
}

////////////

std::string Ti_G12_Descr::label_(int K, bool contracted) {
  ostringstream oss;
  oss << "[T_" << K + 1 << "," << (contracted ? "G12" : "g12") << "]";
  return oss.str();
}

std::string Ti_G12_Descr::symbol_(int K, bool contracted) {
  ostringstream oss;
  oss << "T" << K + 1 << "_" << (contracted ? "G12" : "g12");
  return oss.str();
}

////////////

std::string G12_Ti_G12_Descr::label_(int K, bool contracted) {
  ostringstream oss;
  oss << "[" << (contracted ? "G12" : "g12") << ",[T_" << K + 1 << ","
      << (contracted ? "G12" : "g12") << "]]";
  return oss.str();
}

std::string G12_Ti_G12_Descr::symbol_(int K, bool contracted) {
  ostringstream oss;
  oss << (contracted ? "G12" : "g12") << "_"
      << "T" << K + 1 << "_" << (contracted ? "G12" : "g12");
  return oss.str();
}

////////////

std::string DivG12prime_xTx_Descr::label_(int I) {
  ostringstream oss;
  oss << "(\\nabla_" << I + 1 << " \\cdot g_{12}') (g_{12}' \\cdot \\nabla_"
      << I + 1 << ")";
  return oss.str();
}

std::string DivG12prime_xTx_Descr::symbol_(int I) {
  ostringstream oss;
  oss << "Div" << I + 1 << "G12prime_xTx";
  return oss.str();
}
