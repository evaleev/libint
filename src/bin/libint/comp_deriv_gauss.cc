/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <fstream>
#include <comp_deriv_gauss.h>

using namespace libint2;

CR_DerivGauss_GenericInstantiator CR_DerivGauss_GenericInstantiator::instance_;

CR_DerivGauss_GenericInstantiator::CR_DerivGauss_GenericInstantiator() {}

CR_DerivGauss_GenericInstantiator::~CR_DerivGauss_GenericInstantiator() {
  if (not template_instances_.empty()) {
    std::ofstream ofile("GenericGaussDeriv.cc");

    ofile << "#include \"libint2.h\"" << std::endl;
    ofile << "#include \"GenericGaussDeriv.impl.h\"" << std::endl << std::endl;
    for(auto v=template_instances_.begin();
        v != template_instances_.end();
        ++v) {
      ofile << "template struct libint2::GenericGaussDeriv<" << v->first << "," << (v->second ? "true" : "false")
            << ">;" << std::endl;
    }
  }
}

CR_DerivGauss_GenericInstantiator&
CR_DerivGauss_GenericInstantiator::instance() {
  return instance_;
}

void
CR_DerivGauss_GenericInstantiator::add(unsigned int L, bool vectorize) {
  template_instances_.insert(std::make_pair(L, vectorize));
}

