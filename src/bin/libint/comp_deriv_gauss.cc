
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

