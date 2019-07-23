//
// Created by Eduard Valeyev on 7/3/18.
//

#ifndef LIBINT_FIXTURE_H
#define LIBINT_FIXTURE_H

#include <libint2.hpp>

using std::cout;
using std::cerr;
using std::endl;

using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;
using libint2::BraKet;
using libint2::CartesianShellNormalization;

namespace libint2 {
namespace unit {

class DefaultFixture {
public:
  DefaultFixture() : atoms{ {8, 0.,0.,0.}, {8, 0.,0.,2.}, {1, 0.,-1.,-1.}, {1, 0.,1.,3.}},
    obs("6-31g*", atoms),
    dfbs("aug-cc-pvdz", atoms) {}
protected:
  std::vector<Atom> atoms;
  BasisSet obs, dfbs;
};

}
}

#endif //LIBINT_FIXTURE_H
