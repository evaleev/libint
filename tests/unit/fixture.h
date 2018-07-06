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

static auto atoms = std::vector<Atom>{ {8, 0.,0.,0.}, {8, 0.,0.,2.}, {1, 0.,-1.,-1.}, {1, 0.,1.,3.}};
static auto obs = BasisSet("6-31g*", atoms);
static auto dfbs = BasisSet("aug-cc-pvdz", atoms);

}
}

#endif //LIBINT_FIXTURE_H
