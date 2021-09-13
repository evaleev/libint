#include "catch.hpp"
#include "fixture.h"

TEST_CASE("Basis", "[basis]") {
  std::stringstream sstr;
  sstr << "2\n\nO 0 0 0\nO 0 0 1.5";
  auto atoms = libint2::read_dotxyz(sstr);
  REQUIRE_NOTHROW(libint2::BasisSet(
      "sto-3g", atoms)); // sto-3g.g94 given in old EMSL BSE format
  REQUIRE_NOTHROW(libint2::BasisSet(
      "sto-6g", atoms)); // sto-6g.g94 given in new MolSSI BSE format
}