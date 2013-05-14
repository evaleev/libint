#include <util.h>
#include <cassert>

std::string libint2::to_string(FunctionPosition pos) {
  switch (pos) {
    case InBra:
      return "InBra";
    case InKet:
      return "InKet";
    default:
      assert(false);
      break;
  }
  return ""; // pacify picky compilers
}
