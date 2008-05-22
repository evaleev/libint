#include <util.h>

std::string libint2::to_string(FunctionPosition pos) {
  switch (pos) {
    case InBra:
      return "InBra";
    case InKet:
      return "InKet";
  }
}
