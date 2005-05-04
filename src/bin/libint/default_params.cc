
#include <default_params.h>

const std::string libint2::StaticDefinitions::source_directory("./");
const char libint2::StaticDefinitions::am_letters[StaticDefinitions::num_am_letters] = "spdfghiklmnoqrtuvwxyz";

std::string
libint2::label_to_funcname(const std::string& label)
{
  std::string result("compute");
  result += label;
  return result;
}

