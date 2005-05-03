
#include <string>

#ifndef _libint2_src_bin_libint_defaultparams_h_
#define _libint2_src_bin_libint_defaultparams_h_

/**
  Defaults definitions for various parameters assumed by Libint
*/

namespace libint2 {

  struct StaticDefinitions {
    static const unsigned int num_am_letters = 22;
    static const char am_letters[num_am_letters];
  };

  const int Libint2_DefaultVectorLength = 64;
  
  /// Converts a computation label to the name of the function
  std::string label_to_funcname(const std::string& label);

};

#endif

