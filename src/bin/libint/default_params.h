
#include <string>

#ifndef _libint2_src_bin_libint_defaultparams_h_
#define _libint2_src_bin_libint_defaultparams_h_

/**
  Defaults definitions for various parameters assumed by Libint
*/

namespace libint2 {

  struct StaticDefinitions {
    /// Do not vectorize by default
    static const int default_vector_length = 1;
    /// Produce quartet-level code by default
    static const int unroll_threshold = 1;
    /// Where to put generated library source
    static const std::string source_directory;

    /// De facto am limit
    static const unsigned int num_am_letters = 22;
    /// am -> char conversion
    static const char am_letters[num_am_letters];
  };

  /// Converts a computation label to the name of the function
  std::string label_to_funcname(const std::string& label);

};

#endif

