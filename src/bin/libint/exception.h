
#ifndef _libint2_src_bin_libint_exception_h_
#define _libint2_src_bin_libint_exception_h_

#include <stdexcept>

namespace libint2 {

  class InvalidDecrement : public std::logic_error {
    
  public:
    InvalidDecrement(const std::string& a) :
      logic_error(a) {};
    
  };

  /// This exception used to indicate that some property is not set
  template <class T>
    class NotSet : public std::logic_error {

    public:
    NotSet(const std::string& a) :
      logic_error(a) {};
  };


};

#endif
