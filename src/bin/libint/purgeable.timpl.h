
#ifndef _libint2_src_bin_libint_purgeabletimpl_h_
#define _libint2_src_bin_libint_purgeabletimpl_h_

#include <boost/type_traits.hpp>
#include <dgvertex.h>

namespace libint2 {

  template <typename T>
  bool DefaultPurgingPolicy<T>::purgeable() {

    bool result = false;

    if (boost::is_base_of<DGVertex,T>::value == true) { // can only purge DGVertex objects
      result = true;
    }

    return result;
  }

  template <typename T>
  bool DefaultPurgingPolicy<T>::purge(const T* ref) {

    bool result = false;

    try {
      const DGVertex* dgv_ptr = dynamic_cast<const DGVertex*>(ref);
      if (dgv_ptr->dg() == 0)
        result = true;
    }
    catch(...) {
    }

    return result;
  }

};

#endif
