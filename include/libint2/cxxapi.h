/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_cxxapi_h_
#define _libint2_src_lib_libint_cxxapi_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "Libint2 C++ API requires C++11 support"
#endif

#include <atomic>
#include <cassert>

#include <libint2.h>  // NB this loads libint2/config.h

#ifdef LIBINT_USER_DEFINED_REAL
# error "C++11 API does not support with user-defined real types yet; omit --with-real-type when configuring"
#endif

#include <libint2/util/deprecated.h>
#include <libint2/util/singleton.h>

namespace libint2 {

  namespace detail {
    struct __initializer {
        __initializer() {
          libint2_static_init();
        }
        ~__initializer() {
          libint2_static_cleanup();
        }
    };

    static std::atomic<bool>& verbose_accessor() {
      static std::atomic<bool> value{true};
      return value;
    }
  } // namespace libint2::detail

  inline bool initialized() {
    using namespace detail;
    return managed_singleton<__initializer>::instance_exists();
  }
  /// initializes the libint library.
  /// @param[in] verbose boolean flag that controls the verbosity of messages produced by libint in std::cout . If false, no messages
  //             will be produced. The default is true.
  inline void initialize(bool verbose = true) {
    if (!initialized()) {
      using namespace detail;
      __initializer *x = managed_singleton<__initializer>::instance();
      (void) x;  // to suppress unused variable warning (not guaranteed to work) TODO revise when upgrade to C++17
      assert(x != nullptr);
      verbose_accessor() = verbose;
    }
  }
  /// finalizes the libint library.
  inline void finalize() {
    if (initialized()) {
      using namespace detail;
      managed_singleton<__initializer>::delete_instance();
      verbose_accessor() = true;
    }
  }
  /// Accessor for the verbose flag
  /// @return true if the library is permitted to generate diagnostic messages to std::cout
  inline bool verbose() {
    return detail::verbose_accessor();
  }
}

#include <libint2/chemistry/elements.h>
#include <libint2/atom.h>
#include <libint2/basis.h>
#include <libint2/solidharmonics.h>

#include <libint2/engine.h> // this is the end-user stuff, needs to check if library is initialized

#endif /* _libint2_src_lib_libint_cxxapi_h_ */
