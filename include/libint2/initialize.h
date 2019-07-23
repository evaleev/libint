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

#ifndef _libint2_src_lib_libint_initialize_h_
#define _libint2_src_lib_libint_initialize_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "Libint2 C++ API requires C++11 support"
#endif

#include <atomic>
#include <cassert>
#include <iostream>

#include <libint2.h>
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

    inline std::atomic<bool>& verbose_accessor() {
      static std::atomic<bool> value{false};
      return value;
    }
    inline std::ostream*& verbose_stream_accessor() {
      static std::ostream* value = &std::clog;
      return value;
    }
  } // namespace libint2::detail

  /// checks if the libint has been initialized.
  /// @return true, if libint2::initialize() has been called since the last (if any) call to libint2::finalize()
  inline bool initialized() {
    using namespace detail;
    return managed_singleton<__initializer>::instance_exists();
  }
  /// initializes the libint library if not currently initialized, no-op otherwise
  /// @param[in] verbose boolean flag that controls the verbosity of messages produced by libint in std::clog . If false, no messages
  ///            will be produced. The default is false.
  inline void initialize(bool verbose = false) {
    if (!initialized()) {
      using namespace detail;
      __initializer *x = managed_singleton<__initializer>::instance();
      (void) x;  // to suppress unused variable warning (not guaranteed to work) TODO revise when upgrade to C++17
      assert(x != nullptr);
      verbose_accessor() = verbose;
    }
  }

  /// initializes the libint library if not currently initialized, no-op otherwise
  /// @param[in] os the output stream to which verbose diagnostics will be written (if @c initialize(true) is used, will write to @c std::clog )
  inline void initialize(std::ostream& os) {
    if (!initialized()) {
      initialize(true);
      using namespace detail;
      verbose_stream_accessor() = &os;
    }
  }
  /// finalizes the libint library.
  inline void finalize() {
    if (initialized()) {
      using namespace detail;
      managed_singleton<__initializer>::delete_instance();
      verbose_accessor() = true;
      verbose_stream_accessor() = &std::clog;
    }
  }
  /// Accessor for the disgnostics stream
  /// @return the stream to which the diagnostics will be written if verbose() returns true
  inline std::ostream& verbose_stream() {
    return *detail::verbose_stream_accessor();
  }
  /// Accessor for the verbose flag
  /// @return true if the library is permitted to generate diagnostic messages to the stream returned by verbose_stream(); always returns false
  ///         if @c initialized()==false
  inline bool verbose() {
    if (initialized()) {
      return detail::verbose_accessor();
    } else {
      return false;
    }
  }
}

#endif /* _libint2_src_lib_libint_initialize_h_ */

