/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_timer_h_
#define _libint2_src_lib_libint_timer_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <chrono>

namespace libint2 {

  /// Timers aggregates \c N C++11 "timers"; used to high-resolution profile stages of integral computation
  /// @tparam N the number of timers
  /// @note member functions are not reentrant, use one Timers object per thread
  template <size_t N>
  class Timers {
    public:
      typedef std::chrono::duration<double> dur_t;
      typedef std::chrono::high_resolution_clock clock_t;
      typedef std::chrono::time_point<clock_t> time_point_t;

      static time_point_t now() {
        return clock_t::now();
      }

      /// use this to report the overhead of now() call; if set, the reported timings will be adjusted for this overhead
      /// @note this is clearly compiler and system dependent, please measure carefully
      void set_now_overhead(size_t ns) {
        overhead_ = std::chrono::nanoseconds(ns);
      }

      /// starts timer \c t
      void start(size_t t) {
        tstart_[t] = now();
      }
      /// stops timer \c t
      void stop(size_t t) {
        const auto tstop = now();
        timers_[t] += (tstop - tstart_[t]) - overhead_;
      }
      /// reads timer \c t
      double read(size_t t) const {
        return timers_[t].count();
      }
      void clear() {
        for(auto t=0; t!=ntimers; ++t) {
          timers_[t] = dur_t::zero();
          tstart_[t] = time_point_t();
        }
      }

    private:
      constexpr static size_t ntimers = N;
      dur_t timers_[ntimers];
      time_point_t tstart_[ntimers];
      dur_t overhead_; // the duration of now() call ... use this to automatically adjust reported timings is you need fine-grained timing
  };

} // namespace libint2

#endif // header guard


