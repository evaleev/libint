/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_timer_h_
#define _libint2_src_lib_libint_timer_h_

#ifdef __cplusplus

#include <libint2/util/cxxstd.h>

#if LIBINT2_CPLUSPLUS_STD >= 2011

#include <chrono>

namespace libint2 {

/// Timers aggregates \c N C++11 "timers"; used to high-resolution profile
/// stages of integral computation
/// @tparam N the number of timers
/// @note member functions are not reentrant, use one Timers object per thread
template <size_t N>
class Timers {
 public:
  typedef std::chrono::duration<double> dur_t;
  typedef std::chrono::high_resolution_clock clock_t;
  typedef std::chrono::time_point<clock_t> time_point_t;

  Timers() {
    clear();
    set_now_overhead(0);
  }

  /// returns the current time point
  static time_point_t now() { return clock_t::now(); }

  /// use this to report the overhead of now() call; if set, the reported
  /// timings will be adjusted for this overhead
  /// @note this is clearly compiler and system dependent, please measure
  /// carefully (turn off turboboost, etc.)
  ///       using src/bin/profile/chrono.cc
  void set_now_overhead(size_t ns) { overhead_ = std::chrono::nanoseconds(ns); }

  /// starts timer \c t
  void start(size_t t) { tstart_[t] = now(); }
  /// stops timer \c t
  /// @return the duration, corrected for overhead, elapsed since the last call
  /// to \c start(t)
  dur_t stop(size_t t) {
    const auto tstop = now();
    const dur_t result = (tstop - tstart_[t]) - overhead_;
    timers_[t] += result;
    return result;
  }
  /// reads value (in seconds) of timer \c t , converted to \c double
  double read(size_t t) const { return timers_[t].count(); }
  /// resets timers to zero
  void clear() {
    for (auto t = 0; t != ntimers; ++t) {
      timers_[t] = dur_t::zero();
      tstart_[t] = time_point_t();
    }
  }

 private:
  constexpr static auto ntimers = N;
  dur_t timers_[ntimers];
  time_point_t tstart_[ntimers];
  dur_t overhead_;  // the duration of now() call ... use this to automatically
                    // adjust reported timings is you need fine-grained timing
};

#ifdef LIBINT2_BOYS_TIMING
namespace detail {
/// Thread-local timer accumulating *bare* FmEval->eval(...) time across all
/// engine paths. Slot 0 is started/stopped around each Fm call inside
/// q_gau_gm_eval (boys.h) and inside the Operator::coulomb case in
/// engine.impl.h. Engine::boys_time_seconds() reads slot 0 of this timer;
/// Engine::reset_timers() clears it. set_now_overhead(25) calibrates out
/// the ~27 ns std::chrono::high_resolution_clock::now() cost on Apple M2.
inline ::libint2::Timers<1>& boys_fm_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(25);
    return t;
  }();
  return tm;
}

/// Thread-local timer accumulating time inside the libint-generated kernel
/// (the OS recurrence relations: HRR + VRR + contraction). Wrapped around
/// the kernel call sites in engine.impl.h: the 1-body buildfnptrs call (and
/// the L=0 manual-contraction fallback), and the 2-body compute_directly
/// path and buildfnptrs call. Engine::rr_time_seconds() reads it.
inline ::libint2::Timers<1>& rr_kernel_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(25);
    return t;
  }();
  return tm;
}

/// Thread-local timer accumulating *total* time inside q_gau_gm_eval (the
/// q_gau core evaluator that loops over SAP primitives, calling Fm and doing
/// per-primitive sqrt/pow/multiply/accumulate). Wrapped around the
/// `core_eval_ptr->eval(...)` call inside compute_primdata at engine.impl.h
/// for Operator::q_gau / op_q_gau_op. The benchmark derives the q_gau
/// "post-Fm" cost as this total minus boys_fm_timer.
inline ::libint2::Timers<1>& q_gau_total_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(25);
    return t;
  }();
  return tm;
}

/// Thread-local timer accumulating time inside the solid-harmonic transform
/// stage (compute1: lines ~538–565; compute2: lines ~2019–2200). Engine::
/// tform_time_seconds() reads it.
inline ::libint2::Timers<1>& tform_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(25);
    return t;
  }();
  return tm;
}
}  // namespace detail
#endif

}  // namespace libint2

#endif  // C++11 or later

#endif  // defined(__cplusplus)

#endif  // header guard
