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
/// Active-timer kind. Only the timer with this kind fires; all others are
/// no-ops. Used to measure each component (Boys / RR / tform) in its OWN
/// dedicated run, without nesting overhead. Set to TimerKind::None to
/// disable all inner timers (used for measuring whole-compute() wall time
/// without nested chrono bias).
enum class TimerKind : int { None = 0, Boys = 1, RR = 2, Tform = 3 };

inline TimerKind& active_timer_kind() {
  thread_local TimerKind k = TimerKind::None;
  return k;
}

/// Per-call chrono::high_resolution_clock::now() cost. Measured on the
/// host machine via src/bin/profile/chrono.cc — adjust if you re-measure.
constexpr int kChronoNowOverheadNs = 16;

inline ::libint2::Timers<1>& boys_fm_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(kChronoNowOverheadNs);
    return t;
  }();
  return tm;
}

inline ::libint2::Timers<1>& rr_kernel_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(kChronoNowOverheadNs);
    return t;
  }();
  return tm;
}

inline ::libint2::Timers<1>& tform_timer() {
  thread_local ::libint2::Timers<1> tm = []{
    ::libint2::Timers<1> t;
    t.set_now_overhead(kChronoNowOverheadNs);
    return t;
  }();
  return tm;
}
}  // namespace detail
#endif

}  // namespace libint2

// ---------------------------------------------------------------------------
// Timing-accumulator macros. Each LIBINT2_TIME_*_CALL(EXPR) wraps a single
// statement so that the time spent in EXPR is added to the corresponding
// thread-local timer. Pair-style _BEGIN/_END is provided for multi-statement
// blocks. When LIBINT2_BOYS_TIMING is undefined the wrappers degenerate to
// the bare expression / no-op, with zero runtime cost.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// One-active-timer policy: each macro fires only when the matching
// active_timer_kind is set, so component timings are taken in separate runs
// without nested chrono overhead. Each call pair (start+stop) carries a
// single `now()` overhead on this machine (~16 ns), already calibrated via
// `set_now_overhead` in the timer accessors above.
// ---------------------------------------------------------------------------
#ifdef LIBINT2_BOYS_TIMING
#  define LIBINT2_TIME_BOYS_CALL(EXPR)                                       \
     do {                                                                    \
       if (::libint2::detail::active_timer_kind() ==                         \
           ::libint2::detail::TimerKind::Boys) {                             \
         ::libint2::detail::boys_fm_timer().start(0);                        \
         EXPR;                                                               \
         ::libint2::detail::boys_fm_timer().stop(0);                         \
       } else { EXPR; }                                                      \
     } while (0)
#  define LIBINT2_TIME_RR_BEGIN                                              \
     do {                                                                    \
       if (::libint2::detail::active_timer_kind() ==                         \
           ::libint2::detail::TimerKind::RR)                                 \
         ::libint2::detail::rr_kernel_timer().start(0);                      \
     } while (0)
#  define LIBINT2_TIME_RR_END                                                \
     do {                                                                    \
       if (::libint2::detail::active_timer_kind() ==                         \
           ::libint2::detail::TimerKind::RR)                                 \
         ::libint2::detail::rr_kernel_timer().stop(0);                       \
     } while (0)
#  define LIBINT2_TIME_TFORM_BEGIN                                           \
     do {                                                                    \
       if (::libint2::detail::active_timer_kind() ==                         \
           ::libint2::detail::TimerKind::Tform)                              \
         ::libint2::detail::tform_timer().start(0);                          \
     } while (0)
#  define LIBINT2_TIME_TFORM_END                                             \
     do {                                                                    \
       if (::libint2::detail::active_timer_kind() ==                         \
           ::libint2::detail::TimerKind::Tform)                              \
         ::libint2::detail::tform_timer().stop(0);                           \
     } while (0)
// Deprecated wraps kept as no-ops to avoid touching every call site.
#  define LIBINT2_TIME_RR_CALL(EXPR)   do { EXPR; } while (0)
#  define LIBINT2_TIME_QGAU_CALL(EXPR) do { EXPR; } while (0)
#  define LIBINT2_TIME_PREP_CALL(EXPR) do { EXPR; } while (0)
#  define LIBINT2_TIME_PREP_BEGIN      ((void)0)
#  define LIBINT2_TIME_PREP_END        ((void)0)
#else
#  define LIBINT2_TIME_BOYS_CALL(EXPR) do { EXPR; } while (0)
#  define LIBINT2_TIME_RR_CALL(EXPR)   do { EXPR; } while (0)
#  define LIBINT2_TIME_QGAU_CALL(EXPR) do { EXPR; } while (0)
#  define LIBINT2_TIME_RR_BEGIN        ((void)0)
#  define LIBINT2_TIME_RR_END          ((void)0)
#  define LIBINT2_TIME_TFORM_BEGIN     ((void)0)
#  define LIBINT2_TIME_TFORM_END       ((void)0)
#  define LIBINT2_TIME_PREP_CALL(EXPR) do { EXPR; } while (0)
#  define LIBINT2_TIME_PREP_BEGIN      ((void)0)
#  define LIBINT2_TIME_PREP_END        ((void)0)
#endif

#endif  // C++11 or later

#endif  // defined(__cplusplus)

#endif  // header guard
