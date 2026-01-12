/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_tforms_h_
#define _libint2_src_bin_libint_tforms_h_

#include <memory.h>

#include <string>

namespace libint2 {

/**
  Externally accessible registry of information about a graph. Used to control
  how many algorithms behave (e.g. which transformation rules are allowed, etc.)
*/
class GraphRegistry {
 public:
  GraphRegistry();
  ~GraphRegistry();

  GraphRegistry* clone() const;

  /// Accumulate (true) or assign (false) target vertices? The default is to
  /// assign.
  bool accumulate_targets() const { return accumulate_targets_; }
  void accumulate_targets(bool at) { accumulate_targets_ = at; }
  /// Return pointers to targets via Libint_t::targets? Default is true.
  bool return_targets() const { return return_targets_; }
  void return_targets(bool rt) { return_targets_ = rt; }
  /// Will unroll the integral sets with size <= unroll_threshold. Default is 0
  /// (no unrolling).
  unsigned int unroll_threshold() const { return unroll_threshold_; }
  void unroll_threshold(unsigned int ut) { unroll_threshold_ = ut; }
  /// Minimum size when can unroll
  /// Will uncontract the integral sets if true. Default is no uncontracting
  bool uncontract() const { return uncontract_; }
  void uncontract(bool uc) { uncontract_ = uc; }
  /// Ignore missing prerequisites -- generate the code as is, without
  /// generating code for the prerequisites. This is a hack.
  bool ignore_missing_prereqs() const { return ignore_missing_prereqs_; }
  void ignore_missing_prereqs(bool imp) { ignore_missing_prereqs_ = imp; }
  /// Do Common Subexpression Elimination (CSE)? The default is false.
  bool do_cse() const { return do_cse_; }
  void do_cse(bool dc) { do_cse_ = dc; }
  /// Names the primary scratch space for the intermediates. The default is
  /// "libint->stack".
  const std::string& stack_name() const { return stack_name_; }
  void stack_name(const std::string& stack_name) { stack_name_ = stack_name; }
  /// Condense expressions? The default is false.
  bool condense_expr() const { return condense_expr_; }
  void condense_expr(bool ce) { condense_expr_ = ce; }
  /// if -1, no profiling, otherwise, indicates the current timer
  int current_timer() const { return current_timer_; }
  void current_timer(int ct) { current_timer_ = ct; }

 private:
  bool accumulate_targets_;
  bool return_targets_;
  unsigned int unroll_threshold_;
  bool uncontract_;
  bool ignore_missing_prereqs_;
  bool do_cse_;
  bool condense_expr_;
  std::string stack_name_;
  int current_timer_;
};

/**
  Internal registry of information. Encapsulates info which should not be
  controled by the user,
*/
class InternalGraphRegistry {
 public:
  InternalGraphRegistry();
  ~InternalGraphRegistry();

  /// Are targets computed, then accumulated, or accumulated directly? The
  /// latter possible only if all targets are unrolled.
  bool accumulate_targets_directly() const {
    return accumulate_targets_directly_;
  }
  void accumulate_targets_directly(bool atd) {
    accumulate_targets_directly_ = atd;
  }
  /// The size of the buffer at the beginning of stack allocated to hold
  /// accumulated targets
  MemoryManager::Size size_of_target_accum() const {
    return size_of_target_accum_;
  }
  void size_of_target_accum(const MemoryManager::Size& sota) {
    size_of_target_accum_ = sota;
  }

 private:
  bool accumulate_targets_directly_;
  MemoryManager::Size size_of_target_accum_;
};

}  // namespace libint2

#endif
