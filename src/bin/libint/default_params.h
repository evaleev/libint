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

#ifndef _libint2_src_bin_libint_defaultparams_h_
#define _libint2_src_bin_libint_defaultparams_h_

#include <smart_ptr.h>

#include <climits>
#include <map>
#include <ostream>
#include <string>
#include <vector>

/**
  Defaults definitions for various parameters assumed by Libint
*/

namespace libint2 {

/// These are the parameters received by the compiler
class CompilationParameters {
  std::string default_task_name_;

 public:
  /// Use default parameters
  CompilationParameters();
  ~CompilationParameters();

  /// returns max AM for task \c t and center \c c
  unsigned int max_am(std::string t = "", unsigned int c = 0) const;
  /// returns max AM for which to produce optimal code for task t
  unsigned int max_am_opt(std::string t = "") const;
  /// returns number of basis functions in integrals for task t
  unsigned int num_bf(std::string t = "") const;
  /// returns max vector length
  unsigned int max_vector_length() const { return max_vector_length_; }
  /// returns whether to vectorize line-by-line
  bool vectorize_by_line() const { return vectorize_by_line_; }
  /// returns unroll threshold
  unsigned int unroll_threshold() const { return unroll_threshold_; }
  /// returns alignment size (in units of sizeof(LIBINT_FLOAT))
  unsigned int align_size() const { return align_size_; }
  /// returns directory path for the generated source
  const std::string& source_directory() const { return source_directory_; }
  /// the API prefix
  const std::string& api_prefix() const { return api_prefix_; }
  /// generate single evaluator type (i.e. not specific to each task)?
  bool single_evaltype() const { return single_evaltype_; }
  /// returns whether to use C-style linking
  bool use_C_linking() const { return use_C_linking_; }
  /// whether FLOP counting is enabled
  bool count_flops() const { return count_flops_; }
  /// whether profiling instrumentation is enabled
  bool profile() const { return profile_; }
  /// whether target integrals are accumulated
  bool accumulate_targets() const { return accumulate_targets_; }
  /// name of the floating-point type
  const std::string& realtype() const { return realtype_; }
  /// whether contracted targets are supported
  bool contracted_targets() const { return contracted_targets_; }
  /// default task name
  const std::string& default_task_name() const { return default_task_name_; }

  /// set max AM for task \c t and center \c c
  void max_am(const std::string& t, unsigned int a, unsigned int c = 0);
  /// set max AM for task t
  void max_am_opt(const std::string& t, unsigned int a);
  /// set num of basis functions for task t
  void num_bf(const std::string& t, unsigned int a);
  /// set max vector length
  void max_vector_length(unsigned int a) { max_vector_length_ = a; }
  /// set vectorize_by_line flag
  void vectorize_by_line(bool flag) { vectorize_by_line_ = flag; }
  /// set alignment size (in units of sizeof(LIBINT_FLOAT))
  void align_size(unsigned int a) { align_size_ = a; }
  /// set unroll threshold
  void unroll_threshold(unsigned int a) { unroll_threshold_ = a; }
  /// set generated source directory
  void source_directory(const std::string& a) { source_directory_ = a; }
  /// API prefix
  void api_prefix(const std::string& a) { api_prefix_ = a; }
  /// generate a single evaluator type (i.e. not specific to each task)?
  void single_evaltype(bool se) { single_evaltype_ = se; }
  /// set whether to use C style linking
  void use_C_linking(bool a) { use_C_linking_ = a; }
  /// set whether to count FLOPs
  void count_flops(bool a) { count_flops_ = a; }
  /// set to instrument profiling
  void profile(bool a) { profile_ = a; }
  /// accumulate targets?
  void accumulate_targets(bool a) { accumulate_targets_ = a; }
  /// which floating-point type to use
  void realtype(const std::string& realtype) { realtype_ = realtype; }
  /// support contracted targets?
  void contracted_targets(bool c) { contracted_targets_ = c; }
  /// default task name
  void default_task_name(const std::string& s) { default_task_name_ = s; }

  /// print params out
  void print(std::ostream& os) const;

 private:
  struct Defaults {
    /// By default compile general integrals for p-functions
    static const unsigned int max_am = 1;
    /// By default optimize general integrals for up to p-functions
    static const unsigned int max_am_opt = 1;
    /// By default compute integrals with 4 basis functions
    static const unsigned int num_bf = 4;
    /// Do not vectorize by default
    static const unsigned int max_vector_length = 1;
    /// Vectorize all body by default
    static const bool vectorize_by_line = false;
    /// Use default alignment by default
    static const unsigned int align_size = 0;
    /// Produce quartet-level code by default
    static const unsigned int unroll_threshold = 0;
    /// Where to put generated library source
    static const std::string source_directory;
    /// API prefix
    static const std::string api_prefix;
    /// generate single evaltype?
    static const bool single_evaltype = true;
    /// Use C-style linking convention by default
    static const bool use_C_linking = true;
    /// Do not count FLOPs by default
    static const bool count_flops = false;
    /// Do not profile by default
    static const bool profile = false;
    /// Do not accumulate targets by default
    static const bool accumulate_targets = false;
    /// Use double for computations
    static const std::string realtype;
    /// Do not support contracted targets
    static const bool contracted_targets = false;
    /// task name
    static const std::string task_name;
  };

  struct TaskParameters {
    /// max AM
    std::vector<unsigned int> max_am;
    /// max AM for "optimized" integrals
    unsigned int max_am_opt;
    /// number of basis functions
    unsigned int num_bf;

    TaskParameters() : max_am(1, 0u), max_am_opt(0u), num_bf(0u) {}
  };
  /// Parameters for tasks
  std::map<std::string, TaskParameters> task_params_;
  /// tests if task t is known globally (i.e. to LibraryTaskManager). Throws if
  /// not.
  void task_exists(const std::string& t) const;
  /// adds a task using default task params
  void add_task(const std::string& t);

  /// max vector length
  unsigned int max_vector_length_;
  /// whether to vectorize line-by-line
  bool vectorize_by_line_;
  /** alignment size in units of sizeof(LIBINT2_REALTYPE).
      UINT_MAX => standard compiler/library default for scalar code, veclen for
     vectorized code
    */
  unsigned int align_size_;
  /// unroll threshold
  unsigned int unroll_threshold_;
  /// source directory
  std::string source_directory_;
  /// API prefix
  std::string api_prefix_;
  /// generate single evaluator type?
  bool single_evaltype_;
  /// whether to use C linking
  bool use_C_linking_;
  /// whether to count FLOPs
  bool count_flops_;
  /// whether to turn on profiling instrumentation
  bool profile_;
  /// whether to accumulate targets
  bool accumulate_targets_;
  /// name of the floating-point type
  std::string realtype_;
  /// whether to support contracted targets
  bool contracted_targets_;
};

/** This class maintains various parameters for each task type
    which can only be determined during the source generation
    (max stack size, etc.).
  */
class TaskParameters {
 public:
  TaskParameters();
  ~TaskParameters() {}

  /// returns the max number of targets
  unsigned int max_ntarget() const { return max_ntarget_; }
  /// returns the max quantum number of targets
  unsigned int max_am() const { return max_stack_size_.size() - 1; }
  /// returns max stack size needed for quantum numbers up to am
  unsigned int max_stack_size(unsigned int am) const {
    return max_stack_size_.at(am);
  }
  /** returns max vector stack size.
      vector stack is only used to hold intermediate quantities
      in set-level RR code. This is only needed when doing linewise
     vectorization.
    */
  unsigned int max_vector_stack_size(unsigned int am) const {
    if (am + 1 > max_vector_stack_size_.size())
      return 0;
    else
      return max_vector_stack_size_[am];
  }
  /** returns max rank of high-significance functions in a HRR call.
      This is only needed when doing linewise vectorization.
    */
  unsigned int max_hrr_hsrank(unsigned int am) const {
    if (am + 1 > max_hrr_hsrank_.size())
      return 0;
    else
      return max_hrr_hsrank_[am];
  }
  /** returns max rank of low-significance functions in a HRR call.
      This is only needed when doing linewise vectorization.
    */
  unsigned int max_hrr_lsrank(unsigned int am) const {
    if (am + 1 > max_hrr_lsrank_.size())
      return 0;
    else
      return max_hrr_lsrank_[am];
  }

  /// if max_ntarget_ < ntarget then set max_ntarget_=ntarget
  void max_ntarget(unsigned int ntarget) {
    if (max_ntarget_ < ntarget) max_ntarget_ = ntarget;
  }

  /// if max_stack_size_ < size then set max_stack_size_=size
  void max_stack_size(unsigned int am, unsigned int size) {
    extend_max_size(max_stack_size_, am, size);
  }

  /// if max_vector_stack_size_ < size then set max_vector_stack_size_=size
  void max_vector_stack_size(unsigned int am, unsigned int size) {
    extend_max_size(max_vector_stack_size_, am, size);
  }

  /// if max_hrr_hsrank_ < rank then set max_hrr_hsrank_=rank
  void max_hrr_hsrank(unsigned int am, unsigned int rank) {
    extend_max_size(max_hrr_hsrank_, am, rank);
  }

  /// if max_hrr_lsrank_ < rank then set max_hrr_lsrank_=rank
  void max_hrr_lsrank(unsigned int am, unsigned int rank) {
    extend_max_size(max_hrr_lsrank_, am, rank);
  }

 private:
  unsigned int max_ntarget_;
  std::vector<unsigned int>
      max_stack_size_;  //< keeps track of stack size needed to compute
                        // integrals up to a given quantum number
  std::vector<unsigned int> max_vector_stack_size_;
  std::vector<unsigned int> max_hrr_hsrank_;
  std::vector<unsigned int> max_hrr_lsrank_;

  static void extend_max_size(std::vector<unsigned int>& max_size,
                              unsigned int am, unsigned int size) {
    const int max_am = (int)max_size.size() - 1u;
    if (max_am < (int)am) {
      max_size.resize(am + 1);
      for (int l = std::max(max_am + 1, 1); l <= (int)am; ++l)
        max_size[l] = max_size[l - 1];
    }
    if (max_size[am] < size) max_size[am] = size;
  }
};

/// Converts a label, e.g. name of the target node, to the name of the function
/// to compute it
std::string label_to_funcname(const std::string& label);

/// need to condense expressions? Makes sense if vectorizing the code or the
/// compiler somehow prefers long expressions It does not make sense if there
/// will be only set-level RR calls
bool condense_expr(unsigned int unroll_threshold, bool vectorize);
};  // namespace libint2

#endif
