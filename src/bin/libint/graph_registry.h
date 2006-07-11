
#ifndef _libint2_src_bin_libint_tforms_h_
#define _libint2_src_bin_libint_tforms_h_

#include <string>
#include <memory.h>

using namespace std;


namespace libint2 {
  
  /**
    Externally accessible registry of information about a graph. Used to control how many algorithms
    behave (e.g. which transformation rules are allowed, etc.)
  */
  class GraphRegistry {
    public:
    GraphRegistry();
    ~GraphRegistry();
    
    /// Accumulate (true) or assign (false) target vertices? The default is to assign.
    bool accumulate_targets() const { return accumulate_targets_; }
    void accumulate_targets(bool at) { accumulate_targets_ = at; }
    /// Return pointers to targets via Libint_t::targets? Default is true.
    bool return_targets() const { return return_targets_; }
    void return_targets(bool rt) { return_targets_ = rt; }
    /// Can unroll the integral sets? The default is yes.
    bool can_unroll() const { return can_unroll_; }
    void can_unroll(bool cu) { can_unroll_ = cu; }
    /// Do Common Subexpression Elimination (CSE)? The default is false.
    bool do_cse() const { return do_cse_; }
    void do_cse(bool dc) { do_cse_ = dc; }
    /// Names the primary scratch space for the intermediates. The default is "libint->stack".
    const std::string& stack_name() const { return stack_name_; }
    void stack_name(const std::string& stack_name) { stack_name_ = stack_name; }
    
    private:
    bool accumulate_targets_;
    bool return_targets_;
    bool can_unroll_;
    bool do_cse_;
    std::string stack_name_;
  };
  
  /**
    Internal registry of information. Encapsulates info which should not be controled by the user,
  */
  class InternalGraphRegistry {
    public:
    InternalGraphRegistry();
    ~InternalGraphRegistry();
    
    /// Are targets computed, then accumulated, or accumulated directly? The latter possible only if all targets are unrolled.
    bool accumulate_targets_directly() const { return accumulate_targets_directly_; }
    void accumulate_targets_directly(bool atd) { accumulate_targets_directly_ = atd; }
    /// The size of the buffer at the beginning of stack allocated to hold accumulated targets
    MemoryManager::Size size_of_target_accum() const { return size_of_target_accum_; }
    void size_of_target_accum(const MemoryManager::Size& sota) { size_of_target_accum_ = sota; }
    
    private:
    bool accumulate_targets_directly_;
    MemoryManager::Size size_of_target_accum_;
  };
  
}

#endif

