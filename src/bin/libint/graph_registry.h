
#ifndef _libint2_src_bin_libint_tforms_h_
#define _libint2_src_bin_libint_tforms_h_

#include <string>

using namespace std;


namespace libint2 {
  
  /**
    Registry of information about a graph
    */
  class GraphRegistry {
    public:
    GraphRegistry();
    ~GraphRegistry();
    
    /// Can unroll the integral sets?
    bool can_unroll() const { return can_unroll_; }
    void can_unroll(bool cu) { can_unroll_ = cu; }
    /// Do Common Subexpression Elimination (CSE)?
    bool do_cse() const { return do_cse_; }
    void do_cse(bool dc) { do_cse_ = dc; }
    /// Names the primary scratch space for the intermediates
    const std::string& stack_name() const { return stack_name_; }
    void stack_name(const std::string& stack_name) { stack_name_ = stack_name; }
    
    private:
    bool can_unroll_;
    bool do_cse_;
    std::string stack_name_;
  };
  
}

#endif

