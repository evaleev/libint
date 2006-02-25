
#ifndef _libint2_src_bin_libint_tforms_h_
#define _libint2_src_bin_libint_tforms_h_

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
    
    private:
    bool can_unroll_;
    bool do_cse_;
  };
  
}

#endif

