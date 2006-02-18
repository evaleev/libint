
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
    
    bool can_unroll() const { return can_unroll_; }
    void can_unroll(bool cu) { can_unroll_ = cu; }
    
    private:
    bool can_unroll_;
  };
  
}

#endif

