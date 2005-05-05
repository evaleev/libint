

#include <smart_ptr.h>
#include <rr.h>

#ifndef _libint2_src_bin_libint_tactic_h_
#define _libint2_src_bin_libint_tactic_h_

using namespace std;

namespace libint2 {

  class DirectedGraph;
  
  /** Tactic is used to choose the optimal (in some sense) recurrence relation to reduce
      a vertex.
  */
  class Tactic {
    public:
    typedef SafePtr<RecurrenceRelation> RR;
    typedef vector<RR> rr_stack;
    
    Tactic() {}
    virtual ~Tactic() {}
    
    virtual RR optimal_rr(const rr_stack& stack) const =0;
  };
  
  /** FirstChoiceTactic simply chooses the first RR
    */
  class FirstChoiceTactic : public Tactic {
    public:
    FirstChoiceTactic() : Tactic() {}
    ~FirstChoiceTactic() {}
    
    RR optimal_rr(const rr_stack& stack) const;
  };
  
  /** FewestNewVerticesTactic chooses RR which adds fewest new vertices to
      DirectedGraph dg
    */
  class FewestNewVerticesTactic : public Tactic {
    public:
    FewestNewVerticesTactic(const SafePtr<DirectedGraph>& dg) : Tactic(), dg_(dg) {}
    ~FewestNewVerticesTactic() {}
    
    RR optimal_rr(const rr_stack& stack) const;

    private:
    SafePtr<DirectedGraph> dg_;
  };
  
  /** ZeroNewVerticesTactic chooses first RR which adds no new vertices on
      DirectedGraph dg
    */
  class ZeroNewVerticesTactic : public Tactic {
    public:
    ZeroNewVerticesTactic(const SafePtr<DirectedGraph>& dg) : Tactic(), dg_(dg) {}
    ~ZeroNewVerticesTactic() {}
    
    RR optimal_rr(const rr_stack& stack) const;
    
    private:
    SafePtr<DirectedGraph> dg_;
  };

  /** RandomChoiceTactic chooses randomly among the applicable RRs
    */
  class RandomChoiceTactic : public Tactic {
    public:
    RandomChoiceTactic();
    ~RandomChoiceTactic() {}
    
    RR optimal_rr(const rr_stack& stack) const;
  };
  
  /** NullTactic always returns null RecurrenceRelation
    */
  class NullTactic : public Tactic {
    public:
    NullTactic() : Tactic() {}
    ~NullTactic() {}
    
    RR optimal_rr(const rr_stack& stack) const;
    
  };

};

#endif

