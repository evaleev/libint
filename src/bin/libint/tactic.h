
#include <cmath>
#include <cstdlib>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_tactic_h_
#define _libint2_src_bin_libint_tactic_h_

using namespace std;

namespace libint2 {

  class DirectedGraph;
  class RecurrenceRelation;

  class DummyRandomizePolicy;
  class StdRandomizePolicy;

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
  template <class RandomizePolicy = DummyRandomizePolicy>
    class FirstChoiceTactic : public Tactic {
    public:
    FirstChoiceTactic(const SafePtr<RandomizePolicy>& rpolicy = SafePtr<RandomizePolicy>(new RandomizePolicy)) : Tactic(), rpolicy_(rpolicy) {}
    virtual ~FirstChoiceTactic() {}

    RR optimal_rr(const rr_stack& stack) const {
      if (!stack.empty())
        return stack[0 + rpolicy_->noise(stack.size())];
      else
        return RR();
    }

    private:
    SafePtr<RandomizePolicy> rpolicy_;
  };

  /** FewestNewVerticesTactic chooses RR which adds fewest new vertices to
      DirectedGraph dg
    */
  class FewestNewVerticesTactic : public Tactic {
    public:
    FewestNewVerticesTactic(const SafePtr<DirectedGraph>& dg) : Tactic(), dg_(dg) {}
    virtual ~FewestNewVerticesTactic() {}

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
    virtual ~ZeroNewVerticesTactic() {}

    RR optimal_rr(const rr_stack& stack) const;

    private:
    SafePtr<DirectedGraph> dg_;
  };

  /** RandomChoiceTactic chooses randomly among the applicable RRs
    */
  class RandomChoiceTactic : public Tactic {
    public:
    RandomChoiceTactic();
    virtual ~RandomChoiceTactic() {}

    RR optimal_rr(const rr_stack& stack) const;
  };

  /** NullTactic always returns null RecurrenceRelation
    */
  class NullTactic : public Tactic {
    public:
    NullTactic() : Tactic() {}
    virtual ~NullTactic() {}

    RR optimal_rr(const rr_stack& stack) const;

  };


  /////////////////////////////////

  struct DummyRandomizePolicy {
    unsigned int noise(unsigned int nrrs) const { return 0; }
  };

  /** The shift parameter is computed as follows:
      delta = floor(nrrs*scale*random()/RAND_MAX)
      where nrrs is the number of possibilities, scale
      is the user-specified parameter.
  */
  class StdRandomizePolicy {
    public:
    StdRandomizePolicy(double scale) : scale_(scale) {
      // Initialize state randomly
      time_t crap;
      srandom(time(&crap));
    }

    unsigned int noise(unsigned int nrrs) const {
      unsigned long rand = random();
      const unsigned long range = RAND_MAX;
      const unsigned int result = static_cast<unsigned int>(std::floor(nrrs*scale_*rand/range));
      return result;
    }

    private:
    double scale_;
  };


};

#endif

