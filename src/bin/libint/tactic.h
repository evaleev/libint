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

#include <smart_ptr.h>

#include <cmath>
#include <cstdlib>
#include <ctime>

#ifndef _libint2_src_bin_libint_tactic_h_
#define _libint2_src_bin_libint_tactic_h_

namespace libint2 {

class DirectedGraph;
class RecurrenceRelation;

struct DummyRandomizePolicy;
class StdRandomizePolicy;

/** Tactic is used to choose the optimal (in some sense) recurrence relation to
   reduce a vertex.
*/
class Tactic {
 public:
  typedef std::shared_ptr<RecurrenceRelation> RR;
  typedef std::vector<RR> rr_stack;

  Tactic() {}
  virtual ~Tactic() {}

  virtual RR optimal_rr(const rr_stack& stack) const = 0;

  // make return true if need to debug Tactic classes
  static constexpr bool class_debug() { return false; }
};

/** FirstChoiceTactic simply chooses the first RR
 */
template <class RandomizePolicy = DummyRandomizePolicy>
class FirstChoiceTactic : public Tactic {
 public:
  FirstChoiceTactic(const std::shared_ptr<RandomizePolicy>& rpolicy =
                        std::shared_ptr<RandomizePolicy>(new RandomizePolicy))
      : Tactic(), rpolicy_(rpolicy) {}
  virtual ~FirstChoiceTactic() {}

  RR optimal_rr(const rr_stack& stack) const {
    if (!stack.empty())
      return stack[0 + rpolicy_->noise(stack.size())];
    else
      return RR();
  }

 private:
  std::shared_ptr<RandomizePolicy> rpolicy_;
};

/** FewestNewVerticesTactic chooses RR which adds fewest new vertices to
    DirectedGraph dg
  */
class FewestNewVerticesTactic : public Tactic {
 public:
  FewestNewVerticesTactic(const std::shared_ptr<DirectedGraph>& dg)
      : Tactic(), dg_(dg) {}
  virtual ~FewestNewVerticesTactic() {}

  RR optimal_rr(const rr_stack& stack) const;

 private:
  std::shared_ptr<DirectedGraph> dg_;
};

/** ZeroNewVerticesTactic chooses first RR which adds no new vertices on
    DirectedGraph dg
  */
class ZeroNewVerticesTactic : public Tactic {
 public:
  ZeroNewVerticesTactic(const std::shared_ptr<DirectedGraph>& dg)
      : Tactic(), dg_(dg) {}
  virtual ~ZeroNewVerticesTactic() {}

  RR optimal_rr(const rr_stack& stack) const;

 private:
  std::shared_ptr<DirectedGraph> dg_;
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

/**
 * ParticleDirectionTactic returns the first RR that transfers the quantum
 * numbers between particles in the desired direction.
 */
class ParticleDirectionTactic : public Tactic {
 public:
  /**
   * @param increasing if true, quanta should be transferred from lower to
   * higher particle indices.
   */
  ParticleDirectionTactic(bool increase) : Tactic(), increase_(increase) {}
  virtual ~ParticleDirectionTactic() {}

  RR optimal_rr(const rr_stack& stack) const;

 private:
  bool increase_;
};

/**
 * TwoCenter_OS_Tactic decides graph build for <bra0|ket0>
 */
class TwoCenter_OS_Tactic : public Tactic {
 public:
  /**
   * @param lbra0
   * @param lket0
   */
  TwoCenter_OS_Tactic(unsigned lbra0, unsigned lket0)
      : Tactic(), lbra0_(lbra0), lket0_(lket0) {}
  virtual ~TwoCenter_OS_Tactic() {}

  RR optimal_rr(const rr_stack& stack) const;

 private:
  unsigned lbra0_;
  unsigned lket0_;
};

/**
 * FourCenter_OS_Tactic decides graph build for (bra0 ket0| bra1 ket1) = <bra0
 * bra1|ket0 ket1>
 */
class FourCenter_OS_Tactic : public Tactic {
 public:
  /**
   * @param lbra0
   * @param lbra1
   * @param lket0
   * @param lket1
   */
  FourCenter_OS_Tactic(unsigned lbra0, unsigned lket0, unsigned lbra1,
                       unsigned lket1)
      : Tactic(), lbra0_(lbra0), lket0_(lket0), lbra1_(lbra1), lket1_(lket1) {}
  virtual ~FourCenter_OS_Tactic() {}

  RR optimal_rr(const rr_stack& stack) const;

 private:
  unsigned lbra0_;
  unsigned lket0_;
  unsigned lbra1_;
  unsigned lket1_;
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
    time_t t;
    srandom(time(&t));
  }

  unsigned int noise(unsigned int nrrs) const {
    unsigned long rand = random();
    const unsigned long range = RAND_MAX;
    const unsigned int result =
        static_cast<unsigned int>(std::floor(nrrs * scale_ * rand / range));
    return result;
  }

 private:
  double scale_;
};

};  // namespace libint2

#endif
