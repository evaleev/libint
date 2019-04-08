/*
 *  Copyright (C) 2004-2019 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <rr.h>
#include <integral.h>
#include <iter.h>
#include <algebra.h>

#ifndef _libint2_src_bin_libint_intsettoints_h_
#define _libint2_src_bin_libint_intsettoints_h_

using namespace std;


namespace libint2 {

  /** IntegralSet_to_Integrals_base is dummy class used for dynamic casts only
   */
  class IntegralSet_to_Integrals_base {
  protected:
    virtual ~IntegralSet_to_Integrals_base() {}
  };

  /** IntegralSet_to_Integrals converts I, a set of integrals, to individual integrals. Although this is
  technically not a recurrence relation, it can be treated as one.
  */
  template <class I>
  class IntegralSet_to_Integrals : public RecurrenceRelation,
    public IntegralSet_to_Integrals_base {
  public:
    typedef I TargetType;
    typedef typename I::iter_type ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    IntegralSet_to_Integrals(const SafePtr<I>&);
    virtual ~IntegralSet_to_Integrals() {}

    /// Return an instance if applicable, or a null pointer otherwise
    static SafePtr<IntegralSet_to_Integrals<I>> Instance(const SafePtr<TargetType>& Tint, unsigned int dir) {
      assert(dir == 0);
      // attempt to construct
      SafePtr<IntegralSet_to_Integrals<I>> this_ptr(new IntegralSet_to_Integrals<I>(Tint));
      // if succeeded (nchildren > 0) do post-construction
      assert(this_ptr->num_children() != 0);
      return this_ptr;
    }
    static bool directional() { return false; }


    /// Implementation of RecurrenceRelation::num_children()
    unsigned int num_children() const { return children_.size(); };
    /// target() returns pointer to target
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation's target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation's child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return true;
    }
    /// Reimplementation of RecurrenceRelation::invariant_type()
    bool invariant_type() const {
      // Converts from one BFSet to another!
      return false;
    }

  private:
    SafePtr<TargetType> target_;
    vector< SafePtr<ChildType> > children_;

    /// Implementation of RecurrenceRelation::generate_label()
    std::string generate_label() const {
      return "IntegralSet_to_Integrals";
      //throw std::runtime_error("IntegralSet_to_Integrals::label() -- code for this RR is never generated, so this function should never be used");
    }
    /// Reimplementation of RecurrenceRelation::spfunction_call()
    std::string spfunction_call(const SafePtr<CodeContext>& context,
                                const SafePtr<ImplicitDimensions>& dims) const
    {
      throw logic_error("IntegralSet_to_Integrals::spfunction_call -- should not call this function");
    }

  };


  template <class I>
    IntegralSet_to_Integrals<I>::IntegralSet_to_Integrals(const SafePtr<I>& Tint) :
    target_(Tint)
    {
      target_ = Tint;

      // Construct a subiterator for I
      SubIteratorBase<I> siter(Tint);

      // Set children pointers
      for(siter.init(); siter; ++siter)
        children_.push_back(siter.elem());
    };

  template <class I>
    SafePtr<typename I::iter_type>
    IntegralSet_to_Integrals<I>::child(unsigned int i) const
    {
      return children_.at(i);
    };


};

#endif

