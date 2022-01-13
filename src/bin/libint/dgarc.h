/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
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

//#include <rr.h>
#include <iostream>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_dgarc_h_
#define _libint2_src_bin_libint_dgarc_h_

namespace libint2 {

  class RecurrenceRelation;
  class DGVertex;
  /** Class DGArc describes arcs in a directed graph.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArc {

    SafePtr<DGVertex> orig_;  // Where this Arc leavs
    SafePtr<DGVertex> dest_;  // Where this Arc leads to

  public:
    DGArc(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest);
    virtual ~DGArc() {}

    SafePtr<DGVertex> orig() const { return orig_; }
    SafePtr<DGVertex> dest() const { return dest_; }

    /// Print out the arc
    virtual void print(std::ostream& os) const =0;

  };

  /** Class DGArcDirect describes arcs that does not correspond to any relationship.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArcDirect : public DGArc {

  public:
    DGArcDirect(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest) : DGArc(orig,dest) {}
    virtual ~DGArcDirect() {}

    /// Overload of DGArc::print()
    void print(std::ostream& os) const override
      {
        os << "DGArcDirect: connects " << orig().get() << " to " << dest().get();
      }
  };

  /** Class DGArcRR describes arcs correspond to recurrence relations.
      Each arc connects vertex orig_ to vertex dest_. */
  class DGArcRR : public DGArc {

  public:
    virtual ~DGArcRR() {}

    /// rr() returns pointer to the RecurrenceRelation describing the arc
    virtual SafePtr<RecurrenceRelation> rr() const =0;

  protected:
    DGArcRR(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest);

  };

  /** Class DGArcRel describes arcs in a directed graph which is
      represented by a relationship ArcRel. */
  // NOTE TO SELF (11/24/2004): need to implement checks on ArcRel
  // It obviously must implement some functions
  template <class ArcRel> class DGArcRel : public DGArcRR {

    SafePtr<ArcRel> rel_;     // Relationship described by the arc

  public:
    DGArcRel(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest,
	     const SafePtr<ArcRel>& rel);
    virtual ~DGArcRel();

    /// Implementation of DGArcRR::rr()
    SafePtr<RecurrenceRelation> rr() const override { return dynamic_pointer_cast<RecurrenceRelation,ArcRel>(rel_); }
    /// Overload of DGArc::print()
    void print(std::ostream& os) const override
      {
        os << "DGArcRel<T>: connects " << orig().get() << " to " << dest().get() << std::endl;
      }

  };

  template <class ArcRel>
    DGArcRel<ArcRel>::DGArcRel(const SafePtr<DGVertex>& orig, const SafePtr<DGVertex>& dest,
			       const SafePtr<ArcRel>& rel) :
    DGArcRR(orig,dest), rel_(rel)
    {
    };

  template <class ArcRel>
    DGArcRel<ArcRel>::~DGArcRel()
    {
    };

};

#endif

