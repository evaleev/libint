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

// #include <rr.h>
#include <smart_ptr.h>

#include <iostream>

#ifndef _libint2_src_bin_libint_dgarc_h_
#define _libint2_src_bin_libint_dgarc_h_

namespace libint2 {

class RecurrenceRelation;
class DGVertex;
/** Class DGArc describes arcs in a directed graph.
    Each arc connects vertex orig_ to vertex dest_. */
class DGArc {
  std::shared_ptr<DGVertex> orig_;  // Where this Arc leavs
  std::shared_ptr<DGVertex> dest_;  // Where this Arc leads to

 public:
  DGArc(const std::shared_ptr<DGVertex>& orig,
        const std::shared_ptr<DGVertex>& dest);
  virtual ~DGArc() {}

  std::shared_ptr<DGVertex> orig() const { return orig_; }
  std::shared_ptr<DGVertex> dest() const { return dest_; }

  /// Print out the arc
  virtual void print(std::ostream& os) const = 0;
};

/** Class DGArcDirect describes arcs that does not correspond to any
   relationship. Each arc connects vertex orig_ to vertex dest_. */
class DGArcDirect : public DGArc {
 public:
  DGArcDirect(const std::shared_ptr<DGVertex>& orig,
              const std::shared_ptr<DGVertex>& dest)
      : DGArc(orig, dest) {}
  virtual ~DGArcDirect() {}

  /// Overload of DGArc::print()
  void print(std::ostream& os) const override {
    os << "DGArcDirect: connects " << orig().get() << " to " << dest().get();
  }
};

/** Class DGArcRR describes arcs correspond to recurrence relations.
    Each arc connects vertex orig_ to vertex dest_. */
class DGArcRR : public DGArc {
 public:
  virtual ~DGArcRR() {}

  /// rr() returns pointer to the RecurrenceRelation describing the arc
  virtual std::shared_ptr<RecurrenceRelation> rr() const = 0;

 protected:
  DGArcRR(const std::shared_ptr<DGVertex>& orig,
          const std::shared_ptr<DGVertex>& dest);
};

/** Class DGArcRel describes arcs in a directed graph which is
    represented by a relationship ArcRel. */
// NOTE TO SELF (11/24/2004): need to implement checks on ArcRel
// It obviously must implement some functions
template <class ArcRel>
class DGArcRel : public DGArcRR {
  std::shared_ptr<ArcRel> rel_;  // Relationship described by the arc

 public:
  DGArcRel(const std::shared_ptr<DGVertex>& orig,
           const std::shared_ptr<DGVertex>& dest,
           const std::shared_ptr<ArcRel>& rel);
  virtual ~DGArcRel();

  /// Implementation of DGArcRR::rr()
  std::shared_ptr<RecurrenceRelation> rr() const override {
    return std::dynamic_pointer_cast<RecurrenceRelation, ArcRel>(rel_);
  }
  /// Overload of DGArc::print()
  void print(std::ostream& os) const override {
    os << "DGArcRel<T>: connects " << orig().get() << " to " << dest().get()
       << std::endl;
  }
};

template <class ArcRel>
DGArcRel<ArcRel>::DGArcRel(const std::shared_ptr<DGVertex>& orig,
                           const std::shared_ptr<DGVertex>& dest,
                           const std::shared_ptr<ArcRel>& rel)
    : DGArcRR(orig, dest), rel_(rel){};

template <class ArcRel>
DGArcRel<ArcRel>::~DGArcRel(){};

};  // namespace libint2

#endif
