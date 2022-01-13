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

/**
   Classes here extract information from DirectedGraph
*/

#ifndef _libint2_src_bin_libint_extract_h_
#define _libint2_src_bin_libint_extract_h_

#include <string>
#include <list>
#include <smart_ptr.h>

namespace libint2 {

  class DGVertex;

  /// This class collects labels of all external non-compile-time constants
  class ExtractExternSymbols {
  public:
    typedef SafePtr<DGVertex> VertexPtr;
    typedef std::list<std::string> Symbols;

    ExtractExternSymbols() {}
    ~ExtractExternSymbols() {}

    /// try v
    void operator()(const VertexPtr& v);

    /// return list of sorted symbols
    const Symbols& symbols() const;

  private:
    mutable Symbols symbols_;
    // symbols are stored as a map
    typedef std::map<std::string,bool> LabelMap;
    LabelMap map_;
  };

  /// This class collects all unique RRs. It uses RRStack to get their InstanceID
  class ExtractRR {
  public:
    typedef SafePtr<DGVertex> VertexPtr;
    typedef RRStack::InstanceID RRid;
    typedef std::list<RRid> RRList;

    ExtractRR() {}
    ~ExtractRR() {}

    /// try v
    void operator()(const VertexPtr& v);

    /// return list of sorted RRs
    const RRList& rrlist() const;

  private:
    mutable RRList rrlist_;
    // RRid are stored in a map
    typedef std::map<RRid,bool> RRMap;
    RRMap map_;
  };

};

#endif // header guard
