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

#ifndef _libint2_src_bin_libint_drtree_h_
#define _libint2_src_bin_libint_drtree_h_

#include <smart_ptr.h>

namespace libint2 {

class DGVertex;

/// This is a directed rooted tree
class DRTree : public std::enable_shared_from_this<DRTree> {
 public:
  typedef DRTree this_type;

  /// If v is not on a DRTree, make a new one using v as root
  static std::shared_ptr<DRTree> CreateRootedAt(
      const std::shared_ptr<DGVertex>& v);
  ~DRTree();

  /// number of vertices
  unsigned int nvertices() const { return nvertices_; }
  /// the root of the tree
  const std::shared_ptr<DGVertex>& root() const;
  /// remove all references from vertices to the tree and vice versa
  void detach();
  /// will try to add v to this subtree. Should not be used by the user
  void add_vertex(const std::shared_ptr<DGVertex>& v);
  /// recurively detach v from this
  void detach_from(const std::shared_ptr<DGVertex>& v);

 private:
  /// Create a tree starting at root
  DRTree(const std::shared_ptr<DGVertex>& root);
  /// grows the tree from the root. Must be called after the constructor
  void grow();

  unsigned int nvertices_;
  std::shared_ptr<DGVertex> root_;
};

}  // namespace libint2

#endif  // ifndef
