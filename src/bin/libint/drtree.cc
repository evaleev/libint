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

#include <dgvertex.h>
#include <drtree.h>

#define LOCAL_DEBUG 0

using namespace libint2;

std::shared_ptr<DRTree> DRTree::CreateRootedAt(
    const std::shared_ptr<DGVertex>& v) {
  std::shared_ptr<DRTree> stree = v->subtree();
  if (!stree) {
    std::shared_ptr<DRTree> result(new DRTree(v));
    // Ugly that I have to add the vertex outside the constructor,
    // but enable_shared_from_this requires that a valid shared_ptr to this
    // already exists
    result->grow();
    return result;
  } else
    return stree;
}

DRTree::DRTree(const std::shared_ptr<DGVertex>& r) : nvertices_(0), root_(r) {}

DRTree::~DRTree() {}

void DRTree::grow() { add_vertex(root()); }

const std::shared_ptr<DGVertex>& DRTree::root() const { return root_; }

void DRTree::add_vertex(const std::shared_ptr<DGVertex>& vertex) {
  // If not root and has more than 1 parent -- it is not on the tree
  if ((vertex->num_entry_arcs() <= 1 || vertex == root())) {
    if (vertex->subtree_)
      throw ProgrammingError(
          "DRTree::add_vertex() -- vertex is on a subtree already");
    vertex->subtree_ =
        std::enable_shared_from_this<this_type>::shared_from_this();
    ++nvertices_;
#if LOCAL_DEBUG
    std::cout << "Vertex " << vertex->label()
              << " is on the following subtree:" << std::endl;
    std::cout << "  Root = " << root()->label() << std::endl;
    std::cout << "  nvertices = " << nvertices_ << std::endl;
#endif

    typedef DGVertex::ArcSetType::const_iterator aciter;
    const aciter abegin = vertex->first_exit_arc();
    const aciter aend = vertex->plast_exit_arc();
    for (aciter a = abegin; a != aend; ++a) {
      add_vertex((*a)->dest());
    }
  }
}

void DRTree::detach() { detach_from(root()); }

void DRTree::detach_from(const std::shared_ptr<DGVertex>& v) {
  if (v->subtree_.get() != this)
    return;
  else {
    v->subtree_ = std::shared_ptr<DRTree>();
    typedef DGVertex::ArcSetType::const_iterator aciter;
    const aciter abegin = v->first_exit_arc();
    const aciter aend = v->plast_exit_arc();
    for (aciter a = abegin; a != aend; ++a) {
      detach_from((*a)->dest());
    }
  }
}
