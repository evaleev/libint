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

#include <stdexcept>

#ifndef _libint2_src_bin_libint_exception_h_
#define _libint2_src_bin_libint_exception_h_

namespace libint2 {

class DGVertex;

class InvalidDecrement : public std::logic_error {
 public:
  InvalidDecrement(const std::string& a) : logic_error(a){};
};

class CannotAddArc : public std::logic_error {
 public:
  CannotAddArc(const std::string& a) : logic_error(a){};
};

#if 0
  /** This exception class is used to pass the pointer to the vertex on the graph
   */
  class VertexAlreadyOnStack : public std::logic_error {
    
  public:
    VertexAlreadyOnStack(const std::shared_ptr<DGVertex>& vertex) :
      logic_error("DirectedGraph -- vertex already on stack"), vertex_(vertex) {}
    ~VertexAlreadyOnStack() noexcept {}

    std::shared_ptr<DGVertex> vertex() const { return vertex_; }

  private:
    // Vertex on the stack
    std::shared_ptr<DGVertex> vertex_;
    
  };
#endif

/** This exception class is used to notify that a graph operation cannot be
 * performed
 */
class CannotPerformOperation : public std::logic_error {
 public:
  CannotPerformOperation(const std::string& msg) : logic_error(msg) {}
  virtual ~CannotPerformOperation() noexcept {}
};

/// This exception used to indicate that some property is not set
template <class T>
class NotSet : public std::logic_error {
 public:
  NotSet(const std::string& a) : logic_error(a){};
};

/** This exception used to indicate that some code hasn't been
    developed or generalized yet */
class CodeDoesNotExist : public std::logic_error {
 public:
  CodeDoesNotExist(const std::string& a) : logic_error(a) {}
};

/** This exception used to indicate some programming error */
class ProgrammingError : public std::logic_error {
 public:
  ProgrammingError(const std::string& a) : logic_error(a) {}
};

/** This exception used to indicate some error in the user-provided input */
class InputError : public std::logic_error {
 public:
  InputError(const std::string& a) : logic_error(a) {}
};

};  // namespace libint2

#endif
