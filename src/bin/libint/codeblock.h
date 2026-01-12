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

#include <entity.h>

#include <string>

#ifndef _libint2_src_bin_libint_codeblock_h_
#define _libint2_src_bin_libint_codeblock_h_

namespace libint2 {

class CodeContext;

class CodeBlock {
 public:
  CodeBlock(const std::shared_ptr<CodeContext>& context) : context_(context) {}
  virtual ~CodeBlock() {}

  std::shared_ptr<CodeContext> context() const { return context_; }

  /// Opens a code block
  virtual std::string open() = 0;
  /// Close a code block
  virtual std::string close() = 0;

 private:
  std::shared_ptr<CodeContext> context_;
};

class ForLoop : public CodeBlock {
 public:
  ForLoop(const std::shared_ptr<CodeContext>& context, std::string& varname,
          const std::shared_ptr<Entity>& less_than,
          const std::shared_ptr<Entity>& start_at);
  virtual ~ForLoop();

  /// Implementation of CodeBlock::open()
  std::string open() override;
  /// Implementation of CodeBlock::close()
  std::string close() override;

 private:
  std::string varname_;
  std::shared_ptr<Entity> less_than_;
  std::shared_ptr<Entity> start_at_;

  // checks less_than_ and start_at_ and initializes
  // lt_expr_, sa_expr_, and dummy_loop_
  void init_();
  // string representations of less_than_ and start_at_
  std::string lt_expr_;
  std::string sa_expr_;
  // dummy_loop_ is true if the loop is guaranteed to iterate once
  // i.e. can replace the loop with a constant variable definition
  bool dummy_loop_;
};

};  // namespace libint2

#endif
