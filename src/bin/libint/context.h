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

#include <default_params.h>
#include <entity.h>

#ifndef _libint2_src_bin_libint_codecontext_h_
#define _libint2_src_bin_libint_codecontext_h_

namespace libint2 {

class ForLoop;

/**
   CodeContext provides context for generating code
*/
class CodeContext {
 public:
  virtual ~CodeContext() {}

  /// Returns the CompilationParameters used to create this context
  const std::shared_ptr<CompilationParameters>& cparams() const;

  /// turn comments on and off (the default is on)
  void turn_comments(bool on);
  /// returns true if to print comments
  bool comments_on() const;

  /// this function resets the context to be used for the next source file
  virtual void reset();

  /// produces prefix to function declarations or definitions
  virtual std::string code_prefix() const = 0;
  /// produces postfix to function declarations or definitions
  virtual std::string code_postfix() const = 0;
  /// copyright() returns the copyright statement
  virtual std::string copyright() const = 0;
  /// std_header() returns declarations necessary for every source file
  virtual std::string std_header() const = 0;
  /// std_function_header() returns declarations and definitions necessary for
  /// every function
  virtual std::string std_function_header() const = 0;
  /// label_to_name(label) converts label to a name valid within the context of
  /// the language
  virtual std::string label_to_name(const std::string& label) const = 0;
  /** declare returns a statement which declares variable named 'name' of
      type 'type'
  */
  virtual std::string declare(const std::string& type,
                              const std::string& name) const = 0;
  /** declare_v returns a statement which declares a vector 'name' of 'nelem'
     elements of type 'type'
  */
  virtual std::string declare_v(const std::string& type,
                                const std::string& name,
                                const std::string& nelem) const = 0;
  /** decldef returns a statement which declares variable named 'name' of
      type 'type' and defines its value to be 'value'
  */
  virtual std::string decldef(const std::string& type, const std::string& name,
                              const std::string& value) = 0;
  /** assign returns a statement which assigns variable 'value'
      to variable 'name'
  */
  virtual std::string assign(const std::string& name,
                             const std::string& value) = 0;
  /** accumulate returns a statement which assigns variable 'value'
      to variable 'name'
  */
  virtual std::string accumulate(const std::string& name,
                                 const std::string& value) = 0;
  /** assign_binary_expr returns a statement which assigns binary
      expression 'left oper right' to variable 'name'
  */
  virtual std::string assign_binary_expr(const std::string& name,
                                         const std::string& left,
                                         const std::string& oper,
                                         const std::string& right) = 0;
  /** assign_ternary_expr returns a statement which assigns ternary
      expression 'arg1 oper1 arg2 oper2 arg3' to variable 'name'
  */
  virtual std::string assign_ternary_expr(const std::string& name,
                                          const std::string& arg1,
                                          const std::string& oper1,
                                          const std::string& arg2,
                                          const std::string& oper2,
                                          const std::string& arg3) = 0;
  /** accumulate_binary_expr returns a statement which accumulates binary
      expression 'left oper right' to variable 'name'
  */
  virtual std::string accumulate_binary_expr(const std::string& name,
                                             const std::string& left,
                                             const std::string& oper,
                                             const std::string& right) = 0;
  /** accumulate_ternary_expr returns a statement which accumulates ternary
      expression 'arg1 oper1 arg2 oper2 arg3' to variable 'name'
  */
  virtual std::string accumulate_ternary_expr(const std::string& name,
                                              const std::string& arg1,
                                              const std::string& oper1,
                                              const std::string& arg2,
                                              const std::string& oper2,
                                              const std::string& arg3) = 0;
  /// converts an address on the stack to its string representation
  virtual std::string stack_address(const DGVertex::Address& a) const = 0;

  /// #define a macro
  virtual std::string macro_define(const std::string& name) const = 0;
  /// #define a macro
  virtual std::string macro_define(const std::string& name,
                                   const std::string& value) const = 0;
  /// #if macro
  virtual std::string macro_if(const std::string& name) const = 0;
  /// #ifdef macro
  virtual std::string macro_ifdef(const std::string& name) const = 0;
  /// #endif
  virtual std::string macro_endif() const = 0;

  /// comment(statement) returns commented statement
  virtual std::string comment(const std::string& statement) const = 0;
  /// open a code block
  virtual std::string open_block() const = 0;
  /// close a code block
  virtual std::string close_block() const = 0;
  /// end a statement
  virtual std::string end_of_stat() const = 0;
  /// converts a value to a pointer
  virtual std::string value_to_pointer(const std::string& val) const = 0;

  /** returns a ForLoop object.
   */
  virtual std::shared_ptr<ForLoop> for_loop(
      std::string& varname, const std::shared_ptr<Entity>& less_than,
      const std::shared_ptr<Entity>& start_at =
          std::shared_ptr<Entity>(new CTimeEntity<int>(0))) const = 0;

  /// unique_name<T> returns a unique name for a variable of type T
  template <typename T>
  std::string unique_name() const;
  /// type_name<T> returns name for type T
  template <typename T>
  std::string type_name() const;
  /// returns the name of the evaluator type for task 'task'
  virtual std::string inteval_type_name(const std::string& task) const = 0;
  /// returns the name of the specialized evaluator type for task 'task'
  virtual std::string inteval_spec_type_name(const std::string& task) const = 0;
  /// returns the name of the generic evaluator type (works for any task)
  virtual std::string inteval_gen_type_name() const = 0;
  /// returns the modifier for constant variables
  virtual std::string const_modifier() const = 0;
  /// returns the modifier for mutable variables
  virtual std::string mutable_modifier() const = 0;

 protected:
  /// Lone constructor takes CompilationParams
  CodeContext(const std::shared_ptr<CompilationParameters>& cparams);
  /// generates a unique name for a floating-point variable
  virtual std::string unique_fp_name() const = 0;
  /// generates a unique name for an integer
  virtual std::string unique_int_name() const = 0;

  /// next fp index
  unsigned int next_fp_index() const;
  /// next int index
  unsigned int next_int_index() const;
  /// replaces every appearance of From with To in S
  static std::string replace_chars(const std::string& S,
                                   const std::string& From,
                                   const std::string& To);

  /// returns name of void type
  virtual std::string void_type() const = 0;
  /// returns name of integer type
  virtual std::string int_type() const = 0;
  /// returns name of array size type (e.g. size_t)
  virtual std::string size_type() const = 0;
  /// returns name of floating-point type
  virtual std::string fp_type() const = 0;
  /// returns name of pointer to floating-point type
  virtual std::string ptr_fp_type() const = 0;

 private:
  std::shared_ptr<CompilationParameters> cparams_;
  mutable unsigned int next_index_[EntityTypes::ntypes];
  void zero_out_counters() const;
  bool comments_on_;
};

/**
   CppCodeContext is an implementation of CodeContext for C++
*/
class CppCodeContext : public CodeContext,
                       public std::enable_shared_from_this<CppCodeContext> {
 public:
  CppCodeContext(const std::shared_ptr<CompilationParameters>& cparams,
                 bool vectorize = false);
  virtual ~CppCodeContext();

  /// Implementation of CodeContext::code_prefix()
  std::string code_prefix() const override;
  /// Implementation of CodeContext::code_postfix()
  std::string code_postfix() const override;
  /// Implementation of CodeContext::copyright()
  std::string copyright() const override;
  /// Implementation of CodeContext::std_header()
  std::string std_header() const override;
  /// Implementation of CodeContext::std_function_header()
  std::string std_function_header() const override;
  /// Implementation of CodeContext::label_to_name(label)
  std::string label_to_name(const std::string& label) const override;
  /// Implementation of CodeContext::declare()
  std::string declare(const std::string& type,
                      const std::string& name) const override;
  /// Implementation of CodeContext::declare_v()
  std::string declare_v(const std::string& type, const std::string& name,
                        const std::string& nelem) const override;
  // Implementation of CodeContext::decldef()
  std::string decldef(const std::string& type, const std::string& name,
                      const std::string& value) override;
  /// Implementation of CodeContext::assign()
  std::string assign(const std::string& name,
                     const std::string& value) override;
  /// Implementation of CodeContext::accumulate()
  std::string accumulate(const std::string& name,
                         const std::string& value) override;
  /// Implementation of CodeContext::assign_binary_expr()
  std::string assign_binary_expr(const std::string& name,
                                 const std::string& left,
                                 const std::string& oper,
                                 const std::string& right) override;
  // Implementation of CodeContext::assign_ternary_expr()
  std::string assign_ternary_expr(const std::string& name,
                                  const std::string& arg1,
                                  const std::string& oper1,
                                  const std::string& arg2,
                                  const std::string& oper2,
                                  const std::string& arg3) override;
  /// Implementation of CodeContext::accumulate_binary_expr()
  std::string accumulate_binary_expr(const std::string& name,
                                     const std::string& left,
                                     const std::string& oper,
                                     const std::string& right) override;
  // Implementation of CodeContext::accumulate_ternary_expr()
  std::string accumulate_ternary_expr(const std::string& name,
                                      const std::string& arg1,
                                      const std::string& oper1,
                                      const std::string& arg2,
                                      const std::string& oper2,
                                      const std::string& arg3) override;
  /// Implementation of CodeContext::stack_address()
  std::string stack_address(const DGVertex::Address& a) const override;

  /// Implementation of CodeContext::macro_define()
  std::string macro_define(const std::string& name) const override;
  /// Implementation of CodeContext::macro_define()
  std::string macro_define(const std::string& name,
                           const std::string& value) const override;
  /// Implementation of CodeContext::macro_if()
  std::string macro_if(const std::string& name) const override;
  /// Implementation of CodeContext::macro_ifdef()
  std::string macro_ifdef(const std::string& name) const override;
  /// Implementation of CodeContext::macro_endif()
  std::string macro_endif() const override;

  /// Implementation of CodeContext::comment(statement)
  std::string comment(const std::string& statement) const override;
  /// Implementation of CodeContext::open_block()
  std::string open_block() const override;
  /// Implementation of CodeContext::close_block()
  std::string close_block() const override;
  /// Implementation of CodeContext::end_of_stat()
  std::string end_of_stat() const override;
  /// Implementation of CodeContext::value_to_pointer()
  std::string value_to_pointer(const std::string& val) const override;
  /// Implementation of CodeContext::for_loop()
  std::shared_ptr<ForLoop> for_loop(
      std::string& varname, const std::shared_ptr<Entity>& less_than,
      const std::shared_ptr<Entity>& start_at) const override;

  /// Implementation of CodeContext::inteval_type_name()
  std::string inteval_type_name(const std::string& task) const override;
  /// Implementation of CodeContext::inteval_spec_type_name()
  std::string inteval_spec_type_name(const std::string& task) const override;
  /// Implementation of CodeContext::inteval_spec_type_name()
  std::string inteval_gen_type_name() const override;

 private:
  bool vectorize_;

  /// Implementation of CodeContext::unique_fp_name()
  std::string unique_fp_name() const override;
  /// Implementation of CodeContext::unique_int_name()
  std::string unique_int_name() const override;
  ///
  std::string symbol_to_pointer(const std::string& symbol);

  std::string void_type() const override;
  std::string int_type() const override;
  std::string size_type() const override;
  std::string fp_type() const override;
  std::string ptr_fp_type() const override;
  std::string const_modifier() const override;
  std::string mutable_modifier() const override;

  std::string start_expr() const;
  std::string end_expr() const;

  /// assign/accumulate if accum=false/true
  std::string assign_(const std::string& name, const std::string& value,
                      bool accum);
  /// assign/accumulate if accum=false/true
  std::string assign_binary_expr_(const std::string& name,
                                  const std::string& left,
                                  const std::string& oper,
                                  const std::string& right, bool accum);
  /// assign/accumulate if accum=false/true
  std::string assign_ternary_expr_(const std::string& name,
                                   const std::string& arg1,
                                   const std::string& oper1,
                                   const std::string& arg2,
                                   const std::string& oper2,
                                   const std::string& arg3, bool accum);
};

};  // namespace libint2

#endif
