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

#include <codeblock.h>
#include <context.h>
#include <default_params.h>

#include <cassert>
#include <cstdio>

using namespace std;
using namespace libint2;

namespace libint2 {

template <>
std::string CodeContext::unique_name<EntityTypes::FP>() const {
  return unique_fp_name();
}
template <>
std::string CodeContext::unique_name<EntityTypes::Int>() const {
  return unique_int_name();
}

template <>
std::string CodeContext::type_name<void>() const {
  return void_type();
}
template <>
std::string CodeContext::type_name<int>() const {
  return int_type();
}
template <>
std::string CodeContext::type_name<size_t>() const {
  return size_type();
}
template <>
std::string CodeContext::type_name<const int>() const {
  return const_modifier() + int_type();
}
template <>
std::string CodeContext::type_name<double>() const {
  return fp_type();
}
template <>
std::string CodeContext::type_name<double*>() const {
  return ptr_fp_type();
}
template <>
std::string CodeContext::type_name<const double*>() const {
  return const_modifier() + ptr_fp_type();
}
template <>
std::string CodeContext::type_name<double* const>() const {
  return ptr_fp_type() + const_modifier();
}
};  // namespace libint2

CodeContext::CodeContext(const std::shared_ptr<CompilationParameters>& cparams)
    : cparams_(cparams), comments_on_(false) {
  zero_out_counters();
}

const std::shared_ptr<CompilationParameters>& CodeContext::cparams() const {
  return cparams_;
}

bool CodeContext::comments_on() const { return comments_on_; }

unsigned int CodeContext::next_fp_index() const {
  return next_index_[EntityTypes::FP::type2int()]++;
}

unsigned int CodeContext::next_int_index() const {
  return next_index_[EntityTypes::Int::type2int()]++;
}

void CodeContext::zero_out_counters() const {
  for (unsigned int i = 0; i < EntityTypes::ntypes; i++) next_index_[i] = 0;
}

void CodeContext::reset() { zero_out_counters(); }

std::string CodeContext::replace_chars(const std::string& S,
                                       const std::string& From,
                                       const std::string& To) {
  typedef std::string::size_type size_type;

  const unsigned int max_niter = 1000;
  unsigned int niter = 0;
  std::string curr_str(S);
  size_type curr_pos = curr_str.find(From, 0);
  while (curr_pos != std::string::npos) {
    niter++;
    curr_str.replace(curr_pos, From.length(), To, 0, To.length());
    curr_pos += To.length() - From.length();
    curr_pos = curr_str.find(From, curr_pos);
    if (niter >= max_niter)
      throw std::runtime_error(
          "CodeContext::replace_chars() -- infinite recursion detected");
  }
  return curr_str;
}

//////////////

namespace ForbiddenCppCharacters {
static const unsigned int nchars = 16;
static const char chars[nchars][2] = {"{", "}", "(", ")", " ", "+", "-", "/",
                                      "*", "|", "^", "[", "]", ",", "<", ">"};
static const char subst_chars[nchars][20] = {
    "",        "",  "__",   "__",   "",     "_plus_", "_minus_", "_over_",
    "_times_", "_", "_up_", "_sB_", "_Sb_", "_c_",    "_aB_",    "_Ab_"};
};  // namespace ForbiddenCppCharacters

CppCodeContext::CppCodeContext(
    const std::shared_ptr<CompilationParameters>& cparams, bool vectorize)
    : CodeContext(cparams), vectorize_(vectorize) {}

CppCodeContext::~CppCodeContext() {}

std::string CppCodeContext::code_prefix() const {
  if (cparams()->use_C_linking()) {
    return "#ifdef __cplusplus\nLIBINT_PRAGMA_CLANG(diagnostic "
           "push)\nLIBINT_PRAGMA_CLANG(diagnostic ignored "
           "\"-Wunused-variable\")\nLIBINT_PRAGMA_GCC(diagnostic "
           "push)\nLIBINT_PRAGMA_GCC(diagnostic ignored "
           "\"-Wunused-variable\")\nextern \"C\" {\n#endif\n";
  }
  return "";
}

std::string CppCodeContext::code_postfix() const {
  if (cparams()->use_C_linking()) {
    return "#ifdef __cplusplus\n};\nLIBINT_PRAGMA_CLANG(diagnostic "
           "pop)\nLIBINT_PRAGMA_GCC(diagnostic pop)\n#endif\n";
  }
  return "";
}

std::string CppCodeContext::copyright() const {
  std::ostringstream oss;
  using std::endl;
  oss << "/*" << endl
      << " *  Copyright (C) 2004-2026 Edward F. Valeev" << endl
      << " *" << endl
      << " *  This file is part of Libint library." << endl
      << " *" << endl
      << " *  Libint library is free software: you can redistribute it and/or "
         "modify"
      << endl
      << " *  it under the terms of the GNU Lesser General Public License as "
         "published by"
      << endl
      << " *  the Free Software Foundation, either version 3 of the License, or"
      << endl
      << " *  (at your option) any later version." << endl
      << " *" << endl
      << " *  Libint library is distributed in the hope that it will be useful,"
      << endl
      << " *  but WITHOUT ANY WARRANTY; without even the implied warranty of"
      << endl
      << " *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
      << endl
      << " *  GNU Lesser General Public License for more details." << endl
      << " *" << endl
      << " *  You should have received a copy of the GNU Lesser General Public "
         "License"
      << endl
      << " *  along with Libint library.  If not, see "
         "<http://www.gnu.org/licenses/>."
      << endl
      << " *" << endl
      << " */" << endl
      << endl;
  return oss.str();
}

std::string CppCodeContext::std_header() const {
  std::string result("#include <libint2.h>\n");
  return result;
}

std::string CppCodeContext::std_function_header() const {
  ostringstream oss;
  if (vectorize_) {
    oss << "const int veclen = inteval->veclen;\n";
  }
  return oss.str();
}

std::string CppCodeContext::label_to_name(const std::string& label) const {
  std::string str = label;
  for (unsigned int c = 0; c < ForbiddenCppCharacters::nchars; c++) {
    str = replace_chars(str, ForbiddenCppCharacters::chars[c],
                        ForbiddenCppCharacters::subst_chars[c]);
  }
  return str;
}

std::string CppCodeContext::declare(const std::string& type,
                                    const std::string& name) const {
  ostringstream oss;

  oss << type << " " << name << end_of_stat() << endl;

  return oss.str();
}

std::string CppCodeContext::declare_v(const std::string& type,
                                      const std::string& name,
                                      const std::string& nelem) const {
  ostringstream oss;

  oss << type << " " << name << "[" << nelem << "]" << end_of_stat() << endl;

  return oss.str();
}

std::string CppCodeContext::decldef(const std::string& type,
                                    const std::string& name,
                                    const std::string& value) {
  ostringstream oss;

  oss << type << " " << assign(name, value);

  return oss.str();
}

std::string CppCodeContext::assign(const std::string& name,
                                   const std::string& value) {
  return assign_(name, value, false);
}

std::string CppCodeContext::accumulate(const std::string& name,
                                       const std::string& value) {
  return assign_(name, value, true);
}

std::string CppCodeContext::assign_(const std::string& name,
                                    const std::string& value, bool accum) {
  ostringstream oss;

  if (vectorize_) {
    std::string symb0 = unique_fp_name();
    std::string symb1 = unique_fp_name();
    std::string ptr0 = symbol_to_pointer(name);
    std::string ptr1 = symbol_to_pointer(value);
    bool symb1_is_a_const = (ptr1.length() == 0);
    oss << "LIBINT2_REALTYPE* " << symb0 << " = " << symbol_to_pointer(name)
        << end_of_stat() << endl;
    oss << "__assume_aligned(" << symb0 << ", 16)" << end_of_stat() << endl;
    if (!symb1_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb1 << " = " << symbol_to_pointer(value)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb1 << ", 16)" << end_of_stat() << endl;
    }

    oss << start_expr();
    oss << symb0 << "[v]" << (accum ? " += " : " = ")
        << (symb1_is_a_const ? value : symb1)
        << (symb1_is_a_const ? " " : "[v] ");
  } else {
    oss << start_expr();
    oss << name << (accum ? " += " : " = ") << value;
  }
  oss << end_of_stat() << endl;
  oss << end_expr();

  return oss.str();
}

std::string CppCodeContext::assign_binary_expr(const std::string& name,
                                               const std::string& left,
                                               const std::string& oper,
                                               const std::string& right) {
  return assign_binary_expr_(name, left, oper, right, false);
}

std::string CppCodeContext::accumulate_binary_expr(const std::string& name,
                                                   const std::string& left,
                                                   const std::string& oper,
                                                   const std::string& right) {
  return assign_binary_expr_(name, left, oper, right, true);
}

std::string CppCodeContext::assign_binary_expr_(const std::string& name,
                                                const std::string& left,
                                                const std::string& oper,
                                                const std::string& right,
                                                bool accum) {
  ostringstream oss;

  if (vectorize_) {
    std::string symb0 = unique_fp_name();
    std::string symb1 = unique_fp_name();
    std::string symb2 = unique_fp_name();
    std::string ptr0 = symbol_to_pointer(name);
    std::string ptr1 = symbol_to_pointer(left);
    std::string ptr2 = symbol_to_pointer(right);
    bool symb1_is_a_const = (ptr1.length() == 0);
    bool symb2_is_a_const = (ptr2.length() == 0);
    oss << "LIBINT2_REALTYPE* " << symb0 << " = " << symbol_to_pointer(name)
        << end_of_stat() << endl;
    oss << "__assume_aligned(" << symb0 << ", 16)" << end_of_stat() << endl;
    if (!symb1_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb1 << " = " << symbol_to_pointer(left)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb1 << ", 16)" << end_of_stat() << endl;
    }
    if (!symb2_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb2 << " = " << symbol_to_pointer(right)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb2 << ", 16)" << end_of_stat() << endl;
    }

    oss << start_expr();
    oss << symb0 << "[v]" << (accum ? " += " : " = ")
        << (symb1_is_a_const ? left : symb1)
        << (symb1_is_a_const ? " " : "[v] ") << oper << " "
        << (symb2_is_a_const ? right : symb2)
        << (symb2_is_a_const ? "" : "[v]");
  } else {
    oss << start_expr();
    oss << name << (accum ? " += " : " = ") << left << " " << oper << " "
        << right;
  }
  oss << end_of_stat() << endl;
  oss << end_expr();

  return oss.str();
}

std::string CppCodeContext::assign_ternary_expr(const std::string& name,
                                                const std::string& arg1,
                                                const std::string& oper1,
                                                const std::string& arg2,
                                                const std::string& oper2,
                                                const std::string& arg3) {
  return assign_ternary_expr_(name, arg1, oper1, arg2, oper2, arg3, false);
}

std::string CppCodeContext::assign_ternary_expr_(
    const std::string& name, const std::string& arg1, const std::string& oper1,
    const std::string& arg2, const std::string& oper2, const std::string& arg3,
    bool accum) {
  ostringstream oss;

  // this should only be invoked for FMA, i.e. oper1 = "*" and oper2 = "+" or
  // "-"
  assert(oper1 == "*");
  assert(oper2 == "+" || oper2 == "-");

  if (vectorize_) {
    std::string symb0 = unique_fp_name();
    std::string symb1 = unique_fp_name();
    std::string symb2 = unique_fp_name();
    std::string symb3 = unique_fp_name();
    std::string ptr0 = symbol_to_pointer(name);
    std::string ptr1 = symbol_to_pointer(arg1);
    std::string ptr2 = symbol_to_pointer(arg2);
    std::string ptr3 = symbol_to_pointer(arg3);
    bool symb1_is_a_const = (ptr1.length() == 0);
    bool symb2_is_a_const = (ptr2.length() == 0);
    bool symb3_is_a_const = (ptr3.length() == 0);
    oss << "LIBINT2_REALTYPE* " << symb0 << " = " << symbol_to_pointer(name)
        << end_of_stat() << endl;
    oss << "__assume_aligned(" << symb0 << ", 16)" << end_of_stat() << endl;
    if (!symb1_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb1 << " = " << symbol_to_pointer(arg1)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb1 << ", 16)" << end_of_stat() << endl;
    }
    if (!symb2_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb2 << " = " << symbol_to_pointer(arg2)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb2 << ", 16)" << end_of_stat() << endl;
    }
    if (!symb3_is_a_const) {
      oss << "LIBINT2_REALTYPE* " << symb3 << " = " << symbol_to_pointer(arg3)
          << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb3 << ", 16)" << end_of_stat() << endl;
    }

    oss << start_expr();
    oss << symb0 << "[v]" << (accum ? " += " : " = ");
#if LIBINT_GENERATE_FMA
    oss << "libint2::fma_" << (oper2 == "+" ? "plus" : "minus") << "("
        << (symb1_is_a_const ? arg1 : symb1)
        << (symb1_is_a_const ? " " : "[v] ") << ","
        << (symb2_is_a_const ? arg2 : symb2) << (symb2_is_a_const ? "" : "[v]")
        << "," << (symb3_is_a_const ? arg3 : symb3)
        << (symb3_is_a_const ? "" : "[v]") << ")";
#else
    oss << (symb1_is_a_const ? arg1 : symb1)
        << (symb1_is_a_const ? " " : "[v] ") << "* "
        << (symb2_is_a_const ? arg2 : symb2)
        << (symb2_is_a_const ? " " : "[v] ") << oper2
        << (symb3_is_a_const ? arg3 : symb3) << (symb3_is_a_const ? "" : "[v]");
#endif
  } else {
    oss << start_expr();
    oss << name << (accum ? " += " : " = ");
#if LIBINT_GENERATE_FMA
    oss << "libint2::fma_" << (oper2 == "+" ? "plus" : "minus") << "(" << arg1
        << ", " << arg2 << ", " << arg3 << ")";
#else
    oss << arg1 << " * " << arg2 << " " << oper2 << " " << arg3;
#endif
  }
  oss << end_of_stat() << endl;
  oss << end_expr();

  return oss.str();
}

std::string CppCodeContext::accumulate_ternary_expr(const std::string& name,
                                                    const std::string& arg1,
                                                    const std::string& oper1,
                                                    const std::string& arg2,
                                                    const std::string& oper2,
                                                    const std::string& arg3) {
  return assign_ternary_expr_(name, arg1, oper1, arg2, oper2, arg3, true);
}

std::string CppCodeContext::symbol_to_pointer(const std::string& symbol) {
  std::string::size_type loc = symbol.find("stack");
  // if this quantity is on stack then the symbol is a scalar
  if (loc != std::string::npos) {
    ostringstream oss;
    oss << "(&(" << symbol << "))";
    return oss.str();
  }

  // if this quantity is a part of Libint_t then the symbol is a vector
  // otherwise it's a constant
  loc = symbol.find("inteval");
  if (loc != std::string::npos)
    return symbol;
  else
    return "";
}

std::string CppCodeContext::start_expr() const {
  if (vectorize_)
    return "#ifdef __INTEL_COMPILER\n#pragma ivdep\n#endif\nfor(int v=0; "
           "v<veclen; v++) {\n";
  else
    return "";
}

std::string CppCodeContext::end_expr() const {
  if (vectorize_)
    return "}\n";
  else
    return "";
}

std::string CppCodeContext::stack_address(const DGVertex::Address& a) const {
  ostringstream oss;
  if (vectorize_)
    oss << "(" << a << ")*veclen";
  else
    oss << a;
  return oss.str();
}

std::string CppCodeContext::macro_define(const std::string& name) const {
  ostringstream oss;
  oss << "#define " << name << endl;
  return oss.str();
}

std::string CppCodeContext::macro_define(const std::string& name,
                                         const std::string& value) const {
  ostringstream oss;
  oss << "#define " << name << " " << value << endl;
  return oss.str();
}

std::string CppCodeContext::macro_if(const std::string& name) const {
  ostringstream oss;
  oss << "#if " << name << endl;
  return oss.str();
}

std::string CppCodeContext::macro_ifdef(const std::string& name) const {
  ostringstream oss;
  oss << "#ifdef " << name << endl;
  return oss.str();
}

std::string CppCodeContext::macro_endif() const {
  ostringstream oss;
  oss << "#endif" << endl;
  return oss.str();
}

std::string CppCodeContext::comment(const std::string& statement) const {
  std::string result("/** ");
  result += statement;
  result += " */";
  return result;
}

std::string CppCodeContext::open_block() const { return " {\n"; }

std::string CppCodeContext::close_block() const { return "}\n"; }

std::string CppCodeContext::end_of_stat() const {
  static const std::string ends(";");
  return ends;
}

std::string CppCodeContext::value_to_pointer(const std::string& val) const {
  if (!vectorize_) {
    std::string ptr("&(");
    ptr += val;
    ptr += ")";
    return ptr;
  } else {
    return val;
  }
}

std::shared_ptr<ForLoop> CppCodeContext::for_loop(
    std::string& varname, const std::shared_ptr<Entity>& less_than,
    const std::shared_ptr<Entity>& start_at) const {
  // not implemented
  abort();
}

std::string CppCodeContext::unique_fp_name() const {
  char result[80];
  snprintf(result, 80, "fp%d", next_fp_index());
  return result;
}

std::string CppCodeContext::unique_int_name() const {
  char result[80];
  snprintf(result, 80, "i%d", next_int_index());
  return result;
}

std::string CppCodeContext::void_type() const { return "void"; }
std::string CppCodeContext::int_type() const { return "int"; }
std::string CppCodeContext::size_type() const { return "size_t"; }
std::string CppCodeContext::fp_type() const {
  if (!vectorize_)
    return "LIBINT2_REALTYPE";
  else
    return ptr_fp_type();
}
std::string CppCodeContext::ptr_fp_type() const { return "LIBINT2_REALTYPE*"; }
std::string CppCodeContext::const_modifier() const { return "const "; }
std::string CppCodeContext::mutable_modifier() const {
  return "#ifdef __cplusplus\nmutable \n#endif\n";
}

std::string CppCodeContext::inteval_type_name(const std::string& tlabel) const {
  if (cparams()->single_evaltype())
    return inteval_gen_type_name();
  else
    return inteval_spec_type_name(tlabel);
}

std::string CppCodeContext::inteval_spec_type_name(
    const std::string& tlabel) const {
  ostringstream oss;
  oss << "Libint_" << tlabel << "_t";
  return oss.str();
}

std::string CppCodeContext::inteval_gen_type_name() const { return "Libint_t"; }
