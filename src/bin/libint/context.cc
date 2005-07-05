
#include <context.h>
#include <codeblock.h>
#include <default_params.h>

using namespace libint2;

namespace libint2 {

  template <>
  std::string CodeContext::unique_name<EntityTypes::FP>()
  {
    return unique_fp_name();
  }
  template <>
  std::string CodeContext::unique_name<EntityTypes::Int>()
  {
    return unique_int_name();
  }

  template <>
  std::string CodeContext::type_name<void>()
  {
    return void_type();
  }
  template <>
  std::string CodeContext::type_name<int>()
  {
    return int_type();
  }
  template <>
  std::string CodeContext::type_name<const int>()
  {
    return const_modifier() + int_type();
  }
  template <>
  std::string CodeContext::type_name<double>()
  {
    return fp_type();
  }
  template <>
  std::string CodeContext::type_name<double*>()
  {
    return ptr_fp_type();
  }
};

CodeContext::CodeContext(const SafePtr<CompilationParameters>& cparams) :
  cparams_(cparams),
  comments_on_(false)
{
  zero_out_counters();
}

const SafePtr<CompilationParameters>&
CodeContext::cparams() const
{
  return cparams_;
}

bool
CodeContext::comments_on() const { return comments_on_; }

unsigned int
CodeContext::next_fp_index()
{
  return next_index_[EntityTypes::FP::type2int()]++;
}

unsigned int
CodeContext::next_int_index()
{
  return next_index_[EntityTypes::Int::type2int()]++;
}

void
CodeContext::zero_out_counters()
{
  for(int i=0; i<EntityTypes::ntypes; i++)
    next_index_[i] = 0;
}

void
CodeContext::reset()
{
  zero_out_counters();
}

std::string
CodeContext::replace_chars(const std::string& S, const std::string& From, const std::string& To)
{
  typedef std::string::size_type size_type;

  const unsigned int max_niter = 1000;
  unsigned int niter = 0;
  std::string curr_str(S);
  size_type curr_pos = curr_str.find(From,0);
  while (curr_pos != std::string::npos) {
    niter++;
    curr_str.replace(curr_pos,From.length(),To,0,To.length());
    curr_pos += To.length() - From.length();
    curr_pos = curr_str.find(From,curr_pos);
    if (niter >= max_niter)
      throw std::runtime_error("CodeContext::replace_chars() -- infinite recursion detected");
  }
  return curr_str;
}

//////////////

namespace ForbiddenCppCharacters {
  static const unsigned int nchars = 16;
  static const char chars[nchars][2] = {
    "{",
    "}",
    "(",
    ")",
    " ",
    "+",
    "-",
    "/",
    "*",
    "|",
    "^",
    "[",
    "]",
    ",",
    "<",
    ">"
  };
  static const char subst_chars[nchars][20] = {
    "",
    "",
    "__",
    "__",
    "",
    "_plus_",
    "_minus_",
    "_over_",
    "_times_",
    "_",
    "_up_",
    "_sB_",
    "_Sb_",
    "_c_",
    "_aB_",
    "_Ab_"
  };
};

CppCodeContext::CppCodeContext(const SafePtr<CompilationParameters>& cparams, bool vectorize) :
  CodeContext(cparams), vectorize_(vectorize)
{
}

CppCodeContext::~CppCodeContext()
{
}

std::string
CppCodeContext::code_prefix() const
{
  if (cparams()->use_C_linking()) {
    return "extern \"C\" {\n";
  }
}

std::string
CppCodeContext::code_postfix() const
{
  if (cparams()->use_C_linking()) {
    return "};\n";
  }
}

std::string
CppCodeContext::std_header() const
{
  std::string result("#include <libint2.h>\n");
  return result;
}

std::string
CppCodeContext::std_function_header() const
{
  ostringstream oss;
  if(vectorize_) {
    oss << "const unsigned int veclength = libint->veclength;\n";
  }
  return oss.str();
}

std::string
CppCodeContext::label_to_name(const std::string& label) const
{
  std::string str = label;
  for(int c=0; c<ForbiddenCppCharacters::nchars; c++) {
    str = replace_chars(str,ForbiddenCppCharacters::chars[c],ForbiddenCppCharacters::subst_chars[c]);
  }
  return str;
}

std::string
CppCodeContext::declare(const std::string& type,
                        const std::string& name) const
{
  ostringstream oss;
  
  oss << type << " " << name << end_of_stat() << endl;

  return oss.str();
}

std::string
CppCodeContext::decldef(const std::string& type,
                        const std::string& name,
                        const std::string& value)
{
  ostringstream oss;
  
  oss << type << " " << assign(name,value);

  return oss.str();
}

std::string
CppCodeContext::assign(const std::string& name,
                       const std::string& value)
{
  ostringstream oss;
  
  if (vectorize_) {
    std::string symb0 = unique_fp_name();
    std::string symb1 = unique_fp_name();
    std::string ptr0 = symbol_to_pointer(name);
    std::string ptr1 = symbol_to_pointer(value);
    bool symb1_is_a_const = (ptr1.length() == 0);
    oss << "REALTYPE* " << symb0 << " = "
    << symbol_to_pointer(name) << end_of_stat() << endl;
    oss << "__assume_aligned(" << symb0 << ", 16)" << end_of_stat() << endl;
    if (!symb1_is_a_const) {
      oss << "REALTYPE* " << symb1 << " = "
      <<  symbol_to_pointer(value) << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb1 << ", 16)" << end_of_stat() << endl;
    }
    
    oss << start_expr();
    oss << symb0 << "[v] = "
    << (symb1_is_a_const ? value : symb1)
    << (symb1_is_a_const ? " " : "[v] ");
  }
  else {
    oss << start_expr();
    oss << name << " = " << value;
  }
  oss << end_of_stat() << endl;
  oss << end_expr();

  return oss.str();
}

std::string
CppCodeContext::assign_binary_expr(const std::string& name,
                                   const std::string& left,
                                   const std::string& oper,
                                   const std::string& right)
{
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
    oss << "REALTYPE* " << symb0 << " = "
    << symbol_to_pointer(name) << end_of_stat() << endl;
    oss << "__assume_aligned(" << symb0 << ", 16)" << end_of_stat() << endl;
    if (!symb1_is_a_const) {
      oss << "REALTYPE* " << symb1 << " = "
      <<  symbol_to_pointer(left) << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb1 << ", 16)" << end_of_stat() << endl;
    }
    if (!symb2_is_a_const) {
      oss << "REALTYPE* " << symb2 << " = "
      << symbol_to_pointer(right) << end_of_stat() << endl;
      oss << "__assume_aligned(" << symb2 << ", 16)" << end_of_stat() << endl;
    }
    
    oss << start_expr();
    oss << symb0 << "[v] = "
    << (symb1_is_a_const ? left : symb1)
    << (symb1_is_a_const ? " " : "[v] ")
    << oper << " "
    << (symb2_is_a_const ? right : symb2)
    << (symb2_is_a_const ? "" : "[v]");
  }
  else {
    oss << start_expr();
    oss << name << " = " << left << " "
        << oper << " " << right;
  }
  oss << end_of_stat() << endl;
  oss << end_expr();

  return oss.str();
}

std::string
CppCodeContext::symbol_to_pointer(const std::string& symbol)
{
  std::string::size_type loc = symbol.find("stack");
  // if this quantity is on stack then the symbol is a scalar
  if (loc != std::string::npos) {
    ostringstream oss;
    oss << "(&(" << symbol << "))";
    return oss.str();
  }
  
  // if this quantity is a part of Libint_t then the symbol is a vector
  // otherwise it's a constant
  loc = symbol.find("libint");
  if (loc != std::string::npos)
    return symbol;
  else
    return "";
}

std::string
CppCodeContext::start_expr() const
{
  if (vectorize_)
    return "for(int v=0; v<veclength; v++) {\n";
  else
    return "";
}


std::string
CppCodeContext::end_expr() const
{
  if (vectorize_)
    return "}\n";
  else
    return "";
}


std::string
CppCodeContext::stack_address(const DGVertex::Address& a) const
{
  ostringstream oss;
  if (vectorize_)
    oss << "(" << a << ")*veclength";
  else
    oss << a;
  return oss.str();
}


std::string
CppCodeContext::comment(const std::string& statement) const
{
  const char endl_str[] = "\n";
  const char comment_str1[] = "/// ";
  const char comment_str2[] = "\n/// ";
  std::string result(comment_str1);
  result += statement;
  return replace_chars(result,endl_str,comment_str2);
}

std::string
CppCodeContext::open_block() const
{
  return " {\n";
}

std::string
CppCodeContext::close_block() const
{
  return "}\n";
}

std::string
CppCodeContext::end_of_stat() const
{
  static const std::string ends(";");
  return ends;
}

std::string
CppCodeContext::value_to_pointer(const std::string& val) const
{
  if (!vectorize_) {
    std::string ptr("&(");
    ptr += val; ptr += ")";
    return ptr;
  }
  else {
    return val;
  }
}

SafePtr<ForLoop>
CppCodeContext::for_loop(std::string& varname, const SafePtr<Entity>& less_than,
                         const SafePtr<Entity>& start_at) const
{
  //return SafePtr<ForLoop>(new ForLoop(SafePtr_from_this(), varname, less_than, start_at));
}

std::string
CppCodeContext::unique_fp_name()
{
  char result[80];
  sprintf(result,"fp%d", next_fp_index());
  return result;
}

std::string
CppCodeContext::unique_int_name()
{
  char result[80];
  sprintf(result,"i%d", next_int_index());
  return result;
}

std::string
CppCodeContext::void_type() const { return "void"; }
std::string
CppCodeContext::int_type() const { return "int"; }
std::string
CppCodeContext::fp_type() const
{
  if (!vectorize_)
    return "REALTYPE";
  else
    return ptr_fp_type();
}
std::string
CppCodeContext::ptr_fp_type() const { return "REALTYPE*"; }
std::string
CppCodeContext::const_modifier() const { return "const "; }

