
#include <context.h>

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

CodeContext::CodeContext() :
  comments_on_(true)
{
  zero_out_counters();
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
  static const unsigned int nchars = 11;
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
    "^"
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
    "_up_"
  };
};

CppCodeContext::CppCodeContext() :
  CodeContext()
{
}

CppCodeContext::~CppCodeContext()
{
}

std::string
CppCodeContext::std_header() const
{
  std::string result("#include <libint2.h>\n");
  return result;
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
  std::string ptr("&(");
  ptr += val; ptr += ")";
  return ptr;
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
CppCodeContext::fp_type() const { return "REALTYPE"; }
std::string
CppCodeContext::ptr_fp_type() const { return "REALTYPE*"; }
