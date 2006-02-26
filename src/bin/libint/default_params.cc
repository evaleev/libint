
#include <default_params.h>

using namespace libint2;

const std::string CompilationParameters::Defaults::source_directory("./");
const std::string CompilationParameters::Defaults::realtype("double");

CompilationParameters::CompilationParameters() :
  max_am_eri_(Defaults::max_am_eri), max_am_g12_(Defaults::max_am_g12),
  max_vector_length_(Defaults::max_vector_length),
  vectorize_by_line_(Defaults::vectorize_by_line), unroll_threshold_(Defaults::unroll_threshold),
  source_directory_(Defaults::source_directory), use_C_linking_(Defaults::use_C_linking),
  count_flops_(Defaults::count_flops), realtype_(Defaults::realtype)
{
}

CompilationParameters::~CompilationParameters()
{
}

void
CompilationParameters::print(std::ostream& os) const
{
  using namespace std;
#ifdef INCLUDE_ERI
  os << "ERI_MAX_AM           = " << max_am_eri() << endl;
#endif
#ifdef INCLUDE_G12
  os << "ERI_MAX_G12          = " << max_am_g12() << endl;
#endif
  os << "MAX_VECTOR_LENGTH    = " << max_vector_length() << endl;
  if (max_vector_length() > 1)
    os << "VECTORIZE_BY_LINE    = " << (vectorize_by_line() ? "true" : "false") << endl;
  os << "UNROLL_THRESH        = " << unroll_threshold() << endl;
  os << "SOURCE_DIRECTORY     = " << source_directory() << endl;
  os << "USE_C_LINKING        = " << (use_C_linking() ? "true" : "false") << endl;
  os << "COUNT_FLOPS          = " << (count_flops() ? "true" : "false") << endl;
  os << "REALTYPE             = " << (realtype()) << endl;
  os << endl;
}

//////////

LibraryParameters::LibraryParameters() :
  max_stack_size_(1), max_vector_stack_size_(1),
  max_hrr_hsrank_(1), max_hrr_lsrank_(1)
{
}

LibraryParameters&
LibraryParameters::get_library_params() {
  return LP_obj_;
}

LibraryParameters LibraryParameters::LP_obj_;

//////////

const char libint2::StaticDefinitions::am_letters[StaticDefinitions::num_am_letters] = "spdfghiklmnoqrtuvwxyz";

std::string
libint2::label_to_funcname(const std::string& label)
{
  std::string result("compute");
  result += label;
  return result;
}

