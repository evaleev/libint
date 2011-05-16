
#include <string>
#include <vector>

#ifndef _libint2_src_bin_libint_code_h_
#define _libint2_src_bin_libint_code_h_

using namespace std;

namespace libint2 {

  /** Class CodeSymbols specifies a set of symbols used in a code
    */
  
  class CodeSymbols {
  public:
  CodeSymbols();
  ~CodeSymbols();
  
  void append_symbol(const std::string& s);
  unsigned int n() const;
  const std::string& symbol(unsigned int i) const;
  
  private:
  vector<std::string> symbols_;
  };
  
};

#endif

