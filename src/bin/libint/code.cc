
#include <code.h>

using namespace std;
using namespace libint2;

CodeSymbols::CodeSymbols(): symbols_() {}

CodeSymbols::~CodeSymbols() {}

void
CodeSymbols::append_symbol(const std::string& s)
{
  symbols_.push_back(s);
}

unsigned int
CodeSymbols::n() const
{
  return symbols_.size();
}

const std::string&
CodeSymbols::symbol(unsigned int i) const
{
  return symbols_.at(i);
}

