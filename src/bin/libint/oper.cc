
#include <stdexcept>
#include <cmath>
#include <rr.h>
#include <oper.h>

using namespace std;
using namespace libint2;

////////////

std::string
R12_k_G12_Descr::label_(int K)
{
  ostringstream oss;
  oss << "R12^" << K << " * G12";
  return oss.str();
}

std::string
R12_k_G12_Descr::symbol_(int K)
{
  ostringstream oss;
  oss << "R12_" << (K<0 ? "minus_" : "") << std::abs(K) << "_G12";
  return oss.str();
}

////////////

std::string
Ti_G12_Descr::label_(int K)
{
  ostringstream oss;
  oss << "[T_" << K << ",G12]";
  return oss.str();
}

std::string
Ti_G12_Descr::symbol_(int K)
{
  ostringstream oss;
  oss << "T" << K << "_G12";
  return oss.str();
}
