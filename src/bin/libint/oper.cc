
#include <stdexcept>
#include <cmath>
#include <rr.h>
#include <oper.h>
#include <exception.h>

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

int
R12k_R12l_G12_Descr::K(int xyz) const {
  switch(xyz) {
    case 0: return get<0>(K_);
    case 1: return get<1>(K_);
    case 2: return get<2>(K_);
    default: throw ProgrammingError("R12k_R12l_G12_Descr::K -- argument out of range");
  }
  return 0;
}

int
R12k_R12l_G12_Descr::L(int xyz) const {
  switch(xyz) {
    case 0: return get<0>(L_);
    case 1: return get<1>(L_);
    case 2: return get<2>(L_);
    default: throw ProgrammingError("R12k_R12l_G12_Descr::L -- argument out of range");
  }
  return 0;
}

unsigned int
R12k_R12l_G12_Descr::key() const {
  unsigned int k =  
    (((((get<0>(K_)*kmax
    +get<1>(K_))*kmax
    +get<2>(K_))*kmax
    +get<0>(L_))*kmax
    +get<1>(L_))*kmax
    +get<2>(L_));
}

std::string
R12k_R12l_G12_Descr::label_(const IntVec3& K, const IntVec3& L)
{
  ostringstream oss;
  oss << "(R12x^" << get<0>(K)
      << "*R12y^" << get<1>(K)
      << "*R12z^" << get<2>(K) << ")*"
      << "(R12x^" << get<0>(L)
      << "*R12y^" << get<1>(L)
      << "*R12z^" << get<2>(L) << ") * G12";
  return oss.str();
}

std::string
R12k_R12l_G12_Descr::symbol_(const IntVec3& K, const IntVec3& L)
{
  ostringstream oss;
  oss << "(R12x" << get<0>(K)
      << "_R12y" << get<1>(K)
      << "_R12z" << get<2>(K) << "__"
      << "R12x" << get<0>(L)
      << "_R12y" << get<1>(L)
      << "_R12z" << get<2>(L) << "__G12";
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
