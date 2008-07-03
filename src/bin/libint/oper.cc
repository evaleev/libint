
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

unsigned int
R12k_R12l_G12_Descr::key() const {
  unsigned int k =  
    (((((K_[0]*kmax
    +K_[1])*kmax
    +K_[2])*kmax
    +L_[0])*kmax
    +L_[1])*kmax
    +L_[2]);
}

std::string
R12k_R12l_G12_Descr::label_(const IntVec3& K, const IntVec3& L)
{
  ostringstream oss;
  oss << "(R12x^" << K[0]
      << "*R12y^" << K[1]
      << "*R12z^" << K[2] << ")*"
      << "(R12x^" << L[0]
      << "*R12y^" << L[1]
      << "*R12z^" << L[2] << ") * G12";
  return oss.str();
}

std::string
R12k_R12l_G12_Descr::symbol_(const IntVec3& K, const IntVec3& L)
{
  ostringstream oss;
  oss << "(R12x" << K[0]
      << "_R12y" << K[1]
      << "_R12z" << K[2] << "__"
      << "R12x" << L[0]
      << "_R12y" << L[1]
      << "_R12z" << L[2] << "__G12";
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
