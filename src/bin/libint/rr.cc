
#include <stdexcept>

#include <rr.h>

using namespace libint2;

CartesianGaussian::CartesianGaussian(int nx, int ny, int nz) :
  GaussianShell(nx+ny+nz)
{
    n_[0] = nx;
    n_[1] = ny;
    n_[2] = nz;
}

Operator::Operator(const std::string& descr, const std::string& id, char np,
                   const vector<char>& psymm) :
  descr_(descr), id_(id), np_(np), psymm_(psymm)
{
  int nij = np_*(np_-1)/2;
  if (nij != psymm_.size())
    throw std::runtime_error("Operator::Operator() -- size of psymm is inconsistent with np");
}

Operator::~Operator()
{
}

const std::string&
Operator::descr() const
{
  return descr_;
}

const std::string&
Operator::id() const
{
  return id_;
}

const int
Operator::psymm(int i, int j) const
{
  if (i < 0 || i >= np_)
    throw std::runtime_error("Operator::psymm(i,j) -- index i out of bounds");
  if (j < 0 || j >= np_)
    throw std::runtime_error("Operator::psymm(i,j) -- index j out of bounds");
  if (i == j)
    return 1;
  int ij = (i > j) ? i*(i-1)/2+j : j*(j-1)/2+i;
  int pij = psymm_[ij];
  return pij;
}

namespace libint2 {

static vector<char> symm_2p(1,1);
Operator TwoERep("Two-Electron Repulsion Energy",
                 "2ERep",
                 2,
                 symm_2p);
Operator TwoEDist("Two-Electron Distance",
                  "2EDist",
                  2,
                  symm_2p);
                  
};

IntegralType::IntegralType(const Operator& O, const vector<char>& nbra, const vector<char>& nket) :
  O_(O), nbra_(nbra), nket_(nket)
{
}

IntegralType::~IntegralType()
{
}


RecurrenceRelation::RecurrenceRelation()
{
}

RecurrenceRelation::~RecurrenceRelation()
{
}

