/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

/// This program tests Libint library by computing 2-body repulsion integrals (4, 3, and 2-center varieties)
/// and (optionally) their derivatives using Libint and a dumb but fool-proof reference method

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/chrono.hpp>

#include <rr.h>
#include <iter.h>
#include <util.h>
#include <libint2/intpart_iter.h>
#include <policy_spec.h>
#include <global_macros.h>

#include <libint2.h>
#include <test_eri/eri.h>
#include <test_eri/prep_libint2.h>

// Andrei's integrals package
#ifdef INCLUDE_RYSQ
#include <rysq.hpp>
#endif

// MPQC used to make realistic basis set examples
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/uncontract.h>
#include <chemistry/qc/basis/gaussshell.h>

using namespace std;
using namespace libint2;

#define INCLUDE_RYSQ 1
const double RELATIVE_DEVIATION_THRESHOLD = 1.0E-9; // indicate failure if any integral differs in relative sense by more than this

/// change to true to skip verification and do some timing simulation
const bool do_timing_only = true;

libint2::FmEval_Chebyshev3 fmeval_chebyshev(18);
libint2::FmEval_Taylor<double,6> fmeval_taylor(18, 1e-15);

// profile using a realistic basis set
struct TestBasis {
    TestBasis();
    sc::Ref<sc::GaussianBasisSet> basis;
    std::map<unsigned int, std::vector<unsigned int> > l_to_shells;  // l -> shells of this angular momentum
    sc::Ref<sc::GaussianBasisSet> df_basis;
    std::map<unsigned int, std::vector<unsigned int> > df_l_to_shells;  // l -> shells of this angular momentum
#ifdef INCLUDE_RYSQ
    std::vector< rysq::Shell > rysq_basis;
    std::map<unsigned int, std::vector<unsigned int> > rysq_l_to_shells;
    std::vector< rysq::Shell > df_rysq_basis;
    std::map<unsigned int, std::vector<unsigned int> > df_rysq_l_to_shells;
    std::vector< rysq::Shell > rysqized_basis; // MPQC basis in RYSQ format
    std::vector< rysq::Shell > rysqized_df_basis; // MPQC df_basis in RYSQ format
    static const bool merge_shells_on_same_atom = true; // true corresponds to the RYSQ representation of basis

    static std::tuple<std::vector< rysq::Shell >, std::map<unsigned int, std::vector<unsigned int> > >
           mpqc_to_rysq_basis(const sc::Ref<sc::GaussianBasisSet>& basis,
                               bool merge_shells);
#endif
};
TestBasis* testbs;

// profile 4, 3, and 2-center integrals
#ifdef INCLUDE_ERI
void profile_4eri(unsigned int deriv_order);
#endif
#ifdef INCLUDE_ERI3
void profile_3eri(unsigned int deriv_order);
#endif
#ifdef INCLUDE_ERI2
#endif
#ifdef INCLUDE_RYSQ
void test_rysq(unsigned int deriv_order);
void profile_rysq_4eri(unsigned int deriv_order);
void profile_rysq_3eri(unsigned int deriv_order);
#endif

/// give optional derivative order (default = 0, i.e. regular integrals)
int main(int argc, char** argv) {
  assert(argc == 1 || argc == 2);
  const unsigned int deriv_order = (argc == 2) ? atoi(argv[1]) : 0u;

  testbs = new TestBasis;

  // static initialziation of the library (one needs to happen once per process)
  LIBINT2_PREFIXED_NAME(libint2_static_init)();

#ifdef INCLUDE_RYSQ
  rysq::initialize();
#endif

  // run the tests
#ifdef INCLUDE_ERI
  profile_4eri(deriv_order);
# ifdef INCLUDE_RYSQ
   //test_rysq(deriv_order);
   profile_rysq_4eri(deriv_order);
# endif
#endif
#ifdef INCLUDE_ERI3
   profile_3eri(deriv_order);
# ifdef INCLUDE_RYSQ
   profile_rysq_3eri(deriv_order);
# endif
#endif
#ifdef INCLUDE_ERI2
#endif

  // cleanup static library data (once per process)
  LIBINT2_PREFIXED_NAME(libint2_static_cleanup)();

  return 0;
}

TestBasis::TestBasis()
{

  using namespace sc;

  // benzene, cartesian coordinates in angstroms
  double mol_geom[][4] = {
    {6.0, -1.204339000,  0.542851000, -0.047482000},
    {6.0, -1.205433000, -0.838264000,  0.124329000},
    {6.0, -0.000006000, -1.529539000,  0.208334000},
    {6.0,  1.205441000, -0.838254000,  0.124328000},
    {6.0,  1.204331000,  0.542844000, -0.047481000},
    {6.0,  0.000004000,  1.233142000, -0.133724000},
    {1.0, -2.134105000,  1.075912000, -0.125005000},
    {1.0, -2.136514000, -1.371792000,  0.187422000},
    {1.0,  0.000000000, -2.596462000,  0.339326000},
    {1.0,  2.136514000, -1.371793000,  0.187422000},
    {1.0,  2.134107000,  1.075913000, -0.125006000},
    {1.0, -0.000000000,  2.296090000, -0.286888000}
  };
  const size_t natoms = sizeof(decltype(mol_geom)) / (sizeof(double) * 4);
  const double angstrom_to_bohr = sc::Units("angstrom").to_atomic_units();

  //
  // construct a Molecule object
  //
  Ref<Molecule> mol = new Molecule;
  for(size_t a=0; a<natoms; ++a)
    mol->add_atom((int)mol_geom[a][0],
                  mol_geom[a][1] * angstrom_to_bohr,
                  mol_geom[a][2] * angstrom_to_bohr,
                  mol_geom[a][3] * angstrom_to_bohr);
  const bool use_symmetry = false;
  if (not use_symmetry) {
    Ref<PointGroup> c1_ptgrp = new PointGroup("C1");
    mol->set_point_group(c1_ptgrp);
  }
  else {
    mol->symmetrize(mol->highest_point_group(1e-4));
  }

  ExEnv::out0() << std::endl << indent << "constructed Molecule object:" << std::endl;
  mol->print(ExEnv::out0());
  ExEnv::out0() << std::endl;

  //
  // construct a GaussianBasisSet object
  //
  Ref<AssignedKeyVal> akv = new AssignedKeyVal;
  akv->assign("molecule", mol.pointer());
  //akv->assign("name", "cc-pVDZ");
  //akv->assign("name", "cc-pVTZ");
  //akv->assign("name", "Def2-SVP");
  akv->assign("name", "Def2-TZVPP");
  basis = new GaussianBasisSet(Ref<KeyVal>(akv));
  // get rid of general constractions for simplicity
  if (basis->max_ncontraction() > 1) {
    Ref<GaussianBasisSet> split_basis = new SplitBasisSet(basis, basis->name());
    basis = split_basis;
  }
  const bool uncontract = false;
  if (uncontract) {
    Ref<GaussianBasisSet> unc_basis = new UncontractedBasisSet(basis);
    basis = unc_basis;
  }
  const int nshell = basis->nshell();

  ExEnv::out0() << std::endl << indent << "constructed GaussianBasisSet object:" << std::endl;
  basis->print(ExEnv::out0());
  ExEnv::out0() << std::endl;

  // map l -> shells
  for(unsigned int l=0; l<=basis->max_angular_momentum(); ++l) {
    l_to_shells[l] = std::vector<unsigned int>();
    for(int s=0; s<nshell; ++s) {
      if (basis->shell(s).max_angular_momentum() == l) {
        l_to_shells[l].push_back(s);
      }
    }

    std::cout<< "l_to_shells[" << l << "] = ";
    for(int k=0; k<l_to_shells[l].size(); ++k)
      std::cout << " " << l_to_shells[l][k] << " ";
    std::cout << std::endl;
  }

  //
  // construct a GaussianBasisSet object for density fitting
  //
  akv = new AssignedKeyVal;
  akv->assign("molecule", mol.pointer());
  //akv->assign("name", "cc-pVDZ-RI");
  akv->assign("name", "cc-pVTZ-RI");
  df_basis = new GaussianBasisSet(Ref<KeyVal>(akv));

  ExEnv::out0() << std::endl << indent << "constructed GaussianBasisSet object for density fitting:" << std::endl;
  df_basis->print(ExEnv::out0());
  ExEnv::out0() << std::endl;

  // map l -> shells
  for(unsigned int l=0; l<=df_basis->max_angular_momentum(); ++l) {
    df_l_to_shells[l] = std::vector<unsigned int>();
    for(int s=0; s<nshell; ++s) {
      if (df_basis->shell(s).max_angular_momentum() == l) {
        df_l_to_shells[l].push_back(s);
      }
    }

    std::cout<< "df_l_to_shells[" << l << "] = ";
    for(int k=0; k<df_l_to_shells[l].size(); ++k)
      std::cout << " " << df_l_to_shells[l][k] << " ";
    std::cout << std::endl;
  }

  // create RYSQ bases
#ifdef INCLUDE_RYSQ
  std::tie(rysq_basis, rysq_l_to_shells) = mpqc_to_rysq_basis(basis, merge_shells_on_same_atom);
  std::tie(df_rysq_basis, df_rysq_l_to_shells) = mpqc_to_rysq_basis(df_basis, merge_shells_on_same_atom);
  std::tie(rysqized_basis, std::ignore) = mpqc_to_rysq_basis(basis, false);
  std::tie(rysqized_df_basis, std::ignore) = mpqc_to_rysq_basis(df_basis, false);
#endif
}

#ifdef INCLUDE_RYSQ
std::tuple<std::vector< rysq::Shell >, std::map<unsigned int, std::vector<unsigned int> > >
TestBasis::mpqc_to_rysq_basis(const sc::Ref<sc::GaussianBasisSet>& basis,
                               bool merge_shells) {

  using namespace sc;

  std::vector< rysq::Shell > rysq_basis;
  std::map<unsigned int, std::vector<unsigned int> > rysq_l_to_shells;

  if (merge_shells) {

    const int ncenter = basis->ncenter();
    for(int c=0; c<ncenter; ++c) {

      std::array<double, 3> R;
      const double* rr = basis->molecule()->r(c);
      std::copy(rr, rr+3, R.begin());

      const unsigned int lmax = basis->max_angular_momentum();
      for(unsigned int l=0; l<=lmax; ++l) {

        const int shell_start = basis->shell_on_center(c, 0);
        const int nshell_on_center = basis->nshell_on_center(c);

        std::vector<double> exponents; // unique primitives only
        std::vector<std::vector<double> > ccoefs;// contraction coefficients

        for(int s=0; s<nshell_on_center; ++s) {

          const GaussianShell& sh = basis->shell( basis->shell_on_center(c, s) );
          const int ncontr = sh.ncontraction();
          for(int c=0; c<ncontr; ++c) {
            if (sh.am(c) == l) {
              ccoefs.push_back(std::vector<double>(ccoefs.empty() ? 0 : ccoefs.back().size(), 0.0));

              //
              for (size_t p=0; p<sh.nprimitive(); ++p) {

                const double e = sh.exponent(p);
                auto eiter = std::find(exponents.begin(), exponents.end(), e);
                size_t pp;
                if (eiter != exponents.end()) { // already have this primitive
                  pp = eiter - exponents.begin();
                }
                else { // new primitive
                  exponents.push_back(e);
                  pp = exponents.size() - 1;
                  for(auto& contr: ccoefs) {
                    contr.push_back(0.0);
                  }
                }

                ccoefs.back()[pp] = sh.coefficient_unnorm(c,p);
              }

            }
          }

        }

        // if there is a shell
        if (not ccoefs.empty()) {
          std::vector<std::pair<int, int> > prim_range;
          for(const auto& contr : ccoefs) {
            double v = 0.0;
            int p = 0;
            while(v == 0.0) {
              v = contr[p++];
            }
            const size_t first = p - 1;
            while(v != 0.0) {
              v = p < contr.size() ? contr[p] : 0.0;
              ++p;
            }
            const size_t plast = p - 1;
            prim_range.push_back(std::make_pair(first,(int)plast));
          }

          rysq_basis.push_back(rysq::Shell(false, R,
                      l, exponents,
                      ccoefs, prim_range))
          );

#if 0
          std::cout << "ryqs_shell " << rysq_basis.size() - 1
          << " l = " << l << " nprim = " << exponents.size() << std::endl;
          std::cout<< "  exponents = ";
          for(auto e: exponents)
          std::cout << " " << e << " ";
          std::cout << std::endl;
          assert(ccoefs.size() == prim_range.size());
          for(int c=0; c<ccoefs.size(); ++c) {
            std::cout<< "  contraction " << c << " : [" << prim_range[c].first << "," << prim_range[c].second << "]" << std::endl;
            std::cout << "   coefs = ";
            for(auto cf: ccoefs[c])
            std::cout << " " << cf << " ";
            std::cout << std::endl;
          }
#endif

        }
      }
    }

  }
  else { // no merging shells

    auto nshell = basis->nshell();
    rysq_basis.resize(nshell);

    for(int s=0; s<nshell; ++s) {

      const auto& sh = basis->shell(s);
      const unsigned int l = sh.max_angular_momentum();

      std::array<double, 3> R;
      const double* rr = basis->molecule()->r(basis->shell_to_center(s));
      std::copy(rr, rr+3, R.begin());

      std::vector<double> exponents;
      for (size_t p=0; p<sh.nprimitive(); ++p)
      exponents.push_back(sh.exponent(p));

      std::vector<std::vector<double> > coefs(1);
      for(size_t p=0; p<sh.nprimitive(); ++p)
      coefs[0].push_back(sh.coefficient_unnorm(0,p));

      std::vector<std::pair<int, int> > prim_range(1,std::make_pair(0,(int)sh.nprimitive()));

      rysq_basis[s] = rysq::Shell(false, R,
              l, exponents,
              coefs, prim_range);
    }

  }

  // map l -> shells
  const unsigned int lmax = basis->max_angular_momentum();
  for(unsigned int l=0; l<=lmax; ++l) {
    rysq_l_to_shells[l] = std::vector<unsigned int>();
    for(int s=0; s<rysq_basis.size(); ++s) {
      if (rysq_basis[s]->angular_number() == l) {
        rysq_l_to_shells[l].push_back(s);
      }
    }

    std::cout<< "rysq_l_to_shells[" << l << "] = ";
    for(int k=0; k<rysq_l_to_shells[l].size(); ++k)
    std::cout << " " << rysq_l_to_shells[l][k] << " ";
    std::cout << std::endl;
  }

  return std::tie(rysq_basis, rysq_l_to_shells);
}
#endif

#ifdef INCLUDE_ERI

template<typename LibintEval>
inline void prep_libint2(std::vector<LibintEval>& erievals,
                  const std::array<rysq::Shell,4>& shells,
                  int norm_flag,
                  int deriv_order = 0) {
  using rysq::Shell;
  const Shell& sh0 = shells[0];
  const Shell& sh1 = shells[1];
  const Shell& sh2 = shells[2];
  const Shell& sh3 = shells[3];
  const unsigned int am01 = sh0.angular_number() + sh1.angular_number();
  const unsigned int am23 = sh2.angular_number() + sh3.angular_number();
  const unsigned int amtot = am01 + am23 + deriv_order;

  const double* exps0 = sh0.exponents_pointer();
  const double* exps1 = sh1.exponents_pointer();
  const double* exps2 = sh2.exponents_pointer();
  const double* exps3 = sh3.exponents_pointer();

  const auto& coefs0 = sh0.contractions()[0];
  const auto& coefs1 = sh1.contractions()[0];
  const auto& coefs2 = sh2.contractions()[0];
  const auto& coefs3 = sh3.contractions()[0];

  const auto& A = sh0.position();
  const auto& B = sh1.position();
  const auto& C = sh2.position();
  const auto& D = sh3.position();

  const uint np0 = sh0.num_primitive();
  const uint np1 = sh1.num_primitive();
  const uint np2 = sh2.num_primitive();
  const uint np3 = sh3.num_primitive();
  uint p0123 = 0;
  for (uint p0 = 0; p0 < np0; p0++) {
    for (uint p1 = 0; p1 < np1; p1++) {
      for (uint p2 = 0; p2 < np2; p2++) {
        for (uint p3 = 0; p3 < np3; p3++, p0123++) {

          LibintEval* erieval = &erievals[p0123];
          erieval->veclen = 1;
#if LIBINT2_FLOP_COUNT
          erieval->nflops = erievals[0].nflops;
#endif

          {
            const uint v = 0;

            const double alpha0 = exps0[p0];
            const double alpha1 = exps1[p1];
            const double alpha2 = exps2[p2];
            const double alpha3 = exps3[p3];

            const double c0 = coefs0[p0];
            const double c1 = coefs1[p1];
            const double c2 = coefs2[p2];
            const double c3 = coefs3[p3];

            const double gammap = alpha0 + alpha1;
            const double oogammap = 1.0 / gammap;
            const double rhop = alpha0 * alpha1 * oogammap;
            const double Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
            const double Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
            const double Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
            const double AB_x = A[0] - B[0];
            const double AB_y = A[1] - B[1];
            const double AB_z = A[2] - B[2];
            const double AB2 = AB_x * AB_x + AB_y * AB_y + AB_z * AB_z;

            if (am01) {
            const double PAx = Px - A[0];
            const double PAy = Py - A[1];
            const double PAz = Pz - A[2];
            const double PBx = Px - B[0];
            const double PBy = Py - B[1];
            const double PBz = Pz - B[2];

#if LIBINT2_DEFINED(eri,PA_x)
            erieval->PA_x[v] = PAx;
#endif
#if LIBINT2_DEFINED(eri,PA_y)
            erieval->PA_y[v] = PAy;
#endif
#if LIBINT2_DEFINED(eri,PA_z)
            erieval->PA_z[v] = PAz;
#endif
#if LIBINT2_DEFINED(eri,PB_x)
            erieval->PB_x[v] = PBx;
#endif
#if LIBINT2_DEFINED(eri,PB_y)
            erieval->PB_y[v] = PBy;
#endif
#if LIBINT2_DEFINED(eri,PB_z)
            erieval->PB_z[v] = PBz;
#endif

#if LIBINT2_DEFINED(eri,AB_x)
            erieval->AB_x[v] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
            erieval->AB_y[v] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
            erieval->AB_z[v] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,BA_x)
            erieval->BA_x[v] = -AB_x;
#endif
#if LIBINT2_DEFINED(eri,BA_y)
            erieval->BA_y[v] = -AB_y;
#endif
#if LIBINT2_DEFINED(eri,BA_z)
            erieval->BA_z[v] = -AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
            erieval->oo2z[v] = 0.5*oogammap;
#endif
            }

            const double gammaq = alpha2 + alpha3;
            const double oogammaq = 1.0 / gammaq;
            const double rhoq = alpha2 * alpha3 * oogammaq;
            const double one_o_gammap_plus_gammaq = 1.0 / (gammap + gammaq);
            const double gammapq = gammap * gammaq * one_o_gammap_plus_gammaq;
            const double gammap_o_gammapgammaq = gammapq * oogammaq;
            const double gammaq_o_gammapgammaq = gammapq * oogammap;
            const double Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
            const double Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
            const double Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
            const double CD_x = C[0] - D[0];
            const double CD_y = C[1] - D[1];
            const double CD_z = C[2] - D[2];
            const double CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

            const double PQx = Px - Qx;
            const double PQy = Py - Qy;
            const double PQz = Pz - Qz;
            const double PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

            if (am23) {
            const double QCx = Qx - C[0];
            const double QCy = Qy - C[1];
            const double QCz = Qz - C[2];
            const double QDx = Qx - D[0];
            const double QDy = Qy - D[1];
            const double QDz = Qz - D[2];

#if LIBINT2_DEFINED(eri,QC_x)
            erieval->QC_x[v] = QCx;
#endif
#if LIBINT2_DEFINED(eri,QC_y)
            erieval->QC_y[v] = QCy;
#endif
#if LIBINT2_DEFINED(eri,QC_z)
            erieval->QC_z[v] = QCz;
#endif
#if LIBINT2_DEFINED(eri,QD_x)
            erieval->QD_x[v] = QDx;
#endif
#if LIBINT2_DEFINED(eri,QD_y)
            erieval->QD_y[v] = QDy;
#endif
#if LIBINT2_DEFINED(eri,QD_z)
            erieval->QD_z[v] = QDz;
#endif

#if LIBINT2_DEFINED(eri,CD_x)
            erieval->CD_x[v] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
            erieval->CD_y[v] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
            erieval->CD_z[v] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,DC_x)
            erieval->DC_x[v] = -CD_x;
#endif
#if LIBINT2_DEFINED(eri,DC_y)
            erieval->DC_y[v] = -CD_y;
#endif
#if LIBINT2_DEFINED(eri,DC_z)
            erieval->DC_z[v] = -CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
            erieval->oo2e[v] = 0.5*oogammaq;
#endif
            }

            if (amtot) {
    // Prefactors for interelectron transfer relation
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
            erieval->TwoPRepITR_pfac0_0_0_x[v] = - (alpha1*AB_x + alpha3*CD_x)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
            erieval->TwoPRepITR_pfac0_0_0_y[v] = - (alpha1*AB_y + alpha3*CD_y)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
            erieval->TwoPRepITR_pfac0_0_0_z[v] = - (alpha1*AB_z + alpha3*CD_z)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
            erieval->TwoPRepITR_pfac0_1_0_x[v] = - (alpha1*AB_x + alpha3*CD_x)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
            erieval->TwoPRepITR_pfac0_1_0_y[v] = - (alpha1*AB_y + alpha3*CD_y)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
            erieval->TwoPRepITR_pfac0_1_0_z[v] = - (alpha1*AB_z + alpha3*CD_z)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
            erieval->TwoPRepITR_pfac0_0_1_x[v] = (alpha0*AB_x + alpha2*CD_x)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
            erieval->TwoPRepITR_pfac0_0_1_y[v] = (alpha0*AB_y + alpha2*CD_y)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
            erieval->TwoPRepITR_pfac0_0_1_z[v] = (alpha0*AB_z + alpha2*CD_z)*oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
            erieval->TwoPRepITR_pfac0_1_1_x[v] = (alpha0*AB_x + alpha2*CD_x)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
            erieval->TwoPRepITR_pfac0_1_1_y[v] = (alpha0*AB_y + alpha2*CD_y)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
            erieval->TwoPRepITR_pfac0_1_1_z[v] = (alpha0*AB_z + alpha2*CD_z)*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,eoz)
            erieval->eoz[v] = gammaq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,zoe)
            erieval->zoe[v] = gammap*oogammaq;
#endif
            }

            if (am01 > 1 || am23 > 1) {
            const double Wx = (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
            const double Wy = (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
            const double Wz = (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri,WP_x)
            erieval->WP_x[v] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
            erieval->WP_y[v] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
            erieval->WP_z[v] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
            erieval->WQ_x[v] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
            erieval->WQ_y[v] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
            erieval->WQ_z[v] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
            erieval->oo2ze[v] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
            erieval->roz[v] = gammapq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
            erieval->roe[v] = gammapq*oogammaq;
#endif
            }

            if (deriv_order > 0) {
            // prefactors for derivative ERI relations
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
            erieval->alpha1_rho_over_zeta2[v] = alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
            erieval->alpha2_rho_over_zeta2[v] = alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
            erieval->alpha3_rho_over_eta2[v] = alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
            erieval->alpha4_rho_over_eta2[v] = alpha3 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
            erieval->alpha1_over_zetapluseta[v] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
            erieval->alpha2_over_zetapluseta[v] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
            erieval->alpha3_over_zetapluseta[v] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
            erieval->alpha4_over_zetapluseta[v] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
            erieval->rho12_over_alpha1[v] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
            erieval->rho12_over_alpha2[v] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
            erieval->rho34_over_alpha3[v] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
            erieval->rho34_over_alpha4[v] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
            erieval->two_alpha0_bra[v] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
            erieval->two_alpha0_ket[v] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_bra)
            erieval->two_alpha1_bra[v] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_ket)
            erieval->two_alpha1_ket[v] = 2.0 * alpha3;
#endif
            }

            const double K1 = exp(- rhop * AB2) * oogammap;
            const double K2 = exp(- rhoq * CD2) * oogammaq;
            const double two_times_M_PI_to_25 = 34.986836655249725693;
            double pfac = two_times_M_PI_to_25 * K1 * K2 * sqrt(one_o_gammap_plus_gammaq);
            pfac *= c0 * c1 * c2 * c3;

            // veclen=1, hence use the erieval directly
            {
              fmeval_chebyshev.eval(erieval->LIBINT_T_SS_EREP_SS(0),PQ2*gammapq,amtot);
              //fmeval_taylor.eval(erieval->LIBINT_T_SS_EREP_SS(0),PQ2*gammapq,amtot);
              LIBINT2_REALTYPE* ssss_ptr = erieval->LIBINT_T_SS_EREP_SS(0);
              for(int l=0; l<=amtot; ++l, ++ssss_ptr)
                *ssss_ptr *= pfac;
            }

          }

        }
      }
    }
  } // end of primitive loops
}


void profile_4eri(unsigned int deriv_order) {

  typedef unsigned int uint;

  if (deriv_order > INCLUDE_ERI) return;

  // record start wall time
  const auto start = boost::chrono::high_resolution_clock::now();

  const uint veclen = LIBINT2_MAX_VECLEN;
  assert(veclen == 1);
  auto bs = testbs->basis;
  const auto& bbs = testbs->rysqized_basis;
  const uint max_contrdepth = bs->max_nprimitive_in_shell();
  const uint max_contrdepth4 = max_contrdepth * max_contrdepth * max_contrdepth * max_contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = min(LIBINT2_MAX_AM_ERI, bs->max_angular_momentum());
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) lmax = min(LIBINT2_MAX_AM_ERI1, bs->max_angular_momentum());;
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) lmax = min(LIBINT2_MAX_AM_ERI2, bs->max_angular_momentum());;
#endif
  std::vector<Libint_t> inteval(max_contrdepth4);
  if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_init_eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_init_eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_init_eri2)(&inteval[0], lmax, 0);
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of integrals
          // for example, Libint does not compute (ss|ss) integrals (although it does compute derivatives of (ss|ss)
          // another reason is a given integral type is not unique and can be computed using other functions in Libint
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_eri2)[l0][l1][l2][l3] == 0)
            continue;
#endif

          CartesianDerivIterator<4> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();
          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;

          {
          CGShell sh0(l0);
          CGShell sh1(l1);
          CGShell sh2(l2);
          CGShell sh3(l3);

          cout << "Timing "
               << " (" << sh0.label() << sh1.label() << "|"
               << sh2.label() << sh3.label() << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          }

          const auto start = boost::chrono::high_resolution_clock::now();

          size_t nshellsets = 0;
          size_t total_contrdepth4 = 0;

          // loop over shell quartets of this type
          for (auto s0 = testbs->l_to_shells[l0].cbegin();
               s0 != testbs->l_to_shells[l0].cend();
               ++s0) {
            const auto& sh0 = bbs[*s0];
            for (auto s1 = testbs->l_to_shells[l1].cbegin();
                 s1 != testbs->l_to_shells[l1].cend();
                 ++s1) {
              const auto& sh1 = bbs[*s1];
              for (auto s2 = testbs->l_to_shells[l2].cbegin();
                   s2 != testbs->l_to_shells[l2].cend();
                   ++s2) {
                const auto& sh2 = bbs[*s2];
                for (auto s3 = testbs->l_to_shells[l3].cbegin();
                     s3 != testbs->l_to_shells[l3].cend();
                     ++s3,
                     ++nshellsets) {

                  const auto& sh3 = bbs[*s3];

                  auto contrdepth4 = sh0->num_primitive() *
                      sh1->num_primitive() *
                      sh2->num_primitive() *
                      sh3->num_primitive();
                  total_contrdepth4 += contrdepth4;

                  std::array<std::shared_ptr<const Shell>, 4> sh0123{{sh0, sh1, sh2, sh3}};

                  // this prepares the data
                  prep_libint2(inteval,
                               sh0123,
                               0, deriv_order);

                  //  now use Libint to compute
        #if LIBINT_CONTRACTED_INTS
                  inteval[0].contrdepth = contrdepth4;
        #endif
                  if (deriv_order == 0)
                    LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]](
                        &inteval[0]);
        #if INCLUDE_ERI >= 1
                  else if (deriv_order == 1)
                    LIBINT2_PREFIXED_NAME(libint2_build_eri1)[am[0]][am[1]][am[2]][am[3]](
                        &inteval[0]);
        #endif
        #if INCLUDE_ERI >= 2
                  else if (deriv_order == 2)
                    LIBINT2_PREFIXED_NAME(libint2_build_eri2)[am[0]][am[1]][am[2]][am[3]](
                        &inteval[0]);
        #endif


                }
              }
            }
          }

          cout << " average contrdepth = " << (double)total_contrdepth4 / nshellsets
               << " #(shell sets) = " << nshellsets;
          cout << ": ";

          const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
          std::cout << "wall time = " << elapsed << " seconds" << std::endl;

        }
      }
    }
  }

  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif

  // record end wall time, compute total wall time spent here
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
  std::cout << "wall time = " << elapsed << " seconds" << std::endl;
}

#endif // INCLUDE_ERI

#ifdef INCLUDE_ERI3


void profile_3eri(unsigned int deriv_order) {

  typedef unsigned int uint;

  if (deriv_order > INCLUDE_ERI) return;

  const std::shared_ptr<const Shell> shS(new Shell(false));

  // record start wall time
  const auto start = boost::chrono::high_resolution_clock::now();

  const uint veclen = LIBINT2_MAX_VECLEN;
  assert(veclen == 1);

  auto bs = testbs->basis;
  const auto& bbs = testbs->rysqized_basis;
  const uint max_contrdepth = bs->max_nprimitive_in_shell();

  auto dfbs = testbs->df_basis;
  const auto& bdfbs = testbs->rysqized_df_basis;
  const uint max_contrdepth_df = dfbs->max_nprimitive_in_shell();

  const uint max_contrdepth3 = max_contrdepth_df * max_contrdepth * max_contrdepth;

  unsigned int lmax = max(bs->max_angular_momentum(), dfbs->max_angular_momentum());
  if (deriv_order == 0) {
    assert(LIBINT2_MAX_AM_3ERI >= bs->max_angular_momentum());
    assert(LIBINT2_MAX_AM_3ERI >= dfbs->max_angular_momentum());
  }
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) {
    assert(LIBINT2_MAX_AM_3ERI1 >= bs->max_angular_momentum());
    assert(LIBINT2_MAX_AM_3ERI1 >= dfbs->max_angular_momentum());
  }
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) {
    assert(LIBINT2_MAX_AM_3ERI2 >= bs->max_angular_momentum());
    assert(LIBINT2_MAX_AM_3ERI2 >= dfbs->max_angular_momentum());
  }
#endif
  std::vector<Libint_t> inteval(max_contrdepth3);
  if (deriv_order == 0) LIBINT2_PREFIXED_NAME(libint2_init_3eri)(&inteval[0], lmax, 0);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1) LIBINT2_PREFIXED_NAME(libint2_init_3eri1)(&inteval[0], lmax, 0);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2) LIBINT2_PREFIXED_NAME(libint2_init_3eri2)(&inteval[0], lmax, 0);
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    if (testbs->df_l_to_shells[l0].empty()) continue;
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      if (testbs->l_to_shells[l1].empty()) continue;
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        if (testbs->l_to_shells[l2].empty()) continue;

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of integrals
          // for example, Libint does not compute (ss|ss) integrals (although it does compute derivatives of (ss|ss)
          // another reason is a given integral type is not unique and can be computed using other functions in Libint
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_3eri)[l0][l1][l2] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[l0][l1][l2] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[l0][l1][l2] == 0)
            continue;
#endif

          CartesianDerivIterator<3> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();
          unsigned int am[3];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;

          {
          CGShell sh0(l0);
          CGShell sh1(l1);
          CGShell sh2(l2);

          cout << "Timing "
               << " (" << sh0.label() << "|" << sh1.label()
               << sh2.label() << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          }

          const auto start = boost::chrono::high_resolution_clock::now();

          size_t nshellsets = 0;
          size_t total_contrdepth = 0;

          // loop over shell quartets of this type
          for (auto s0 = testbs->df_l_to_shells[l0].cbegin();
               s0 != testbs->df_l_to_shells[l0].cend();
               ++s0) {
            const auto& sh0 = bdfbs[*s0];
            for (auto s1 = testbs->l_to_shells[l1].cbegin();
                 s1 != testbs->l_to_shells[l1].cend();
                 ++s1) {
              const auto& sh1 = bbs[*s1];
              for (auto s2 = testbs->l_to_shells[l2].cbegin();
                   s2 != testbs->l_to_shells[l2].cend();
                   ++s2,
                   ++nshellsets) {
                const auto& sh2 = bbs[*s2];

                  auto contrdepth = sh0->num_primitive() *
                      sh1->num_primitive() *
                      sh2->num_primitive();
                  total_contrdepth += contrdepth;

                  std::array<std::shared_ptr<const Shell>, 4> shells{{sh0, shS, sh1, sh2}};

                  // this prepares the data
                  prep_libint2(inteval,
                               shells,
                               0, deriv_order);

                  //  now use Libint to compute
        #if LIBINT_CONTRACTED_INTS
                  inteval[0].contrdepth = contrdepth;
        #endif
                  if (deriv_order == 0)
                    LIBINT2_PREFIXED_NAME(libint2_build_3eri)[am[0]][am[1]][am[2]](
                        &inteval[0]);
        #if INCLUDE_ERI >= 1
                  else if (deriv_order == 1)
                    LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[am[0]][am[1]][am[2]](
                        &inteval[0]);
        #endif
        #if INCLUDE_ERI >= 2
                  else if (deriv_order == 2)
                    LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[am[0]][am[1]][am[2]](
                        &inteval[0]);
        #endif


              }
            }
          }

          cout << " average contrdepth = " << (double)total_contrdepth / nshellsets
               << " #(shell sets) = " << nshellsets;
          cout << ": ";

          const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
          std::cout << "wall time = " << elapsed << " seconds" << std::endl;

      }
    }
  }

  if (deriv_order == 0)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#if INCLUDE_ERI >= 1
  if (deriv_order == 1)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif
#if INCLUDE_ERI >= 2
  if (deriv_order == 2)
    LIBINT2_PREFIXED_NAME(libint2_cleanup_eri)(&inteval[0]);
#endif

  // record end wall time, compute total wall time spent here
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
  std::cout << "wall time = " << elapsed << " seconds" << std::endl;
}

#endif // INCLUDE_ERI3

#ifdef INCLUDE_ERI2
#endif // INCLUDE_ERI2


#ifdef INCLUDE_RYSQ
void test_rysq(unsigned int deriv_order) {

  if (deriv_order > INCLUDE_RYSQ) return;

  // record start wall time
  struct timeval tod;
  gettimeofday(&tod,0);
  const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

  typedef unsigned int uint;

  const uint veclen = 1;
  const uint max_contrdepth = 3;
  const uint max_contrdepth4 = max_contrdepth * max_contrdepth * max_contrdepth * max_contrdepth;

  unsigned int lmax;
  if (deriv_order == 0) lmax = LIBINT2_MAX_AM_ERI;
#if INCLUDE_RYSQ >= 1
  if (deriv_order == 1) lmax = LIBINT2_MAX_AM_ERI1;
#endif
#if INCLUDE_RYSQ >= 2
  if (deriv_order == 2) lmax = LIBINT2_MAX_AM_ERI2;
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {

          // record start wall time
          struct timeval tod;
          gettimeofday(&tod,0);
          const double start_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;

#if LIBINT_CONTRACTED_INTS
          const uint contrdepth = do_timing_only ? std::min((4*lmax+4) / (l0+l1+l2+l3+4), max_contrdepth) : max_contrdepth;
#else
          const uint contrdepth = 1;
#endif
          const uint contrdepth4 = contrdepth * contrdepth * contrdepth * contrdepth;

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of integrals
          // for example, Libint does not compute (ss|ss) integrals (although it does compute derivatives of (ss|ss)
          // another reason is a given integral type is not unique and can be computed using other functions in Libint
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_eri2)[l0][l1][l2][l3] == 0)
            continue;
#endif

          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;
          RandomShellSet<4u> rsqset(am, veclen, contrdepth);

          CartesianDerivIterator<4> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();

          CGShell sh0(am[0]);
          CGShell sh1(am[1]);
          CGShell sh2(am[2]);
          CGShell sh3(am[3]);

          const double* A = &(rsqset.R[0][0]);
          const double* B = &(rsqset.R[1][0]);
          const double* C = &(rsqset.R[2][0]);
          const double* D = &(rsqset.R[3][0]);
          LIBINT2_REF_REALTYPE Aref[3]; for(int i=0; i<3; ++i) Aref[i] = A[i];
          LIBINT2_REF_REALTYPE Bref[3]; for(int i=0; i<3; ++i) Bref[i] = B[i];
          LIBINT2_REF_REALTYPE Cref[3]; for(int i=0; i<3; ++i) Cref[i] = C[i];
          LIBINT2_REF_REALTYPE Dref[3]; for(int i=0; i<3; ++i) Dref[i] = D[i];

          const int nrepeats = do_timing_only ? 200*(lmax-l0+1)*(lmax-l1+1)*(lmax-l2+1)*(lmax-l3+1) : 1;

          cout << (do_timing_only ? "Timing " : "Testing ")
               << " (" << sh0.label() << sh1.label() << "|"
               << sh2.label() << sh3.label() << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          if (do_timing_only) {
            cout << " contrdepth = " << contrdepth
                 << " #(repeats) = " << nrepeats;
          }
          cout << ": ";

          std::array<double, 3> R0; std::copy(A, A+3, R0.begin());
          std::array<double, 3> R1; std::copy(B, B+3, R1.begin());
          std::array<double, 3> R2; std::copy(C, C+3, R2.begin());
          std::array<double, 3> R3; std::copy(D, D+3, R3.begin());
          std::vector<std::pair<int, int> > prim_range(1,std::make_pair(0,(int)contrdepth));
          using rysq::Shell;
          Shell s0(false, R0,
                   am[0], rsqset.exp[0][0],
                   rsqset.coef[0], prim_range);
          std::shared_ptr<Shell> s1(new Shell(false, R1,
                                              am[1], rsqset.exp[1][0],
                                              rsqset.coef[1], prim_range));
          std::shared_ptr<Shell> s2(new Shell(false, R2,
                                              am[2], rsqset.exp[2][0],
                                              rsqset.coef[2], prim_range));
          std::shared_ptr<Shell> s3(new Shell(false, R3,
                                              am[3], rsqset.exp[3][0],
                                              rsqset.coef[3], prim_range));
          std::array<std::shared_ptr<const Shell>, 4> s0123{{s0, s1, s2, s3}};

          for(int k=0; k<nrepeats; ++k) {

            //  now use bagel to compute
            ERIBatch batch(s0123, 0.0);
            //SlaterBatch batch(s0123, 0.0, 1.0, false);
            batch.compute();

            // print out integrals
            if (not do_timing_only) {
              const double* buffer = batch.data();

              const int nbf0 = ((l0+1)*(l0+2)/2);
              const int nbf1 = ((l1+1)*(l1+2)/2);
              const int nbf2 = ((l2+1)*(l2+2)/2);
              const int nbf3 = ((l3+1)*(l3+2)/2);
              for(int i0=0; i0<nbf0; ++i0) {
                for(int i1=0; i1<nbf1; ++i1) {
                  for(int i2=0; i2<nbf2; ++i2) {
                    for(int i3=0; i3<nbf3; ++i3) {

                      int i = i0;
                      int j = i1;
                      int k = i2;
                      int l = i3;
                      int ni = nbf0;
                      int nj = nbf1;
                      int nk = nbf2;
                      int nl = nbf3;

                      if (batch.swap01()) {
                        std::swap(i,j);
                        std::swap(ni,nj);
                      }
                      if (batch.swap23()) {
                        std::swap(k,l);
                        std::swap(nk,nl);
                      }
                      if (batch.swap0123()) {
                        std::swap(i,k);
                        std::swap(j,l);
                        std::swap(ni,nk);
                        std::swap(nj,nl);
                      }

                      const int lkji = ((l*nk+k)*nj+j)*ni+i;
                      std::cout << i << " " << j << " " << k << " " << l << " "
                                << std::setprecision(15) << buffer[lkji] << std::endl;
                    }
                  }
                }
              }
            }

          } // end of nrepeats

          if (do_timing_only) {
            // record end wall time, compute total wall time spent here
            gettimeofday(&tod,0);
            const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
            std::cout << "wall time = " << (end_wall_time - start_wall_time) << " seconds" << std::endl;
          }

        }
      }
    }
  }

  // record end wall time, compute total wall time spent here
  gettimeofday(&tod,0);
  const double end_wall_time = tod.tv_sec + 0.000001 * tod.tv_usec;
  std::cout << "wall time = " << (end_wall_time - start_wall_time) << " seconds" << std::endl;
}

void profile_bagel_4eri(unsigned int deriv_order) {

  typedef unsigned int uint;

  if (deriv_order > INCLUDE_RYSQ) return;

  resources__ = new Resources(1);

  // record start wall time
  const auto start = boost::chrono::high_resolution_clock::now();

  const uint veclen = 1;
  const auto& bs = testbs->rysq_basis;

  unsigned int lmax;
  if (deriv_order == 0) lmax = min(LIBINT2_MAX_AM_ERI, testbs->basis->max_angular_momentum());
#if INCLUDE_RYSQ >= 1
  if (deriv_order == 1) lmax = min(LIBINT2_MAX_AM_ERI1, testbs->basis->max_angular_momentum());;
#endif
#if INCLUDE_RYSQ >= 2
  if (deriv_order == 2) lmax = min(LIBINT2_MAX_AM_ERI2, testbs->basis->max_angular_momentum());;
#endif

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        for (unsigned int l3 = 0; l3 <= lmax; ++l3) {

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of integrals
          // for example, Libint does not compute (ss|ss) integrals (although it does compute derivatives of (ss|ss)
          // another reason is a given integral type is not unique and can be computed using other functions in Libint
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_eri)[l0][l1][l2][l3] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_eri1)[l0][l1][l2][l3] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_eri2)[l0][l1][l2][l3] == 0)
            continue;
#endif
          //if (not (l0 >= l1 && l2 >= l3 && l0+l1<=l2+l3))
          //  continue;
          //if (not (l0 <= l1 && l2 <= l3 && l0+l1>=l2+l3))
          //  continue;

          CartesianDerivIterator<4> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();
          unsigned int am[4];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;
          am[3] = l3;

          {
          CGShell sh0(l0);
          CGShell sh1(l1);
          CGShell sh2(l2);
          CGShell sh3(l3);

          cout << "Timing "
               << " (" << sh0.label() << sh1.label() << "|"
               << sh2.label() << sh3.label() << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          }

          const auto start = boost::chrono::high_resolution_clock::now();

          size_t nshellsets = 0;
          size_t total_contrdepth4 = 0;

          // loop over shell quartets of this type
          for (auto s0 = testbs->rysq_l_to_shells[l0].cbegin();
               s0 != testbs->rysq_l_to_shells[l0].cend();
               ++s0) {
            const auto& sh0 = bs[*s0];
            for (auto s1 = testbs->rysq_l_to_shells[l1].cbegin();
                 s1 != testbs->rysq_l_to_shells[l1].cend();
                 ++s1) {
              const auto& sh1 = bs[*s1];
              for (auto s2 = testbs->rysq_l_to_shells[l2].cbegin();
                   s2 != testbs->rysq_l_to_shells[l2].cend();
                   ++s2) {
                const auto& sh2 = bs[*s2];
                for (auto s3 = testbs->rysq_l_to_shells[l3].cbegin();
                     s3 != testbs->rysq_l_to_shells[l3].cend();
                     ++s3,
                     ++nshellsets) {

                  const auto& sh3 = bs[*s3];

                  auto contrdepth4 = sh0->num_primitive() *
                      sh1->num_primitive() *
                      sh2->num_primitive() *
                      sh3->num_primitive();
                  total_contrdepth4 += contrdepth4;

                  std::array<std::shared_ptr<const Shell>, 4> sh0123{{sh0, sh1, sh2, sh3}};

                  ERIBatch batch(sh0123, 0.0);
                  //SlaterBatch batch(sh0123, 0.0, 1.0, false);
                  batch.compute();

                }
              }
            }
          }

          cout << " average contrdepth = " << (double)total_contrdepth4 / nshellsets
               << " #(shell sets) = " << nshellsets;
          cout << ": ";

          const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
          std::cout << "wall time = " << elapsed << " seconds" << std::endl;

        }
      }
    }
  }

  // record end wall time, compute total wall time spent here
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
  std::cout << "wall time = " << elapsed << " seconds" << std::endl;
}

void profile_bagel_3eri(unsigned int deriv_order) {

  typedef unsigned int uint;

  if (deriv_order > INCLUDE_RYSQ) return;

  const std::shared_ptr<const Shell> shS(new Shell(false));;

  resources__ = new Resources(1);

  // record start wall time
  const auto start = boost::chrono::high_resolution_clock::now();

  const uint veclen = 1;
  auto bs = testbs->rysq_basis;
  auto dfbs = testbs->df_rysq_basis;

  unsigned int lmax = max(testbs->basis->max_angular_momentum(), testbs->df_basis->max_angular_momentum());

  for (unsigned int l0 = 0; l0 <= lmax; ++l0) {
    if (testbs->df_rysq_l_to_shells[l0].empty()) continue;
    for (unsigned int l1 = 0; l1 <= lmax; ++l1) {
      if (testbs->rysq_l_to_shells[l1].empty()) continue;
      for (unsigned int l2 = 0; l2 <= lmax; ++l2) {
        if (testbs->rysq_l_to_shells[l2].empty()) continue;

          // can compute this? skip, if not.
          // there are many reason why Libint could not compute a type of integrals
          // for example, Libint does not compute (ss|ss) integrals (although it does compute derivatives of (ss|ss)
          // another reason is a given integral type is not unique and can be computed using other functions in Libint
          if (deriv_order == 0 && LIBINT2_PREFIXED_NAME(libint2_build_3eri)[l0][l1][l2] == 0)
            continue;
#if INCLUDE_ERI >= 1
          if (deriv_order == 1 && LIBINT2_PREFIXED_NAME(libint2_build_3eri1)[l0][l1][l2] == 0)
            continue;
#endif
#if INCLUDE_ERI >= 2
          if (deriv_order == 2 && LIBINT2_PREFIXED_NAME(libint2_build_3eri2)[l0][l1][l2] == 0)
            continue;
#endif
          //if (not (l0 >= l1 && l2 >= l3 && l0+l1<=l2+l3))
          //  continue;
          //if (not (l0 <= l1 && l2 <= l3 && l0+l1>=l2+l3))
          //  continue;

          CartesianDerivIterator<3> diter(deriv_order);
          const unsigned int nderiv = diter.range_size();
          unsigned int am[3];
          am[0] = l0;
          am[1] = l1;
          am[2] = l2;

          {
          CGShell sh0(l0);
          CGShell sh1(l1);
          CGShell sh2(l2);

          cout << "Timing "
               << " (" << sh0.label() << "|" << sh1.label()
               << sh2.label() << ") ";
          if (deriv_order > 0) {
            cout << " deriv order = " << deriv_order;
          }
          }

          const auto start = boost::chrono::high_resolution_clock::now();

          size_t nshellsets = 0;
          size_t total_contrdepth = 0;

          // loop over shell quartets of this type
          for (auto s0 = testbs->df_rysq_l_to_shells[l0].cbegin();
               s0 != testbs->df_rysq_l_to_shells[l0].cend();
               ++s0) {
            const auto& sh0 = dfbs[*s0];
            for (auto s1 = testbs->rysq_l_to_shells[l1].cbegin();
                 s1 != testbs->rysq_l_to_shells[l1].cend();
                 ++s1) {
              const auto& sh1 = bs[*s1];
              for (auto s2 = testbs->rysq_l_to_shells[l2].cbegin();
                   s2 != testbs->rysq_l_to_shells[l2].cend();
                   ++s2,
                   ++nshellsets) {
                const auto& sh2 = bs[*s2];

                  auto contrdepth = sh0->num_primitive() *
                      sh1->num_primitive() *
                      sh2->num_primitive();
                  total_contrdepth += contrdepth;

                  std::array<std::shared_ptr<const Shell>, 4> shells{{shS, sh0, sh1, sh2}};

                  ERIBatch batch(shells, 0.0);
                  //SlaterBatch batch(sh0123, 0.0, 1.0, false);
                  batch.compute();

              }
            }
          }

          cout << " average contrdepth = " << (double)total_contrdepth / nshellsets
               << " #(shell sets) = " << nshellsets;
          cout << ": ";

          const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
          std::cout << "wall time = " << elapsed << " seconds" << std::endl;

      }
    }
  }

  // record end wall time, compute total wall time spent here
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
  std::cout << "wall time = " << elapsed << " seconds" << std::endl;
}

#endif // INCLUDE_RYSQ
