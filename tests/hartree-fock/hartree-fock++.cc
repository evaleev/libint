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

#if __cplusplus <= 199711L
# error "Hartree-Fock test requires C++11 support"
#endif

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.h>
#include <libint2/cxxapi.h>
#include <libint2/chemistry/elements.h>
#include <libint2/basis.h>
#include <libint2/diis.h>

#if defined(_OPENMP)
# include <omp.h>
#endif

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library

using libint2::Shell;
using libint2::Atom;
using libint2::BasisSet;

std::vector<Atom> read_geometry(const std::string& filename);
Matrix compute_soad(const std::vector<Atom>& atoms);

Matrix compute_1body_ints(const BasisSet& obs,
                          libint2::OneBodyEngine::integral_type t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());
Matrix compute_2body_fock(const BasisSet& obs,
                          const Matrix& D);
// an Fock builder that can accept densities expressed a separate basis
Matrix compute_2body_fock_general(const BasisSet& obs,
                                  const Matrix& D,
                                  const BasisSet& D_bs);


int main(int argc, char *argv[]) {

  using std::cout;
  using std::cerr;
  using std::endl;

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    std::vector<Atom> atoms = read_geometry(filename);

    // announce OpenMP
#if defined(_OPENMP)
    std::cout << "Will use OpenMP to scale up to " << omp_get_max_threads() << " threads" << std::endl;
#endif

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij*xij + yij*yij + zij*zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc << endl;

    BasisSet obs("aug-cc-pVDZ", atoms);

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::init();

    // compute one-body integrals
    auto S = compute_1body_ints(obs, libint2::OneBodyEngine::overlap);
    auto T = compute_1body_ints(obs, libint2::OneBodyEngine::kinetic);
    auto V = compute_1body_ints(obs, libint2::OneBodyEngine::nuclear, atoms);
    Matrix H = T + V;
    T.resize(0,0);
    V.resize(0,0);

    Matrix D;
    {  // use SOAD as the guess density
      auto D_minbs = compute_soad(atoms); // compute guess in minimal basis
      BasisSet minbs("sto-3g", atoms);
      if (minbs == obs)
        D = D_minbs;
      else { // if basis != minimal basis, map non-representable SOAD guess into the AO basis
             // by diagonalizing a Fock matrix
        std::cout << "projecting SOAD into AO basis" << std::endl;
        auto F = H;
        F += compute_2body_fock_general(obs, D_minbs, minbs);

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();
      }
    }

    /*** =========================== ***/
    /***          SCF loop           ***/
    /*** =========================== ***/

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rms_error = 0.0;
    auto ediff_rel = 0.0;
    auto ehf = 0.0;
    auto n2 = D.cols() * D.rows();
    libint2::DIIS<Matrix> diis;
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Last iteration's energy
      auto ehf_last = ehf;

      // build a new Fock matrix
      auto F = H;
      F += compute_2body_fock(obs, D);

      // compute HF energy with the non-extrapolated Fock matrix
      ehf = D.cwiseProduct(H+F).sum();
      ediff_rel = std::abs((ehf - ehf_last)/ehf);

      // DIIS extrapolate F
      Matrix FD_comm = F*D*S - S*D*F;
      diis.extrapolate(F,FD_comm);
      // compute SCF error
      rms_error = FD_comm.norm()/n2;

      // solve F C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        std::cout <<
        "\n\nIter         E(HF)                 D(E)/E         RMS([F,D])/nn       Time(s)\n";
      printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
             ediff_rel, rms_error, time_elapsed.count());

    } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    libint2::cleanup(); // done with libint

  } // end of try block; if any exceptions occurred, report them and exit cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

std::vector<Atom> read_geometry(const std::string& filename) {

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return libint2::read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

// computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
// in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<Atom>& atoms) {

  // compute number of atomic orbitals
  size_t nao = 0;
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) // H, He
      nao += 1;
    else if (Z <= 10) // Li - Ne
      nao += 5;
    else
      throw "SOAD with Z > 10 is not yet supported";
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0; // first AO of this atom
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) { // H, He
      D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
      ao_offset += 1;
    }
    else if (Z <= 10) {
      D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
      D(ao_offset+1, ao_offset+1) = (Z == 3) ? 1 : 2; // Li? only 1 electron in 2s, else 2 electrons
      // smear the remaining electrons in 2p orbitals
      const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4)/3 : 0;
      for(auto xyz=0; xyz!=3; ++xyz)
        D(ao_offset+2+xyz, ao_offset+2+xyz) = num_electrons_per_2p;
      ao_offset += 5;
    }
  }

  return D * 0.5; // we use densities normalized to # of electrons/2
}

Matrix compute_1body_ints(const BasisSet& obs,
                          libint2::OneBodyEngine::integral_type obtype,
                          const std::vector<Atom>& atoms)
{
  const auto n = obs.nbf();
  const auto nshells = obs.size();
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  std::vector<Matrix> results(nthreads, Matrix::Zero(n,n));

  // construct the overlap integrals engine
  std::vector<libint2::OneBodyEngine> engines(nthreads);
  engines[0] = libint2::OneBodyEngine(obtype, obs.max_nprim(), obs.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == libint2::OneBodyEngine::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engines[0].set_q(q);
  }
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {

      auto bf1 = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();

      for(auto s2=0; s2<=s1; ++s2) {

        if (s12 % nthreads != thread_id)
          continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();

        // compute shell pair; return is the pointer to the buffer
        const auto* buf = engines[thread_id].compute(obs[s1], obs[s2]);

        auto& r = results[thread_id];

        // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
        Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
        r.block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
        r.block(bf2, bf1, n2, n1) = buf_mat.transpose();

      }
    }
  } // omp parallel

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    results[0] += results[i];
  }

  return results[0];
}

Matrix compute_2body_fock(const BasisSet& obs,
                          const Matrix& D) {

  const auto n = obs.nbf();
  const auto nshells = obs.size();
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  std::vector<Matrix> G(nthreads, Matrix::Zero(n,n));

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(obs.max_nprim(), obs.max_l(), 0);
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

#ifndef _OPENMP
  std::chrono::duration<double> time_elapsed = std::chrono::duration<double>::zero();
#endif // not defined _OPENMP

  auto shell2bf = obs.shell2bf();

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s1234=0l; s1!=nshells; ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();// number of basis functions in this shell

      for(auto s2=0; s2<=s1; ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        for(auto s3=0; s3<=s1; ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for(auto s4=0; s4<=s4_max; ++s4, ++s1234) {

            if (s1234 % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#ifndef _OPENMP
            const auto tstart = std::chrono::high_resolution_clock::now();
#endif

            const auto* buf = engines[thread_id].compute(obs[s1], obs[s2], obs[s3], obs[s4]);

#ifndef _OPENMP
            const auto tstop = std::chrono::high_resolution_clock::now();
            time_elapsed += tstop - tstart;
#endif

            // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    auto& g = G[thread_id];

                    g(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                    g(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                    g(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                    g(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                    g(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                    g(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                  }
                }
              }
            }

          }
        }
      }
    }
  } // omp parallel

#ifndef _OPENMP
  std::cout << "time for integrals = " << time_elapsed.count() << std::endl;
  engines[0].print_timers();
#endif // not defined _OPENMP

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}

Matrix compute_2body_fock_general(const BasisSet& obs,
                                  const Matrix& D,
                                  const BasisSet& D_bs) {

  const auto n = obs.nbf();
  const auto n_D = D_bs.nbf();
  assert(D.cols() == D.rows() && D.cols() == n_D);

#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  std::vector<Matrix> G(nthreads, Matrix::Zero(n,n));

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(std::max(obs.max_nprim(),D_bs.max_nprim()),
                                   std::max(obs.max_l(), D_bs.max_l()), 0);
  engines[0].set_precision(1e-8);
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }
  auto shell2bf = obs.shell2bf();
  auto shell2bf_D = D_bs.shell2bf();

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    auto s1234 = 0ul;
    // loop over permutationally-unique set of shells
    for(auto s1=0; s1!=obs.size(); ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();   // number of basis functions in this shell

      for(auto s2=0; s2<=s1; ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        for(auto s3=0; s3<D_bs.size(); ++s3) {

          auto bf3_first = shell2bf_D[s3];
          auto n3 = D_bs[s3].size();

          for(auto s4=0; s4<D_bs.size(); ++s4, ++s1234) {

            if (s1234 % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf_D[s4];
            auto n4 = D_bs[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

            if (s3 >= s4) {
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s1234_deg = s12_deg * s34_deg;
              //auto s1234_deg = s12_deg;
              const auto* buf_J = engines[thread_id].compute(obs[s1], obs[s2], D_bs[s3], D_bs[s4]);

              for(auto f1=0, f1234=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for(auto f3=0; f3!=n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      auto& g = G[thread_id];

                      const auto value = buf_J[f1234];
                      const auto value_scal_by_deg = value * s1234_deg;
                      g(bf1,bf2) += 2.0 * D(bf3,bf4) * value_scal_by_deg;
                    }
                  }
                }
              }
            }

            const auto* buf_K = engines[thread_id].compute(obs[s1], D_bs[s3], obs[s2], D_bs[s4]);

            for(auto f1=0, f1324=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f2=0; f2!=n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                    const auto bf4 = f4 + bf4_first;

                    auto& g = G[thread_id];

                    const auto value = buf_K[f1324];
                    const auto value_scal_by_deg = value * s12_deg;
                    g(bf1,bf2) -= D(bf3,bf4) * value_scal_by_deg;
                  }
                }
              }
            }

          }
        }
      }
    }

  } // omp parallel

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}
