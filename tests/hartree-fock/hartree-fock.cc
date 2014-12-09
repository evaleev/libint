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

// BTAS tensor algebra library: not yet ready
//#include <btas/btas.h>

// Libint Gaussian integrals library
#include <libint2.h>
#include <libint2/cxxapi.h>
#include <libint2/chemistry/elements.h>

#if defined(_OPENMP)
# include <omp.h>
#endif

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library

struct Atom {
    int atomic_number;
    double x, y, z;
};

std::vector<Atom> read_geometry(const std::string& filename);
std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_basis(const std::vector<Atom>& atoms, std::string file_dot_g94);
size_t nbasis(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
Matrix compute_soad(const std::vector<Atom>& atoms);
Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::OneBodyEngine::type t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());

// simple-to-read, but inefficient Fock builder; computes ~16 times as many ints as possible
Matrix compute_2body_fock_simple(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);
// an efficient Fock builder; *integral-driven* hence computes permutationally-unique ints once
Matrix compute_2body_fock(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);
// an efficient Fock builder that can accept densities expressed a separate basis
Matrix compute_2body_fock_general(const std::vector<libint2::Shell>& shells,
                                  const Matrix& D,
                                  const std::vector<libint2::Shell>& shells_D);
// threaded versions of the above
#if defined(_OPENMP)
Matrix compute_1body_ints_openmp(const std::vector<libint2::Shell>& shells,
                                 libint2::OneBodyEngine::type t,
                                 const std::vector<Atom>& atoms = std::vector<Atom>());
Matrix compute_2body_fock_openmp(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D);
#endif


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
    cout << "\tNuclear repulsion energy = " << enuc << endl;

    /*** =========================== ***/
    /*** create basis set            ***/
    /*** =========================== ***/

    auto shells = make_sto3g_basis(atoms);
    //auto shells = make_augccpvdz_basis(atoms);
    size_t nao = 0;
    for (auto s=0; s<shells.size(); ++s)
      nao += shells[s].size();

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::init();

    // compute overlap integrals
#if defined(_OPENMP)
      auto S = compute_1body_ints_openmp(shells, libint2::OneBodyEngine::overlap);
      { // validate overlap integrals; this also tests 1-thread function
        auto S_1thread = compute_1body_ints(shells, libint2::OneBodyEngine::overlap);
        if ((S - S_1thread).norm() > 1e-12)
          throw "overlap sanity check: threaded result != sequential result";
      }
#else
      auto S = compute_1body_ints(shells, libint2::OneBodyEngine::overlap);
#endif
    cout << "\n\tOverlap Integrals:\n";
    cout << S << endl;

    // compute kinetic-energy integrals
#if defined(_OPENMP)
    auto T = compute_1body_ints_openmp(shells, libint2::OneBodyEngine::kinetic);
#else
    auto T = compute_1body_ints(shells, libint2::OneBodyEngine::kinetic);
#endif
    cout << "\n\tKinetic-Energy Integrals:\n";
    cout << T << endl;

    // compute nuclear-attraction integrals
#if defined(_OPENMP)
    Matrix V = compute_1body_ints_openmp(shells, libint2::OneBodyEngine::nuclear, atoms);
#else
    Matrix V = compute_1body_ints(shells, libint2::OneBodyEngine::nuclear, atoms);
#endif
    cout << "\n\tNuclear Attraction Integrals:\n";
    cout << V << endl;

    // Core Hamiltonian = T + V
    Matrix H = T + V;
    cout << "\n\tCore Hamiltonian:\n";
    cout << H << endl;

    // T and V no longer needed, free up the memory
    T.resize(0,0);
    V.resize(0,0);

    /*** =========================== ***/
    /*** build initial-guess density ***/
    /*** =========================== ***/

    const auto use_hcore_guess = false ;  // use core Hamiltonian eigenstates to guess density?
                                         // set to true to match the result of versions 0, 1, and 2 of the code
                                         // HOWEVER !!! even for medium-size molecules hcore will usually fail !!!
                                         // thus set to false to use Superposition-Of-Atomic-Densities (SOAD) guess
    Matrix D;
    if (use_hcore_guess) { // hcore guess
      // solve H C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();
      cout << "\n\tInitial C Matrix:\n";
      cout << C << endl;

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
    }
    else {  // use SOAD as the guess density
      auto D_minbs = compute_soad(atoms); // compute guess in minimal basis
      auto minbs_shells = make_sto3g_basis(atoms);
      if (minbs_shells == shells)
        D = D_minbs;
      else { // if basis != minimal basis, map non-representable SOAD guess into the AO basis
             // by diagonalizing a Fock matrix
        std::cout << "projecting SOAD into AO basis" << std::endl;
        auto F = H;
        F += compute_2body_fock_general(shells, D_minbs, minbs_shells);
        cout << "\n\tSOAD Fock Matrix:\n";
        cout << F << endl;

        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();
      }
    }

    cout << "\n\tInitial Density Matrix:\n";
    cout << D << endl;

    /*** =========================== ***/
    /*** main iterative loop         ***/
    /*** =========================== ***/

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rmsd = 0.0;
    auto ediff = 0.0;
    auto ehf = 0.0;
    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Save a copy of the energy and the density
      auto ehf_last = ehf;
      auto D_last = D;

      // build a new Fock matrix
      auto F = H;
      //F += compute_2body_fock_simple(shells, D);
#if defined(_OPENMP)
      F += compute_2body_fock_openmp(shells, D);
      if (iter == 1) {
        auto G_1thread = compute_2body_fock(shells, D);
        if ((F - H - G_1thread).norm() > 1e-12)
          throw "Fock sanity check: threaded result != sequential result";
      }
#else
      F += compute_2body_fock(shells, D);
#endif

      if (iter == 1) {
        cout << "\n\tFock Matrix:\n";
        cout << F << endl;
      }

      // solve F C = e S C
      Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
      auto eps = gen_eig_solver.eigenvalues();
      auto C = gen_eig_solver.eigenvectors();

      // compute density, D = C(occ) . C(occ)T
      auto C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();

      // compute HF energy
      ehf = 0.0;
      for (auto i = 0; i < nao; i++)
        for (auto j = 0; j < nao; j++)
          ehf += D(i,j) * (H(i,j) + F(i,j));

      // compute difference with last iteration
      ediff = ehf - ehf_last;
      rmsd = (D - D_last).norm();

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        std::cout <<
        "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n";
      printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iter, ehf, ehf + enuc,
             ediff, rmsd, time_elapsed.count());

    } while (((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iter < maxiter));

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

namespace {
  bool strcaseequal(const std::string& a, const std::string& b) {
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
                                              [](char a, char b) {return ::tolower(a) == ::tolower(b);}
                                             );
  }
}

// this reads the geometry in the standard xyz format supported by most chemistry software
std::vector<Atom> read_dotxyz(std::istream& is) {
  size_t natom;
  is >> natom;

  std::string comment;
  std::getline(is, comment);

  std::vector<Atom> atoms(natom);
  for (auto i = 0; i < natom; i++) {
    std::string element_symbol;
    double x, y, z;
    is >> element_symbol >> x >> y >> z;

    // .xyz files report element labels, hence convert to atomic numbers
    int Z = -1;
    using libint2::chemistry::element_info;
    for(const auto& e: element_info) {
      if (strcaseequal(e.symbol, element_symbol)) {
        Z = e.Z;
        break;
      }
    }
    if (Z == -1) {
      std::cerr << "read_dotxyz: element symbol \"" << element_symbol << "\" is not recognized" << std::endl;
      throw "Did not recognize element symbol in .xyz file";
    }

    atoms[i].atomic_number = Z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 1 / 0.52917721092; // 2010 CODATA value
    atoms[i].x = x * angstrom_to_bohr;
    atoms[i].y = y * angstrom_to_bohr;
    atoms[i].z = z * angstrom_to_bohr;
  }

  return atoms;
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
    return read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms) {

  std::vector<libint2::Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {71.616837000, 13.045096000, 3.530512200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
              {99.106169000, 18.052312000, 4.885660200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.780455900, 0.878496600, 0.285714400},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
          {3.780455900, 0.878496600, 0.285714400},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {130.709320000, 23.808861000, 6.443608300},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  // technical step: rescale contraction coefficients to include primitive normalization coefficients
  for(auto& s: shells) {
    s.renorm();
  }

  return shells;
}

unsigned short to_l(std::string amlabel) {
  assert(amlabel.size() == 1);
  const char amchar = ::toupper(amlabel.c_str()[0]);
  switch (amchar) {
    case 'S': return 0;
    case 'P': return 1;
    case 'D': return 2;
    case 'F': return 3;
    case 'G': return 4;
    case 'H': return 5;
    case 'I': return 6;
    case 'K': return 7;
    case 'L': return 8;
    case 'M': return 9;
    case 'N': return 10;
    case 'O': return 11;
    case 'Q': return 12;
    case 'R': return 13;
    case 'T': return 14;
    case 'U': return 15;
    case 'V': return 16;
    case 'W': return 17;
    case 'X': return 18;
    case 'Y': return 19;
    case 'Z': return 20;
    default: throw "invalid angular momentum label";
  }
}

std::vector<std::vector<libint2::Shell>> read_g94_basis_library(std::string file_dot_g94) {

  std::cout << "Will read basis set from " << file_dot_g94 << std::endl;
  std::ifstream is(file_dot_g94);
  assert(is.good());
  std::vector<std::vector<libint2::Shell>> ref_shells(118); // 118 = number of chemical elements

  std::string comment, rest;

  // skip till first basis
  while(std::getline(is, comment) && comment != "****") {
  }

  size_t Z;
  auto nextbasis = true, nextshell = false;
  while(std::getline(is, comment) && comment != "") {
    if (comment == "****") {
      nextbasis = true;
      nextshell = false;
      continue;
    }
    if (nextbasis) {
      nextbasis = false;
      std::istringstream iss(comment);
      std::string elemsymbol;
      iss >> elemsymbol >> rest;

      bool found = false;
      using libint2::chemistry::element_info;
      for(const auto& e: element_info) {
        if (strcaseequal(e.symbol, elemsymbol)) {
          Z = e.Z;
          found = true;
          break;
        }
      }
      if (not found) {
        std::ostringstream oss;
        oss << "in file " << file_dot_g94
            << " found G94 basis set for element symbol \""
            << elemsymbol << "\", not found in Periodic Table.";
        throw std::runtime_error(oss.str());
      }

      nextshell = true;
      continue;
    }
    if (nextshell) {
      std::istringstream iss(comment);
      std::string amlabel;
      unsigned nprim;
      iss >> amlabel >> nprim >> rest;
      if (amlabel != "SP" && amlabel != "sp") {
        auto l = to_l(amlabel);
        std::vector<double> exps;
        std::vector<double> coeffs;
        for(auto p = 0; p!=nprim; ++p) {
          std::getline(is, comment);
          std::istringstream iss(comment);
          double e, c;
          iss >> e >> c;
          exps.emplace_back(e);
          coeffs.emplace_back(c);
        }
        auto pure = l>1;
        ref_shells[Z].push_back(
            libint2::Shell{
              exps,
              {
                {l, pure, coeffs}
              },
              {{0.0, 0.0, 0.0}}
            }
        );
      }
      else { // split the SP shells
        std::vector<double> exps;
        std::vector<double> coeffs_s, coeffs_p;
        for(auto p = 0; p!=nprim; ++p) {
          std::getline(is, comment);
          std::istringstream iss(comment);
          double e, c1, c2;
          iss >> e >> c1 >> c2;
          exps.emplace_back(e);
          coeffs_s.emplace_back(c1);
          coeffs_p.emplace_back(c2);
        }
        ref_shells[Z].push_back(
            libint2::Shell{exps,
          {
           {0, false, coeffs_s}
          },
          {{0.0, 0.0, 0.0}}
        }
        );
        ref_shells[Z].push_back(
            libint2::Shell{ exps,
          {
           {1, false, coeffs_p}
          },
          {{0.0, 0.0, 0.0}}
        }
        );
      }
    }
  }

  return ref_shells;
}

std::vector<libint2::Shell> make_basis(const std::vector<Atom>& atoms,
                                       std::string file_dot_g94) {

  // read in the file contents
  std::vector<std::vector<libint2::Shell>> ref_shells = read_g94_basis_library(file_dot_g94);

  // for each atom find the corresponding basis
  std::vector<libint2::Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    auto Z = atoms[a].atomic_number;
    if (ref_shells[Z].empty())
      throw std::string("did not find the basis for this Z in ") + file_dot_g94;

    for(auto s: ref_shells[Z]) {
      shells.push_back(std::move(s));
      shells.back().move({{atoms[a].x, atoms[a].y, atoms[a].z}});
    }

  }

  // technical step: rescale contraction coefficients to include primitive normalization coefficients
  for(auto& s: shells) {
    s.renorm();
  }

  return shells;
}

size_t nbasis(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
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

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::OneBodyEngine::type obtype,
                          const std::vector<Atom>& atoms)
{
  const auto n = nbasis(shells);
  Matrix result(n,n);

  // construct the overlap integrals engine
  libint2::OneBodyEngine engine(obtype, max_nprim(shells), max_l(shells), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == libint2::OneBodyEngine::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engine.set_q(q);
  }

  auto shell2bf = map_shell_to_basis_function(shells);

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      const auto* buf = engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;
}

#if defined(_OPENMP)
Matrix compute_1body_ints_openmp(const std::vector<libint2::Shell>& shells,
                                 libint2::OneBodyEngine::type obtype,
                                 const std::vector<Atom>& atoms)
{
  const auto n = nbasis(shells);
  const auto nthreads = omp_get_max_threads();
  std::vector<Matrix> results(nthreads, Matrix::Zero(n,n));

  // construct the overlap integrals engine
  std::vector<libint2::OneBodyEngine> engines(nthreads);
  engines[0] = libint2::OneBodyEngine(obtype, max_nprim(shells), max_l(shells), 0);
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

  auto shell2bf = map_shell_to_basis_function(shells);

#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0l, s12=0l; s1!=shells.size(); ++s1) {

      auto bf1 = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();

      for(auto s2=0; s2<=s1; ++s2) {

        if (s12 % nthreads != thread_id)
          continue;

        auto bf2 = shell2bf[s2];
        auto n2 = shells[s2].size();

        // compute shell pair; return is the pointer to the buffer
        const auto* buf = engines[thread_id].compute(shells[s1], shells[s2]);

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
#endif // defined(_OPENMP)

Matrix compute_2body_fock_simple(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D) {

  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the electron repulsion integrals engine
  libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  // loop over shell pairs of the Fock matrix, {s1,s2}
  // Fock matrix is symmetric, but skipping it here for simplicity (see compute_2body_fock)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2!=shells.size(); ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      // loop over shell pairs of the density matrix, {s3,s4}
      // again symmetry is not used for simplicity
      for(auto s3=0; s3!=shells.size(); ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        for(auto s4=0; s4!=shells.size(); ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
          const auto* buf_1234 = engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

          // we don't have an analog of Eigen for tensors (yet ... see github.com/BTAS/BTAS, under development)
          // hence some manual labor here:
          // 1) loop over every integral in the shell set (= nested loops over basis functions in each shell)
          // and 2) add contribution from each integral
          for(auto f1=0, f1234=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f2=0; f2!=n2; ++f2) {
              const auto bf2 = f2 + bf2_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) += D(bf3,bf4) * 2.0 * buf_1234[f1234];
                }
              }
            }
          }

          // exchange contribution to the Fock matrix is from {s1,s3,s2,s4} integrals
          const auto* buf_1324 = engine.compute(shells[s1], shells[s3], shells[s2], shells[s4]);

          for(auto f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;
                  G(bf1,bf2) -= D(bf3,bf4) * buf_1324[f1324];
                }
              }
            }
          }

        }
      }
    }
  }

  return G;
}

Matrix compute_2body_fock(const std::vector<libint2::Shell>& shells,
                          const Matrix& D) {

  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  // construct the 2-electron repulsion integrals engine
  libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim(shells), max_l(shells), 0);

  auto shell2bf = map_shell_to_basis_function(shells);

  // The problem with the simple Fock builder is that permutational symmetries of the Fock,
  // density, and two-electron integrals are not taken into account to reduce the cost.
  // To make the simple Fock builder efficient we must rearrange our computation.
  // The most expensive step in Fock matrix construction is the evaluation of 2-e integrals;
  // hence we must minimize the number of computed integrals by taking advantage of their permutational
  // symmetry. Due to the multiplicative and Hermitian nature of the Coulomb kernel (and realness
  // of the Gaussians) the permutational symmetry of the 2-e ints is given by the following relations:
  //
  // (12|34) = (21|34) = (12|43) = (21|43) = (34|12) = (43|12) = (34|21) = (43|21)
  //
  // (here we use chemists' notation for the integrals, i.e in (ab|cd) a and b correspond to
  // electron 1, and c and d -- to electron 2).
  //
  // It is easy to verify that the following set of nested loops produces a permutationally-unique
  // set of integrals:
  // foreach a = 0 .. n-1
  //   foreach b = 0 .. a
  //     foreach c = 0 .. a
  //       foreach d = 0 .. (a == c ? b : c)
  //         compute (ab|cd)
  //
  // The only complication is that we must compute integrals over shells. But it's not that complicated ...
  //
  // The real trick is figuring out to which matrix elements of the Fock matrix each permutationally-unique
  // (ab|cd) contributes. STOP READING and try to figure it out yourself. (to check your answer see below)

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3<=s1; ++s3) {

        auto bf3_first = shell2bf[s3];
        auto n3 = shells[s3].size();

        const auto s4_max = (s1 == s3) ? s2 : s3;
        for(auto s4=0; s4<=s4_max; ++s4) {

          auto bf4_first = shell2bf[s4];
          auto n4 = shells[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
          auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
          auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

          const auto* buf = engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

          // ANSWER
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

                  G(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                  G(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                  G(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                  G(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                  G(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                  G(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}

#if defined(_OPENMP)
Matrix compute_2body_fock_openmp(const std::vector<libint2::Shell>& shells,
                                 const Matrix& D) {

  const auto n = nbasis(shells);
  const auto nthreads = omp_get_max_threads();
  std::vector<Matrix> G(nthreads, Matrix::Zero(n,n));

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(max_nprim(shells), max_l(shells), 0);
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = map_shell_to_basis_function(shells);

#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s1234=0l; s1!=shells.size(); ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();// number of basis functions in this shell

      for(auto s2=0; s2<=s1; ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        for(auto s3=0; s3<=s1; ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for(auto s4=0; s4<=s4_max; ++s4, ++s1234) {

            if (s1234 % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            const auto* buf = engines[thread_id].compute(shells[s1], shells[s2], shells[s3], shells[s4]);

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

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}
#endif // defined(_OPENMP)

Matrix compute_2body_fock_general(const std::vector<libint2::Shell>& shells,
                                  const Matrix& D,
                                  const std::vector<libint2::Shell>& shells_D) {

  const auto n = nbasis(shells);
  Matrix G = Matrix::Zero(n,n);

  const auto n_D = nbasis(shells_D);
  assert(D.cols() == D.rows() && D.cols() == n_D);

  // construct the 2-electron repulsion integrals engine
  libint2::TwoBodyEngine<libint2::Coulomb> engine(std::max(max_nprim(shells),max_nprim(shells_D)),
                                                  std::max(max_l(shells), max_l(shells_D)), 0);

  auto shell2bf = map_shell_to_basis_function(shells);
  auto shell2bf_D = map_shell_to_basis_function(shells_D);

  // loop over permutationally-unique set of shells
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1_first = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();   // number of basis functions in this shell

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2_first = shell2bf[s2];
      auto n2 = shells[s2].size();

      for(auto s3=0; s3<shells_D.size(); ++s3) {

        auto bf3_first = shell2bf_D[s3];
        auto n3 = shells_D[s3].size();

        for(auto s4=0; s4<shells_D.size(); ++s4) {

          auto bf4_first = shell2bf_D[s4];
          auto n4 = shells_D[s4].size();

          // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
          auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

          if (s3 >= s4) {
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s1234_deg = s12_deg * s34_deg;
            //auto s1234_deg = s12_deg;
            const auto* buf_J = engine.compute(shells[s1], shells[s2], shells_D[s3], shells_D[s4]);

            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_J[f1234];
                    const auto value_scal_by_deg = value * s1234_deg;
                    G(bf1,bf2) += 2.0 * D(bf3,bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }

          const auto* buf_K = engine.compute(shells[s1], shells_D[s3], shells[s2], shells_D[s4]);

          for(auto f1=0, f1324=0; f1!=n1; ++f1) {
            const auto bf1 = f1 + bf1_first;
            for(auto f3=0; f3!=n3; ++f3) {
              const auto bf3 = f3 + bf3_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                  const auto bf4 = f4 + bf4_first;

                  const auto value = buf_K[f1324];
                  const auto value_scal_by_deg = value * s12_deg;
                  G(bf1,bf2) -= D(bf3,bf4) * value_scal_by_deg;
                }
              }
            }
          }

        }
      }
    }
  }

  // symmetrize the result and return
  Matrix Gt = G.transpose();
  return 0.5 * (G + Gt);
}
