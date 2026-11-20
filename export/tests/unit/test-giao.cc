/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// Tests of the zero-field London-orbital (GIAO) derivative integrals
// (Operator::giao_overlap, giao_kinetic, giao_nuclear) and of the angular
// momentum integrals (Operator::angular_momentum)

#include "catch.hpp"
#include "fixture.h"

namespace {

using libint2::BasisSet;
using libint2::Engine;
using libint2::Operator;
using libint2::Shell;

/// matrices of all components of a 1-body operator, [comp][i * n + j]
using Matrices = std::vector<std::vector<double>>;

/// the 1-body operators tested here, all have 3 (x,y,z) components
const std::vector<Operator> magnetic_opers = {
    Operator::giao_overlap, Operator::giao_kinetic, Operator::giao_nuclear,
    Operator::angular_momentum};

/// point charges of the atoms
std::vector<std::pair<double, std::array<double, 3>>> charges(
    const std::vector<libint2::Atom>& atoms) {
  std::vector<std::pair<double, std::array<double, 3>>> q;
  for (const auto& a : atoms)
    q.push_back({static_cast<double>(a.atomic_number), {{a.x, a.y, a.z}}});
  return q;
}

/// computes all components of the integrals of \p engine over \p shells
Matrices compute(Engine& engine, const std::vector<Shell>& shells) {
  const auto shell2bf = BasisSet::compute_shell2bf(shells);
  const auto n = libint2::nbf(shells);
  const auto& res = engine.results();
  Matrices M;
  for (size_t s1 = 0; s1 != shells.size(); ++s1) {
    for (size_t s2 = 0; s2 != shells.size(); ++s2) {
      engine.compute(shells[s1], shells[s2]);
      if (M.empty()) M.resize(res.size(), std::vector<double>(n * n));
      const auto n1 = shells[s1].size();
      const auto n2 = shells[s2].size();
      for (size_t c = 0; c != res.size(); ++c)
        for (size_t i = 0; i != n1; ++i)
          for (size_t j = 0; j != n2; ++j)
            M[c][(shell2bf[s1] + i) * n + shell2bf[s2] + j] =
                res[c] ? res[c][i * n2 + j] : 0.0;
    }
  }
  return M;
}

/// computes the integrals of \p oper over \p shells, with \p atoms providing
/// point charges and \p origin the origin for origin-dependent operators
Matrices compute(Operator oper, const std::vector<Shell>& shells,
                 const std::vector<libint2::Atom>& atoms,
                 const std::array<double, 3>& origin = {{0, 0, 0}},
                 int deriv_order = 0) {
  Engine engine(oper, libint2::max_nprim(shells), libint2::max_l(shells),
                deriv_order);
  if (libint2::is_nuclear_type(oper)) engine.set_params(charges(atoms));
  if (oper == Operator::emultipole1 || oper == Operator::angular_momentum)
    engine.set_params(origin);
  return compute(engine, shells);
}

/// \f$ (\vec{A} - \vec{B}) \f$ for each pair of basis functions
std::vector<std::array<double, 3>> center_differences(
    const std::vector<Shell>& shells) {
  const auto n = libint2::nbf(shells);
  std::vector<std::array<double, 3>> AB(n * n);
  size_t i = 0;
  for (const auto& s1 : shells)
    for (size_t f1 = 0; f1 != s1.size(); ++f1, ++i) {
      size_t j = 0;
      for (const auto& s2 : shells)
        for (size_t f2 = 0; f2 != s2.size(); ++f2, ++j)
          for (int xyz = 0; xyz != 3; ++xyz)
            AB[i * n + j][xyz] = s1.O[xyz] - s2.O[xyz];
    }
  return AB;
}

/// requires \p A and \p B to agree elementwise to within \p tol (absolute)
/// or 1e-10 (relative); reports the worst element only
void require_equal(const Matrices& A, const Matrices& B, double tol,
                   const std::string& label) {
  REQUIRE(A.size() == B.size());
  double worst_excess = 0;  // worst |A-B| - max(tol, 1e-10 |B|)
  size_t worst_c = 0, worst_ij = 0;
  for (size_t c = 0; c != A.size(); ++c) {
    REQUIRE(A[c].size() == B[c].size());
    for (size_t ij = 0; ij != A[c].size(); ++ij) {
      const auto excess = std::abs(A[c][ij] - B[c][ij]) -
                          std::max(tol, 1e-10 * std::abs(B[c][ij]));
      if (excess > worst_excess) {
        worst_excess = excess;
        worst_c = c;
        worst_ij = ij;
      }
    }
  }
  INFO(label << ": component " << worst_c << ", element " << worst_ij << ": "
             << A[worst_c][worst_ij] << " vs " << B[worst_c][worst_ij]);
  REQUIRE(worst_excess <= 0);
}

/// off-axis 3-center system with Cartesian s, p, d shells used to compare to
/// libcint
const std::vector<libint2::Atom>& reference_atoms() {
  static const std::vector<libint2::Atom> atoms = {
      {8, 0.1, -0.2, 0.05}, {1, 1.5, 0.3, 1.1}, {1, -1.2, 1.3, -0.6}};
  return atoms;
}
const std::vector<Shell>& reference_shells() {
  static const std::vector<Shell> shells = {
      Shell{{5.0, 1.2}, {{0, false, {0.4, 0.7}}}, {{0.1, -0.2, 0.05}}},
      Shell{{3.0, 0.8}, {{1, false, {0.5, 0.6}}}, {{0.1, -0.2, 0.05}}},
      Shell{{1.1}, {{2, false, {1.0}}}, {{0.1, -0.2, 0.05}}},
      Shell{{1.3, 0.3}, {{0, false, {0.5, 0.6}}}, {{1.5, 0.3, 1.1}}},
      Shell{{0.9}, {{1, false, {1.0}}}, {{1.5, 0.3, 1.1}}},
      Shell{{0.7}, {{0, false, {1.0}}}, {{-1.2, 1.3, -0.6}}},
      Shell{{0.6}, {{2, false, {1.0}}}, {{-1.2, 1.3, -0.6}}}};
  return shells;
}

/// computes all components of the 2-body integrals of \p engine over all
/// quartets of \p shells, [comp][((i * n + j) * n + k) * n + l]
Matrices compute4(Engine& engine, const std::vector<Shell>& shells) {
  const auto shell2bf = BasisSet::compute_shell2bf(shells);
  const auto n = libint2::nbf(shells);
  const auto& res = engine.results();
  Matrices M(res.size(), std::vector<double>(n * n * n * n, 0.0));
  for (size_t s1 = 0; s1 != shells.size(); ++s1)
    for (size_t s2 = 0; s2 != shells.size(); ++s2)
      for (size_t s3 = 0; s3 != shells.size(); ++s3)
        for (size_t s4 = 0; s4 != shells.size(); ++s4) {
          engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
          if (res[0] == nullptr) continue;  // screened out
          const auto n1 = shells[s1].size(), n2 = shells[s2].size(),
                     n3 = shells[s3].size(), n4 = shells[s4].size();
          for (size_t c = 0; c != res.size(); ++c)
            for (size_t i = 0, ijkl = 0; i != n1; ++i)
              for (size_t j = 0; j != n2; ++j)
                for (size_t k = 0; k != n3; ++k)
                  for (size_t l = 0; l != n4; ++l, ++ijkl)
                    M[c][((((shell2bf[s1] + i) * n + shell2bf[s2] + j) * n +
                           shell2bf[s3] + k) *
                              n +
                          shell2bf[s4] + l)] = res[c][ijkl];
        }
  return M;
}

}  // namespace

#if defined(LIBINT2_SUPPORT_ONEBODY)

TEST_CASE("GIAO and angular momentum integrals match libcint",
          "[engine][1-body][giao]") {
  const auto& atoms = reference_atoms();
  const auto& shells = reference_shells();
  const std::array<double, 3> gauge_origin = {{0.3, -0.4, 0.7}};

  // reference values computed with libcint 6.2 (int1e_igovlp, int1e_igkin,
  // int1e_ignuc, int1e_giao_irjxp, int1e_cg_irxp), converted to the libint
  // conventions and normalization:
  //   giao_overlap     = -int1e_igovlp
  //   giao_kinetic     = -int1e_igkin - int1e_giao_irjxp / 2
  //   giao_nuclear     = -int1e_ignuc
  //   angular_momentum = int1e_cg_irxp (common origin = gauge_origin)
  // these are also the conventions of the PySCF NMR module, where
  // s10 = giao_overlap and h10 = giao_kinetic + giao_nuclear
  struct Reference {
    Operator oper;
    size_t s1, s2;                     // shell indices
    std::vector<double> values;        // [comp][i][j]
  };
  const std::vector<Reference> refs = {
      {Operator::giao_overlap,
       1,
       6,
       {
        3.021122020186132e-04, 3.988397005132536e-03, 6.384317859364704e-03,
        -4.863004225484325e-02, -3.633088560636573e-02, 3.827905676833900e-02,
        -3.868701241705084e-02, 1.257695951545168e-02, -6.062975557788412e-02,
        -3.684583247072864e-02, 2.679126365995700e-02, -6.676810577165847e-02,
        -1.083306560041157e-01, 9.307138051414647e-02, -1.266257209914329e-02,
        -1.531225442873199e-01, 4.160131753503353e-03, -9.871540046441271e-02,
        3.971862939387043e-02, -1.481206338879153e-02, 1.233796783458727e-02,
        5.581659336928924e-02, -5.898384660701712e-02, 5.836058168606938e-02,
        4.126225193009704e-02, -8.884539690966840e-03, -3.094668894757280e-02,
        2.956957069480363e-19, 1.776907938193368e-02, -4.126225193009705e-02,
        -1.091563820230047e-01, 8.120194169020445e-02, -2.386593242726590e-03,
        -1.116331867385785e-01, 1.481206338879154e-02, -8.700204013687485e-02,
        9.105415112027916e-02, -4.215847875363014e-02, 1.570359774570276e-02,
        2.260676076695847e-01, -6.345479788038494e-02, 5.812015189271289e-02,
        1.725946062112487e-01, -4.565670293313453e-02, 4.984407512290794e-02,
        7.369166494145728e-02, -1.257695951545168e-02, 3.831563016616990e-02,
        -3.523803112177950e-02, 1.246335179871137e-03, 1.981762133045598e-02,
        4.863004225484325e-02, 2.586142123635838e-02, -3.343137848578098e-03}},
      {Operator::giao_overlap,
       2,
       3,
       {
        -4.292091392581461e-02, 3.711978227866305e-02, -3.345367270487419e-02,
        2.928779037285300e-03, 2.048538625414737e-02, -6.630723250126254e-02,
        -9.288081389117771e-02, -1.635648860988657e-02, 2.331026548526158e-02,
        -4.164242342799797e-03, 1.990835428255348e-02, 8.121323555890211e-02,
        1.014568442302184e-01, -4.170423893827141e-02, 3.350477051827912e-02,
        -1.922066172190021e-03, -3.679401704484100e-02, 4.973667402125381e-02}},
      {Operator::giao_nuclear,
       1,
       6,
       {
        2.460442029331187e-03, -4.366826681587838e-02, -5.759891736993161e-02,
        3.377611372954764e-01, 2.362062902718692e-01, -2.554252228641896e-01,
        3.366955365562792e-01, -1.762565168981901e-01, 4.540408423826872e-01,
        3.311974058010283e-01, -2.484648268297157e-01, 4.693330549244857e-01,
        9.202122892227845e-01, -8.643847592180640e-01, 1.909038098475657e-01,
        1.278456749678579e+00, -1.545792588752561e-01, 6.557448521222560e-01,
        -3.294829209560933e-01, 2.006579846240895e-01, -1.442876873410087e-01,
        -4.812926677195740e-01, 4.456572929277943e-01, -4.061141040453767e-01,
        -2.787455270698468e-01, 6.856532544092653e-02, 2.116950480793609e-01,
        4.044658477962266e-04, -1.391580084441719e-01, 2.828459374490686e-01,
        8.939642671803126e-01, -7.465862205213493e-01, 1.163207145484744e-01,
        9.615785874168974e-01, -2.029823461965428e-01, 5.799881045504384e-01,
        -7.652660862650316e-01, 5.503934212258096e-01, -2.177737514316954e-01,
        -1.786197661636123e+00, 5.560273262127098e-01, -4.263359482224903e-01,
        -1.316649981735281e+00, 5.107407078908260e-01, -4.195546507360799e-01,
        -6.614614288763730e-01, 1.757957880190346e-01, -2.859447157357363e-01,
        2.225698842782292e-01, 5.878240309937316e-03, -1.133752015063445e-01,
        -3.378859899335487e-01, -1.592622811645867e-01, 2.694438317957669e-02}},
      {Operator::giao_nuclear,
       2,
       3,
       {
        3.267722662547983e-01, -2.441137397605864e-01, 2.368597027880620e-01,
        -6.308444494430653e-03, -1.322461273386365e-01, 4.793674898470501e-01,
        6.319791835914814e-01, 1.162486359701446e-01, -1.550126308611172e-01,
        3.086427037374170e-02, -1.401469309294065e-01, -5.462201692012262e-01,
        -7.366388233832938e-01, 2.701284930283320e-01, -2.419973985454553e-01,
        -6.286012280636120e-03, 2.430648035607565e-01, -3.790518106526258e-01}},
      {Operator::giao_kinetic,
       1,
       6,
       {
        1.037122738312315e-03, 1.380889459930640e-02, 1.805635176994818e-02,
        -1.150959040211149e-02, 4.010815444454655e-03, -1.167965893897036e-05,
        6.112552923653641e-03, 2.551263407430333e-02, -6.815347839297765e-02,
        7.777442626611892e-03, 9.970718659299174e-02, -3.849383176613848e-02,
        -9.855229449389378e-03, 1.395522710307952e-01, -3.061488085167377e-02,
        -1.665995855336069e-01, 7.250709263599402e-02, 8.651828913569444e-02,
        -1.491192021885934e-02, -3.275319321677727e-02, 8.924619837516913e-02,
        1.218729086832224e-02, -8.406113260269307e-02, 3.942555863110796e-02,
        1.059519617797541e-03, -8.922651002972426e-03, -7.946397133482045e-04,
        4.381097227530276e-19, 1.784530200594488e-02, -1.059519617797559e-03,
        -1.207415806924580e-01, 1.331909224278590e-01, -6.556833090270556e-02,
        -2.437458173664463e-02, 3.275319321677725e-02, 7.171430386796093e-02,
        -3.648636905860770e-02, -1.607141995131127e-01, 3.711978899714818e-02,
        1.838639711367739e-01, -7.762614263595982e-02, -4.171485524675719e-02,
        1.052441759760995e-01, -1.379761377044467e-01, 4.820499541840010e-02,
        -1.555488525322378e-02, -2.551263407430332e-02, -4.048161829112986e-02,
        6.516587804413372e-03, -9.855414842789682e-03, -2.372163467699243e-02,
        1.150959040211151e-02, -1.191777495748812e-02, -7.542030883786703e-03}},
      {Operator::giao_kinetic,
       2,
       3,
       {
        -2.442898526223603e-02, -1.665595730651388e-03, -2.522406920069620e-02,
        6.265610796321857e-03, -2.683466547548116e-03, -1.808151843426069e-02,
        -2.050537176196199e-02, -4.516088027403439e-03, 1.760889978999806e-03,
        3.599457313195011e-05, 2.887864133770545e-03, 1.507120165050615e-02,
        4.233644309343944e-02, 4.371312415822540e-03, 3.279357322997598e-02,
        -8.371288001349122e-03, 2.202781999697229e-03, 1.693192855496372e-02}},
      {Operator::angular_momentum,
       1,
       6,
       {
        7.143628858920698e-03, 1.109658786952543e-02, -1.682033628879182e-02,
        -2.443076838421388e-03, 4.549749415969957e-02, -8.994320740519685e-02,
        -5.717791463284815e-02, 1.992551534740832e-01, 1.236608137832722e-01,
        -3.202037857058914e-01, -4.559677382016371e-02, 1.103250364167103e-01,
        -3.566164634584098e-01, 1.267739986447860e-01, -5.837228599112469e-02,
        -8.613745552590392e-02, 4.892217921807013e-02, -7.771592852972897e-02,
        2.575266150427902e-01, -2.320856296594157e-01, -1.944140182639518e-02,
        1.218985525490967e-01, 1.213711866672886e-01, -8.613999432072017e-02,
        -3.314798102414900e-02, -6.692048419704509e-03, 3.118274094778163e-02,
        2.453200280535523e-02, -3.616922794421403e-02, 9.184105596124971e-02,
        -8.351911226299238e-03, 6.858977690209898e-02, -3.619616226382969e-02,
        -3.821486033955018e-01, 6.303071622441520e-02, -6.719473263475798e-02,
        2.961034936737166e-01, -1.299283700480787e-01, 1.083074249372067e-01,
        4.073806517282538e-02, -4.504279169236408e-02, 1.450421162596349e-01,
        -3.729400712145510e-02, -7.460197362957013e-02, 1.180517613038327e-02,
        3.327737910725738e-01, -1.099131620933877e-01, 1.214953145452036e-01,
        2.600435216522820e-02, -4.404539449820905e-03, -1.436240465898984e-02,
        -2.208892596693390e-02, -9.328266215485517e-03, -1.897848556052905e-03}},
      {Operator::angular_momentum,
       2,
       3,
       {
        1.445054896484761e-01, 6.857508713122527e-02, -1.652942425193157e-02,
        1.086327330044552e-01, 2.537009664265523e-02, 2.123878011814298e-03,
        -2.724685031368239e-01, -5.005847271566102e-02, 1.997105708477781e-02,
        -1.031662395112292e-01, 1.830161411846225e-02, 4.164694738576923e-02,
        4.330341154401356e-02, -1.181229341412690e-01, 1.463892285743353e-02,
        -1.453572798687144e-01, -1.081381146352746e-01, -7.925379196053901e-02}}};

  for (const auto& ref : refs) {
    if (std::max(shells[ref.s1].contr[0].l, shells[ref.s2].contr[0].l) >
        LIBINT2_MAX_AM_giao_overlap)
      continue;
    Engine engine(ref.oper, 2, 2);
    if (ref.oper == Operator::giao_nuclear) engine.set_params(charges(atoms));
    if (ref.oper == Operator::angular_momentum)
      engine.set_params(gauge_origin);
    engine.compute(shells[ref.s1], shells[ref.s2]);
    const auto n12 = shells[ref.s1].size() * shells[ref.s2].size();
    REQUIRE(engine.results().size() == 3);
    REQUIRE(ref.values.size() == 3 * n12);
    for (size_t c = 0; c != 3; ++c) {
      for (size_t ij = 0; ij != n12; ++ij) {
        INFO("oper " << static_cast<int>(ref.oper) << ", shells (" << ref.s1
                     << "," << ref.s2 << "), component " << c << ", element "
                     << ij);
        REQUIRE(engine.results()[c][ij] ==
                Approx(ref.values[c * n12 + ij]).epsilon(1e-10).margin(1e-12));
      }
    }
  }
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "GIAO and angular momentum integral identities",
                 "[engine][1-body][giao]") {
  if (obs.max_l() > LIBINT2_MAX_AM_giao_overlap) return;
  const std::vector<Shell>& shells = obs.shells();
  const auto n = libint2::nbf(shells);
  const auto AB = center_differences(shells);

  SECTION("matrices are antisymmetric") {
    for (auto oper : magnetic_opers) {
      const auto X = compute(oper, shells, atoms, {{0.3, -0.4, 0.7}});
      REQUIRE(X.size() == 3);
      auto minus_Xt = X;
      for (size_t c = 0; c != 3; ++c)
        for (size_t i = 0; i != n; ++i)
          for (size_t j = 0; j != n; ++j)
            minus_Xt[c][i * n + j] = -X[c][j * n + i];
      require_equal(X, minus_Xt, 1e-12,
                    "oper " + std::to_string(static_cast<int>(oper)));
    }
  }

  SECTION("giao_overlap = 1/2 R_AB x dipole") {
    // emultipole1 components: overlap, x, y, z
    const auto D = compute(Operator::emultipole1, shells, atoms);
    Matrices X_ref(3, std::vector<double>(n * n));
    for (size_t ij = 0; ij != n * n; ++ij)
      for (int k = 0; k != 3; ++k) {
        const auto l = (k + 1) % 3, m = (k + 2) % 3;
        X_ref[k][ij] =
            0.5 * (AB[ij][l] * D[1 + m][ij] - AB[ij][m] * D[1 + l][ij]);
      }
    require_equal(compute(Operator::giao_overlap, shells, atoms), X_ref,
                  1e-12, "giao_overlap");
  }

  SECTION("translation shifts GIAO integrals by 1/2 (R_AB x d) O") {
    // the London phase factors are defined relative to the coordinate
    // origin, hence translating the system by d changes
    // X_k by 1/2 (R_AB x d)_k times the integral of the parent operator
    const std::array<double, 3> d = {{0.7, -1.1, 0.4}};
    auto shifted_shells = shells;
    for (auto& s : shifted_shells)
      s.move({{s.O[0] + d[0], s.O[1] + d[1], s.O[2] + d[2]}});
    auto shifted_atoms = atoms;
    for (auto& a : shifted_atoms) {
      a.x += d[0];
      a.y += d[1];
      a.z += d[2];
    }
    const std::vector<std::pair<Operator, Operator>> giao_and_parent = {
        {Operator::giao_overlap, Operator::overlap},
        {Operator::giao_kinetic, Operator::kinetic},
        {Operator::giao_nuclear, Operator::nuclear}};
    for (const auto& op : giao_and_parent) {
      const auto X = compute(op.first, shells, atoms);
      const auto O = compute(op.second, shells, atoms)[0];
      auto X_ref = X;
      for (size_t ij = 0; ij != n * n; ++ij)
        for (int k = 0; k != 3; ++k) {
          const auto l = (k + 1) % 3, m = (k + 2) % 3;
          X_ref[k][ij] +=
              0.5 * (AB[ij][l] * d[m] - AB[ij][m] * d[l]) * O[ij];
        }
      require_equal(compute(op.first, shifted_shells, shifted_atoms), X_ref,
                    1e-11, "translated GIAO");
    }
  }

#if LIBINT2_MAX_DERIV_ORDER > 0 && LIBINT2_TASK_EXISTS_2overlap1
  SECTION("angular_momentum origin dependence") {
    // (r - O') x nabla = (r - O) x nabla + (O - O') x nabla, and
    // (a| nabla_m |b) = - d(a|b)/dB_m
    const std::array<double, 3> O1 = {{0.3, -0.4, 0.7}};
    const std::array<double, 3> O2 = {{-0.5, 0.9, 0.2}};
    const auto L1 = compute(Operator::angular_momentum, shells, atoms, O1);
    const auto L2 = compute(Operator::angular_momentum, shells, atoms, O2);
    // overlap derivatives: d/dA_{x,y,z}, d/dB_{x,y,z}
    const auto dS = compute(Operator::overlap, shells, atoms, {{0, 0, 0}}, 1);
    auto L2_ref = L1;
    for (size_t ij = 0; ij != n * n; ++ij)
      for (int k = 0; k != 3; ++k) {
        const auto l = (k + 1) % 3, m = (k + 2) % 3;
        const auto nabla_l = -dS[3 + l][ij];
        const auto nabla_m = -dS[3 + m][ij];
        L2_ref[k][ij] +=
            (O1[l] - O2[l]) * nabla_m - (O1[m] - O2[m]) * nabla_l;
      }
    require_equal(L2, L2_ref, 1e-12, "angular_momentum");
  }

  SECTION("geometric derivatives are not supported") {
    for (auto oper : magnetic_opers)
      REQUIRE_THROWS_AS(Engine(oper, 1, 1, 1), std::invalid_argument);
  }
#endif
}

#endif  // defined(LIBINT2_SUPPORT_ONEBODY)

#if defined(LIBINT2_SUPPORT_GIAO_ERI)

TEST_CASE("GIAO electron repulsion integrals match libcint",
          "[engine][2-body][giao]") {
  const auto& shells = reference_shells();

  // reference values computed with libcint 6.2, converted to the libint
  // conventions and normalization:
  //   X(1)(ab|cd) = -int2e_ig1(ab|cd), X(2)(ab|cd) = -int2e_ig1(cd|ab)
  // (these are also the conventions of the PySCF NMR module)
  // the quartets exercise all shell permutations done by the Engine:
  // none, bra and ket function swaps, bra-ket swap
  struct Reference {
    std::array<size_t, 4> s;     // shell indices
    std::vector<double> values;  // [comp][i][j][k][l]
  };
  const std::vector<Reference> refs = {
      {{{1, 0, 4, 3}},
       {
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, 0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, 0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, -0.000000000000000e+00, -0.000000000000000e+00,
        -0.000000000000000e+00, 0.000000000000000e+00, -0.000000000000000e+00}},
      {{{0, 4, 5, 1}},
       {
        -2.065954639780780e-03, 1.748107884968873e-03, -5.779538547816009e-04,
        -3.004752832249725e-03, 3.784505391540075e-03, -1.326620812183039e-03,
        -3.712142151877517e-04, -2.070666028592435e-04, 1.859900113005574e-05,
        2.891816662586581e-03, -2.900133009039144e-03, -1.067742973781234e-04,
        8.724632918739669e-05, 1.049526578608104e-04, -2.770093994813214e-04,
        -2.745570390948598e-03, 4.486663821807764e-03, -2.336566552615843e-03,
        1.377550632761716e-03, -9.497947947017615e-04, 8.214500432221932e-04,
        3.964791238624681e-03, -5.095984644844295e-03, 1.900736987425634e-03,
        1.802366758797287e-03, -1.860417778000897e-03, 1.087852071167470e-03,
        5.081627480860406e-05, -4.704384959912405e-03, -1.011515287968272e-02,
        -1.451973201262578e-04, -1.010581110967298e-03, -3.073322118072156e-03,
        -2.707626785150762e-04, -2.976859204637611e-03, -7.479382737047251e-03,
        4.602984449351487e-03, -4.928982936260262e-04, -8.614888412861848e-03,
        1.334656529797825e-03, -2.904258233299913e-05, -2.593053645962293e-03,
        3.105690721371375e-03, 2.222013558388264e-04, -6.443723508243654e-03,
        1.052063925657853e-02, 8.271312319149370e-03, 3.497940373765485e-04,
        3.370371247478266e-03, 1.954140878089213e-03, 1.626742839236361e-04,
        7.708503944810248e-03, 6.466490768903284e-03, 8.863430122453294e-05}},
      {{{1, 4, 0, 3}},
       {
        4.710258148921315e-03, 1.347407271372066e-02, 1.460091585942131e-03,
        -1.061092429024923e-02, -3.919680458138358e-03, -9.729744651068968e-03,
        1.292913777830875e-02, 1.243009102775367e-02, 2.126272673334751e-03,
        3.953079956810819e-03, 5.336876676684140e-03, 2.509640425875196e-02,
        -3.264907967590822e-03, -4.108272944356591e-04, 5.301856545353345e-03,
        -2.508067283285199e-02, -6.167189647013270e-03, -3.107828140059753e-03,
        -8.162763225614522e-03, -2.050680013100096e-02, -1.389745747589996e-02,
        1.570261713347079e-02, 5.421872655820505e-03, 1.044829927506656e-02,
        -5.295672831624999e-03, -1.363669772890333e-02, -1.355112069179786e-03,
        4.710258148921313e-03, 4.771279749388771e-03, 5.604278711814460e-03,
        3.464264492530932e-04, -3.114481252720347e-03, -3.166476577353912e-05,
        7.711351711879068e-03, 4.039501280779865e-03, 1.321073467916741e-03,
        1.698522181640389e-03, 1.423537126861882e-03, 5.274667251682687e-03,
        -1.567676231865755e-04, -4.108272944356587e-04, 1.157669419481015e-03,
        -3.943777123599978e-03, -9.494035805835935e-04, -8.532703648893272e-04,
        -7.089164285057173e-03, -7.039581154833542e-03, -9.984117926077700e-03,
        -3.872506832009930e-04, 4.348273715263157e-03, -5.090514644357644e-04,
        -8.403813176029241e-03, -4.933904764571442e-03, -1.355112069179784e-03}},
      {{{2, 3, 5, 0}},
       {
        -2.087681400185191e-03, 1.996066337152117e-03, -1.715855584329187e-03,
        1.550928393364600e-03, 8.959297546302851e-04, -3.524481105488151e-03,
        -4.913814898192811e-03, -1.531244390940009e-03, 1.173127372182342e-03,
        -1.418948718667479e-04, 1.950923307859170e-03, 4.758089576636779e-03,
        5.123487056529212e-03, -1.932257787183771e-03, 1.729175363780658e-03,
        -2.000335537882920e-03, -2.123584105154270e-03, 2.433551199395259e-03,
        -1.392507377629879e-04, 7.965900088086791e-05, 1.510328078046807e-04,
        2.510549622653471e-04, 4.493852984378240e-04, 2.323447153609037e-04,
        -2.238757891564961e-05, -1.320680149216821e-04, 1.055796288942382e-04,
        1.083728716016340e-04, 3.336683461278232e-04, 4.142483944608074e-04,
        2.268378318744764e-04, -4.640903438886944e-04, -5.842031816111950e-05,
        -2.520186823730769e-04, -1.287667211960560e-04, 4.912684026492867e-04}}};

  Engine engine(Operator::giao_coulomb, 2, 2);
  for (const auto& ref : refs) {
    engine.compute(shells[ref.s[0]], shells[ref.s[1]], shells[ref.s[2]],
                   shells[ref.s[3]]);
    size_t n1234 = 1;
    for (auto s : ref.s) n1234 *= shells[s].size();
    REQUIRE(engine.results().size() == 6);
    REQUIRE(ref.values.size() == 6 * n1234);
    for (size_t c = 0; c != 6; ++c) {
      for (size_t ijkl = 0; ijkl != n1234; ++ijkl) {
        INFO("shells (" << ref.s[0] << "," << ref.s[1] << "|" << ref.s[2]
                        << "," << ref.s[3] << "), component " << c
                        << ", element " << ijkl);
        REQUIRE(engine.results()[c][ijkl] ==
                Approx(ref.values[c * n1234 + ijkl])
                    .epsilon(1e-10)
                    .margin(1e-12));
      }
    }
  }
}

TEST_CASE_METHOD(libint2::unit::DefaultFixture,
                 "GIAO electron repulsion integral identities",
                 "[engine][2-body][giao]") {
  if (obs.max_l() > LIBINT2_MAX_AM_giao_eri) return;
  const std::vector<Shell>& shells = obs.shells();
  const auto n = libint2::nbf(shells);
  auto idx = [n](size_t i, size_t j, size_t k, size_t l) {
    return ((i * n + j) * n + k) * n + l;
  };
  Engine engine(Operator::giao_coulomb, obs.max_nprim(), obs.max_l());
  const auto X = compute4(engine, shells);
  REQUIRE(X.size() == 6);

  SECTION("permutational symmetry") {
    // X(1) is antisymmetric w.r.t. a<->b, symmetric w.r.t. c<->d,
    // vice versa for X(2); X(1)(ab|cd) = X(2)(cd|ab)
    auto X_ab = X, X_cd = X, X_braket = X;  // X with indices permuted
    for (size_t c = 0; c != 6; ++c) {
      const auto sign_ab = c < 3 ? -1 : 1;
      const auto sign_cd = c < 3 ? 1 : -1;
      for (size_t i = 0; i != n; ++i)
        for (size_t j = 0; j != n; ++j)
          for (size_t k = 0; k != n; ++k)
            for (size_t l = 0; l != n; ++l) {
              X_ab[c][idx(i, j, k, l)] = sign_ab * X[c][idx(j, i, k, l)];
              X_cd[c][idx(i, j, k, l)] = sign_cd * X[c][idx(i, j, l, k)];
              X_braket[c][idx(i, j, k, l)] = X[(c + 3) % 6][idx(k, l, i, j)];
            }
    }
    require_equal(X, X_ab, 1e-12, "a<->b");
    require_equal(X, X_cd, 1e-12, "c<->d");
    require_equal(X, X_braket, 1e-12, "(ab|<->|cd)");
  }

  SECTION("translation shifts GIAO ERIs by 1/2 (R x d) (ab|cd)") {
    const std::array<double, 3> d = {{0.7, -1.1, 0.4}};
    auto shifted_shells = shells;
    for (auto& s : shifted_shells)
      s.move({{s.O[0] + d[0], s.O[1] + d[1], s.O[2] + d[2]}});
    Engine eri_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    const auto eri = compute4(eri_engine, shells)[0];
    const auto X_shifted = compute4(engine, shifted_shells);
    const auto AB = center_differences(shells);
    auto X_ref = X;
    for (size_t i = 0; i != n; ++i)
      for (size_t j = 0; j != n; ++j)
        for (size_t k = 0; k != n; ++k)
          for (size_t l = 0; l != n; ++l)
            for (int x = 0; x != 3; ++x) {
              const auto y = (x + 1) % 3, z = (x + 2) % 3;
              const auto& Rab = AB[i * n + j];
              const auto& Rcd = AB[k * n + l];
              const auto ijkl = idx(i, j, k, l);
              X_ref[x][ijkl] +=
                  0.5 * (Rab[y] * d[z] - Rab[z] * d[y]) * eri[ijkl];
              X_ref[3 + x][ijkl] +=
                  0.5 * (Rcd[y] * d[z] - Rcd[z] * d[y]) * eri[ijkl];
            }
    require_equal(X_shifted, X_ref, 1e-11, "translated giao_coulomb");
  }

  SECTION("unsupported uses throw") {
    Engine e(Operator::giao_coulomb, 1, 1);
    REQUIRE_THROWS_AS(e.set(libint2::BraKet::xs_xx), std::invalid_argument);
#if LIBINT2_MAX_DERIV_ORDER > 0
    REQUIRE_THROWS_AS(Engine(Operator::giao_coulomb, 1, 1, 1),
                      std::invalid_argument);
#endif
  }
}

#endif  // defined(LIBINT2_SUPPORT_GIAO_ERI)
