#include "catch.hpp"

#include <libint2.hpp>

using std::cout;
using std::cerr;
using std::endl;

using libint2::Atom;
using libint2::BasisSet;
using libint2::Shell;
using libint2::Engine;
using libint2::Operator;

std::tuple<int, int, int, int> permute4_1234(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i1,i2,i3,i4);
};
std::tuple<int, int, int, int> permute4_2134(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i2,i1,i3,i4);
};
std::tuple<int, int, int, int> permute4_1243(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i1,i2,i4,i3);
};
std::tuple<int, int, int, int> permute4_2143(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i2,i1,i4,i3);
};
std::tuple<int, int, int, int> permute4_3412(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i3,i4,i1,i2);
};
std::tuple<int, int, int, int> permute4_3421(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i3,i4,i2,i1);
};
std::tuple<int, int, int, int> permute4_4312(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i4,i3,i1,i2);
};
std::tuple<int, int, int, int> permute4_4321(int i1, int i2, int i3, int i4) {
  return std::make_tuple(i4,i3,i2,i1);
};

template <unsigned int deriv_order>
void validate4(const BasisSet& obs, const std::vector<Atom>& atoms) {
  if (deriv_order > 1 || deriv_order > LIBINT2_MAX_DERIV_ORDER)
    return;

  Engine engine_ref(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  Engine engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);

  const auto nderiv_shellset =
      libint2::num_geometrical_derivatives(4, deriv_order); // # of derivs for each shell quartet

  auto shell2bf = obs.shell2bf(); // maps shell index to basis function index
  auto shell2atom = obs.shell2atom(atoms);

  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    for (auto s2 = 0; s2 <= s1; ++s2) {
      for (auto s3 = 0; s3 <= s1; ++s3) {
        auto s4_max = (s1 == s3) ? s2 : s3;
        for (auto s4 = 0; s4 <= s4_max; ++s4) {

          engine_ref.compute(obs[s1], obs[s2], obs[s3], obs[s4]);
          const auto& shellset_ref = engine_ref.results();  // location of the computed integrals
          if (shellset_ref[0] == nullptr)
            continue;  // nullptr returned if the entire shell-set was screened out

          auto n1 = obs[s1].size(); // number of basis functions in first shell
          auto n2 = obs[s2].size(); // number of basis functions in second shell
          auto n3 = obs[s3].size(); // number of basis functions in third shell
          auto n4 = obs[s4].size(); // number of basis functions in fourth shell

          auto validate = [&](std::tuple<int, int, int, int> (&permuter)(int,int,int,int)) {
            int ss1, ss2, ss3, ss4;
            std::tie(ss1, ss2, ss3, ss4) = permuter(s1, s2, s3, s4);

            engine.compute(obs[ss1], obs[ss2], obs[ss3], obs[ss4]);
            const auto& shellset = engine.results();  // location of the computed integrals

            for(int d=0; d!=nderiv_shellset; ++d) {
              REQUIRE(shellset[d] != nullptr);

              auto compare = [&](int dd) {
                int nn1, nn2, nn3, nn4;
                std::tie(nn1, nn2, nn3, nn4) = permuter(n1, n2, n3, n4);
                for (auto f1 = 0; f1 != n1; ++f1) {
                  for (auto f2 = 0; f2 != n2; ++f2) {
                    for (auto f3 = 0; f3 != n3; ++f3) {
                      for (auto f4 = 0; f4 != n4; ++f4) {
                        int ff1, ff2, ff3, ff4;
                        std::tie(ff1, ff2, ff3, ff4) = permuter(f1, f2, f3, f4);

                        auto f1234 = ((f1 * n2 + f2) * n3 + f3) * n4 + f4;
                        auto ff1234 = ((ff1 * nn2 + ff2) * nn3 + ff3) * nn4 + ff4;

                        REQUIRE(shellset[dd][ff1234] == Approx(shellset_ref[d][f1234]).margin(deriv_level == 1 ? 1e-12 : 1e-14));
                      }
                    }
                  }
                }
              };

              // permuting geometric deriv index is order-specific
              int dd;
              switch (deriv_order) {
                case 0:
                  dd = 0;
                  break;

                case 1: {
                  const auto xyz = d % 3;
                  const auto atom = d / 3;
                  int a[4] = {0, 0, 0, 0};
                  a[atom] = 1;
                  int aa[4];
                  std::tie(aa[0], aa[1], aa[2], aa[3]) = permuter(a[0], a[1], a[2], a[3]);
                  if (aa[0]) dd = xyz;
                  else if (aa[1]) dd = 3 + xyz;
                  else if (aa[2]) dd = 6 + xyz;
                  else if (aa[3]) dd = 9 + xyz;
                  else
                    assert(false);
                }
                break;

                default:
                  assert(false);  // can't compare yet deriv_order > 1
              }
              compare(dd);

            }  // d

          };

          validate(permute4_2134);
          validate(permute4_1243);
          validate(permute4_2143);
          validate(permute4_3412);
          validate(permute4_3421);
          validate(permute4_4312);
          validate(permute4_4321);

        }
      }
    }
  }
}

TEST_CASE( "2-e 4-c integrals permute correctly", "[permute-2e-4c]" ) {

  auto atoms = std::vector<Atom>{ {8, 0.,0.,0.}, {8, 0.,0.,2.}, {1, 0.,-1.,-1.}, {1, 0.,1.,3.}};
  auto obs = BasisSet("6-31g**", atoms);

  SECTION( "deriv_order=0" ) {
    validate4<0>(obs, atoms);
  }  // section
  SECTION( "deriv_order=1" ) {
    validate4<1>(obs, atoms);
  }  // section
}


