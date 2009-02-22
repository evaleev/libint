
#ifndef _libint2_src_bin_libint_cgshellinfo_h_
#define _libint2_src_bin_libint_cgshellinfo_h_

#include <libint2_config.h>
#include <exception.h>
#include <utility>
#include <algorithm>

namespace libint2 {

  namespace detail {
    int notxyz(int a, int b) {
      if (a == b)
        throw libint2::ProgrammingError("notxyz(a,b) -- a equals b");
      int amax = std::max(a,b);
      int amin = std::min(a,b);
      if (amin == 0 && amax == 1)
        return 2;
      if (amin == 0 && amax == 2)
        return 1;
      if (amin == 1 && amax == 2)
        return 0;
    }

    std::pair<int,int> notxyz(int a) {
      switch(a) {
      case 0: return std::make_pair(1,2); break;
      case 1: return std::make_pair(0,2); break;
      case 2: return std::make_pair(0,1); break;
      }
    }
  }


  template <CGShellOrdering Ord, unsigned int lmax> struct CGShellOrderingGenerator {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]);
  };
  template <unsigned int lmax> struct CGShellOrderingGenerator<CGShellOrdering_GAMESS,lmax> {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]) {

      for (unsigned int am = 0; am <= lmax; ++am) {

        if (am == 0) {
          cartindex[0][0][0] = 0;
          continue;
        }

        unsigned int current_index = 0;
        unsigned int qn[3] = { 0, 0, 0 };

        const int ammin = ((int) am + 2) / 3;
        for (int am1 = am; am1 >= ammin; --am1) {

          for (int xyz1 = 0; xyz1 < 3; ++xyz1) {

            qn[xyz1] = am1;

            // distribute the remaining quanta according to the following rules

            // "nothing to distribute" is a special case
            if (am - am1 == 0) {
              std::pair<int, int> xyz(detail::notxyz( xyz1));
              qn[xyz.first] = 0;
              qn[xyz.second] = 0;
              cartindex[am][qn[0]][qn[1]] = current_index;
              ++current_index;
            } else {
              int am23 = (int) am - qn[xyz1];
              const int maxam23 = std::min((int) qn[xyz1], am23);
              const int minam23 = (am23 + 1) / 2;
              for (int am2 = maxam23; am2 >= minam23; --am2) {
                const int xyz2min = (am2 == qn[xyz1]) ? xyz1 + 1 : 0;
                for (int xyz2 = xyz2min; xyz2 < 3; ++xyz2) {
                  if (xyz1 == xyz2)
                    continue;
                  qn[xyz2] = am2;
                  const int xyz3 = detail::notxyz(xyz1, xyz2);
                  qn[xyz3] = am23 - am2;
                  if (qn[xyz3] == qn[xyz1] && xyz3 < xyz1 || qn[xyz3] == qn[xyz2]
                      && xyz3 < xyz2)
                    continue;
                  {
                    cartindex[am][qn[0]][qn[1]] = current_index;
                    ++current_index;
                  }
                }
              }
            }

          }
        }
      }
    }
  };
  template <unsigned int lmax> struct CGShellOrderingGenerator<CGShellOrdering_ORCA,lmax> {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]) {
      // same as GAMESS
      CGShellOrderingGenerator<CGShellOrdering_GAMESS,lmax>::compute(cartindex);
      // except for p-type shells: z, x, y
      if (lmax > 0) {
        cartindex[1][0][0] = 0;
        cartindex[1][1][0] = 1;
        cartindex[1][0][1] = 2;
      }

      for(int l=0; l<=lmax; ++l) {
        for(int i=0; i<=l; ++i) {
          for(int j=0; j<=l-i; ++j) {
            std::cout << "CGShellOrderingGenerator<CGShellOrdering_ORCA>: cartindex[" << l << "]["
                      << i << "][" << j << "] = " << cartindex[l][i][j] << std::endl;
          }
        }
      }
    }
  };

  template <CGShellOrdering Ord, unsigned int lmax> struct CGShellOrderingData {

    struct Triple {
      Triple() : i(0), j(0), k(0) {}
      Triple(int ii, int jj, int kk) : i(ii), j(jj), k(kk) {}
      int i, j, k;
    };

    CGShellOrderingData() {
      // compute cartindex
      CGShellOrderingGenerator<Ord,lmax>::compute(cartindex);
      // then use it to compute cartindex_to_ijk
      for(unsigned int l=0; l<=lmax; ++l) {
        for(int i=0; i<=l; ++i) {
          for(int j=0; j<=l-i; ++j) {
            const int c = cartindex[l][i][j];
            cartindex_to_ijk[l][c] = Triple(i,j,l-i-j);
          }
        }
      }
    }

    int cartindex[lmax+1][lmax+1][lmax+1];
    Triple cartindex_to_ijk[lmax+1][(lmax+1)*(lmax+2)/2];
  };

  /// provides ordering maps for up to angular momentum lmax and ordering specified by CGShellOrderingSpec
  template <typename OrderingData> struct CGShellInfo {
    // computes cartindex xyz from am i j
    static int cartindex(unsigned int am, int i, int j) {
      return data_.cartindex[am][i][j];
    }
    // computes i j k from cartindex xyz
    static int cartindex_to_ijk(unsigned int am,
                                int xyz,
                                int& i,
                                int& j,
                                int& k) {
      const typename OrderingData::Triple& ijk = data_.cartindex_to_ijk[am][xyz];
      i = ijk.i;
      j = ijk.j;
      k = ijk.k;
    }

    private:
      static OrderingData data_;
  };

  template <typename OrderingData> OrderingData CGShellInfo<OrderingData>::data_;



};

#endif // header guard
