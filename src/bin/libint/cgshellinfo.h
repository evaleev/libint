
#ifndef _libint2_src_bin_libint_cgshellinfo_h_
#define _libint2_src_bin_libint_cgshellinfo_h_

#include <cassert>
#include <utility>
#include <algorithm>

namespace libint2 {

  namespace detail {
    inline int notxyz(int a, int b) {
      assert(a != b);
      int amax = std::max(a,b);
      int amin = std::min(a,b);
      if (amin == 0 && amax == 1)
        return 2;
      if (amin == 0 && amax == 2)
        return 1;
      if (amin == 1 && amax == 2)
        return 0;
      abort(); // unreachable
    }

   inline std::pair<int,int> notxyz(int a) {
      switch(a) {
      case 0: return std::make_pair(1,2); break;
      case 1: return std::make_pair(0,2); break;
      case 2: return std::make_pair(0,1); break;
      }
      abort(); // unreachable
    }
  }


  template <CGShellOrdering Ord, unsigned int lmax> struct CGShellOrderingGenerator {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]);
  };
  template <unsigned int lmax> struct CGShellOrderingGenerator<CGShellOrdering_Standard,lmax> {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]) {

      // see cgshell_ordering.h
      for (unsigned int am = 0; am <= lmax; ++am) {
        int count = 0;
        for(int i=(am);(i)>=0;(i)--) {
          for(int j=(am)-(i);(j)>=0;(j)--, ++count) {
            cartindex[am][i][j] = count;
          }
        }
      }

    }
  };
  template <unsigned int lmax> struct CGShellOrderingGenerator<CGShellOrdering_GAMESS,lmax> {
    static void compute(int (&cartindex)[lmax+1][lmax+1][lmax+1]) {

      for (unsigned int am = 0; am <= lmax; ++am) {

        if (am == 0) {
          cartindex[0][0][0] = 0;
          continue;
        }
        /// 12/15/2011 GAMESS ordering for g functions different from what I had assumed before (thanks to Luke Roskop)
        if (am == 4) {
          cartindex[4][4][0] = 0;
          cartindex[4][0][4] = 1;
          cartindex[4][0][0] = 2;
          cartindex[4][3][1] = 3;
          cartindex[4][3][0] = 4;
          cartindex[4][1][3] = 5;
          cartindex[4][0][3] = 6;
          cartindex[4][1][0] = 7;
          cartindex[4][0][1] = 8;
          cartindex[4][2][2] = 9;
          cartindex[4][2][0] = 10;
          cartindex[4][0][2] = 11;
          cartindex[4][2][1] = 12;
          cartindex[4][1][2] = 13;
          cartindex[4][1][1] = 14;
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
                  if ((qn[xyz3] == qn[xyz1] && xyz3 < xyz1) ||
                      (qn[xyz3] == qn[xyz2] && xyz3 < xyz2)
                     )
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

      // identical to GAMESS for s through g functions
      // identical to standard for h and higher
      CGShellOrderingGenerator<CGShellOrdering_Standard,lmax>::compute(cartindex);

      for (unsigned int am = 0; am <= std::min(4u,lmax); ++am) {

        if (am == 0) {
          cartindex[0][0][0] = 0;
          continue;
        }
        if (am == 1) {
          cartindex[1][0][0] = 0;
          cartindex[1][1][0] = 1;
          cartindex[1][0][1] = 2;
          continue;
        }
        if (am == 2) {
          cartindex[2][2][0] = 0;
          cartindex[2][0][2] = 1;
          cartindex[2][0][0] = 2;
          cartindex[2][1][1] = 3;
          cartindex[2][1][0] = 4;
          cartindex[2][0][1] = 5;
          continue;
        }
        if (am == 3) {
          cartindex[3][3][0] = 0;
          cartindex[3][0][3] = 1;
          cartindex[3][0][0] = 2;
          cartindex[3][2][1] = 3;
          cartindex[3][2][0] = 4;
          cartindex[3][1][2] = 5;
          cartindex[3][0][2] = 6;
          cartindex[3][1][0] = 7;
          cartindex[3][0][1] = 8;
          cartindex[3][1][1] = 9;
          continue;
        }
        if (am == 4) {
          cartindex[4][4][0] = 0;
          cartindex[4][0][4] = 1;
          cartindex[4][0][0] = 2;
          cartindex[4][3][1] = 3;
          cartindex[4][3][0] = 4;
          cartindex[4][1][3] = 5;
          cartindex[4][0][3] = 6;
          cartindex[4][1][0] = 7;
          cartindex[4][0][1] = 8;
          cartindex[4][2][2] = 9;
          cartindex[4][2][0] = 10;
          cartindex[4][0][2] = 11;
          cartindex[4][2][1] = 12;
          cartindex[4][1][2] = 13;
          cartindex[4][1][1] = 14;
          continue;
        }
        /*
        if (am == 5) {
          cartindex[5][5][0] = 0;
          cartindex[5][4][1] = 1;
          cartindex[5][4][0] = 2;
          cartindex[5][3][2] = 3;
          cartindex[5][3][1] = 4;
          cartindex[5][3][0] = 5;
          cartindex[5][2][3] = 6;
          cartindex[5][2][2] = 7;
          cartindex[5][2][1] = 8;
          cartindex[5][2][0] = 9;
          cartindex[5][1][4] = 10;
          cartindex[5][1][3] = 11;
          cartindex[5][1][2] = 12;
          cartindex[5][1][1] = 13;
          cartindex[5][1][0] = 14;
          cartindex[5][0][5] = 15;
          cartindex[5][0][4] = 16;
          cartindex[5][0][3] = 17;
          cartindex[5][0][2] = 18;
          cartindex[5][0][1] = 19;
          cartindex[5][0][0] = 20;
          continue;
        }
        if (am == 6) {
          cartindex[6][6][0] = 0;
          cartindex[6][5][1] = 1;
          cartindex[6][5][0] = 2;
          cartindex[6][4][2] = 3;
          cartindex[6][4][1] = 4;
          cartindex[6][4][0] = 5;
          cartindex[6][3][3] = 6;
          cartindex[6][3][2] = 7;
          cartindex[6][3][1] = 8;
          cartindex[6][3][0] = 9;
          cartindex[6][2][4] = 10;
          cartindex[6][2][3] = 11;
          cartindex[6][2][2] = 12;
          cartindex[6][2][1] = 13;
          cartindex[6][2][0] = 14;
          cartindex[6][1][5] = 15;
          cartindex[6][1][4] = 16;
          cartindex[6][1][3] = 17;
          cartindex[6][1][2] = 18;
          cartindex[6][1][1] = 19;
          cartindex[6][1][0] = 20;
          cartindex[6][0][6] = 21;
          cartindex[6][0][5] = 22;
          cartindex[6][0][4] = 23;
          cartindex[6][0][3] = 24;
          cartindex[6][0][2] = 25;
          cartindex[6][0][1] = 26;
          cartindex[6][0][0] = 27;
          continue;
        }
        if (am == 7) {
          cartindex[7][7][0] = 0;
          cartindex[7][6][1] = 1;
          cartindex[7][6][0] = 2;
          cartindex[7][5][2] = 3;
          cartindex[7][5][1] = 4;
          cartindex[7][5][0] = 5;
          cartindex[7][4][3] = 6;
          cartindex[7][4][2] = 7;
          cartindex[7][4][1] = 8;
          cartindex[7][4][0] = 9;
          cartindex[7][3][4] = 10;
          cartindex[7][3][3] = 11;
          cartindex[7][3][2] = 12;
          cartindex[7][3][1] = 13;
          cartindex[7][3][0] = 14;
          cartindex[7][2][5] = 15;
          cartindex[7][2][4] = 16;
          cartindex[7][2][3] = 17;
          cartindex[7][2][2] = 18;
          cartindex[7][2][1] = 19;
          cartindex[7][2][0] = 20;
          cartindex[7][1][6] = 21;
          cartindex[7][1][5] = 22;
          cartindex[7][1][4] = 23;
          cartindex[7][1][3] = 24;
          cartindex[7][1][2] = 25;
          cartindex[7][1][1] = 26;
          cartindex[7][1][0] = 27;
          cartindex[7][0][7] = 28;
          cartindex[7][0][6] = 29;
          cartindex[7][0][5] = 30;
          cartindex[7][0][4] = 31;
          cartindex[7][0][3] = 32;
          cartindex[7][0][2] = 33;
          cartindex[7][0][1] = 34;
          cartindex[7][0][0] = 35;
          continue;
        }
        if (am == 8) {
          cartindex[8][8][0] = 0;
          cartindex[8][7][1] = 1;
          cartindex[8][7][0] = 2;
          cartindex[8][6][2] = 3;
          cartindex[8][6][1] = 4;
          cartindex[8][6][0] = 5;
          cartindex[8][5][3] = 6;
          cartindex[8][5][2] = 7;
          cartindex[8][5][1] = 8;
          cartindex[8][5][0] = 9;
          cartindex[8][4][4] = 10;
          cartindex[8][4][3] = 11;
          cartindex[8][4][2] = 12;
          cartindex[8][4][1] = 13;
          cartindex[8][4][0] = 14;
          cartindex[8][3][5] = 15;
          cartindex[8][3][4] = 16;
          cartindex[8][3][3] = 17;
          cartindex[8][3][2] = 18;
          cartindex[8][3][1] = 19;
          cartindex[8][3][0] = 20;
          cartindex[8][2][6] = 21;
          cartindex[8][2][5] = 22;
          cartindex[8][2][4] = 23;
          cartindex[8][2][3] = 24;
          cartindex[8][2][2] = 25;
          cartindex[8][2][1] = 26;
          cartindex[8][2][0] = 27;
          cartindex[8][1][7] = 28;
          cartindex[8][1][6] = 29;
          cartindex[8][1][5] = 30;
          cartindex[8][1][4] = 31;
          cartindex[8][1][3] = 32;
          cartindex[8][1][2] = 33;
          cartindex[8][1][1] = 34;
          cartindex[8][1][0] = 35;
          cartindex[8][0][8] = 36;
          cartindex[8][0][7] = 37;
          cartindex[8][0][6] = 38;
          cartindex[8][0][5] = 39;
          cartindex[8][0][4] = 40;
          cartindex[8][0][3] = 41;
          cartindex[8][0][2] = 42;
          cartindex[8][0][1] = 43;
          cartindex[8][0][0] = 44;
          continue;
        }
        */
      }

#if 0
      for(int l=0; l<=lmax; ++l) {
        for(int i=0; i<=l; ++i) {
          for(int j=0; j<=l-i; ++j) {
            std::cout << "CGShellOrderingGenerator<CGShellOrdering_ORCA>: cartindex[" << l << "]["
                      << i << "][" << j << "] = " << cartindex[l][i][j] << std::endl;
          }
        }
      }
#endif

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
    static void cartindex_to_ijk(unsigned int am,
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
