#
# check if CXX provides C++11 support with CXXFLAGS, possibly by appending -std=c++11 to CXXFLAGS
#
AC_DEFUN([ACX_CHECK_CPP11_GENERAL], [
  AC_LANG_SAVE
  AC_LANG([C++])

  acx_have_cxx11=no
  AC_MSG_CHECKING([CXX for general C++11 support])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
      [[#include <vector>
        #include <iostream>
        #include "$srcdir/include/libint2/util/cxxstd.h"
        #if LIBINT2_CPLUSPLUS_STD < 2011
        # error "no C++11 support"
        #endif
        auto f(int x, double y) -> decltype(x*y) {
          return x*y;
        }
      ]],
      [[
       std::vector<double> reals(10, 0.0);
       for (const auto &i: reals) {
         std::cout << i << "\n";
       }
      ]]
     )
    ],
    [acx_have_cxx11=yes
     AC_DEFINE(LIBINT_HAS_CXX11)
     AC_MSG_RESULT([yes])
    ]
  )
  # if not found, try adding -std flag
  if test "X$acx_have_cxx11" = "Xno"; then
    old_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS -std=c++11"
    CXXFLAGS_ADDED_DASHSTD=1
    AC_COMPILE_IFELSE(
     [AC_LANG_PROGRAM(
       [[#include <vector>
         #include <iostream>
         #include "$srcdir/include/libint2/util/cxxstd.h"
         #if LIBINT2_CPLUSPLUS_STD < 2011
         # error "no C++11 support"
         #endif
         auto f(int x, double y) -> decltype(x*y) {
           return x*y;
         }
       ]],
       [[
        std::vector<double> reals(10, 0.0);
        for (const auto &i: reals) {
          std::cout << i << "\n";
        }
       ]]
      )
     ],
     [acx_have_cxx11=yes
      AC_MSG_RESULT([yes (with -std=c++11)])
      AC_DEFINE(LIBINT_HAS_CXX11)
     ],
     [
      AC_MSG_RESULT([no])
      CXXFLAGS=$old_CXXFLAGS
     ]
   )
  fi

  AC_LANG_RESTORE  
])

#
# check if CXXGEN provides C++11 support with CXXGENFLAGS, possibly by appending -std=c++11 to CXXGENFLAGS
#
AC_DEFUN([ACX_CHECK_CPP11_CXXGEN], [
  AC_LANG_SAVE
  AC_LANG([C++])

  AC_MSG_CHECKING([CXXGEN for general C++11 support])
  ref_CXX=$CXX
  CXX=$CXXGEN
  ref_CXXFLAGS=$CXXFLAGS
  CXXFLAGS=$CXXGENFLAGS
  
  CXXGEN_SUPPORTS_CPP11=no
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
      [[#include <vector>
        #include <iostream>
        #include "$srcdir/include/libint2/util/cxxstd.h"
        #if LIBINT2_CPLUSPLUS_STD < 2011
        # error "no C++11 support"
        #endif
        auto f(int x, double y) -> decltype(x*y) {
          return x*y;
        }
      ]],
      [[
       std::vector<double> reals(10, 0.0);
       for (const auto& i: reals) {
         std::cout << i << "\n";
       }
      ]]
     )
    ],
    [
     CXXGEN_SUPPORTS_CPP11=yes
     AC_DEFINE(LIBINT_HAS_CXX11_CXXGEN)
     AC_MSG_RESULT([yes])
    ]
  )
  # if not found, try adding -std flag
  if test "X$CXXGEN_SUPPORTS_CPP11" = "Xno"; then
    CXXFLAGS="$CXXFLAGS -std=c++11"
    CXXGENFLAGS="$CXXGENFLAGS -std=c++11"
    CXXGENFLAGS_ADDED_DASHSTD=1
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM(
       [[#include <vector>
         #include <iostream>
         #include "$srcdir/include/libint2/util/cxxstd.h"
         #if LIBINT2_CPLUSPLUS_STD < 2011
         # error "no C++11 support"
         #endif
          auto f(int x, double y) -> decltype(x*y) {
            return x*y;
          }
       ]],
       [[
         std::vector<double> reals(10, 0.0);
         for (const auto& i: reals) {
           std::cout << i << "\n";
         }
       ]]
      )],
      [
        CXXGEN_SUPPORTS_CPP11=yes
        AC_DEFINE(LIBINT_HAS_CXX11_CXXGEN)
        AC_MSG_RESULT([yes])
      ],
      [
        AC_MSG_RESULT([no])
      ]
    )
  fi

  AC_SUBST(CXXGEN_SUPPORTS_CPP11)
  # revert CXX and CXXFLAGS
  CXX=$ref_CXX
  CXXFLAGS=$ref_CXXFLAGS

  AC_LANG_RESTORE  
])

AC_DEFUN([ACX_CHECK_CPP11],
[
  ACX_CHECK_CPP11_GENERAL
])
