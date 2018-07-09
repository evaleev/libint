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
# check if CXXGEN provides C++11 support with CXXGENFLAGS, without any additional flags
#
AC_DEFUN([ACX_CHECK_CPP11_CXXGEN], [
  AC_LANG_SAVE
  AC_LANG([C++])

  AC_MSG_CHECKING([CXXGEN for general C++11 support])
  ref_CXX=$CXX
  CXX=$CXXGEN
  ref_CXXFLAGS=$CXXFLAGS
  CXXFLAGS=$CXXGENFLAGS
  
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
    ],
    [
     CXXGEN_SUPPORTS_CPP11=no
     AC_MSG_RESULT([no])
    ]
  )
  AC_SUBST(CXXGEN_SUPPORTS_CPP11)
  # revert CXX and CXXFLAGS
  CXX=$ref_CXX
  CXXFLAGS=$ref_CXXFLAGS

  AC_LANG_RESTORE  
])

# borrowed shamelessly from MADNESS (GPLv3)
# written by Justus Calvin (justus.c79@gmail.com)

AC_DEFUN([ACX_CHECK_SHARED_PTR], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for shared_ptr in std namespace
  AC_MSG_CHECKING([CXX for shared_ptr])
  acx_shared_ptr=no
  
  # Check for std::shared_ptr in <memory>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <memory>]],
        [[std::shared_ptr<int> p;]]
      )
    ],
    [
      AC_DEFINE([LIBINT_USE_MEMORY],[1],[define if Libint is using <memory>.])
      AC_DEFINE([LIBINT_HAS_STD_SHARED_PTR],[1],[define if std::shared_ptr is available.])
      acx_shared_ptr=yes
    ]
  )
  
  # Check for std::tr1::shared_ptr in <memory> unless 
  # disabled by user
  if test "$enable_cpptr1" = yes; then
    if test "$acx_shared_ptr" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <memory>]],
            [[std::tr1::shared_ptr<int> p;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_MEMORY],[1],[define if Libint is using <memory>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
          acx_shared_ptr=yes
        ]
      )
    fi

    # Check for std::tr1::shared_ptr in <tr1/memory>
    if test "$acx_shared_ptr" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <tr1/memory>]],
            [[std::tr1::shared_ptr<int> p;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_TR1_MEMORY],[1],[define if Libint is using <tr1/memory>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
          acx_shared_ptr=yes
        ]
      )
    fi
  fi
  
  # Check if we should use boost tr1 memory
  if test "$acx_with_boost$acx_shared_ptr" = yesno; then
    AC_DEFINE([LIBINT_USE_BOOST_TR1_MEMORY_HPP],[1],[define if Libint is using <boost/tr1/memory.hpp>.])
    AC_DEFINE([LIBINT_HAS_STD_TR1_SHARED_PTR],[1],[define if std::tr1::shared_ptr is available.])
    acx_shared_ptr=yes
  fi
  
  # post shared_ptr results
  AC_MSG_RESULT([$acx_shared_ptr])
  if test "$acx_shared_ptr" = no; then
    AC_MSG_WARN([std::shared_ptr and std::tr1::shared_ptr not found ... using default implementation])
  fi
  
  #Check for std::make_shared and std::allocate_shared
  acx_std_make_shared=no
  AC_MSG_CHECKING([CXX for std::make_shared and std::allocate_shared])
  
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <memory>]],
        [[using ::std::make_shared; using ::std::allocate_shared;]]
      )
    ],
    [
      AC_DEFINE([LIBINT_HAS_STD_MAKE_SHARED],[1],
        [define if std::make_shared and std::allocate_shared are available.])
      acx_std_make_shared=yes
    ]
  )
    
  # post make_shared results
  AC_MSG_RESULT([$acx_std_make_shared])
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_TYPE_TRAITS], [
  AC_LANG_SAVE
  AC_LANG([C++])

  # Check for type traits in <type_traits> and std namespace
  AC_MSG_CHECKING([CXX for type_traits])
  acx_type_traits=no

  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <type_traits>]],
        [[typedef std::is_same<int, double> sameT;]]
      )
    ],
    [
      AC_DEFINE([LIBINT_USE_TYPE_TRAITS],[1],[define if Libint is using <type_traits>.])
      AC_DEFINE([LIBINT_HAS_STD_TYPE_TRAITS],[1],[define if std type traits are available.])
      acx_type_traits=yes
    ]
  )

  # look in tr1 unless disabled by the user
  if test "$enable_cpptr1" = yes; then
    if test "$acx_type_traits" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <type_traits>]],
            [[typedef std::tr1::is_same<int, double> sameT;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_TYPE_TRAITS],[1],[define if Libint is using <type_traits>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
          acx_type_traits=yes
        ]
      )
    fi

    if test "$acx_type_traits" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <tr1/type_traits>]],
            [[typedef std::tr1::is_same<int, double> sameT;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_TR1_TYPE_TRAITS],[1],[define if Libint is using <tr1/type_traits>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
          acx_type_traits=yes
        ]
      )
    fi
  fi
  
  # Check if we should use boost tr1 type_traits
  if test "$acx_with_boost$acx_type_traits" = yesno; then
    AC_DEFINE([LIBINT_USE_BOOST_TR1_TYPE_TRAITS_HPP],[1],[define if Libint is using <boost/tr1/type_traits.hpp>.])
    AC_DEFINE([LIBINT_HAS_STD_TR1_TYPE_TRAITS],[1],[define if std::tr1 type traits are available.])
    acx_type_traits=yes
  fi
  
  # post type traits results
  AC_MSG_RESULT([$acx_type_traits])
  if test "$acx_type_traits" = no; then
    AC_MSG_WARN([std or std::tr1 type traits not found ... using default implentation])
  fi
  
  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_ARRAY], [
  AC_LANG_SAVE
  AC_LANG([C++])
  
  # Check for array in std namespace
  AC_MSG_CHECKING([CXX for array])
  acx_array=no
  
  # Check for std::array in <array>
  AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM(
        [[#include <array>]],
        [[std::array<int,10> a;]]
      )
    ],
    [
      AC_DEFINE([LIBINT_USE_ARRAY],[1],[define if Libint is using <array>.])
      AC_DEFINE([LIBINT_HAS_STD_ARRAY],[1],[define if std::array is available.])
      AC_DEFINE([LIBINT_ARRAY_HAS_FILL],[1],[define if array has fill member function.])
      acx_array=yes
    ]
  )
  
  # Check for std::tr1::array in <array>
  if test "$enable_cpptr1" = yes; then
    if test "$acx_array" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <array>]],
            [[std::tr1::array<int,10> a;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_ARRAY],[1],[define if Libint is using <array>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])

          # Check to see if array has fill function
          AC_COMPILE_IFELSE(
            [
              AC_LANG_PROGRAM(
                [[#include <array>]],
                [[std::tr1::array<int,10> a; a.fill(0);]]
              )
            ],
            [AC_DEFINE([LIBINT_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
          )
          acx_array=yes
        ]
      )
    fi

    # Check for std::tr1::array in <tr1/array>
    if test "$acx_array" = no; then
      AC_COMPILE_IFELSE(
        [
          AC_LANG_PROGRAM(
            [[#include <tr1/array>]],
            [[std::tr1::array<int,10> a;]]
          )
        ],
        [
          AC_DEFINE([LIBINT_USE_TR1_ARRAY],[1],[define if Libint is using <tr1/array>.])
          AC_DEFINE([LIBINT_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])        
          # Check to see if array has fill function
          AC_COMPILE_IFELSE(
            [
              AC_LANG_PROGRAM(
                [[#include <tr1/array>]],
                [[std::tr1::array<int,10> a; a.fill(0);]]
              )
            ],
            [AC_DEFINE([LIBINT_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
          )
          acx_array=yes
        ]
      )
    fi
  fi
  
  # Check if we should use boost tr1 array
  if test "$acx_with_boost$acx_array" = yesno; then
    AC_DEFINE([LIBINT_USE_BOOST_TR1_ARRAY_HPP],[1],[define if Libint is using <boost/tr1/array.hpp>.])
    AC_DEFINE([LIBINT_HAS_STD_TR1_ARRAY],[1],[define if std::tr1::array is available.])
    # Check to see if array has fill function
    AC_COMPILE_IFELSE(
      [
        AC_LANG_PROGRAM(
          [[#include <tr1/array>]],
          [[std::tr1::array<int,10> a; a.fill(0);]]
        )
      ],
      [AC_DEFINE([LIBINT_ARRAY_HAS_FILL],[1],[define if array has fill member function.])]
    )
    acx_array=yes
  fi
  
  # post array results
  AC_MSG_RESULT([$acx_array])
  if test "$acx_array" = no; then
    AC_MSG_WARN([std::array or std::tr1::array not supported ... using default implementation])
  fi

  AC_LANG_RESTORE
])

AC_DEFUN([ACX_CHECK_CPP11],
[
  ACX_CHECK_CPP11_GENERAL
  ACX_CHECK_SHARED_PTR
  ACX_CHECK_TYPE_TRAITS
  ACX_CHECK_ARRAY
])
