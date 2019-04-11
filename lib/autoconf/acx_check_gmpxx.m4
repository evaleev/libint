#
# check gmp c++ interface
#
AC_DEFUN([ACX_CHECK_GMPXX], [
 AC_LANG_SAVE
 AC_LANG([C++])

 acx_have_gmpxx=no
 AC_CHECK_HEADER(gmpxx.h,[
  AC_MSG_CHECKING([if GMP C++ library is usable])
  LIBS="-lgmpxx -lgmp $LIBS"
  AC_LINK_IFELSE(
   [AC_LANG_PROGRAM(
    [[#include <cstddef>
      #include <gmpxx.h>
    ]],
    [[
     mpf_class a = 1.0;
    ]]
    )
   ],
   [AC_MSG_RESULT([yes])
    acx_have_gmpxx=yes
   ],
   [
    AC_MSG_RESULT([no])
   ]
  )

  ]
 )
 
 AC_LANG_RESTORE
]
) # ACX_CHECK_GMPXX
