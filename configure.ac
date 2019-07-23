define([libint_mmm_version],[2.6.0])
define([libint_buildid],[])
define([libint_so_version],[2:3:0])

dnl --------- Begin ---------
dnl Process this file with autoconf to produce a configure script.
define([AC_CACHE_LOAD], )dnl for debugging configure.in
define([AC_CACHE_SAVE], )dnl for debugging configure.in
AC_PREREQ([2.68])
AC_INIT([Libint2],[libint_mmm_version],[libint@valeyev.net],[libint2])
AC_CONFIG_SRCDIR(src/bin/libint/rr.h)
AC_CONFIG_HEADERS(include/libint2/config.h)
AC_CONFIG_AUX_DIR(bin)
AC_CONFIG_MACRO_DIR([lib/autoconf])

dnl --- Compiling C++ only ---
AC_LANG(C++)

AC_CANONICAL_TARGET
AC_DEFINE_UNQUOTED(LIBINT_HOST_ARCH, "$host")
AC_DEFINE_UNQUOTED(LIBINT_TARGET_ARCH, "$target")

define([default_prefix_dash],ifelse(libint_buildid, ,[],[-]))
define([default_prefix],"/usr/local")

LIBINT_MMM_VERSION=libint_mmm_version
LIBINT_BUILDID=libint_buildid
AC_SUBST(LIBINT_MMM_VERSION)
AC_SUBST(LIBINT_BUILDID)
changequote(<<, >>)dnl
LIBINT_MAJOR_VERSION=`echo $LIBINT_MMM_VERSION|sed 's/\([0-9]*\)\.[0-9]*\.[0-9]*/\1/'`
LIBINT_MINOR_VERSION=`echo $LIBINT_MMM_VERSION|sed 's/[0-9]*\.\([0-9]*\)\.[0-9]*/\1/'`
LIBINT_MICRO_VERSION=`echo $LIBINT_MMM_VERSION|sed 's/[0-9]*\.[0-9]*\.\([0-9]*\)/\1/'`
changequote([, ])dnl
AC_DEFINE_UNQUOTED(LIBINT_MAJOR_VERSION,$LIBINT_MAJOR_VERSION)
AC_DEFINE_UNQUOTED(LIBINT_MINOR_VERSION,$LIBINT_MINOR_VERSION)
AC_DEFINE_UNQUOTED(LIBINT_MICRO_VERSION,$LIBINT_MICRO_VERSION)
AC_SUBST(LIBINT_MAJOR_VERSION)
AC_SUBST(LIBINT_MINOR_VERSION)
AC_SUBST(LIBINT_MICRO_VERSION)


LIBINT_SO_VERSION=libint_so_version
AC_SUBST(LIBINT_SO_VERSION)

dnl --------- Features ---------

AC_ARG_WITH(api-prefix,
AS_HELP_STRING([--with-api-prefix],[Prepend this string to every name in the library API (except for the types).]),
[
AC_DEFINE_UNQUOTED(LIBINT_API_PREFIX,"$withval")
echo Using API prefix: $withval
])

LIBINT_CGSHELL_ORDERING_STANDARD=1
LIBINT_CGSHELL_ORDERING_INTV3=2
LIBINT_CGSHELL_ORDERING_GAMESS=3
LIBINT_CGSHELL_ORDERING_ORCA=4
LIBINT_CGSHELL_ORDERING_BAGEL=5
libint_cgshell_ordering=standard
AC_ARG_WITH(cartgauss-ordering,
[  --with-cartgauss-ordering
                         Use one of the following known orderings for shells of cartesian Gaussians:
                           standard -- standard ordering (xxx, xxy, xxz, xyy, xyz, xzz, yyy, ...)
                           intv3  -- intv3 ordering (yyy, yyz, yzz, zzz, xyy, xyz, xzz, xxy, xxz, xxx)
                           gamess -- GAMESS ordering (xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz)
                           orca -- ORCA ordering (hydrid between GAMESS and standard)
                           bagel -- axis-permuted version of intv3 (xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz)
                         The default is standard.],
[
libint_cgshell_ordering=$withval
])
AC_MSG_RESULT([Using $libint_cgshell_ordering ordering of cartesian Gaussians])
case $libint_cgshell_ordering in
  standard)
    libint_cgshell_ordering=$LIBINT_CGSHELL_ORDERING_STANDARD
  ;;
  intv3)
    libint_cgshell_ordering=$LIBINT_CGSHELL_ORDERING_INTV3
  ;;
  gamess)
    libint_cgshell_ordering=$LIBINT_CGSHELL_ORDERING_GAMESS
  ;;
  orca)
    libint_cgshell_ordering=$LIBINT_CGSHELL_ORDERING_ORCA
  ;;
  bagel)
    libint_cgshell_ordering=$LIBINT_CGSHELL_ORDERING_BAGEL
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --with-cartgauss-ordering ($withval)])
  ;;
esac
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING,$libint_cgshell_ordering)
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING_STANDARD,$LIBINT_CGSHELL_ORDERING_STANDARD)
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING_INTV3,$LIBINT_CGSHELL_ORDERING_INTV3)
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING_GAMESS,$LIBINT_CGSHELL_ORDERING_GAMESS)
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING_ORCA,$LIBINT_CGSHELL_ORDERING_ORCA)
AC_DEFINE_UNQUOTED(LIBINT_CGSHELL_ORDERING_BAGEL,$LIBINT_CGSHELL_ORDERING_BAGEL)

LIBINT_SHGSHELL_ORDERING_STANDARD=1
LIBINT_SHGSHELL_ORDERING_GAUSSIAN=2
libint_shgshell_ordering=standard
AC_ARG_WITH(shgauss-ordering,
[  --with-shgauss-ordering
                         Use one of the following known orderings for shells of solid harmonic Gaussians:
                           standard -- standard ordering (-l, -l+1 ... l)
                           gaussian -- the Gaussian ordering (0, 1, -1, 2, -2, ... l, -l)
                         The default is standard.],
[
libint_shgshell_ordering=$withval
])
AC_MSG_RESULT([Using $libint_shgshell_ordering ordering of solid harmonic Gaussians])
case $libint_shgshell_ordering in
  standard)
    libint_shgshell_ordering=$LIBINT_SHGSHELL_ORDERING_STANDARD
  ;;
  gaussian)
    libint_shgshell_ordering=$LIBINT_SHGSHELL_ORDERING_GAUSSIAN
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --with-solidharmgauss-ordering ($withval)])
  ;;
esac
AC_DEFINE_UNQUOTED(LIBINT_SHGSHELL_ORDERING,$libint_shgshell_ordering)
AC_DEFINE_UNQUOTED(LIBINT_SHGSHELL_ORDERING_STANDARD,$LIBINT_SHGSHELL_ORDERING_STANDARD)
AC_DEFINE_UNQUOTED(LIBINT_SHGSHELL_ORDERING_GAUSSIAN,$LIBINT_SHGSHELL_ORDERING_GAUSSIAN)

LIBINT_SHELL_SET_STANDARD=1   AC_SUBST(LIBINT_SHELL_SET_STANDARD)
LIBINT_SHELL_SET_ORCA=2      AC_SUBST(LIBINT_SHELL_SET_ORCA)
libint_shell_set=standard
AC_ARG_WITH(shell-set,
[  --with-shell-set       The library will support computation of shell sets sets subject to these restrictions:
                           standard -- standard ordering:
                             for (ab|cd):
                               l(a) >= l(b),
                               l(c) >= l(d),
                               l(a)+l(b) <= l(c)+l(d)
                             for (b|cd):
                               l(c) >= l(d)
                           orca -- ORCA ordering:
                             for (ab|cd):
                               l(a) <= l(b),
                               l(c) <= l(d),
                               l(a) < l(c) || (l(a) == l(c) && l(b) < l(d))
                             for (b|cd):
                               l(c) <= l(d)
                          The default is standard.],
[
libint_shell_set=$withval
])
AC_MSG_RESULT([Will generate $libint_shell_set shell quartet sets])
case $libint_shell_set in
  standard)
    libint_shell_set=$LIBINT_SHELL_SET_STANDARD
  ;;
  orca)
    libint_shell_set=$LIBINT_SHELL_SET_ORCA
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --with-shell-set ($withval)])
  ;;
esac
LIBINT_SHELL_SET=$libint_shell_set AC_SUBST(LIBINT_SHELL_SET)
AC_DEFINE_UNQUOTED(LIBINT_SHELL_SET,$libint_shell_set)
AC_DEFINE_UNQUOTED(LIBINT_SHELL_SET_STANDARD,$LIBINT_SHELL_SET_STANDARD)
AC_DEFINE_UNQUOTED(LIBINT_SHELL_SET_ORCA,$LIBINT_SHELL_SET_ORCA)

AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],[Compile with debugging options]),
[
case $enableval in
  yes)
    DEBUG=yes
  ;;
  no)
    DEBUG=no
  ;;
  opt)
    DEBUG=opt
  ;;
  *)
    AC_MSG_ERROR([Invalid value for --enable-debug ($enableval)])
  ;;
esac
],[
    DEBUG=no
]
)

# can be used in the future to exclude source directories
EXCLUDED_DIRS=

AC_ARG_ENABLE(1body,
AS_HELP_STRING([--enable-1body=N],[Compile with support for up to N-th derivatives of 1-body integrals]),
[
case $enableval in
  yes)
    INCLUDE_ONEBODY=0
  ;;
  no)
    INCLUDE_ONEBODY=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-1body ($enableval)])
    fi
    INCLUDE_ONEBODY=$enableval
  ;;
esac
],[
    INCLUDE_ONEBODY=0
]
)

DISABLE_ONEBODY_PROPERTY_DERIVS=1
AC_ARG_ENABLE(1body-property-derivs,
AS_HELP_STRING([--enable-1body-property-derivs],[Enable geometric derivatives of 1-body property integrals (all but overlap, kinetic, elecpot). These derivatives are disabled by default to save compile time.]),
[
case $enableval in
  yes)
    DISABLE_ONEBODY_PROPERTY_DERIVS=0
  ;;
esac
])
if test X$DISABLE_ONEBODY_PROPERTY_DERIVS != 0; then
  AC_DEFINE(DISABLE_ONEBODY_PROPERTY_DERIVS)
fi

AC_ARG_ENABLE(eri,
AS_HELP_STRING([--enable-eri=N],[Compile with support for up to N-th derivatives of electron repulsion integrals]),
[
case $enableval in
  yes)
    INCLUDE_ERI=0
  ;;
  no)
    INCLUDE_ERI=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-eri ($enableval)])
    fi
    INCLUDE_ERI=$enableval
  ;;
esac
],[
    INCLUDE_ERI=0
]
)

AC_ARG_ENABLE(eri3,
AS_HELP_STRING([--enable-eri3=N],[Compile with support for up to N-th derivatives of 3-center electron repulsion integrals]),
[
case $enableval in
  yes)
    INCLUDE_ERI3=0
  ;;
  no)
    INCLUDE_ERI3=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-eri3 ($enableval)])
    fi
    INCLUDE_ERI3=$enableval
  ;;
esac
],[
    INCLUDE_ERI3=no
]
)

AC_ARG_ENABLE(eri2,
AS_HELP_STRING([--enable-eri2=N],[Compile with support for up to N-th derivatives of 2-center electron repulsion integrals]),
[
case $enableval in
  yes)
    INCLUDE_ERI2=0
  ;;
  no)
    INCLUDE_ERI2=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-eri2 ($enableval)])
    fi
    INCLUDE_ERI2=$enableval
  ;;
esac
],[
    INCLUDE_ERI2=no
]
)

AC_ARG_ENABLE(g12,
AS_HELP_STRING([--enable-g12=N],[Compile with support for N-th derivatives of MP2-F12 energies with Gaussian factors]),
[
case $enableval in
  yes)
    INCLUDE_G12=0
  ;;
  no)
    INCLUDE_G12=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-g12 ($enableval)])
    fi
    INCLUDE_G12=$enableval
  ;;
esac
],[
    INCLUDE_G12=no
]
)

AC_ARG_ENABLE(g12dkh,
AS_HELP_STRING([--enable-g12dkh=N],[Compile with support for N-th derivatives of DKH-MP2-F12 energies with Gaussian factors]),
[
case $enableval in
  yes)
    INCLUDE_G12DKH=0
  ;;
  no)
    INCLUDE_G12DKH=no
  ;;
  *)
    if test $enableval -lt 0; then
      AC_MSG_ERROR([Invalid value for --enable-g12dkh ($enableval)])
    fi
    INCLUDE_G12DKH=$enableval
  ;;
esac
],[
    INCLUDE_G12DKH=no
]
)

LIBINT_MAX_AM=4
AC_ARG_WITH(max-am,
AS_HELP_STRING([--with-max-am=N],[Support Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      LIBINT_MAX_AM_LIST="$withval"
      AC_SUBST(LIBINT_MAX_AM_LIST)
      AC_DEFINE_UNQUOTED(LIBINT_MAX_AM_LIST, "$LIBINT_MAX_AM_LIST")
      LIBINT_MAX_AM=`echo $LIBINT_MAX_AM_LIST | tr , "\n" | sort -n | tail -n1`
    ;;
    *)
      if test $withval -le 0; then
        AC_MSG_ERROR([Invalid value for --with-max-am ($withval)])
      else
        LIBINT_MAX_AM=$withval
      fi
      if test $LIBINT_MAX_AM -ge 8; then
        AC_MSG_ERROR([Value for --with-max-am too high ($withval). Are you sure you know what you are doing?])
      fi
    ;;
  esac
]
)
AC_DEFINE_UNQUOTED(LIBINT_MAX_AM,$LIBINT_MAX_AM)

LIBINT_OPT_AM=$((($LIBINT_MAX_AM/2)+1))
AC_ARG_WITH(opt-am,
AS_HELP_STRING([--with-opt-am=N],[Optimize maximally for up to angular momentum N (N <= max-am).]),
[
  case $withval in 
    *,*)
      LIBINT_OPT_AM_LIST="$withval"
      AC_SUBST(LIBINT_OPT_AM_LIST)
      AC_DEFINE_UNQUOTED(LIBINT_OPT_AM_LIST, "$LIBINT_OPT_AM_LIST")
      LIBINT_OPT_AM=`echo $LIBINT_OPT_AM_LIST | tr , "\n" | sort -n | tail -n1`
    ;;
    *)
      if test $withval -le $LIBINT_MAX_AM; then
        LIBINT_OPT_AM=$withval
      fi
    ;;
  esac
]
)
AC_DEFINE_UNQUOTED(LIBINT_OPT_AM,$LIBINT_OPT_AM)

ERI_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(eri-max-am,
AS_HELP_STRING([--with-eri-max-am=N],[Support ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI_MAX_AM_LIST="$withval"
      AC_SUBST(ERI_MAX_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI_MAX_AM_LIST, "$ERI_MAX_AM_LIST")
    ;;
    *) 
     if test $withval -le 0; then
      AC_MSG_ERROR([Invalid value for --with-eri-max-am ($withval)])
     else
      ERI_MAX_AM=$withval
     fi
     if test $ERI_MAX_AM -ge 8; then
      AC_MSG_ERROR([Value for --with-eri-max-am too high ($withval). Are you sure you know what you are doing?])
     fi
     AC_SUBST(ERI_MAX_AM)
     AC_DEFINE_UNQUOTED(ERI_MAX_AM,$ERI_MAX_AM)
    ;; 
  esac
]
)

ERI_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(eri-opt-am,
AS_HELP_STRING([--with-eri-opt-am=N],[Optimize ERIs maximally for up to angular momentum N (N <= max-am).  Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI_OPT_AM_LIST="$withval"
      AC_SUBST(ERI_OPT_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI_OPT_AM_LIST, "$ERI_OPT_AM_LIST")  
    ;;
    *) 
     if test $withval -le $ERI_MAX_AM; then
       ERI_OPT_AM=$withval
       AC_SUBST(ERI_OPT_AM)
       AC_DEFINE_UNQUOTED(ERI_OPT_AM,$ERI_OPT_AM)  
     fi
    ;;
  esac
]
)

if test X$INCLUDE_ERI != Xno; then
  AC_DEFINE_UNQUOTED(INCLUDE_ERI,$INCLUDE_ERI)
  LIBINT_SUPPORTS_ERI=yes
  AC_SUBST(LIBINT_SUPPORTS_ERI)
  LIBINT_ERI_DERIV=$INCLUDE_ERI
  AC_SUBST(LIBINT_ERI_DERIV)  
fi

ONEBODY_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(1body-max-am,
AS_HELP_STRING([--with-1body-max-am=N],[Support 1-e ints for Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ONEBODY_MAX_AM_LIST="$withval"
      AC_SUBST(ONEBODY_MAX_AM_LIST)
      AC_DEFINE_UNQUOTED(ONEBODY_MAX_AM_LIST, "$ONEBODY_MAX_AM_LIST")
    ;;
    *) 
     if test $withval -le 0; then
      AC_MSG_ERROR([Invalid value for --with-1body-max-am ($withval)])
     else
      ONEBODY_MAX_AM=$withval
     fi
     if test $ONEBODY_MAX_AM -ge 8; then
      AC_MSG_ERROR([Value for --with-1body-max-am too high ($withval). Are you sure you know what you are doing?])
     fi
     AC_SUBST(ONEBODY_MAX_AM)
     AC_DEFINE_UNQUOTED(ONEBODY_MAX_AM,$ONEBODY_MAX_AM)
    ;; 
  esac
]
)

ONEBODY_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(1body-opt-am,
AS_HELP_STRING([--with-1body-opt-am=N],[Optimize 1-e ints maximally for up to angular momentum N (N <= max-am).  Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ONEBODY_OPT_AM_LIST="$withval"
      AC_SUBST(ONEBODY_OPT_AM_LIST)
      AC_DEFINE_UNQUOTED(ONEBODY_OPT_AM_LIST, "$ONEBODY_OPT_AM_LIST")  
    ;;
    *) 
     if test $withval -le $ONEBODY_MAX_AM; then
       ONEBODY_OPT_AM=$withval
       AC_SUBST(ONEBODY_OPT_AM)
       AC_DEFINE_UNQUOTED(ONEBODY_OPT_AM,$ONEBODY_OPT_AM)  
     fi
    ;;
  esac
]
)

MULTIPOLE_MAX_ORDER=4
AC_ARG_WITH(multipole-max-order,
AS_HELP_STRING([--with-multipole-max-order=N],[Maximum order of spherical multipole integrals. The default is 4. There is no maximum]),
[
  case $withval in 
    *) 
      MULTIPOLE_MAX_ORDER=$withval
    ;;
  esac
]
)
AC_DEFINE_UNQUOTED(MULTIPOLE_MAX_ORDER,$MULTIPOLE_MAX_ORDER)  

if test X$INCLUDE_ONEBODY != Xno; then
  AC_DEFINE_UNQUOTED(INCLUDE_ONEBODY,$INCLUDE_ONEBODY)
  LIBINT_SUPPORTS_ONEBODY=yes
  AC_SUBST(LIBINT_SUPPORTS_ONEBODY)
  LIBINT_ONEBODY_DERIV=$INCLUDE_ONEBODY
  AC_SUBST(LIBINT_ONEBODY_DERIV)
fi

# check max am for 3-center ERIs, if needed
ERI3_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(eri3-max-am,
AS_HELP_STRING([--with-eri3-max-am=N],[Support 3-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI3_MAX_AM_LIST="$withval"
      AC_SUBST(ERI3_MAX_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI3_MAX_AM_LIST, "$ERI3_MAX_AM_LIST")
    ;;
    *) 
     if test $withval -le 0; then
      AC_MSG_ERROR([Invalid value for --with-eri3-max-am ($withval)])
     else
      ERI3_MAX_AM=$withval
     fi
     if test $ERI3_MAX_AM -ge 8; then
      AC_MSG_ERROR([Value for --with-eri3-max-am too high ($withval). Are you sure you know what you are doing?])
     fi
     AC_SUBST(ERI3_MAX_AM)
     AC_DEFINE_UNQUOTED(ERI3_MAX_AM,$ERI3_MAX_AM)
    ;; 
  esac
]
)

ERI3_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(eri3-opt-am,
AS_HELP_STRING([--with-eri3-opt-am=N],[Optimize 3-center ERIs maximally for up to angular momentum N (N <= max-am).  Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI3_OPT_AM_LIST="$withval"
      AC_SUBST(ERI3_OPT_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI3_OPT_AM_LIST, "$ERI3_OPT_AM_LIST")  
    ;;
    *) 
     if test $withval -le $ERI3_MAX_AM; then
       ERI3_OPT_AM=$withval
       AC_SUBST(ERI3_OPT_AM)
       AC_DEFINE_UNQUOTED(ERI3_OPT_AM,$ERI3_OPT_AM)  
     fi
    ;;
  esac
]
)

AC_ARG_ENABLE(eri3-pure-sh,
AS_HELP_STRING([--enable-eri3-pure-sh],[Assume the "unpaired" center of 3-center electron repulsion integrals will be transformed to pure solid harmonics]),
[
case $enableval in
  yes)
    ERI3_PURE_SH=yes
  ;;
  no)
    ERI3_PURE_SH=no
  ;;
esac
],[
    ERI3_PURE_SH=no
]
)

if test X$INCLUDE_ERI3 != Xno; then
  AC_DEFINE_UNQUOTED(INCLUDE_ERI3,$INCLUDE_ERI3)
  if test X$ERI3_PURE_SH != Xno; then
    AC_DEFINE(ERI3_PURE_SH)
  fi
fi

# check max am for 2-center ERIs, if needed
ERI2_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(eri2-max-am,
AS_HELP_STRING([--with-eri2-max-am=N],[Support 2-center ERIs for Gaussians of angular momentum up to N. Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI2_MAX_AM_LIST="$withval"
      AC_SUBST(ERI2_MAX_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI2_MAX_AM_LIST, "$ERI2_MAX_AM_LIST")
    ;;
    *) 
     if test $withval -le 0; then
      AC_MSG_ERROR([Invalid value for --with-eri2-max-am ($withval)])
     else
      ERI2_MAX_AM=$withval
     fi
     if test $ERI2_MAX_AM -ge 8; then
      AC_MSG_ERROR([Value for --with-eri2-max-am too high ($withval). Are you sure you know what you are doing?])
     fi
     AC_SUBST(ERI2_MAX_AM)
     AC_DEFINE_UNQUOTED(ERI2_MAX_AM,$ERI2_MAX_AM)
    ;; 
  esac
]
)

ERI2_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(eri2-opt-am,
AS_HELP_STRING([--with-eri2-opt-am=N],[Optimize 2-center ERIs maximally for up to angular momentum N (N <= max-am).  Can specify values for each derivative levels as a list N0,N1,N2...]),
[
  case $withval in 
    *,*)
      ERI2_OPT_AM_LIST="$withval"
      AC_SUBST(ERI2_OPT_AM_LIST)
      AC_DEFINE_UNQUOTED(ERI2_OPT_AM_LIST, "$ERI2_OPT_AM_LIST")  
    ;;
    *) 
     if test $withval -le $ERI2_MAX_AM; then
       ERI2_OPT_AM=$withval
       AC_SUBST(ERI2_OPT_AM)
       AC_DEFINE_UNQUOTED(ERI2_OPT_AM,$ERI2_OPT_AM)  
     fi
    ;;
  esac
]
)

AC_ARG_ENABLE(eri2-pure-sh,
AS_HELP_STRING([--enable-eri2-pure-sh],[Assume the 2-center electron repulsion integrals will be transformed to pure solid harmonics]),
[
case $enableval in
  yes)
    ERI2_PURE_SH=yes
  ;;
  no)
    ERI2_PURE_SH=no
  ;;
esac
],[
    ERI2_PURE_SH=no
]
)

if test X$INCLUDE_ERI2 != Xno; then
  AC_DEFINE_UNQUOTED(INCLUDE_ERI2,$INCLUDE_ERI2)
  if test X$ERI2_PURE_SH != Xno; then
    AC_DEFINE(ERI2_PURE_SH)
  fi
fi

# check max am for G12 integrals, if needed
G12_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(g12-max-am,
AS_HELP_STRING([--with-g12-max-am=N],[Support integrals for G12 methods of angular momentum up to N.]),
if test $withval -le 0; then
  AC_MSG_ERROR([Invalid value for --with-g12-max-am ($withval)])
else
  G12_MAX_AM=$withval
fi
)

G12_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(g12-opt-am,
AS_HELP_STRING([--with-g12-opt-am=N],[Optimize G12 integrals for up to angular momentum N (N <= max-am).]),
if test $withval -le $G12_MAX_AM; then
  G12_OPT_AM=$withval
fi
)

AC_ARG_ENABLE(t1g12-support,
AS_HELP_STRING([--disable-t1g12-support],[Libint will not support [Ti,G12] integrals]),
[
case $enableval in
  yes)
    SUPPORT_T1G12=1
  ;;
  no)
    SUPPORT_T1G12=0
  ;;
  *)
    SUPPORT_T1G12=$enableval
  ;;
esac
],[
    SUPPORT_T1G12=1
]
)

AC_ARG_ENABLE(composite-evaluators,
AS_HELP_STRING([--disable-composite-evaluators],[Libint will not use composite evaluators (i.e. every evaluator wil compute one integral type only)]),
[
case $enableval in
  yes)
    LIBINT_USE_COMPOSITE_EVALUATORS=1
    AC_MSG_RESULT([Will generate composite evaluators for multiple integral types])
  ;;
  no)
    LIBINT_USE_COMPOSITE_EVALUATORS=0
    AC_MSG_RESULT([Will generate separate evaluators for each integral type])
  ;;
  *)
    LIBINT_USE_COMPOSITE_EVALUATORS=$enableval
  ;;
esac
],[
    LIBINT_USE_COMPOSITE_EVALUATORS=1
    AC_MSG_RESULT([Will generate composite evaluators for multiple integral types])
]
)
AC_DEFINE_UNQUOTED(LIBINT_USE_COMPOSITE_EVALUATORS,$LIBINT_USE_COMPOSITE_EVALUATORS)

if test "X$INCLUDE_G12" != "Xno"; then
  AC_DEFINE_UNQUOTED(INCLUDE_G12,$INCLUDE_G12)
  AC_SUBST(G12_MAX_AM)
  AC_SUBST(G12_OPT_AM)
  AC_SUBST(SUPPORT_T1G12)
  AC_DEFINE_UNQUOTED(G12_MAX_AM,$G12_MAX_AM)
  AC_DEFINE_UNQUOTED(G12_OPT_AM,$G12_OPT_AM)
  AC_DEFINE_UNQUOTED(SUPPORT_T1G12,$SUPPORT_T1G12)
fi

# check max am for G12DKH integrals, if needed
G12DKH_MAX_AM=$LIBINT_MAX_AM
AC_ARG_WITH(g12dkh-max-am,
AS_HELP_STRING([--with-g12dkh-max-am=N],[Support integrals for relativistic G12 methods of angular momentum up to N.]),
if test $withval -le 0; then
  AC_MSG_ERROR([Invalid value for --with-g12dkh-max-am ($withval)])
else
  G12DKH_MAX_AM=$withval
fi
)

G12DKH_OPT_AM=$LIBINT_OPT_AM
AC_ARG_WITH(g12dkh-opt-am,
AS_HELP_STRING([--with-g12dkh-opt-am=N],[Optimize G12DKH integrals for up to angular momentum N (N <= max-am).]),
if test $withval -le $G12DKH_MAX_AM; then
  G12DKH_OPT_AM=$withval
fi
)

if test "X$INCLUDE_G12DKH" != "Xno"; then
  AC_DEFINE_UNQUOTED(INCLUDE_G12DKH,$INCLUDE_G12DKH)
  AC_SUBST(G12DKH_MAX_AM)
  AC_SUBST(G12DKH_OPT_AM)
  AC_DEFINE_UNQUOTED(G12DKH_MAX_AM,$G12DKH_MAX_AM)
  AC_DEFINE_UNQUOTED(G12DKH_OPT_AM,$G12DKH_OPT_AM)
fi

LIBINT_UNROLLING_THRESHOLD=100
AC_ARG_ENABLE(unrolling,
AS_HELP_STRING([--enable-unrolling],[Unroll shell sets into integrals (--enable-unrolling=S will unroll shell sets larger than S).]),
[
case $enableval in
  yes)
  LIBINT_UNROLLING_THRESHOLD=1000000000
  ;;
  no)
  LIBINT_UNROLLING_THRESHOLD=
  AC_MSG_RESULT([Shell sets will not be unrolled])
  ;;
  *)
  LIBINT_UNROLLING_THRESHOLD=$enableval
  ;;
esac
])
if test X$LIBINT_UNROLLING_THRESHOLD != X; then
  if test X$LIBINT_UNROLLING_THRESHOLD = X1000000000; then
    AC_MSG_RESULT([All shell sets will be unrolled])
  else
    AC_MSG_RESULT([Shell sets will be unrolled up to size $LIBINT_UNROLLING_THRESHOLD])
  fi
  AC_DEFINE_UNQUOTED(LIBINT_ENABLE_UNROLLING, $LIBINT_UNROLLING_THRESHOLD)
fi

AC_ARG_ENABLE(generic-code,
AS_HELP_STRING([--enable-generic-code],[Use manually-written generic code.]),
[
case $enableval in
  yes)
  AC_DEFINE(LIBINT_ENABLE_GENERIC_CODE)
  AC_MSG_RESULT([Generic RR code will be used])
  ;;
  no)
  AC_MSG_RESULT([Generic RR code will not be used])
  ;;
esac
])

LIBINT_VECTOR_LENGTH=1
AC_ARG_WITH(vector-length,
AS_HELP_STRING([--with-vector-length=N],[Compute integrals in vectors of length N.]),
[
if test $withval -gt 0; then
  LIBINT_VECTOR_LENGTH=$withval
  AC_DEFINE_UNQUOTED(LIBINT_VECTOR_LENGTH,$withval)
  AC_MSG_RESULT([Using vector length: $withval])
fi
])

if test $LIBINT_VECTOR_LENGTH -gt 1; then
  LIBINT_VECTOR_METHOD=block
  AC_ARG_WITH(vector-method,
  AS_HELP_STRING([--with-vector-method],[Specifies how to vectorize integrals. Allowed values are 'block' (default), and 'line'.]),
  [
  LIBINT_VECTOR_METHOD=$withval
  ])
  case $LIBINT_VECTOR_METHOD in
    block)
    ;;
    line)
    ;;
    *)
    AC_MSG_ERROR([Unrecognized vector method: $LIBINT_VECTOR_METHOD])
    ;;
  esac
  AC_DEFINE_UNQUOTED(LIBINT_VECTOR_METHOD,"$LIBINT_VECTOR_METHOD")
  AC_MSG_RESULT([Using vector method: $LIBINT_VECTOR_METHOD])
fi

AC_ARG_ENABLE(single-evaltype,
 AS_HELP_STRING([--enable-single-evaltype],[Generate single evaluator type (i.e. all tasks use the same evaluator)]),
 [
  if test X$enableval = Xyes; then
    enable_single_evaltype=1
  else
    enable_single_evaltype=0
  fi
 ],
 [
  enable_single_evaltype=1
 ]
)
if test X$enable_single_evaltype = X1; then
  AC_DEFINE(LIBINT_SINGLE_EVALTYPE)
  AC_MSG_RESULT([Will generate single evaluator type])
else
  AC_MSG_RESULT([Will generate specialized evaluator types for each task])
  AC_MSG_ERROR([Cannot generate specialized evaluator types yet])
fi

libint_generate_gma=0
AC_ARG_ENABLE(fma,
AS_HELP_STRING([--enable-fma],[Generate FMA (fused multiply-add) instructions (CXXGEN must support C++11; to benefit must have FMA-capable hardware and compiler)]),
[
case $enableval in
  yes)
    libint_generate_gma=1
    AC_DEFINE(LIBINT_GENERATE_FMA)
    AC_MSG_RESULT([Will generate FMA instructions])
  ;;
  no)
    AC_MSG_RESULT([Will not generate FMA instructions])
  ;;
esac
],
[
    AC_MSG_RESULT([Will not generate FMA instructions])
])

AC_ARG_ENABLE(accum-ints,
AS_HELP_STRING([--enable-accum-ints],[Accumulate integrals to the buffer, rather than copy.]),
[
case $enableval in
  yes)
    AC_DEFINE(LIBINT_ACCUM_INTS)
    AC_MSG_RESULT([Will accumulate integrals to the buffer])
  ;;
  no)
    AC_MSG_RESULT([Will copy integrals to the buffer])
  ;;
esac
],
[
AC_MSG_RESULT([Will copy integrals to the buffer])
])

AC_ARG_ENABLE(flop-counter,
AS_HELP_STRING([--enable-flop-counter],[Support (approximate) FLOP counting by the library. CXXGEN must support C++11!]),
[
case $enableval in
  yes)
    AC_DEFINE(LIBINT_FLOP_COUNT)
    AC_MSG_RESULT([FLOP counter will be supported])
  ;;
  no)
    AC_MSG_RESULT([FLOP counter will not be supported])
  ;;
esac
],
[
    AC_MSG_RESULT([FLOP counter will not be supported])
])

AC_ARG_ENABLE(profile,
AS_HELP_STRING([--enable-profile],[Turn on profiling instrumentation of the library. CXXGEN must support C++11!]),
[
case $enableval in
  yes)
    AC_DEFINE(LIBINT_PROFILE)
    AC_MSG_RESULT([Profiling instrumentation enabled (generated code will require C++11)])
  ;;
  no)
    AC_MSG_RESULT([Profiling instrumentation disabled])
  ;;
esac
],
[
    AC_MSG_RESULT([Profiling instrumentation disabled])
])

acx_contracted_ints=yes
AC_ARG_ENABLE(contracted-ints,
AS_HELP_STRING([--enable-contracted-ints],[Turn on support for contracted integrals.]),
[
case $enableval in
  no)
    acx_contracted_ints=no
  ;;
esac
])
if test "X$acx_contracted_ints" = Xyes; then
  AC_DEFINE(LIBINT_CONTRACTED_INTS)
  AC_MSG_RESULT([Will support contracted integrals])
else
  AC_MSG_RESULT([Will not support contracted integrals])
fi
AC_SUBST(LIBINT_CONTRACTED_INTS, $acx_contracted_ints)


dnl
dnl Preferred strategies for evaluation
dnl
AC_ARG_WITH(eri-strategy,
AS_HELP_STRING([--with-eri-strategy],[Compute ERIs using the following strategy. This option is for experts ONLY.]),
[
  case $withval in
    OS)
      AC_DEFINE_UNQUOTED(LIBINT_ERI_STRATEGY,0)
    ;;
    HGP)
      AC_DEFINE_UNQUOTED(LIBINT_ERI_STRATEGY,1)
    ;;
    HL)
      AC_DEFINE_UNQUOTED(LIBINT_ERI_STRATEGY,2)
    ;;
  esac
AC_MSG_RESULT([Using ERI strategy: $withval])
],
[AC_DEFINE(LIBINT_ERI_STRATEGY,1)]
)


BUILDID="libint_buildid"
AC_ARG_WITH(build-id,
AS_HELP_STRING([--with-build-id],[Gives an identifier for the build.]),
BUILDID=$withval
)
AC_SUBST(BUILDID)
AC_DEFINE_UNQUOTED(LIBINT_BUILDID,"$BUILDID")

if test "$BUILDID"; then
  LIBINT_VERSION="$LIBINT_MMM_VERSION-$BUILDID"
else
  LIBINT_VERSION="$LIBINT_MMM_VERSION"
fi

AC_DEFINE_UNQUOTED(LIBINT_VERSION, "$LIBINT_VERSION")
AC_SUBST(LIBINT_VERSION)

ac_default_prefix="/usr/local/libint/$LIBINT_VERSION"

AC_ARG_WITH(cxx,
AS_HELP_STRING([--with-cxx],[Gives the name of the C++ compiler to use.]),
CXX=$withval
)

CXXFLAGS_OPT=
AC_ARG_WITH(cxx-optflags,
AS_HELP_STRING([--with-cxx-optflags],[Gives optimization flags for the C++ compiler.]),
if test "X$withval" != "Xyes" -a "X$withval" != "Xno"; then
  CXXFLAGS_OPT=$withval
fi
)

AC_ARG_WITH(cxxgen,
AS_HELP_STRING([--with-cxxgen],[Gives the name of the C++ compiler to compile generated code.]),
CXXGEN=$withval
)

CXXGENFLAGS_OPT=
AC_ARG_WITH(cxxgen-optflags,
AS_HELP_STRING([--with-cxxgen-optflags],[Gives optimization flags for the C++ compiler for generated source.]),
if test "X$withval" != "Xyes" -a "X$withval" != "Xno"; then
  CXXGENFLAGS_OPT=$withval
fi
)

AC_ARG_WITH(cxxdepend,
AS_HELP_STRING([--with-cxxdepend],[Gives the name of the C++ compiler with which to extract dependency information.]),
CXXDEPEND=$withval
)

AC_ARG_WITH(ranlib,
AS_HELP_STRING([--with-ranlib],[Gives the name of the ranlib program.]),
RANLIB=$withval
)

AC_ARG_WITH(ar,
AS_HELP_STRING([--with-ar],[Names the archive creator.]),
AR=$withval
)

ARFLAGS=r
AC_ARG_WITH(ar-flags,
AS_HELP_STRING([--with-ar-flags],[Flags for the the archive creator.]),
ARFLAGS=$withval
)
AC_SUBST(ARFLAGS)

AC_ARG_WITH(ld,
AS_HELP_STRING([--with-ld],[Names the object linker.]),
LD=$withval
)

libintincludedir=$includedir
AC_ARG_WITH(libint-includedir,
AS_HELP_STRING([--with-libint-includedir],[Specifies include file install subdir.]),
libintincludedir=$withval
)
AC_SUBST(libintincludedir)

EXPORTDIR=libint-$LIBINT_VERSION
AC_ARG_WITH(libint-exportdir,
AS_HELP_STRING([--with-libint-exportdir],[Specifies export directory name. Default is libint-<version>.]),
EXPORTDIR=$withval
)
AC_SUBST(EXPORTDIR)

AC_ARG_WITH(incdirs,
AS_HELP_STRING([--with-incdirs],[Specifies include directories (-Idir1 -Idir2).]),
EXTRAINCDIRS=$withval
CPPFLAGS=$withval
echo Using extra include directories: $withval
)

AC_ARG_WITH(libs,
AS_HELP_STRING([--with-libs],[Specifies libraries (-llib1 -llib2).]),
LIBS=$withval
echo Using extra libraries: $withval
)

AC_ARG_WITH(libdirs,
AS_HELP_STRING([--with-libdirs],[Specifies library directories (-Ldir1 -Ldir2).]),
LIBDIRS=$withval
LDFLAGS=$withval
echo Using extra library directories: $withval
)

SCRATCHDIR=tmp
AC_ARG_WITH(scratchdir,
AS_HELP_STRING([--with-scratchdir],[Specifies the location for the machine-generated library source and object files.]),
SCRATCHDIR=$withval
)
AC_SUBST(SCRATCHDIR)

PKGCONFIGDIR='$(libdir)/pkgconfig'
AC_ARG_WITH(pkgconfigdir,
AS_HELP_STRING([--with-pkgconfigdir],[Specifies the location to put pkg-config's data file. Default is $(libdir)/pkgconfig.]),
PKGCONFIGDIR=$withval
)
AC_SUBST(PKGCONFIGDIR)

CMAKEDIR='$(libdir)/cmake/libint2'
AC_ARG_WITH(cmakedir,
AS_HELP_STRING([--with-cmakedir],[Specifies the location to put FindLibint2.cmake file. Default is $(libdir)/cmake/libint2.]),
CMAKEDIR=$withval
)
AC_SUBST(CMAKEDIR)

dnl --------- Want absolute path for srcdir, not relative. ---------

if test X$ac_srcdir_defaulted = Xyes; then
  case $srcdir in
    ..*)
      srcdir=`(cd $srcdir; pwd)`
      #srcdir=`(cd $srcdir; echo $PWD)`
      ;;
  esac
fi

dnl --------- Checks for programs. ---------
AC_PROG_MAKE_SET
AC_PROG_LN_S
AC_PROG_INSTALL

# if not given a C++ compiler
if test "X$CXX" = "X"; then
  AC_CHECK_PROGS(CXX, clang++ icpc xlC_r g++ c++ CC, c++)
fi

# if not given a C++ compiler for generated source, use CXX
if test "X$CXXGEN" = X; then
  CXXGEN=$CXX
fi
AC_SUBST(CXXGEN)

# discard CXXFLAGS provided by the C++ detection macros
initial_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
AC_PROG_CXXCPP
CXXFLAGS="$initial_CXXFLAGS"

AC_CHECK_PROG(AR,ar,ar)
AC_CHECK_PROG(PERL,perl,perl)
AC_CHECK_PROGS(TAR,gnutar tar,tar)
AC_CHECK_PROGS(PYTHON,python python3 python2 python2.7,python)
AC_PATH_PROG(DOXYGEN,doxygen,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(PDFLATEX,pdflatex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(PSLATEX,pslatex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(DVIPS,dvips,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(LATEX,latex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(LATEX2HTML,latex2html,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(BIBTEX,bibtex,,$PATH:/usr/local/bin:/opt/local/bin)

dnl The lack of certain tools is a show stopper

changequote(<<, >>)dnl

if [ X$PERL = X ]; then
  echo "Could not find program: perl"
  exit 1
fi
if [ X$PYTHON = X ]; then
  echo "Could not find program: python"
  exit 1
fi

changequote([, ])dnl

dnl -------- Checks for compiler/linker options. ----------

# set the optimization default if neither CXXFLAGS or --with-*-optflags were given
if test "X$CXXFLAGS_OPT" = X -a "X$CXXFLAGS" = X; then
  CXXFLAGS_OPT=-O2
fi
if test "X$CXXGENFLAGS_OPT" = X -a "X$CXXGENFLAGS" = X; then
  AC_MSG_WARN([C++ optimization options are not given! For optimum performance give CXXGENFLAGS or use --with-cxxgen-optflags configure option])
  CXXGENFLAGS_OPT="-O2 -DNDEBUG"
fi
# options needed only for debugging
CXXFLAGS_DBG=-g
# options that are always needed
CXXFLAGS_MISC=

# now form CXXGENFLAGS ...
if test "X$CXXGENFLAGS" = X; then
  CXXGENFLAGS=$CXXFLAGS
fi
if test $DEBUG = yes; then
 CXXGENFLAGS="$CXXGENFLAGS $CXXFLAGS_DBG $CXXFLAGS_MISC"
 LDFLAGS="$LDFLAGS -g"
elif test $DEBUG = opt; then
 CXXGENFLAGS="$CXXGENFLAGS $CXXFLAGS_DBG $CXXGENFLAGS_OPT $CXXFLAGS_MISC"
 LDFLAGS="$LDFLAGS -g"
else
 CXXGENFLAGS="$CXXGENFLAGS $CXXGENFLAGS_OPT $CXXFLAGS_MISC"
fi

# ... and CXXFLAGS
if test $DEBUG = yes; then
 CXXFLAGS="$CXXFLAGS $CXXFLAGS_DBG $CXXFLAGS_MISC"
 LDFLAGS="$LDFLAGS -g"
elif test $DEBUG = opt; then
 CXXFLAGS="$CXXFLAGS $CXXFLAGS_DBG $CXXFLAGS_OPT $CXXFLAGS_MISC"
 LDFLAGS="$LDFLAGS -g"
else
 CXXFLAGS="$CXXFLAGS $CXXFLAGS_OPT $CXXFLAGS_MISC"
fi

OBJSUF=o
LIBSUF=a

dnl -- check how dependency information is built --

# The GNU compilers work with:
CXXDEPENDSUF=none
CXXDEPENDFLAGS=-M

/bin/rm -f depcheck.u depcheck.c depcheck.o

# Check for an IBM visual age C++ compiler
echo "#include <iostream>" > depcheck.cc
$CXX $CPPFLAGS $CXXFLAGS -M -E depcheck.cc > /dev/null 2>&1
if test -f depcheck.u; then
  CXXDEPENDSUF=u
  CXXDEPENDFLAGS="-M -E"
fi
/bin/rm -f depcheck.u depcheck.c depcheck.o

if test "X$CXXDEPEND" = X; then
  CXXDEPEND=$CXX
fi
AC_SUBST(CXXDEPEND)
AC_SUBST(CXXDEPENDSUF)
AC_SUBST(CXXDEPENDFLAGS)

dnl -- special architecture options --

# nothing so far

AC_SUBST(EXTRAINCLUDE)
AC_SUBST(EXTRADEFINES)

AC_SUBST(CXXFLAGS)
AC_SUBST(CXXGENFLAGS)
AC_SUBST(LIBDIRS)
AC_SUBST(LD)
AC_SUBST(LDFLAGS)
AC_SUBST(OBJSUF)
AC_SUBST(LIBSUF)


dnl ---------- Checks for "standard" header files. -----------

AC_HEADER_STDC
AC_CHECK_HEADERS(limits.h stdint.h)

dnl ---------- Check for posix_memalign ----------------------
AC_CHECK_FUNC(posix_memalign,[have_posix_memalign=1 AC_DEFINE(HAVE_POSIX_MEMALIGN)],[have_posix_memalign=0])
if test "X$have_posix_memalign" != "X0"; then
 # control memory alignment if posix_memalign is found 
 # if vectorization is off (default) then let the default alignment be determined by the system
 # if vectorization is on then the default alignment is veclen * sizeof(LIBINT2_REALTYPE)
 # see default_params.cc to see this logic in action
 LIBINT_ALIGN_SIZE=default  # this will be interpreted as "use the default alignment"
 AC_ARG_WITH(align-size,
 AS_HELP_STRING([--with-align-size],[(EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics.]),
  [
   if test "X$withval" != "Xdefault"; then
    LIBINT_ALIGN_SIZE=$withval
    if test $LIBINT_ALIGN_SIZE -lt 1; then
     AC_MSG_ERROR([Invalid alignment: $LIBINT_ALIGN_SIZE])
    fi
   fi
  ]
 )
 if test "X$LIBINT_ALIGN_SIZE" != "Xdefault"; then
  AC_MSG_RESULT([will use custom memory alignment: $LIBINT_ALIGN_SIZE * sizeof(LIBINT2_REALTYPE)])
 else
  LIBINT_ALIGN_SIZE=0
 fi
 AC_DEFINE_UNQUOTED(LIBINT_ALIGN_SIZE, $LIBINT_ALIGN_SIZE)
else # no posix_memalign = death
 echo "did not find posix_memalign ... this SHOULD NOT happen. Cannot proceed."
 exit 1
fi

dnl ------------------ Check for support for C++11 standard features --------------

enable_cpptr1=no   # TR1's array seems to be broken on OS X
acx_with_boost=yes
ACX_CHECK_CPP11

if test "X$acx_have_cxx11" = "Xno"; then
  AC_MSG_ERROR([C++11 support is required to compile the Libint compiler. Obtain a recent C++ compiler that supports C++11, or consider using a pregenerated library])
fi

# to support FMA CXXGEN must be able to compile C++11
# NB -std=c++11 is not added to CXXGENFLAGS automatically!
ACX_CHECK_CPP11_CXXGEN
if test "X$libint_generate_gma" = "X1" -a "X$CXXGEN_SUPPORTS_CPP11" = "Xno"; then
  AC_MSG_ERROR([CXXGEN must support C++11 to provide FMA support. Add -std=c++11 to CXXGENFLAGS.])
fi

dnl ------------------ Check for gmp C++ interface -----------------------
ACX_CHECK_GMPXX

if test "X$acx_have_gmpxx" = "Xno"; then
  AC_MSG_ERROR([cannot link against GMP C++ library. Please get the latest GMP with support for C++. It can be obtained at gmplib.org .])
fi

dnl ---------- Check for optional MPFR library (needed for high-precision testing) --------
AC_ARG_ENABLE(mpfr,
 AS_HELP_STRING([--enable-mpfr],[Use MPFR library to test Libint integrals in high precision.]),
 [
  case $enableval in
   yes)
    AC_CHECK_HEADER(mpfr.h,[
     AC_MSG_CHECKING([if MPFR library is usable])
      LIBS="-lmpfr $LIBS"
      AC_LINK_IFELSE(
       [AC_LANG_PROGRAM(
       [[#include <mpfr.h>
       ]],
       [[
        mpfr_t x;  mpfr_init(x);
       ]]
      )
     ],
     [AC_MSG_RESULT([yes])
      AC_DEFINE(LIBINT_HAS_MPFR)
     ],
     [
      AC_MSG_RESULT([no])
     ]
    )
    ]
   )
   ;;
  esac
 ]
)

dnl ------------------ Check for Boost library -----------------------

dnl Boost.Preprocessor requires 1.30 (which introduced MPL; Preprocessor was introduced in 1.29)
AX_BOOST_BASE([1.30],, [AC_MSG_ERROR([libint compiler needs Boost 1.30 or later (needs Boost 1.57 if don't want to use bundled Boost.Preprocessor and be able to use clang), download from https://www.boost.org/users/download/ , unpack (no need to install), and add -I/path/to/boost to CPPFLAGS])])
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
# no need to update LDFLAGS: we are not using precompoiled boost

AC_MSG_CHECKING([for variadic macro support in Boost.Preprocessor])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include <boost/preprocessor.hpp>
      #if not BOOST_PP_VARIADICS  // no variadic macros? your compiler is out of date! (should not be possible since variadic macros are part of C++11)
      #  error "your compiler does not provide variadic macros (but does support C++11), something is seriously broken, please create an issue at https://github.com/evaleev/libint/issues"
      #endif
    ]],
    [[
    ]]
   )
  ],
  [AC_MSG_RESULT([yes])
   AC_DEFINE(LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)
   HAVE_SYSTEM_BOOST_PREPROCESSOR_VARIADICS=1
  ],
  [
    AC_MSG_RESULT([no])
    HAVE_SYSTEM_BOOST_PREPROCESSOR_VARIADICS=0
    AC_MSG_WARN([Boost.Preprocessor with variadic macro support not detected. Will use the bundled Boost.Preprocessor. Download latest boost from https://www.boost.org/users/download/ , unpack (no need to install), and add -I/path/to/boost to CPPFLAGS])
  ]
)
AC_SUBST(HAVE_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)

AC_MSG_CHECKING([for Boost::mpl::for_each])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include <boost/mpl/list.hpp>
      #include <boost/mpl/for_each.hpp>
      #include <iostream>
      struct P {
      template <class T> void operator()(T x) const { std::cout << x << std::endl; }
      };
    ]],
    [[
      typedef boost::mpl::list<int,double,char> tlist;
      boost::mpl::for_each<tlist>(P());
    ]]
   )
  ],
  [AC_MSG_RESULT([yes])
  ],
  [
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Boost::mpl is required. Download latest boost from https://www.boost.org/users/download/ , unpack (no need to install), and add -I/path/to/boost to CPPFLAGS])
  ]
)

dnl ------------------ Check for shared_ptr -----------------------

HAVE_SHARED_PTR=no
AC_MSG_CHECKING([for std::shared_ptr])
AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
      [[#include <memory>
      ]],
      [[
        std::shared_ptr<int> a(new int);
        std::shared_ptr<int> b = a;
      ]]
     )
    ],
    [AC_MSG_RESULT([yes])
     HAVE_SHARED_PTR=std
    ],
    [
      AC_MSG_RESULT([no])
    ]
  )

if test "X$HAVE_SHARED_PTR" = "Xno"; then
  AC_MSG_CHECKING([for Boost::shared_ptr])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
      [[#include <boost/shared_ptr.hpp>
      ]],
      [[
       boost::shared_ptr<int> a(new int);
       boost::shared_ptr<int> b = a;
      ]]
     )
    ],
    [AC_MSG_RESULT([yes])
     HAVE_SHARED_PTR=boost
     AC_DEFINE(HAVE_SHARED_PTR_IN_BOOST)
    ],
    [
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([shared_ptr is required. Please download Boost from boost.org or use C++ compiler that supports std::shared_ptr.])
    ]
  )
fi

dnl ------------------ Check for Eigen library -----------------------

ACX_CHECK_EIGEN
if test "X$acx_have_eigen" = "Xyes"; then
  AC_SUBST(LIBINT_HAS_EIGEN,$acx_have_eigen)
fi

dnl -- Checks for typedefs, structures, and compiler characteristics. --
AC_TYPE_SIZE_T
``
CPPFLAGS="$CPPFLAGS" CXX="$CXX" CXXFLAGS="$CXXFLAGS" \
LD="$LD" LDFLAGS="$LDFLAGS" LIBS="-lm $LIBS $FLIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB"

dnl --------- expand some variables -----------
DATADIR_ABSOLUTE=`eval echo $datarootdir`/libint/$LIBINT_VERSION
AC_SUBST(DATADIR_ABSOLUTE)

dnl Check if sources compiled with CXX and CXXGEN can be linked together
CXX_COMPATIBLE_WITH_CXXGEN=yes
if test "X$CXX" != "X$CXXGEN"; then
AC_MSG_CHECKING([[for compatible cxx and cxxgen compilers]])
echo "extern void func1(); int main(int argc, char**argv) { func1(); } " > linkcheck1.cc
echo "void func1() {}" > linkcheck2.cc
$CXX $CPPFLAGS $CXXFLAGS -c linkcheck1.cc > /dev/null 2>&1
$CXXGEN $CPPFLAGS $CXXGENFLAGS -c linkcheck2.cc > /dev/null 2>&1
$LD $LDFLAGS -o linkcheck1 linkcheck1.o linkcheck2.o > /dev/null 2>&1
$CXX $LDFLAGS -o linkcheck2 linkcheck1.o linkcheck2.o > /dev/null 2>&1
$CXXGEN $LDFLAGS -o linkcheck3 linkcheck1.o linkcheck2.o > /dev/null 2>&1
if test -f linkcheck1 -a -f linkcheck2 -a -f linkcheck3; then
  AC_MSG_RESULT([[yes]])
else
  AC_MSG_RESULT([[no]])
  AC_MSG_WARN([[Compilers cxx and cxxgen produce incompatible code! You will not be able to test generated code.]])
  CXX_COMPATIBLE_WITH_CXXGEN=no
fi
/bin/rm -f linkcheck1 linkcheck2 linkcheck3  linkcheck1.cc linkcheck2.cc linkcheck1.o linkcheck2.o
fi

AC_SUBST(CXX_COMPATIBLE_WITH_CXXGEN)

dnl ------------------ Check if real type is given and is usable with CXXGEN -----------------------

AC_ARG_WITH(real-type,
AS_HELP_STRING([--with-real-type],[Compute all integrals using the given real type (default: double).]),
[
AC_DEFINE_UNQUOTED(LIBINT_USER_DEFINED_REAL,"$withval")
libint_user_defined_real=$withval
AC_MSG_RESULT([Using real type: $withval])
])

# extra includes may be needed in that case!
if test "X$libint_user_defined_real" != "X"; then
 copy_CXX=$CXX
 copy_CXXFLAGS=$CXXFLAGS
 CXX=$CXXGEN
 CXXFLAGS=$CXXGENFLAGS
 
 AC_ARG_WITH(real-type-includes,
  AS_HELP_STRING([--with-real-type-includes],[Additional includes necessary to use the real type (default: none).]),
  [
   AC_DEFINE_UNQUOTED(LIBINT_USER_DEFINED_REAL_INCLUDES,"$withval")
   libint_user_defined_real_includes=$withval
  ]
 )

 AC_MSG_CHECKING([if $libint_user_defined_real is usable])
 # by default should look for SIMD vector types in src/lib/libint
 copy2_CXXFLAGS=$CXXFLAGS
 CXXFLAGS="$CXXFLAGS -I$srcdir/include -Iinclude -I$srcdir/src/lib/libint"
 AC_COMPILE_IFELSE(
  [AC_LANG_PROGRAM(
    [[#include <libint2/util/vector.h>
      $libint_user_defined_real_includes
      
    ]],
    [[
     $libint_user_defined_real x1;
     $libint_user_defined_real x2 = 2.0;
     $libint_user_defined_real x3 = x2;
     $libint_user_defined_real x4 = x2 + x3;
     $libint_user_defined_real x5 = x2 - x3;
     $libint_user_defined_real x6 = x2 * x3;
     $libint_user_defined_real x7 = 2 * x2;
     $libint_user_defined_real x8 = x2 * 3;
     x6 += x2 * x3;
     x7 -= 3 * x3;
    ]]
   )
  ],
  [AC_MSG_RESULT([yes])
   CXXFLAGS=$copy2_CXXFLAGS
  ],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([[real type $libint_user_defined_real is not usable, perhaps extra -I directories or extra #include's are needed?]])
  ]
 )
 
 # CXXFLAGS may have changed
 CXXGENFLAGS=$CXXFLAGS
 CXX=$copy_CXX
 CXXFLAGS=$copy_CXXFLAGS
fi

dnl --------- Shared library configuration. ---------

dnl by default make only static libraries
LT_INIT([disable-shared])
LT_LANG([C++])
OBJSUF=lo
LIBSUF=la
ENABLESHARED=yes
AC_SUBST(ENABLESHARED)

dnl -- configure tester --
TESTER_ARGS="(\"1 0 1 0 0\", \"1 0 1 0 10000\", \"1 1 1 1 0\", \"1 1 1 1 10000\", \"2 0 2 0 0\", \"2 1 2 0 10000\", \"1 1 1 1 10000 16 0\", \"1 1 1 1 10000 16 1\", \"1 1 1 1 10000 1 1\", \"1 1 1 1 0 16 0\", \"1 1 1 1 0 16 1\", \"1 1 1 1 10000 16 0 1\")"
TESTER_ARGS_LONG="(\"3 0 3 0 10000\", \"4 0 4 0 10000\", \"2 2 2 2 0\", \"2 2 2 2 10000\")"
AC_SUBST(TESTER_ARGS)
AC_SUBST(TESTER_ARGS_LONG)

dnl -- make Makefiles in object tree --
$PERL $srcdir/bin/objectdir.pl $EXCLUDED_DIRS $srcdir

AC_CONFIG_FILES([src/bin/test_eri/test_eri.pl],[chmod +x src/bin/test_eri/test_eri.pl])
AC_CONFIG_FILES([src/bin/test_eri/time_eri.pl],[chmod +x src/bin/test_eri/time_eri.pl])
AC_CONFIG_FILES([src/bin/test_eri/run_time_eri.pl],[chmod +x src/bin/test_eri/run_time_eri.pl])
AC_CONFIG_FILES([src/bin/test_eri/stdtests.pl],[chmod +x src/bin/test_eri/stdtests.pl])
AC_CONFIG_FILES([src/bin/test_eri/run_timing_suite.pl],[chmod +x src/bin/test_eri/run_timing_suite.pl])
AC_CONFIG_FILES(src/bin/MakeVars src/lib/MakeVars src/lib/MakeRules src/lib/libint/MakeVars.features
                tests/MakeVars
                doc/MakeVars doc/MakeRules doc/progman/macros.tex doc/classdoc/doxygen.cfg
                include/libint2/basis.h libint2.pc)

dnl --------- Done! ---------
AC_OUTPUT

