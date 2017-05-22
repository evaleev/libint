define([libint_mmm_version],[2.4.0])
define([libint_buildid],[beta.2])
define([libint_so_version],[2:2:0])

dnl --------- Begin ---------
dnl Process this file with autoconf to produce a configure script.
define([AC_CACHE_LOAD], )dnl for debugging configure.in
define([AC_CACHE_SAVE], )dnl for debugging configure.in
AC_PREREQ(2.57)
AC_INIT([Libint2],[libint_mmm_version],[libint@valeyev.net],[libint2])
AC_CONFIG_SRCDIR(include/libint2.h)
AC_CONFIG_HEADER(include/libint2/config.h)
AC_CONFIG_AUX_DIR(bin)
AC_CONFIG_MACRO_DIR([lib/autoconf])

dnl --- Can only compile in-place!!! ---
abs_builddir=`pwd -P`
rel_srcdir=`dirname $0`
abs_srcdir=`cd $rel_srcdir; pwd -P`
if test "X$abs_builddir" != "X$abs_srcdir"; then
  AC_MSG_ERROR([** Only in-place build is supported for exported libraries!])
fi

dnl --- Compiling C++ only ---
AC_LANG(C++)

AC_CANONICAL_TARGET

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
AC_SUBST(LIBINT_MAJOR_VERSION)
AC_SUBST(LIBINT_MINOR_VERSION)
AC_SUBST(LIBINT_MICRO_VERSION)


LIBINT_SO_VERSION=libint_so_version
AC_SUBST(LIBINT_SO_VERSION)

dnl --------- Features ---------

AC_ARG_ENABLE(debug,
[  --enable-debug          Compile with debugging options],
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

BUILDID="libint_buildid"
AC_ARG_WITH(build-id,
[  --with-build-id         Gives an identifier for the build.],
BUILDID=$withval
)
AC_SUBST(BUILDID)
AC_DEFINE_UNQUOTED(LIBINT2_BUILDID,"$BUILDID")

if test "$BUILDID"; then
  LIBINT_VERSION="$LIBINT_MMM_VERSION-$BUILDID"
else
  LIBINT_VERSION="$LIBINT_MMM_VERSION"
fi

AC_DEFINE_UNQUOTED(LIBINT2_VERSION, "$LIBINT_VERSION")
AC_SUBST(LIBINT_VERSION)

ac_default_prefix="/usr/local/libint/$LIBINT_VERSION"

AC_ARG_WITH(cxx,
[  --with-cxx              Gives the name of the C++ compiler to use.],
CXX=$withval
)

CXXFLAGS_OPT=
AC_ARG_WITH(cxx-optflags,
[  --with-cxx-optflags     Gives optimization flags for the C++ compiler.],
if test "X$withval" != "Xyes" -a "X$withval" = "Xno"; then
  CXXFLAGS_OPT=$withval
fi
)

AC_ARG_WITH(ranlib,
[  --with-ranlib           Gives the name of the ranlib program.],
RANLIB=$withval
)

AC_ARG_WITH(ar,
[  --with-ar               Names the archive creator.],
AR=$withval
)

ARFLAGS=r
AC_ARG_WITH(ar-flags,
[  --with-ar-flags         Flags for the the archive creator.],
ARFLAGS=$withval
)
AC_SUBST(ARFLAGS)

AC_ARG_WITH(ld,
[  --with-ld               Names the object linker.],
LD=$withval
)

libintincludedir=$includedir
AC_ARG_WITH(libint-includedir,
[  --with-libint-includedir       Specifies include file install subdir.],
libintincludedir=$withval
)
AC_SUBST(libintincludedir)

AC_ARG_WITH(incdirs,
[  --with-incdirs          Specifies include directories (-Idir1 -Idir2).],
EXTRAINCDIRS=$withval
CPPFLAGS=$withval
echo Using extra include directories: $withval
)

AC_ARG_WITH(libs,
[  --with-libs             Specifies libraries (-llib1 -llib2).],
LIBS=$withval
echo Using extra libraries: $withval
)

AC_ARG_WITH(libdirs,
[  --with-libdirs          Specifies library directories (-Ldir1 -Ldir2).],
LIBDIRS=$withval
LDFLAGS=$withval
echo Using extra library directories: $withval
)

PKGCONFIGDIR='$(libdir)/pkgconfig'
AC_ARG_WITH(pkgconfigdir,
[  --with-pkgconfigdir     Specifies the location to put pkg-config's data file. Default is $(libdir)/pkgconfig.],
PKGCONFIGDIR=$withval
)
AC_SUBST(PKGCONFIGDIR)

# these definitions must match configure.ac
LIBINT_SHGSHELL_ORDERING_STANDARD=1
LIBINT_SHGSHELL_ORDERING_GAUSSIAN=2
libint_shgshell_ordering=standard
AC_ARG_WITH(shgauss-ordering,
[  --with-shgauss-ordering          Use one of the following known orderings for shells of solid harmonic Gaussians:
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

dnl --------- Want absolute path for srcdir, not relative. ---------

if test X$ac_srcdir_defaulted = Xyes; then
  case $srcdir in
    ../*)
      srcdir=`(cd $srcdir; pwd)`
      #srcdir=`(cd $srcdir; echo $PWD)`
      ;;
  esac
fi

dnl --------- Checks for programs. ---------
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_RANLIB

# if not given a C++ compiler
if test "X$CXX" = "X"; then
  AC_MSG_WARN([No C++ compiler was given! Will try to find one, but not guaranteed to be the best choice...])
  AC_CHECK_PROGS(CXX, clang++ icpc xlC_r g++ c++ CC, g++)
fi
# discard CXXFLAGS provided by the C++ detection macros
initial_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
AC_PROG_CXXCPP
CXXFLAGS="$initial_CXXFLAGS"

AC_CHECK_PROG(AR,ar,ar)

dnl -------- Checks for compiler/linker options. ----------

# set the optimization default if neither CXXFLAGS or --with-*-optflags were given
if test "X$CXXFLAGS_OPT" = X -a "X$CXXFLAGS" = X; then
  AC_MSG_WARN([C++ optimization options are not given! For optimum performance give CXXFLAGS or use --with-cxx-optflags configure option])
  CXXFLAGS_OPT="-O2 -DNDEBUG"
fi

# options needed only for debugging
CXXFLAGS_DBG=-g
# options that are always needed
CXXFLAGS_MISC=

# now form CXXFLAGS ...
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

dnl -- special architecture options --

if test "X$CXXFLAGS_OPT" = X; then
  AC_MSG_WARN([C++ optimization options are not given! For optimum performance use --with-cxx-optflags configure option])
fi

dnl -- language options --

ACX_CHECK_CPP11_GENERAL
AC_SUBST(LIBINT_HAS_CXX11,$acx_have_cxx11)

dnl ------------------ Check for Eigen library -----------------------

ACX_CHECK_EIGEN
if test "X$acx_have_eigen" = "Xyes"; then
  AC_SUBST(LIBINT_HAS_EIGEN,$acx_have_eigen)
fi

dnl -- other options --

AC_MSG_CHECKING([GNU c++ library/template flags])
if test X$GXX = Xyes; then
  AC_MSG_RESULT([GNU])
else
  case $target in
    *)
      AC_MSG_RESULT([generic non GNU])
      ;;
    esac
fi

AC_SUBST(EXTRAINCLUDE)

AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBDIRS)
AC_SUBST(LD)
AC_SUBST(LDFLAGS)
AC_SUBST(OBJSUF)
AC_SUBST(LIBSUF)
AC_SUBST(CXXDEPENDSUF)
AC_SUBST(CXXDEPENDFLAGS)


dnl -------- Checks for architecture specific features. --------
# nothing so far

dnl -------- extra define for internal includes --------
EXTRADEFINES=-D__COMPILING_LIBINT2=1
AC_SUBST(EXTRADEFINES)

dnl ---------- Checks for header files. -----------
AC_HEADER_STDC
AC_CHECK_HEADERS(limits.h stdint.h)

dnl -- Checks for typedefs, structures, and compiler characteristics. --
AC_TYPE_SIZE_T

dnl ---------- Check for posix_memalign ----------------------
AC_CHECK_FUNC(posix_memalign,[have_posix_memalign=1 AC_DEFINE(HAVE_POSIX_MEMALIGN)],[have_posix_memalign=0])
if test "X$have_posix_memalign" != "X0"; then
 # control memory alignment if posix_memalign is found 
 # if vectorization is off (default) then let the default alignment be determined by the system
 # if vectorization is on then the default alignment is veclen * sizeof(LIBINT2_REALTYPE)
 # see default_params.cc to see this logic in action
 LIBINT_ALIGN_SIZE=default  # this will be interpreted as "use the default alignment"
 AC_ARG_WITH(align-size,
  [  --with-align-size       (EXPERT) if posix_memalign is available, this will specify alignment of Libint data, in units of sizeof(LIBINT2_REALTYPE). Default is to use built-in heuristics.],
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
 AC_DEFINE_UNQUOTED(LIBINT2_ALIGN_SIZE,[$LIBINT_ALIGN_SIZE])
else # no posix_memalign = death
 echo "did not find posix_memalign ... this SHOULD NOT happen. Cannot proceed."
 exit 1
fi

dnl ------------------ Check for real type. If given, check if it is usable -----------------------

AC_ARG_WITH(real-type,
 [  --with-real-type       Compute all integrals using the given real type.],
 [
  AC_DEFINE_UNQUOTED(LIBINT2_REALTYPE,[$withval])
  libint_user_defined_real=$withval
  AC_MSG_RESULT([Using real type: $withval])
 ],
 [
   AC_DEFINE_UNQUOTED(LIBINT2_REALTYPE,[double])
 ]
)

# extra includes may be needed in that case!
if test "X$libint_user_defined_real" != "X"; then
 copy_CXX=$CXX
 copy_CXXFLAGS=$CXXFLAGS
 CXXFLAGS="$CXXFLAGS -Iinclude"
 
 AC_ARG_WITH(real-type-includes,
  [  --with-real-type-includes Additional includes necessary to use the real type (default: none).],
  [
   AC_DEFINE_UNQUOTED(LIBINT2_REALTYPE_INCLUDES,"$withval")
   libint_user_defined_real_includes=$withval
  ]
 )

 AC_MSG_CHECKING([if Libint real type is usable])
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
  ],
  [AC_MSG_RESULT([no])
   AC_MSG_ERROR([[Libint real type is not usable, perhaps extra -I directories or extra #include's are needed?]])
  ]
 )
 CXXFLAGS="$copy_CXXFLAGS" 
fi

CPPFLAGS="$CPPFLAGS" \
LD="$LD" LDFLAGS="$LDFLAGS" LIBS="-lm $LIBS $FLIBS" \
LN_S="$LN_S" NM="$NM" RANLIB="$RANLIB"

dnl --------- expand some variables -----------
DATADIR_ABSOLUTE=`eval echo $datarootdir`/libint/$LIBINT_VERSION
echo $DATADIR_ABSOLUTE
AC_SUBST(DATADIR_ABSOLUTE)

dnl --------- Check for documentation tools ---------

AC_PATH_PROG(DOXYGEN,doxygen,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(PDFLATEX,pdflatex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(PSLATEX,pslatex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(DVIPS,dvips,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(LATEX,latex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(LATEX2HTML,latex2html,,$PATH:/usr/local/bin:/opt/local/bin)
AC_PATH_PROG(BIBTEX,bibtex,,$PATH:/usr/local/bin:/opt/local/bin)
AC_CHECK_PROGS(PYTHON,python python3 python2 python2.7,python)

dnl --------- Shared library configuration. ---------

dnl by default make only static libraries
LT_INIT([disable-shared])
LT_LANG([C++])
OBJSUF=lo
LIBSUF=la
ENABLESHARED=yes
AC_SUBST(ENABLESHARED)

AC_CONFIG_FILES(MakeVars MakeRules MakeSuffixRules doc/Makefile include/libint2/basis.h libint2.pc)

dnl --------- Done! ---------
AC_OUTPUT

