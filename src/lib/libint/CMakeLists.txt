cmake_minimum_required(VERSION 3.8)

project(export
        LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH ${Libint2_BINARY_DIR}/cmake)

################################# Main Project #################################
include(GNUInstallDirs)
include(int_computed)  # for macros.tex

# <<<  Install  >>>

#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.F90 $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.cc $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/c_to_f.py $(TOPDIR)/$(EXPORTDIR)/fortran/c_to_f.py
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/make_defs.py $(TOPDIR)/$(EXPORTDIR)/fortran/make_defs.py
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/hftest.cmake $(TOPDIR)/$(EXPORTDIR)/cmake/hftest.cmake

# <<<  Install General  >>>

install(FILES ${Libint2_ROOT}/export/LICENSE.export
        RENAME LICENSE
        DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES ${Libint2_ROOT}/export/INSTALL.export
        RENAME INSTALL
        DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES ${Libint2_ROOT}/COPYING
              ${Libint2_ROOT}/COPYING.LESSER
              ${Libint2_ROOT}/CITATION
              ${Libint2_ROOT}/README.md
        DESTINATION ${CMAKE_INSTALL_PREFIX})

install(FILES ${Libint2_ROOT}/doc/progman/progman.tex
              ${Libint2_ROOT}/doc/progman/refs.bib
              ${Libint2_ROOT}/doc/Libint_Logo3_alt.pdf
              ${Libint2_ROOT}/doc/Libint_Logo3_alt.eps
        DESTINATION ${CMAKE_INSTALL_PREFIX}/doc)

configure_file(${Libint2_ROOT}/doc/progman/macros.tex.in ${CMAKE_INSTALL_PREFIX}/doc/macros.tex @ONLY)

install(DIRECTORY ${Libint2_ROOT}/doc/progman/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/doc
        FILES_MATCHING PATTERN "*.cc")

install(DIRECTORY ${Libint2_ROOT}/tests/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/tests
        FILES_MATCHING PATTERN "*.c"
                       PATTERN "*.cc"
                       PATTERN "*.h"
                       PATTERN "*.hpp"
                       PATTERN "*.py"
                       PATTERN "*.xyz"
                       PATTERN "CMakeLists.txt")

install(FILES ${Libint2_ROOT}/src/bin/test_eri/eri.h
              ${Libint2_ROOT}/src/bin/test_eri/prep_libint2.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/eri)

install(DIRECTORY ${Libint2_ROOT}/lib/basis
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION}
        #DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/basis
        FILES_MATCHING PATTERN "*.g94")

install(FILES ${Libint2_ROOT}/external/boost.tar.gz
        DESTINATION ${CMAKE_INSTALL_PREFIX}/external)

# <<<  Install CMake  >>>

install(FILES ${Libint2_ROOT}/cmake/autocmake_safeguards.cmake
              ${Libint2_ROOT}/cmake/int_orderings.cmake
              ${Libint2_ROOT}/cmake/psi4OptionsTools.cmake
              ${Libint2_ROOT}/cmake/Libint2Config.cmake.in
              ${Libint2_ROOT}/cmake/libint2.pc.cmake.in
              ${Libint2_ROOT}/export/cmake/modules/FindMPFR.cmake
              ${Libint2_BINARY_DIR}/cmake/int_computed.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)

install(FILES CMakeLists.txt.export
        RENAME CMakeLists.txt
        DESTINATION ${CMAKE_INSTALL_PREFIX})

# <<<  Install Headers  >>>

install(FILES ${Libint2_ROOT}/include/libint2.h
              ${Libint2_ROOT}/include/libint2.hpp
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

install(DIRECTORY ${Libint2_ROOT}/include/libint2
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        FILES_MATCHING PATTERN "*.h")

install(FILES ${Libint2_ROOT}/src/bin/libint/util_types.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/libint2)

install(DIRECTORY ./
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/libint2
        FILES_MATCHING PATTERN "*.h")

# this isn't stamping the final install dir, but there's plenty else wrong too
set(DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(${Libint2_ROOT}/include/libint2/basis.h.in ${CMAKE_INSTALL_PREFIX}/include/libint2/basis.h @ONLY)
