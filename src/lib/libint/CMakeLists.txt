cmake_minimum_required(VERSION 3.8)

project(export
        LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH ${Libint2_BINARY_DIR}/cmake)

################################# Main Project #################################
include(GNUInstallDirs)
include(int_computed)  # for macros.tex

# <<<  Install  >>>

#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/lib/basis
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/lib/basis/* $(TOPDIR)/$(EXPORTDIR)/lib/basis
#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests
#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/eri
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/eri/test.cc $(TOPDIR)/$(EXPORTDIR)/tests/eri/test.cc
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/eri.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/eri.h
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/src/bin/test_eri/prep_libint2.h $(TOPDIR)/$(EXPORTDIR)/tests/eri/prep_libint2.h
#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/unit
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/catch.hpp $(TOPDIR)/$(EXPORTDIR)/tests/unit/catch.hpp
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/fixture.h $(TOPDIR)/$(EXPORTDIR)/tests/unit/fixture.h
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/*.cc $(TOPDIR)/$(EXPORTDIR)/tests/unit
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/unit/*.c $(TOPDIR)/$(EXPORTDIR)/tests/unit
#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock.cc $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock.cc
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock++.cc $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock++.cc
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/h2o.xyz $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/h2o.xyz
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/h2o_rotated.xyz $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/h2o_rotated.xyz
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock-validate.py $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock-validate.py
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCTOPDIR)/tests/hartree-fock/hartree-fock++-validate.py $(TOPDIR)/$(EXPORTDIR)/tests/hartree-fock/hartree-fock++-validate.py
#    $(INSTALL) $(INSTALLDIROPT) $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.F90 $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/*.cc $(TOPDIR)/$(EXPORTDIR)/fortran
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/c_to_f.py $(TOPDIR)/$(EXPORTDIR)/fortran/c_to_f.py
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/fortran/make_defs.py $(TOPDIR)/$(EXPORTDIR)/fortran/make_defs.py
#    -$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/cmake/CMakeLists.txt.include.export $(TOPDIR)/$(EXPORTDIR)/include/CMakeLists.txt
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
