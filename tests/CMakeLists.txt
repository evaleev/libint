include(AddCustomTargetSubproject)
add_custom_target_subproject(
  libint2 check  # target check-libint2
  USES_TERMINAL
  COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/"
  )

# Note: if forming compile line by hand rather than using targets, you'll need
#   to include the header toggle if compiling within the build tree (as opposed
#   to linking to installed targets):
#   `target_compile_definitions(... PRIVATE -D__COMPILING_LIBINT2=1)`

## Tests ========================================================================

add_executable(
  eritest-libint2
  EXCLUDE_FROM_ALL
  eri/test.cc
  )
target_include_directories(
  eritest-libint2
  PRIVATE
    eri  # src/bin/test_eri/eri.h
    # see note on eri.h at unit_tests-libint2
  )
target_compile_definitions(
  eritest-libint2
  PRIVATE
    $<$<TARGET_EXISTS:Multiprecision::mpfr>:LIBINT_HAS_MPFR=1>
  )
target_link_libraries(
  eritest-libint2
  PRIVATE
    Libint2::int2  # N.B. plain library
    $<$<TARGET_EXISTS:Multiprecision::mpfr>:Multiprecision::mpfr>
    $<$<TARGET_EXISTS:Multiprecision::mpfr>:Multiprecision::gmpxx>
  )
add_test(  # Test #1
  NAME libint2/eritest/build
  COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target eritest-libint2
  )
set_tests_properties(
  libint2/eritest/build
  PROPERTIES
    FIXTURES_SETUP LIBINT2_ERITEST_EXEC
  )
add_test(  # Test #2
  NAME libint2/eritest/run0
  COMMAND $<TARGET_FILE:eritest-libint2> 0 2
  )
add_test(  # Test #3
  NAME libint2/eritest/run1
  COMMAND $<TARGET_FILE:eritest-libint2> 1 1
  )
add_test(  # Test #4
  NAME libint2/eritest/run2
  COMMAND $<TARGET_FILE:eritest-libint2> 2 1
  )
add_test(  # Test #5
  NAME libint2/eritest/run3
  COMMAND $<TARGET_FILE:eritest-libint2> 3 1
  )
add_test(  # Test #6
  NAME libint2/eritest/run4
  COMMAND $<TARGET_FILE:eritest-libint2> 4 1
  )
set_tests_properties(
  libint2/eritest/run0
  libint2/eritest/run1
  libint2/eritest/run2
  libint2/eritest/run3
  libint2/eritest/run4
  PROPERTIES
    FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC
  )

if (LIBINT2_REQUIRE_CXX_API)

    set(utests_src
            unit/catch.hpp
            unit/fixture.h
            unit/test-1body.cc
            unit/test-2body.cc
            unit/test-basis.cc
            unit/test-core.cc
            unit/test-core-ints.cc
            unit/test-permute.cc
            unit/test-precision.cc
            unit/test-shell-order.cc
            unit/test-util.cc
            )
    # if have C, also test C bindings
    include(CheckLanguage)
    check_language(C)
    if (CMAKE_C_COMPILER)
        enable_language(C)
        list(APPEND unit_tests-src
            unit/c-api.c
            unit/c-api-util.cc
            unit/test-c-api.cc
          )
    endif()
    add_executable(
      unit_tests-libint2
      EXCLUDE_FROM_ALL
      unit/test.cc
      ${utests_src}
      )
    target_compile_definitions(
      unit_tests-libint2
      PRIVATE
        NO_LIBINT_COMPILER_CODE  # src/bin/test_eri/eri.h
        # note there's two ways to use eri.h, depending on the test src. eritest-libint2 uses the other.
        $<$<TARGET_EXISTS:Multiprecision::mpfr>:LIBINT_HAS_MPFR=1>
      )
    target_link_libraries(
      unit_tests-libint2
      PRIVATE
        $<IF:$<TARGET_EXISTS:Libint2::int2-cxx>,Libint2::int2-cxx,Libint2::cxx>
        # N.B. cxx compiled library if LIBINT2_REQUIRE_CXX_API_COMPILED=ON else header-only library
        $<$<TARGET_EXISTS:Multiprecision::mpfr>:Multiprecision::mpfr>
        $<$<TARGET_EXISTS:Multiprecision::mpfr>:Multiprecision::gmpxx>
        $<$<TARGET_EXISTS:Multiprecision::mpfr>:Multiprecision::gmp>
        # MPFR: numeric.h:175:23: error: cannot convert 'const char*' to 'double' in initialization: To to(ss.str().c_str());
      )
    add_test(  # Test #7
      NAME libint2/unit/build
      COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target unit_tests-libint2
      )
    set_tests_properties(
      libint2/unit/build
      PROPERTIES
        FIXTURES_SETUP LIBINT2_UNIT_TESTS_EXEC
      )
    add_test(  # Test #8
      NAME libint2/unit/run
      COMMAND $<TARGET_FILE:unit_tests-libint2>
      )
    add_test(  # Test #9
      NAME libint2/unit/sho=gaussian/run
      COMMAND $<TARGET_FILE:unit_tests-libint2> --shgshell-order=gaussian
      )
    set_tests_properties(
      libint2/unit/run
      libint2/unit/sho=gaussian/run
      PROPERTIES
        FIXTURES_REQUIRED LIBINT2_UNIT_TESTS_EXEC
      )

    add_executable(
      hf-libint2
      EXCLUDE_FROM_ALL
      hartree-fock/hartree-fock.cc
      )
    target_link_libraries(
      hf-libint2
      PRIVATE
        $<IF:$<TARGET_EXISTS:Libint2::int2-cxx>,Libint2::int2-cxx,Libint2::cxx>
        # N.B. cxx compiled library if LIBINT2_REQUIRE_CXX_API_COMPILED=ON else header-only library
      )
    add_test(  # Test #10
      NAME libint2/hf/build
      COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf-libint2
      )
    set_tests_properties(
      libint2/hf/build
      PROPERTIES
        FIXTURES_SETUP LIBINT2_HFTEST_EXEC
      )
    if (Python_Interpreter_FOUND)
        add_test(  # Test #11
          NAME libint2/hf/run
          COMMAND ${CMAKE_COMMAND}
                  -DexecName=hf-libint2
                  -DtestName=hartree-fock
                  -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
                  -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                  -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                  -DpythonExec=${Python_EXECUTABLE}
                  -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake
          )
    else()
        add_test(  # Test #11 alt.
          NAME libint2/hf/run
          COMMAND $<TARGET_FILE:hf-libint2> ${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
          )
    endif()
    set_tests_properties(
      libint2/hf/run
      PROPERTIES
        FIXTURES_REQUIRED LIBINT2_HFTEST_EXEC
      )

    add_executable(
      hf++-libint2
      EXCLUDE_FROM_ALL
      hartree-fock/hartree-fock++.cc
      )
    find_package(Threads)  # for some reason clang does not link in threading support even though we are using C++ threads
    target_link_libraries(
      hf++-libint2
      PRIVATE
        Libint2::cxx
        Threads::Threads
        # mostly errors for int2-cxx (excpet some Linux configurations)
        # $<IF:$<TARGET_EXISTS:Libint2::int2-cxx>,Libint2::int2-cxx,Libint2::cxx>
      )
    add_test(  # Test #12
      NAME libint2/hf++/build
      COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf++-libint2
      )
    set_tests_properties(
      libint2/hf++/build
      PROPERTIES
        FIXTURES_SETUP LIBINT2_HFXXTEST_EXEC
      )
    if (Python_Interpreter_FOUND)
        add_test(  # Test #13
          NAME libint2/hf++/run
          COMMAND ${CMAKE_COMMAND}
                  -DexecName=hf++-libint2
                  -DtestName=hartree-fock++
                  -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
                  -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                  -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                  -DpythonExec=${Python_EXECUTABLE}
                  -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake
          )
    else()
        add_test(  # Test #13 alt.
          NAME libint2/hf++/run
          COMMAND $<TARGET_FILE:hf++-libint2> ${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
          )
    endif()
    set_tests_properties(
      libint2/hf++/run
      PROPERTIES
        FIXTURES_REQUIRED LIBINT2_HFXXTEST_EXEC
      )

endif (LIBINT2_REQUIRE_CXX_API)

if (LIBINT2_ENABLE_FORTRAN)

    # Note: if forming compile line by hand rather than using targets, you'll
    #   need to include the Fortran module file directory:
    #   `target_include_directories(... PRIVATE $<TARGET_PROPERTY:libint_f,Fortran_MODULE_DIRECTORY>)`

    add_executable(
      fortran_example-libint2
      EXCLUDE_FROM_ALL
      fortran/fortran_example.F90
      )
    target_link_libraries(
      fortran_example-libint2
      PRIVATE
        Libint2::int2
        libint_f
      )
    add_test(  # Test #14
      NAME libint2/fortran_example/build
      COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_example-libint2
      )
    set_tests_properties(
      libint2/fortran_example/build
      PROPERTIES
        FIXTURES_SETUP LIBINT2_FORTRAN_EXAMPLE_EXEC
      )
    add_test(  # Test #15
      NAME libint2/fortran_example/run
      COMMAND $<TARGET_FILE:fortran_example-libint2>
      )
    set_tests_properties(
      libint2/fortran_example/run
      PROPERTIES
        FIXTURES_REQUIRED LIBINT2_FORTRAN_EXAMPLE_EXEC
      )

    if (LIBINT2_REQUIRE_CXX_API)

        add_executable(
          fortran_test-libint2
          EXCLUDE_FROM_ALL
          fortran/test.cc
          fortran/test-eri.cc
          $<TARGET_OBJECTS:libint_f>
          )
        target_link_libraries(
          fortran_test-libint2
          PRIVATE
            $<IF:$<TARGET_EXISTS:Libint2::int2-cxx>,Libint2::int2-cxx,Libint2::cxx>
            # N.B. cxx compiled library if LIBINT2_REQUIRE_CXX_API_COMPILED=ON else header-only library
            libint_f
          )
        add_test(  # Test #16
          NAME libint2/fortran_test/build
          COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_test-libint2
          )
        set_tests_properties(
          libint2/fortran_test/build
          PROPERTIES
            FIXTURES_SETUP LIBINT2_FORTRAN_TEST_EXEC
          )
        add_test(  # Test #17
          NAME libint2/fortran_test/run
          COMMAND $<TARGET_FILE:fortran_test-libint2>
          )
        set_tests_properties(
          libint2/fortran_test/run
          PROPERTIES
            FIXTURES_REQUIRED LIBINT2_FORTRAN_TEST_EXEC
          )

    endif (LIBINT2_REQUIRE_CXX_API)
endif (LIBINT2_ENABLE_FORTRAN)
