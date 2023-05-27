include(AddCustomTargetSubproject)
add_custom_target_subproject(libint2 check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/")
add_dependencies(check-libint2 Libint2::int2)

# Note: if forming compile line by hand rather than using targets, you'll need to include the header toggle
#   if compiling within the build tree (as opposed to linking to installed targets)
#   `target_compile_definitions(... PRIVATE -D__COMPILING_LIBINT2=1)`

add_executable(eritest-libint2 EXCLUDE_FROM_ALL eri/test.cc)
target_link_libraries(eritest-libint2 Libint2::int2)  # N.B. plain library
target_include_directories(eritest-libint2 PRIVATE eri
        ${CMAKE_CURRENT_BINARY_DIR}/../include)
if (TARGET Multiprecision::mpfr)
    target_compile_definitions(eritest-libint2 PRIVATE -DLIBINT_HAS_MPFR=1)
    target_link_libraries(eritest-libint2 Multiprecision::mpfr)
endif()

# Test 1
add_test(libint2/eritest/build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target eritest-libint2)
set_tests_properties(libint2/eritest/build PROPERTIES FIXTURES_SETUP LIBINT2_ERITEST_EXEC)
# Test 2
add_test(NAME libint2/eritest/run0
        COMMAND $<TARGET_FILE:eritest-libint2> 0 2)
set_tests_properties(libint2/eritest/run0
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
# Test 3
add_test(NAME libint2/eritest/run1
        COMMAND $<TARGET_FILE:eritest-libint2> 1 1)
set_tests_properties(libint2/eritest/run1
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
# Test 4
add_test(NAME libint2/eritest/run2
        COMMAND $<TARGET_FILE:eritest-libint2> 2 1)
set_tests_properties(libint2/eritest/run2
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
# Test 5
add_test(NAME libint2/eritest/run3
        COMMAND $<TARGET_FILE:eritest-libint2> 3 1)
set_tests_properties(libint2/eritest/run3
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
# Test 6
add_test(NAME libint2/eritest/run4
        COMMAND $<TARGET_FILE:eritest-libint2> 4 1)
set_tests_properties(libint2/eritest/run4
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)

if (LIBINT_HAS_CXX_API)
    set(unit_tests-src
            unit/test-1body.cc
            unit/test-2body.cc
            unit/test-basis.cc
            unit/test-core.cc
            unit/test-core-ints.cc
            unit/test-permute.cc
            unit/test-precision.cc
            unit/test-shell-order.cc
            unit/test-util.cc)

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
    set(unit_tests_s-src ${unit_tests-src} unit/test_s.cc)
    add_executable(unit_tests_s-libint2 EXCLUDE_FROM_ALL
        ${unit_tests_s-src})
    target_link_libraries(unit_tests_s-libint2 Libint2::cxx) # N.B. cxx header-only library
    target_include_directories(unit_tests_s-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include)
    target_compile_definitions(unit_tests_s-libint2 PRIVATE -DNO_LIBINT_COMPILER_CODE=1)
    if (TARGET Multiprecision::mpfr)
        target_compile_definitions(unit_tests_s-libint2 PRIVATE -DLIBINT_HAS_MPFR=1)
        target_link_libraries(unit_tests_s-libint2 Multiprecision::mpfr)
    endif()
    set(unit_tests_g-src ${unit_tests-src} unit/test_g.cc)
    add_executable(unit_tests_g-libint2 EXCLUDE_FROM_ALL
        ${unit_tests_g-src})
    target_link_libraries(unit_tests_g-libint2 Libint2::cxx) # N.B. cxx header-only library
    target_include_directories(unit_tests_g-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include)
    target_compile_definitions(unit_tests_g-libint2 PRIVATE -DNO_LIBINT_COMPILER_CODE=1)
    if (TARGET Multiprecision::mpfr)
        target_compile_definitions(unit_tests_g-libint2 PRIVATE -DLIBINT_HAS_MPFR=1)
        target_link_libraries(unit_tests_g-libint2 Multiprecision::mpfr)
    endif()

    # Test 7
    add_test(libint2/unit_s/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target unit_tests_s-libint2)
    set_tests_properties(libint2/unit_s/build PROPERTIES FIXTURES_SETUP LIBINT2_UNIT_TESTS_EXEC)
    # Test 8
    add_test(NAME libint2/unit_s/run
            COMMAND $<TARGET_FILE:unit_tests_s-libint2>)
    set_tests_properties(libint2/unit_s/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_UNIT_TESTS_EXEC)

    # Test 9
    add_test(libint2/unit_g/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target unit_tests_g-libint2)
    set_tests_properties(libint2/unit_g/build PROPERTIES FIXTURES_SETUP LIBINT2_UNIT_TESTS_EXEC)
    # Test 10
    add_test(NAME libint2/unit_g/run
            COMMAND $<TARGET_FILE:unit_tests_g-libint2>)
    set_tests_properties(libint2/unit_g/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_UNIT_TESTS_EXEC)

    add_executable(hf-libint2 EXCLUDE_FROM_ALL hartree-fock/hartree-fock.cc)
    target_link_libraries(hf-libint2 Libint2::cxx) # N.B. cxx library
    target_include_directories(hf-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include)

    # Test 11
    add_test(libint2/hf/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf-libint2)
    set_tests_properties(libint2/hf/build PROPERTIES FIXTURES_SETUP LIBINT2_HFTEST_EXEC)
    if (Python_Interpreter_FOUND)
        # Test 12
        add_test(NAME libint2/hf/run
                 COMMAND ${CMAKE_COMMAND}
                 -DexecName=hf-libint2
                 -DtestName=hartree-fock
                 -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
                 -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                 -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                 -DpythonExec=${Python_EXECUTABLE}
                 -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake)
    else()
        # Test 12 alt.
        add_test(NAME libint2/hf/run
                COMMAND $<TARGET_FILE:hf-libint2> hartree-fock/h2o.xyz)
    endif()
    set_tests_properties(libint2/hf/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_HFTEST_EXEC)

    find_package(Threads)  # for some reason clang does not link in threading support even though we are using C++ threads
    add_executable(hf++-libint2 EXCLUDE_FROM_ALL hartree-fock/hartree-fock++.cc)
    target_link_libraries(hf++-libint2 Libint2::cxx Threads::Threads ${CMAKE_THREAD_LIBS_INIT})  # N.B. cxx library
    target_include_directories(hf++-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include)

    # Test 13
    add_test(libint2/hf++/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf++-libint2)
    set_tests_properties(libint2/hf++/build PROPERTIES FIXTURES_SETUP LIBINT2_HFXXTEST_EXEC)
    if (Python_Interpreter_FOUND)
        # Test 14
        add_test(NAME libint2/hf++/run
                COMMAND ${CMAKE_COMMAND}
                -DexecName=hf++-libint2
                -DtestName=hartree-fock++
                -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
                -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                -DpythonExec=${Python_EXECUTABLE}
                -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake)
    else()
        # Test 14 alt.
        add_test(NAME libint2/hf++/run
                COMMAND $<TARGET_FILE:hf++-libint2> hartree-fock/h2o_rotated.xyz)
    endif()
    set_tests_properties(libint2/hf++/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_HFXXTEST_EXEC)

endif()

if (LIBINT2_ENABLE_FORTRAN)

    # Note: if forming compile line by hand rather than using targets, you'll need to include the Fortran module file directory:
    #   `target_include_directories(... PRIVATE $<TARGET_PROPERTY:libint_f,Fortran_MODULE_DIRECTORY>)`

    # tests
    add_executable(fortran_example-libint2 EXCLUDE_FROM_ALL fortran/fortran_example.F90)
    target_link_libraries(fortran_example-libint2 Libint2::int2 libint_f)

    # Test 15
    add_test(libint2/fortran_example/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_example-libint2)
    set_tests_properties(libint2/fortran_example/build PROPERTIES FIXTURES_SETUP LIBINT2_FORTRAN_EXAMPLE_EXEC)
    # Test 16
    add_test(NAME libint2/fortran_example/run
            COMMAND $<TARGET_FILE:fortran_example-libint2>)
    set_tests_properties(libint2/fortran_example/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_FORTRAN_EXAMPLE_EXEC)

    if (LIBINT_HAS_CXX_API)
        add_executable(fortran_test-libint2 EXCLUDE_FROM_ALL fortran/test.cc fortran/test-eri.cc $<TARGET_OBJECTS:libint_f>)
        target_link_libraries(fortran_test-libint2 Libint2::cxx libint_f)

        # Test 17
        add_test(libint2/fortran_test/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target fortran_test-libint2)
        set_tests_properties(libint2/fortran_test/build PROPERTIES FIXTURES_SETUP LIBINT2_FORTRAN_TEST_EXEC)
        # Test 18
        add_test(NAME libint2/fortran_test/run
                COMMAND $<TARGET_FILE:fortran_test-libint2>)
        set_tests_properties(libint2/fortran_test/run
                PROPERTIES FIXTURES_REQUIRED LIBINT2_FORTRAN_TEST_EXEC)
    endif()
endif()
