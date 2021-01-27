add_custom_target(check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V)

add_executable(eritest EXCLUDE_FROM_ALL eri/test.cc)
target_link_libraries(eritest int-library)  # plain int2 target, actually
target_include_directories(eritest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/eri
                                           ${CMAKE_CURRENT_BINARY_DIR}/../include
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
target_link_libraries(eritest MPFR::MPFRXX)
target_compile_definitions(eritest PRIVATE -D__COMPILING_LIBINT2=1)

add_test(eritest_build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target eritest)
set_tests_properties(eritest_build PROPERTIES FIXTURES_SETUP ERITEST_EXEC)
add_test(NAME eritest_run0
         COMMAND $<TARGET_FILE:eritest> 0 2)
set_tests_properties(eritest_run0
                     PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)
add_test(NAME eritest_run1
         COMMAND $<TARGET_FILE:eritest> 1 1)
set_tests_properties(eritest_run1
                     PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)
add_test(NAME eritest_run2
         COMMAND $<TARGET_FILE:eritest> 2 1)
set_tests_properties(eritest_run2
                     PROPERTIES FIXTURES_REQUIRED ERITEST_EXEC)

if (ENABLE_CXX11API)
    enable_language(C)
    add_executable(unit_tests EXCLUDE_FROM_ALL
        unit/test.cc
        unit/c-api.c
        unit/c-api-util.cc
        unit/test-1body.cc
        unit/test-2body.cc
        unit/test-core-ints.cc
        unit/test-c-api.cc
        unit/test-core.cc
        unit/test-permute.cc
        unit/test-shell-order.cc
        unit/test-util.cc)
    target_link_libraries(unit_tests int-library) #cxx)
    target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(unit_tests MPFR::MPFRXX)
    target_compile_definitions(unit_tests PRIVATE -D__COMPILING_LIBINT2=1)
    set_source_files_properties(unit_c-api.c PROPERTIES LANGUAGE C)

    add_test(unit_build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target unit_tests)
    set_tests_properties(unit_build PROPERTIES FIXTURES_SETUP UNIT_TESTS_EXEC)
    add_test(NAME unit_run
            COMMAND $<TARGET_FILE:unit_tests>)
    set_tests_properties(unit_run
            PROPERTIES FIXTURES_REQUIRED UNIT_TESTS_EXEC)

    add_executable(hartree-fock EXCLUDE_FROM_ALL hartree-fock/hartree-fock.cc)
    target_link_libraries(hartree-fock int-library) #cxx)
    target_include_directories(hartree-fock PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(hartree-fock MPFR::MPFRXX)
    target_compile_definitions(hartree-fock PRIVATE -D__COMPILING_LIBINT2=1)

    add_test(hftest_build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target hartree-fock)
    set_tests_properties(hftest_build PROPERTIES FIXTURES_SETUP HFTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME hftest_run
                COMMAND COMMAND ${CMAKE_COMMAND}
                -DtestName=hartree-fock
                -DtestArgs=${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
                -DsrcDir=${CMAKE_SOURCE_DIR}
                -DpythonExec=${PYTHON_EXECUTABLE}
                -P ${CMAKE_SOURCE_DIR}/cmake/hftest.cmake)
    else()
        add_test(NAME hftest_run
                COMMAND $<TARGET_FILE:hartree-fock> ${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o.xyz)
    endif()
    set_tests_properties(hftest_run
            PROPERTIES FIXTURES_REQUIRED HFTEST_EXEC)

    find_package(Threads)  # for some reason clang does not link in threading support even though we are using C++ threads
    add_executable(hartree-fock++ EXCLUDE_FROM_ALL hartree-fock/hartree-fock++.cc)
    target_link_libraries(hartree-fock++ int-library ${CMAKE_THREAD_LIBS_INIT})  #cxx
    target_include_directories(hartree-fock++ PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(hartree-fock++ MPFR::MPFRXX)
    target_compile_definitions(hartree-fock++ PRIVATE -D__COMPILING_LIBINT2=1)

    add_test(hf++test_build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target hartree-fock++)
    set_tests_properties(hf++test_build PROPERTIES FIXTURES_SETUP HFXXTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME hf++test_run
                COMMAND COMMAND ${CMAKE_COMMAND}
                -DtestName=hartree-fock++
                -DtestArgs=${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
                -DsrcDir=${CMAKE_SOURCE_DIR}
                -DpythonExec=${PYTHON_EXECUTABLE}
                -P ${CMAKE_SOURCE_DIR}/cmake/hftest.cmake)
    else()
        add_test(NAME hf++test_run
                COMMAND $<TARGET_FILE:hartree-fock++> ${CMAKE_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz)
    endif()
    set_tests_properties(hf++test_run
            PROPERTIES FIXTURES_REQUIRED HFXXTEST_EXEC)
endif ()
