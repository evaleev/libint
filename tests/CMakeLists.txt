include(AddCustomTargetSubproject)
add_custom_target_subproject(libint2 check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/")

add_executable(eritest-libint2 EXCLUDE_FROM_ALL eri/test.cc)
target_link_libraries(eritest-libint2 int-library)  # plain int2 target, actually
target_include_directories(eritest-libint2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/eri
                                           ${CMAKE_CURRENT_BINARY_DIR}/../include
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                           ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
target_link_libraries(eritest-libint2 MPFR::MPFRXX)
target_compile_definitions(eritest-libint2 PRIVATE -D__COMPILING_LIBINT2=1)

add_test(libint2/eritest/build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target eritest-libint2)
set_tests_properties(libint2/eritest/build PROPERTIES FIXTURES_SETUP LIBINT2_ERITEST_EXEC)
add_test(NAME libint2/eritest/run0
        COMMAND $<TARGET_FILE:eritest-libint2> 0 2)
set_tests_properties(libint2/eritest/run0
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
add_test(NAME libint2/eritest/run1
        COMMAND $<TARGET_FILE:eritest-libint2> 1 1)
set_tests_properties(libint2/eritest/run1
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)
add_test(NAME libint2/eritest/run2
        COMMAND $<TARGET_FILE:eritest-libint2> 2 1)
set_tests_properties(libint2/eritest/run2
        PROPERTIES FIXTURES_REQUIRED LIBINT2_ERITEST_EXEC)

if (ENABLE_CXX11API)
    set(unit_tests-src unit/test.cc
            unit/test-1body.cc
            unit/test-2body.cc
            unit/test-core-ints.cc
            unit/test-core.cc
            unit/test-permute.cc
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
    endif(CMAKE_C_COMPILER)
    add_executable(unit_tests-libint2 EXCLUDE_FROM_ALL
        ${unit_tests-src})
    target_link_libraries(unit_tests-libint2 int-library) #cxx)
    target_include_directories(unit_tests-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                  ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(unit_tests-libint2 MPFR::MPFRXX)
    target_compile_definitions(unit_tests-libint2 PRIVATE -D__COMPILING_LIBINT2=1)

    add_test(libint2/unit/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target unit_tests-libint2)
    set_tests_properties(libint2/unit/build PROPERTIES FIXTURES_SETUP LIBINT2_UNIT_TESTS_EXEC)
    add_test(NAME libint2/unit/run
            COMMAND $<TARGET_FILE:unit_tests-libint2>)
    set_tests_properties(libint2/unit/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_UNIT_TESTS_EXEC)

    add_executable(hf-libint2 EXCLUDE_FROM_ALL hartree-fock/hartree-fock.cc)
    target_link_libraries(hf-libint2 int-library) #cxx)
    target_include_directories(hf-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                    ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(hf-libint2 MPFR::MPFRXX)
    target_compile_definitions(hf-libint2 PRIVATE -D__COMPILING_LIBINT2=1)

    add_test(libint2/hf/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf-libint2)
    set_tests_properties(libint2/hf/build PROPERTIES FIXTURES_SETUP LIBINT2_HFTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME libint2/hf/run
                 COMMAND ${CMAKE_COMMAND}
                 -DexecName=hf-libint2
                 -DtestName=hartree-fock
                 -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o.xyz
                 -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                 -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                 -DpythonExec=${PYTHON_EXECUTABLE}
                 -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake)
    else()
        add_test(NAME libint2/hf/run
                COMMAND $<TARGET_FILE:hartree-fock> ${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o.xyz)
    endif()
    set_tests_properties(libint2/hf/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_HFTEST_EXEC)

    find_package(Threads)  # for some reason clang does not link in threading support even though we are using C++ threads
    add_executable(hf++-libint2 EXCLUDE_FROM_ALL hartree-fock/hartree-fock++.cc)
    target_link_libraries(hf++-libint2 int-library ${CMAKE_THREAD_LIBS_INIT})  #cxx
    target_include_directories(hf++-libint2 PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../include
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../src
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                                      ${CMAKE_CURRENT_SOURCE_DIR}/../include/libint2)
    target_link_libraries(hf++-libint2 MPFR::MPFRXX)
    target_compile_definitions(hf++-libint2 PRIVATE -D__COMPILING_LIBINT2=1)

    add_test(libint2/hf++/build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target hf++-libint2)
    set_tests_properties(libint2/hf++/build PROPERTIES FIXTURES_SETUP LIBINT2_HFXXTEST_EXEC)
    if (PYTHONINTERP_FOUND)
        add_test(NAME libint2/hf++/run
                COMMAND ${CMAKE_COMMAND}
                -DexecName=hf++-libint2
                -DtestName=hartree-fock++
                -DtestArgs=${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz
                -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
                -DpythonExec=${PYTHON_EXECUTABLE}
                -P ${PROJECT_SOURCE_DIR}/tests/hartree-fock/hftest.cmake)
    else()
        add_test(NAME libint2/hf++/run
                COMMAND $<TARGET_FILE:hf++-libint2> ${PROJECT_SOURCE_DIR}/tests/hartree-fock/h2o_rotated.xyz)
    endif()
    set_tests_properties(libint2/hf++/run
            PROPERTIES FIXTURES_REQUIRED LIBINT2_HFXXTEST_EXEC)

endif ()
