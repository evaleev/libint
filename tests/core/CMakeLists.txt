add_executable(libint-libcompiler-test EXCLUDE_FROM_ALL test.cc)
target_link_libraries(libint-libcompiler-test libint-libcompiler)

add_executable(libint-kernel-profile EXCLUDE_FROM_ALL profile.cc)
target_link_libraries(libint-kernel-profile libint-libcompiler)

add_executable(libint-chrono-util EXCLUDE_FROM_ALL chrono.cc)
target_link_libraries(libint-chrono-util libint-libcompiler)

set_target_properties(libint-libcompiler-test libint-kernel-profile libint-chrono-util
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

add_test(libint2/compiler/coretest/build ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target libint-libcompiler-test)
set_tests_properties(libint2/compiler/coretest/build PROPERTIES FIXTURES_SETUP LIBINT2_COMPILER_CORETEST_EXEC)
add_test(NAME libint2/compiler/coretest/run
        COMMAND $<TARGET_FILE:libint-libcompiler-test>)
set_tests_properties(libint2/compiler/coretest/run
        PROPERTIES FIXTURES_REQUIRED LIBINT2_COMPILER_CORETEST_EXEC)
