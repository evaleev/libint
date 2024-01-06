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
  )
target_link_libraries(
  eritest-libint2 
  Libint2::int2  # N.B. plain library
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
    set(utests_src
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
    add_executable(
      unit_tests-libint2 
      EXCLUDE_FROM_ALL 
      unit/test.cc 
      ${utests_src}
      )


