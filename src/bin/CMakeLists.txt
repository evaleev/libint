add_subdirectory(libint)

include(AddCustomTargetSubproject)
add_custom_target_subproject(
  libint2compiler
  check
  USES_TERMINAL
  COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/compiler"
  )

if (BUILD_TESTING)
    add_subdirectory(test_eri)
    # borrows main libint/test.cc from libtool layout
    add_subdirectory(profile)
endif()
