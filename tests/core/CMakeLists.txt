add_executable(libint-libcompiler-test EXCLUDE_FROM_ALL test.cc)
target_link_libraries(libint-libcompiler-test libint-libcompiler)

add_executable(libint-kernel-profile EXCLUDE_FROM_ALL profile.cc)
target_link_libraries(libint-kernel-profile libint-libcompiler)

add_executable(libint-chrono-util EXCLUDE_FROM_ALL chrono.cc)
target_link_libraries(libint-chrono-util libint-libcompiler)

set_target_properties(libint-libcompiler-test libint-kernel-profile libint-chrono-util
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

