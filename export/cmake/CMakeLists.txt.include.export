
install(DIRECTORY libint2
        COMPONENT int2
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h"
                       PATTERN "*.h.in" EXCLUDE
        )
if (NOT LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/include/libint2/boost
            COMPONENT int2
            DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}/libint2"
            )
endif()
install(FILES libint2.h
        COMPONENT int2
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
        )
install(FILES libint2_params.h libint2_iface.h libint2_types.h
        COMPONENT int2
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}/libint2"
        )
install(FILES libint2.hpp
        COMPONENT cxx
        DESTINATION "${LIBINT2_INSTALL_INCLUDEDIR}"
        )
