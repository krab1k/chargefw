set_target_properties(chargefw_core chargefw_cli
        PROPERTIES
        INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
)

# Prefer the installed dependency copies over same-named libraries from
# LD_LIBRARY_PATH. This is important for fetched dependencies such as Gemmi,
# whose unversioned SONAME can collide with another installation.
target_link_options(chargefw_core PRIVATE "LINKER:--disable-new-dtags")
target_link_options(chargefw_cli PRIVATE "LINKER:--disable-new-dtags")

install(
        TARGETS chargefw_core
        FILE_SET HEADERS
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME
        DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
        DIRECTORY
        ${PROJECT_SOURCE_DIR}/data/parameters
        DESTINATION
        ${CMAKE_INSTALL_DATADIR}/chargefw
        FILES_MATCHING
        PATTERN "*.json"
)

install(
        FILES
        ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
        DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/chargefw
)

if(CHARGEFW_TBB_FETCHED AND TARGET tbb)
    install(TARGETS tbb
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()

if(CHARGEFW_GEMMI_FETCHED AND TARGET gemmi_cpp)
    install(TARGETS gemmi_cpp
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
