set_target_properties(chargefw_core
        PROPERTIES
        INSTALL_RPATH "$ORIGIN"
)

# Prefer the installed dependency copies over same-named libraries from
# LD_LIBRARY_PATH. This is important for fetched dependencies such as Gemmi,
# whose unversioned SONAME can collide with another installation.
target_link_options(chargefw_core PRIVATE "LINKER:--disable-new-dtags")

if(TARGET chargefw_cli)
    set_target_properties(chargefw_cli
            PROPERTIES
            INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}"
    )
    target_link_options(chargefw_cli PRIVATE "LINKER:--disable-new-dtags")
endif()

install(
        TARGETS chargefw_core
        EXPORT chargefwTargets
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

configure_package_config_file(
        ${PROJECT_SOURCE_DIR}/cmake/chargefwConfig.cmake.in
        ${PROJECT_BINARY_DIR}/chargefwConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)
write_basic_package_version_file(
        ${PROJECT_BINARY_DIR}/chargefwConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)
install(
        EXPORT chargefwTargets
        FILE chargefwTargets.cmake
        NAMESPACE chargefw::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)
install(
        FILES
        ${PROJECT_BINARY_DIR}/chargefwConfig.cmake
        ${PROJECT_BINARY_DIR}/chargefwConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/chargefw
)

if(CHARGEFW_TBB_FETCHED AND TARGET tbb)
    set_target_properties(tbb PROPERTIES EXPORT_NAME _tbb_runtime)
    install(TARGETS tbb
            EXPORT chargefwTargets
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
