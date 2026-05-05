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