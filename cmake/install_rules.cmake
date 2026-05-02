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
