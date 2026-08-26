if(NOT DEFINED CHARGEFW_BUILD_DIR OR NOT DEFINED CHARGEFW_INSTALL_BINDIR OR
   NOT DEFINED CHARGEFW_INSTALL_PREFIX)
    message(FATAL_ERROR "Installed CLI test setup requires build directory, bindir, and prefix")
endif()

set(install_source_prefix "${CHARGEFW_INSTALL_PREFIX}_source")
file(REMOVE_RECURSE "${install_source_prefix}" "${CHARGEFW_INSTALL_PREFIX}")

execute_process(
        COMMAND "${CMAKE_COMMAND}" --install "${CHARGEFW_BUILD_DIR}" --prefix "${install_source_prefix}"
        RESULT_VARIABLE install_result
        ERROR_VARIABLE install_error
)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "ChargeFW installation failed: ${install_error}")
endif()

file(RENAME "${install_source_prefix}" "${CHARGEFW_INSTALL_PREFIX}")
set(CHARGEFW_CLI "${CHARGEFW_INSTALL_PREFIX}/${CHARGEFW_INSTALL_BINDIR}/chargefw")
