set(CHARGEFW_INSTALL_PREFIX
    "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_prefix")
include("${CMAKE_CURRENT_LIST_DIR}/prepare_installed_cli.cmake")

set(consumer_build_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_build")
file(REMOVE_RECURSE "${consumer_build_directory}")

execute_process(
        COMMAND "${CMAKE_COMMAND}"
                -S "${CMAKE_CURRENT_LIST_DIR}/consumer"
                -B "${consumer_build_directory}"
                -DCMAKE_BUILD_TYPE=Debug
                -DCMAKE_PREFIX_PATH=${CHARGEFW_INSTALL_PREFIX}
        RESULT_VARIABLE configure_result
        ERROR_VARIABLE configure_error
)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "Downstream ChargeFW consumer configuration failed: ${configure_error}")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" --build "${consumer_build_directory}"
        RESULT_VARIABLE build_result
        ERROR_VARIABLE build_error
)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "Downstream ChargeFW consumer build failed: ${build_error}")
endif()

execute_process(
        COMMAND "${consumer_build_directory}/chargefw_consumer"
        RESULT_VARIABLE run_result
        ERROR_VARIABLE run_error
)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "Downstream ChargeFW consumer failed: ${run_error}")
endif()
