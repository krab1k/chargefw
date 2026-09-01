set(CHARGEFW_INSTALL_PREFIX
    "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_prefix")
if(NOT DEFINED CHARGEFW_CXX_COMPILER)
    message(FATAL_ERROR "Downstream consumer test requires the configured C++ compiler")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/prepare_moved_install.cmake")

file(GLOB_RECURSE nlohmann_configurations
     "${CHARGEFW_INSTALL_PREFIX}/*nlohmann_json*Config.cmake"
     "${CHARGEFW_INSTALL_PREFIX}/*nlohmann_json*-config.cmake")
if(NOT nlohmann_configurations)
    message(FATAL_ERROR "The installed SDK does not contain the nlohmann/json CMake package")
endif()
list(GET nlohmann_configurations 0 nlohmann_configuration)
get_filename_component(nlohmann_json_directory "${nlohmann_configuration}" DIRECTORY)

file(GLOB_RECURSE gemmi_configurations
     "${CHARGEFW_INSTALL_PREFIX}/*gemmi*Config.cmake"
     "${CHARGEFW_INSTALL_PREFIX}/*gemmi*-config.cmake")
set(gemmi_arguments)
if(gemmi_configurations)
    list(GET gemmi_configurations 0 gemmi_configuration)
    get_filename_component(gemmi_directory "${gemmi_configuration}" DIRECTORY)
    list(APPEND gemmi_arguments -Dgemmi_DIR=${gemmi_directory})
endif()

set(consumer_build_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cmake_consumer_build")
file(REMOVE_RECURSE "${consumer_build_directory}")

execute_process(
        COMMAND "${CMAKE_COMMAND}"
                -S "${CMAKE_CURRENT_LIST_DIR}/consumer"
                -B "${consumer_build_directory}"
                -DCMAKE_BUILD_TYPE=Debug
                 -DCMAKE_CXX_COMPILER=${CHARGEFW_CXX_COMPILER}
                 -DCMAKE_PREFIX_PATH=${CHARGEFW_INSTALL_PREFIX}
                 -Dnlohmann_json_DIR=${nlohmann_json_directory}
                 ${gemmi_arguments}
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
