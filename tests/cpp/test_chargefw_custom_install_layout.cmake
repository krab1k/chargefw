if(NOT DEFINED CHARGEFW_SOURCE_DIR OR NOT DEFINED CHARGEFW_BINARY_DIR OR
   NOT DEFINED CHARGEFW_CXX_COMPILER)
    message(FATAL_ERROR "Custom-layout test requires source, binary, and compiler paths")
endif()

set(test_root "${CMAKE_CURRENT_BINARY_DIR}/chargefw_custom_install_layout")
set(build_directory "${test_root}/build")
set(source_prefix "${test_root}/source-prefix")
set(moved_prefix "${test_root}/moved-prefix")
file(REMOVE_RECURSE "${test_root}")

execute_process(
        COMMAND "${CMAKE_COMMAND}"
                -S "${CHARGEFW_SOURCE_DIR}"
                -B "${build_directory}"
                -G Ninja
                -DCMAKE_BUILD_TYPE=Debug
                -DCMAKE_CXX_COMPILER=${CHARGEFW_CXX_COMPILER}
                -DCMAKE_INSTALL_LIBDIR=lib/chargefw
                -DCMAKE_INSTALL_DATADIR=resources
                -DCHARGEFW_BUILD_TESTS=OFF
                -DCHARGEFW_BUILD_CLI=ON
                -DFETCHCONTENT_SOURCE_DIR_CLI11=${CHARGEFW_BINARY_DIR}/_deps/cli11-src
                -DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON=${CHARGEFW_BINARY_DIR}/_deps/nlohmann_json-src
                -DFETCHCONTENT_SOURCE_DIR_EIGEN=${CHARGEFW_BINARY_DIR}/_deps/eigen-src
                -DFETCHCONTENT_SOURCE_DIR_NANOFLANN=${CHARGEFW_BINARY_DIR}/_deps/nanoflann-src
                -DFETCHCONTENT_SOURCE_DIR_ONETBB=${CHARGEFW_BINARY_DIR}/_deps/onetbb-src
                -DFETCHCONTENT_SOURCE_DIR_GEMMI=${CHARGEFW_BINARY_DIR}/_deps/gemmi-src
        RESULT_VARIABLE configure_result
        ERROR_VARIABLE configure_error
)
if(NOT configure_result EQUAL 0)
    message(FATAL_ERROR "Custom-layout configuration failed: ${configure_error}")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" --build "${build_directory}" --target chargefw_cli
        RESULT_VARIABLE build_result
        ERROR_VARIABLE build_error
)
if(NOT build_result EQUAL 0)
    message(FATAL_ERROR "Custom-layout build failed: ${build_error}")
endif()

execute_process(
        COMMAND "${CMAKE_COMMAND}" --install "${build_directory}" --prefix "${source_prefix}"
        RESULT_VARIABLE install_result
        ERROR_VARIABLE install_error
)
if(NOT install_result EQUAL 0)
    message(FATAL_ERROR "Custom-layout installation failed: ${install_error}")
endif()

file(RENAME "${source_prefix}" "${moved_prefix}")
execute_process(
        COMMAND "${moved_prefix}/bin/chargefw" parameters eem
        RESULT_VARIABLE run_result
        ERROR_VARIABLE run_error
)
if(NOT run_result EQUAL 0)
    message(FATAL_ERROR "Moved custom-layout CLI failed: ${run_error}")
endif()
