include("${CMAKE_CURRENT_LIST_DIR}/support/test_cli.cmake")

if(NOT DEFINED CHARGEFW_SOURCE_DIR OR NOT DEFINED CHARGEFW_BINARY_DIR OR
   NOT DEFINED CHARGEFW_PARENT_CACHE)
    message(FATAL_ERROR "Custom-layout test requires source, binary, and parent settings paths")
endif()

set(test_root "${CMAKE_CURRENT_BINARY_DIR}/chargefw_custom_install_layout")
set(build_directory "${test_root}/build")
set(source_prefix "${test_root}/source-prefix")
set(moved_prefix "${test_root}/moved-prefix")
file(REMOVE_RECURSE "${test_root}")

set(gemmi_source_argument)
if(EXISTS "${CHARGEFW_BINARY_DIR}/_deps/gemmi-src")
    set(gemmi_source_argument
        -DFETCHCONTENT_SOURCE_DIR_GEMMI=${CHARGEFW_BINARY_DIR}/_deps/gemmi-src)
endif()

run_checked("Custom-layout configuration"
            "${CMAKE_COMMAND}"
            -S "${CHARGEFW_SOURCE_DIR}"
            -B "${build_directory}"
            -G Ninja
            -C "${CHARGEFW_PARENT_CACHE}"
            -DCMAKE_INSTALL_LIBDIR=lib/chargefw
            -DCMAKE_INSTALL_DATADIR=resources
            -DCHARGEFW_BUILD_TESTS=OFF
            -DCHARGEFW_BUILD_CLI=ON
            -DCHARGEFW_BUILD_PYTHON=OFF
            -DFETCHCONTENT_SOURCE_DIR_CLI11=${CHARGEFW_BINARY_DIR}/_deps/cli11-src
            -DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON=${CHARGEFW_BINARY_DIR}/_deps/nlohmann_json-src
            -DFETCHCONTENT_SOURCE_DIR_EIGEN=${CHARGEFW_BINARY_DIR}/_deps/eigen-src
            -DFETCHCONTENT_SOURCE_DIR_NANOFLANN=${CHARGEFW_BINARY_DIR}/_deps/nanoflann-src
            -DFETCHCONTENT_SOURCE_DIR_ONETBB=${CHARGEFW_BINARY_DIR}/_deps/onetbb-src
            ${gemmi_source_argument})
run_checked("Custom-layout build"
            "${CMAKE_COMMAND}" --build "${build_directory}" --target chargefw_cli)
run_checked("Custom-layout installation"
            "${CMAKE_COMMAND}" --install "${build_directory}" --prefix "${source_prefix}")

file(RENAME "${source_prefix}" "${moved_prefix}")
run_checked("Moved custom-layout CLI" "${moved_prefix}/bin/chargefw" parameters EEM_Baek1991)
