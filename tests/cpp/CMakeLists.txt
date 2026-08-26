add_library(chargefw_test_support
        support/test_molecules.cpp
        support/test_molecules.h
)

target_link_libraries(chargefw_test_support
        PUBLIC
        chargefw::core
)

target_include_directories(chargefw_test_support
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

function(chargefw_add_test name)
    add_executable(${name}
            ${ARGN}
    )

    target_link_libraries(${name}
            PRIVATE
            chargefw::core
            chargefw_test_support
    )

    target_link_options(${name} PRIVATE "LINKER:--disable-new-dtags")

    add_test(
            NAME ${name}
            COMMAND ${name}
    )
endfunction()

function(chargefw_add_snitch_test name)
    add_executable(${name}
            ${ARGN}
    )

    target_link_libraries(${name}
            PRIVATE
            chargefw::core
            chargefw_test_support
            snitch::snitch
    )

    if(CHARGEFW_ENABLE_CLANG_TIDY)
        # Snitch's expression decomposition intentionally forms an intermediate chained comparison.
        # Suppress the resulting false positive only for Snitch test targets; project targets retain
        # this clang-tidy check.
        set_property(TARGET ${name} APPEND PROPERTY CXX_CLANG_TIDY
                "--checks=-bugprone-chained-comparison")
    endif()

    target_link_options(${name} PRIVATE "LINKER:--disable-new-dtags")

    add_test(
            NAME ${name}
            COMMAND ${name}
    )
endfunction()

if(CHARGEFW_BUILD_CLI)
    add_test(
            NAME test_chargefw_cli_outputs
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_BUILD_DIR=${PROJECT_BINARY_DIR}
                    -DCHARGEFW_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                    -DCHARGEFW_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/tests/cpp/chargefw_cli_outputs_prefix
                    -DCHARGEFW_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/water.sdf
                    -DCHARGEFW_MOL2_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/mol2/aromatic.mol2
                    -DCHARGEFW_MIXED_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/mixed_v2000_v3000.sdf
                    -DCHARGEFW_MALFORMED_INPUT=${PROJECT_SOURCE_DIR}/tests/fixtures/synthetic/sdf/malformed_then_water.sdf
                    -DCHARGEFW_VERSION=${PROJECT_VERSION}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_outputs.cmake
    )
    add_test(
            NAME test_chargefw_cli_structural_outputs
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_BUILD_DIR=${PROJECT_BINARY_DIR}
                    -DCHARGEFW_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
                    -DCHARGEFW_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/tests/cpp/chargefw_cli_structural_outputs_prefix
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_structural_outputs.cmake
    )
endif()

add_test(
        NAME test_chargefw_cmake_consumer
        COMMAND ${CMAKE_COMMAND}
                -DCHARGEFW_BUILD_DIR=${PROJECT_BINARY_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cmake_consumer.cmake
)
add_test(
        NAME test_chargefw_custom_install_layout
        COMMAND ${CMAKE_COMMAND}
                -DCHARGEFW_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                -DCHARGEFW_BINARY_DIR=${PROJECT_BINARY_DIR}
                -DCHARGEFW_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_custom_install_layout.cmake
)

add_subdirectory(core)
add_subdirectory(charges)
add_subdirectory(calculation)
add_subdirectory(adapters)
add_subdirectory(features)
add_subdirectory(methods)
add_subdirectory(parameters)
