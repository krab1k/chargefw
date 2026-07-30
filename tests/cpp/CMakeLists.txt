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

    add_test(
            NAME ${name}
            COMMAND ${name}
    )
endfunction()

if(CHARGEFW_BUILD_CLI)
    add_test(
            NAME test_chargefw_cli_outputs
            COMMAND ${CMAKE_COMMAND}
                    -DCHARGEFW_CLI=$<TARGET_FILE:chargefw_cli>
                    -DCHARGEFW_PARAMETER_DIR=${PROJECT_SOURCE_DIR}/data/parameters
                    -DCHARGEFW_INPUT=${PROJECT_SOURCE_DIR}/tests/water.sdf
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_chargefw_cli_outputs.cmake
    )
endif()

add_subdirectory(core)
add_subdirectory(charges)
add_subdirectory(calculation)
add_subdirectory(adapters)
add_subdirectory(features)
add_subdirectory(methods)
add_subdirectory(parameters)
