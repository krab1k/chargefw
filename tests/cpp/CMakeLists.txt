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

add_subdirectory(core)
add_subdirectory(features)
add_subdirectory(methods)