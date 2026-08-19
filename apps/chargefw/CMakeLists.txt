add_executable(chargefw_cli
        main.cpp
        input.cpp
        reports.cpp
        output.cpp
)

set_target_properties(chargefw_cli
        PROPERTIES
        OUTPUT_NAME chargefw
)

target_link_libraries(chargefw_cli
        PRIVATE
        chargefw::core
        CLI11::CLI11
)

install(TARGETS chargefw_cli
        RUNTIME DESTINATION bin
)
