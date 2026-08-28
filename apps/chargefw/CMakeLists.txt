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

chargefw_enable_optimizations(chargefw_cli)

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_cli)
    # The CLI deliberately uses partial designated initialization for public aggregate result types;
    # omitted fields retain their default member initializers.
    target_compile_options(chargefw_cli PRIVATE -Wno-missing-field-initializers)
endif()

target_link_libraries(chargefw_cli
        PRIVATE
        chargefw::core
        CLI11::CLI11
)

install(TARGETS chargefw_cli
        RUNTIME DESTINATION bin
)
