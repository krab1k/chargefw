target_sources(chargefw_core
        PRIVATE
        method_options.cpp
        method_registry.cpp
        builtin_methods.h
        builtin_methods.cpp

        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/calculation_input.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/feature_requirements.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/method.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/method_metadata.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/method_options.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/methods/method_registry.h
)

add_subdirectory(builtin)