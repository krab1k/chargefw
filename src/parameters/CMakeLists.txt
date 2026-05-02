add_subdirectory(io)
add_subdirectory(models)
add_subdirectory(classification)

target_sources(chargefw_core
        PRIVATE
        parameter_helpers.h
        parameter_helpers.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
)


