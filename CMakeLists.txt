cmake_minimum_required(VERSION 3.27)

project(
        chargefw
        VERSION 0.0.1
        DESCRIPTION "Empirical partial atomic charge calculation framework"
        LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CHARGEFW_ENABLE_CCACHE "Use ccache if available" ON)
option(CHARGEFW_BUILD_TESTS "Build tests" ON)
option(CHARGEFW_BUILD_CLI "Build the chargefw command-line application" ON)
option(CHARGEFW_ENABLE_WARNINGS "Enable project warnings" ON)
option(CHARGEFW_ENABLE_SANITIZERS "Enable sanitizers" OFF)

include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

if(CHARGEFW_ENABLE_CCACHE)
  include(cmake/ccache.cmake)
endif()

add_library(chargefw_core)
add_library(chargefw::core ALIAS chargefw_core)

target_compile_features(chargefw_core
        PUBLIC
        cxx_std_23
)

target_include_directories(chargefw_core
        PRIVATE
        src
)

add_subdirectory(src)

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core)
endif()

if(CHARGEFW_ENABLE_SANITIZERS)
  chargefw_enable_sanitizers(chargefw_core)
endif()

if(CHARGEFW_BUILD_CLI)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()