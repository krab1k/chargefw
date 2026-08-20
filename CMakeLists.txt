cmake_minimum_required(VERSION 3.27)

project(
  chargefw
  VERSION 0.0.1
  DESCRIPTION "Empirical partial atomic charge calculation framework"
  LANGUAGES CXX
)

include(GNUInstallDirs)

option(CHARGEFW_ENABLE_CCACHE "Use ccache if available" ON)
option(CHARGEFW_BUILD_TESTS "Build tests" ON)
option(CHARGEFW_BUILD_CLI "Build the chargefw command-line application" ON)
option(CHARGEFW_ENABLE_WARNINGS "Enable project warnings" ON)
option(CHARGEFW_ENABLE_SANITIZERS "Enable sanitizers" OFF)
option(CHARGEFW_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

include(cmake/dependencies.cmake)
include(cmake/developer_tools.cmake)
include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

chargefw_setup_dependencies()
chargefw_setup_developer_tools()

add_library(chargefw_core)
add_library(chargefw::core ALIAS chargefw_core)

target_compile_features(chargefw_core
  PUBLIC
    cxx_std_23
)

target_include_directories(chargefw_core
  PRIVATE
    src

  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
)

target_link_libraries(chargefw_core
  PUBLIC
    nlohmann_json::nlohmann_json
    gemmi::gemmi_cpp
  PRIVATE
    Eigen3::Eigen
    nanoflann::nanoflann
    TBB::tbb
)

configure_file(
  ${PROJECT_SOURCE_DIR}/include/chargefw/config.h.in
  ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
  @ONLY
)

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core)
endif()

if(CHARGEFW_ENABLE_SANITIZERS)
  chargefw_enable_sanitizers(chargefw_core)
endif()

add_subdirectory(src)

if(CHARGEFW_BUILD_CLI)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

include(cmake/install_rules.cmake)
