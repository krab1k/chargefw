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
option(CHARGEFW_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(CHARGEFW_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(CHARGEFW_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
option(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS "Optimize for the build host CPU" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

include(cmake/dependencies.cmake)
include(cmake/developer_tools.cmake)
include(cmake/optimization.cmake)
include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

chargefw_setup_dependencies()
chargefw_setup_developer_tools()

add_library(chargefw_core SHARED)
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

target_compile_definitions(chargefw_core
  PRIVATE
    EIGEN_DONT_PARALLELIZE
)

chargefw_enable_optimizations(chargefw_core)

configure_file(
  ${PROJECT_SOURCE_DIR}/include/chargefw/config.h.in
  ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
  @ONLY
)

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core)
endif()

if(CHARGEFW_ENABLE_ASAN)
  chargefw_enable_sanitizer(chargefw_core address)
endif()

if(CHARGEFW_ENABLE_UBSAN)
  chargefw_enable_sanitizer(chargefw_core undefined)
endif()

add_subdirectory(src)

if(CHARGEFW_BUILD_CLI)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_BUILD_TESTS)
  chargefw_setup_test_dependencies()
  enable_testing()
  add_subdirectory(tests)
endif()

include(cmake/install_rules.cmake)
