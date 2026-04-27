cmake_minimum_required(VERSION 3.27)

project(
  chargefw
  VERSION 0.0.1
  DESCRIPTION "Empirical partial atomic charge calculation framework"
  LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 23)

option(CHARGEFW_ENABLE_CCACHE "Use ccache if available" ON)
option(CHARGEFW_BUILD_TESTS "Build tests" ON)
option(CHARGEFW_BUILD_CLI "Build the chargefw command-line application" ON)
option(CHARGEFW_ENABLE_WARNINGS "Enable project warnings" ON)

include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

if(CHARGEFW_ENABLE_CCACHE)
  include(cmake/ccache.cmake)
endif()

add_library(chargefw_core)
add_library(chargefw::core ALIAS chargefw_core)

target_sources(chargefw_core
  PRIVATE
        src/core/atom.cpp
        src/core/bond.cpp
        src/core/molecule.cpp
        src/core/conformer.cpp
        src/core/molecule_collection.cpp
        src/features/prepared_molecule.cpp
        src/features/conformer_features.cpp
        src/charges/atomic_charges.cpp
        src/charges/charge_collection.cpp
  PUBLIC
  FILE_SET HEADERS
  BASE_DIRS include
  FILES
        include/chargefw/core/atom.h
        include/chargefw/core/bond.h
        include/chargefw/core/molecule.h
        include/chargefw/core/position.h
        include/chargefw/core/conformer.h
        include/chargefw/core/molecule_collection.h
        include/chargefw/features/prepared_molecule.h
        include/chargefw/features/conformer_features.h
        include/chargefw/charges/atomic_charges.h
        include/chargefw/charges/charge_collection.h
)

target_include_directories(chargefw_core
        PRIVATE
        src
)

if(CHARGEFW_BUILD_CLI)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core)
endif()

if(CHARGEFW_ENABLE_SANITIZERS)
  chargefw_enable_sanitizers(chargefw_core)
endif()

if(CHARGEFW_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
