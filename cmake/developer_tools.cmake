function(chargefw_setup_developer_tools)
  if(CHARGEFW_ENABLE_CCACHE)
    find_program(CCACHE_PROGRAM ccache)

    if(CCACHE_PROGRAM)
      message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
      set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}" PARENT_SCOPE)
    else()
      message(STATUS "ccache not found")
    endif()
  endif()

  if(CHARGEFW_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE clang-tidy REQUIRED)

    set(
      CMAKE_CXX_CLANG_TIDY
      "${CLANG_TIDY_EXE}"
      "--header-filter=^${PROJECT_SOURCE_DIR}/(include|src|apps|tests)/"
      "--exclude-header-filter=^(${PROJECT_BINARY_DIR}|${PROJECT_SOURCE_DIR}/build)/"
      PARENT_SCOPE
    )
  endif()
endfunction()
