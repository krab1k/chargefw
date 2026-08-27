if(CHARGEFW_PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE "${CHARGEFW_PYTHON_EXECUTABLE}")
endif()
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
set(CHARGEFW_RESOLVED_PYTHON_EXECUTABLE "${Python_EXECUTABLE}" PARENT_SCOPE)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE nanobind_ROOT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
list(PREPEND CMAKE_PREFIX_PATH "${nanobind_ROOT}")
find_package(nanobind CONFIG REQUIRED)

set(CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY "${PROJECT_BINARY_DIR}/python/chargefw")
file(MAKE_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}")
file(COPY "${PROJECT_SOURCE_DIR}/data/parameters"
     DESTINATION "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/_data")
file(COPY "${PROJECT_SOURCE_DIR}/python/chargefw/_chargefw"
     DESTINATION "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}")
set(CHARGEFW_PYTHON_PACKAGE_SOURCES
    __init__.py
    _arrays.py
    _resources.py
    calculation.py
    charges.py
    core.py
    methods.py
    parameters.py
)
foreach(package_source IN LISTS CHARGEFW_PYTHON_PACKAGE_SOURCES)
    configure_file("chargefw/${package_source}"
                   "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/${package_source}" COPYONLY)
endforeach()
configure_file(chargefw/py.typed "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}/py.typed" COPYONLY)

nanobind_add_module(chargefw_python NB_STATIC
    src/calculation.cpp
    src/core.cpp
    src/methods.cpp
    src/module.cpp
    src/parameters.cpp
)
set_target_properties(chargefw_python PROPERTIES
    OUTPUT_NAME _chargefw
    LIBRARY_OUTPUT_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY "${CHARGEFW_PYTHON_BUILD_PACKAGE_DIRECTORY}"
    BUILD_RPATH "$ORIGIN/../.."
    INSTALL_RPATH "$ORIGIN"
)
target_link_libraries(chargefw_python PRIVATE chargefw::core)
target_link_options(chargefw_python PRIVATE "LINKER:--disable-new-dtags")

if(CHARGEFW_ENABLE_WARNINGS)
    chargefw_enable_warnings(chargefw_python)
endif()

if(CHARGEFW_ENABLE_ASAN)
    chargefw_enable_sanitizer(chargefw_python address)
endif()

if(CHARGEFW_ENABLE_UBSAN)
    chargefw_enable_sanitizer(chargefw_python undefined)
endif()

# A wheel needs the library beside the private extension. The ordinary native installation remains in
# ${CMAKE_INSTALL_LIBDIR}; this additional copy exists only when Python support is requested.
install(TARGETS chargefw_python
    LIBRARY DESTINATION chargefw COMPONENT python
    RUNTIME DESTINATION chargefw COMPONENT python
)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/python/chargefw/_chargefw"
        DESTINATION chargefw
        COMPONENT python
        FILES_MATCHING PATTERN "*.pyi")
install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/parameters/"
        DESTINATION chargefw/_data/parameters
        COMPONENT python
        FILES_MATCHING PATTERN "*.json")
install(TARGETS chargefw_core
    LIBRARY DESTINATION chargefw COMPONENT python
    RUNTIME DESTINATION chargefw COMPONENT python
)

if(TARGET gemmi_cpp)
    install(TARGETS gemmi_cpp
        LIBRARY DESTINATION chargefw COMPONENT python
        RUNTIME DESTINATION chargefw COMPONENT python
    )
endif()

if(TARGET tbb)
    install(TARGETS tbb
        LIBRARY DESTINATION chargefw COMPONENT python
        RUNTIME DESTINATION chargefw COMPONENT python
    )
endif()
