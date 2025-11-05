# Star-Magic UQFF CMake configuration
cmake_minimum_required(VERSION 3.16)

project(StarMagicUQFF LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Place build outputs in a dedicated folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Helper target: Verify enhanced module behaviour (self-contained mock tests)
add_executable(verify_enhancements
    verify_enhancements.cpp
)

# Lightweight smoke-test target that exercises one enhanced module
add_executable(test_enhanced_modules
    test_enhanced_modules.cpp
)

target_include_directories(test_enhanced_modules PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_compile_definitions(test_enhanced_modules PRIVATE _USE_MATH_DEFINES)

target_include_directories(verify_enhancements PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
