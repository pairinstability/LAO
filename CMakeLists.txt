cmake_minimum_required(VERSION 3.22)
project(lao LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library(lao INTERFACE)
target_include_directories(lao INTERFACE
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_compile_options(lao INTERFACE
# Compiler is GNU:
$<$<CXX_COMPILER_ID:GNU>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wpedantic>>
# Compiler is Clang:
$<$<CXX_COMPILER_ID:Clang>:$<BUILD_INTERFACE:-Wall;-Wextra;-Wpedantic;-Wno-c99-extensions>>
)

option(BUILD_TESTS "Build unit tests" ON)
if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

#option(BUILD_BENCH "Build benchmarking" ON)
#if (BUILD_BENCH)
#    add_subdirectory(bench)
#endif()

option(BUILD_EXAMPLES "Build examples" ON)
if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()