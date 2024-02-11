# google test cmake for building the dependency

include(FetchContent)

find_program(CLANGXX_COMPILER NAMES clang++)
find_program(CLANG_COMPILER NAMES clang)

if(CLANGXX_COMPILER AND CLANG_COMPILER)
    set(CMAKE_CXX_COMPILER "${CLANGXX_COMPILER}")
    set(CMAKE_C_COMPILER "${CLANG_COMPILER}")
endif()

find_program(MOLD_LINKER NAMES mold)
if(MOLD_LINKER)
    set(CMAKE_LINKER "${MOLD_LINKER}")
endif()

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)