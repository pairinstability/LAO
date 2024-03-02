# building unit testing with googletest

# google test 
include(../cmake/GoogleTest.cmake)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

enable_testing()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

function(create_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE lao GTest::GTest)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

create_test(test_dense_matrix)
create_test(test_matrix_math)
create_test(test_astro_body)