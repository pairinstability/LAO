set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_executable(runexamples main.cpp)
target_link_libraries(runexamples PRIVATE lao)