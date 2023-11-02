cmake_minimum_required(VERSION 3.22)
project(templateISD CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
# why is the -lpng or png needed at the end (it is in the main repo, but not in the cum4ri repo..??? libm4ri.a
# alternate names: m4ri or libm4ri.a, png or -lpng
set(LINK_TEST_FLAGS libm4ri.a png)
link_directories("deps/m4ri/.libs")

include_directories("challenges/*")
include_directories("deps/m4ri")
include_directories("algorithms")
include_directories("misc")

link_directories("deps/m4ri/.libs")

add_executable(main main.cpp)
target_link_libraries(main ${LINK_TEST_FLAGS})