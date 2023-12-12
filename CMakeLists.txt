cmake_minimum_required(VERSION 3.22)
project(templateISD CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_VERBOSE_MAKEFILE ON)

# why is the -lpng or png needed at the end (it is in the main repo, but not in the cum4ri repo..??? libm4ri.a
# alternate names: m4ri or libm4ri.a, png or -lpng
set(LINK_TEST_FLAGS libm4ri.a png)

# automatically check if current cpu supports certain avx instruction sets
include("cmake/OptimizeForArchitecture.cmake")
OptimizeForArchitecture()

# list ceated in OptimizeForArchitecture
if("avx2" IN_LIST _available_vector_units_list)
    message("-- AVX2 available, setting appropriate flags.")
	set(AVX2_FLAGS "-DUSE_AVX2 -mavx2 -mavx -mbmi2 ")
endif()

# not needed for now
# if("avx512f" IN_LIST _available_vector_units_list)
#     message("-- AVX512 available, setting appropriate flags.")
#     set(AVX512_FLAGS "-DUSE_AVX512 -mavx512f -mavx512dq  -mavx512bf16 -mavx512cd  -mavx512bw -mavx512vl -mavx512vbmi -mavx512ifma -mavx512vbmi2 -mavx512vnni -mavx512bitalg -mavx512vpopcntdq")
# endif()


set(OPTIMIZE_FLAGS "-g -march=native -O3 -DNDEBUG -ffast-math -ftree-vectorize -funroll-loops -fno-exceptions ${AVX2_FLAGS} ${AVX512_FLAGS}")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPTIMIZE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMIZE_FLAGS}")


link_directories("deps/m4ri/.libs")
include_directories("challenges/*")
include_directories("deps/m4ri")
include_directories("algorithms")
include_directories("misc")

# original, single threaded (working) template prange
add_executable(main main.cpp)
target_link_libraries(main ${LINK_TEST_FLAGS})

# this line is needed together with target_link_libraries(main OpenMP::OpenMP_CXX), but no manual -fopenmp is required...
find_package(OpenMP)

# multithreaded template prange
if(OpenMP_CXX_FOUND)
    add_executable(prange prange.cpp)
    target_link_libraries(prange ${LINK_TEST_FLAGS})
    target_link_libraries(prange OpenMP::OpenMP_CXX)

    add_executable(dumer dumer.cpp)
    target_link_libraries(dumer ${LINK_TEST_FLAGS})
    target_link_libraries(dumer OpenMP::OpenMP_CXX)
endif()