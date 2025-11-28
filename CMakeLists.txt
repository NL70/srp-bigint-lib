cmake_minimum_required(VERSION 3.15)
project(BigIntProject LANGUAGES CXX)

# Compilation options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

# Include directories 
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/BigInt.cpp
    src/BigIntUtils.cpp
)

# Try to find GMP automatically
find_path(GMP_INCLUDE_DIR gmpxx.h)
find_library(GMP_LIBRARY gmp)
find_library(GMPXX_LIBRARY gmpxx)

if (NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "GMP or GMPXX not found. Please install libgmp-dev and libgmpxx4ldbl.")
else()
    message(STATUS "Found GMP include dir: ${GMP_INCLUDE_DIR}")
    message(STATUS "Found GMP libraries: ${GMP_LIBRARY}, ${GMPXX_LIBRARY}")
endif()

include_directories(${GMP_INCLUDE_DIR})

# Google Benchmark
# Assume benchmark was cloned at ${PROJECT_SOURCE_DIR}/benchmark
set(BENCHMARK_ROOT ${PROJECT_SOURCE_DIR}/benchmark)
set(BENCHMARK_INCLUDE_DIR ${BENCHMARK_ROOT}/include)
set(BENCHMARK_LIBRARY ${BENCHMARK_ROOT}/build/src/libbenchmark.a)

if (EXISTS ${BENCHMARK_LIBRARY})
    message(STATUS "Using local Google Benchmark at ${BENCHMARK_LIBRARY}")
    include_directories(${BENCHMARK_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Google Benchmark not found. Please build it in ${BENCHMARK_ROOT}/build")
endif()

# Executable target
add_executable(testBigInt
    tests/testBigInt.cpp
    ${SOURCES}
)

# Link against GMP and Google Benchmark
target_link_libraries(testBigInt
    ${GMPXX_LIBRARY}
    ${GMP_LIBRARY}
    ${BENCHMARK_LIBRARY}
    pthread
)

# Output directory
set_target_properties(testBigInt PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
