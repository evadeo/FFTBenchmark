cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/)
# find_package(Boost COMPONENTS graph program_options REQUIRED)
find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})
# include_directories(${Boost_INCLUDE_DIR})

#find_package(OpenMP REQUIRED)
find_package(CUDA REQUIRED)

cmake_policy(SET CMP0042 NEW)
set(PROJECT_NAME fft)
set(${PROJECT_NAME} CXX CUDA)
set(CMAKE_VERBOSE_MAKEFILE ON)

#CXX
set(CMAKE_CUDA_COMPILER "nvcc")
set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -std=c++14 -pedantic -Ofast -mavx -fopenmp")


set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
enable_language("CUDA")

#INCLUDE
set(PROJECT_INCLUDE src/)

include_directories(${PROJECT_INCLUDE})

set(PROJECT_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterative_fft.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterative_fft_parallel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/recursive_fft.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parallel_fft.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fftw.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/stockham_simple.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/recursive_fft_sse.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/cudafft.cu
)

add_library(_${PROJECT_NAME} ${PROJECT_SRC_FILES})
# add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc)
cuda_add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc)

target_link_libraries(${PROJECT_NAME}
  ${TBB_LIBRARIES}
  _${PROJECT_NAME})

target_link_libraries(_${PROJECT_NAME}
  #   ${Boost_PROGRAM_OPTIONS_LIBRARY}
  fftw3
  ${TBB_LIBRARIES})
