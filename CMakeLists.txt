cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/)
# find_package(Boost COMPONENTS graph program_options REQUIRED)
find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})
# include_directories(${Boost_INCLUDE_DIR})

cmake_policy(SET CMP0042 NEW)
set(PROJECT_NAME fft)
set(${PROJECT_NAME} CXX)

#CXX
set(CMAKE_CXX_COMPILER "clang++")

#MAC OS SETTINGS
#set(CMAKE_CXX_COMPILER "g++-8")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -std=c++14 -fopenmp -pedantic -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")

#INCLUDE
set(PROJECT_INCLUDE src/)

include_directories(${PROJECT_INCLUDE})

set(PROJECT_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterative_fft.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterative_fft_parallel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/iterative_fft_openmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/stockham_openmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/recursive_fft.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parallel_fft.cc
)

add_library(_${PROJECT_NAME} SHARED ${PROJECT_SRC_FILES})
add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cc)

target_link_libraries(${PROJECT_NAME}
  _${PROJECT_NAME})

target_link_libraries(_${PROJECT_NAME}
#   ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${TBB_LIBRARIES})
