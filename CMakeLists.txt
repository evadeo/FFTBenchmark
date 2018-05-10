cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/)
# find_package(Boost COMPONENTS graph program_options REQUIRED)
find_package(TBB)
include_directories(${TBB_INCLUDE_DIRS})
add_definitions(${TBB_DEFINITIONS})
# include_directories(${Boost_INCLUDE_DIR})

set(PROJECT_NAME fft)
set(${PROJECT_NAME} CXX)

#CXX
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -std=c++17 -pedantic -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

#INCLUDE
set(PROJECT_INCLUDE src/)

include_directories(${PROJECT_INCLUDE})
set(PROJECT_SRC_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/recursive_fft.cc
)

add_library(_${PROJECT_NAME} SHARED ${PROJECT_SRC_FILES})
add_executable(${PROJECT_NAME} ${PROJECT_SRC_FILES})

target_link_libraries(${PROJECT_NAME}
  _${PROJECT_NAME})

target_link_libraries(_${PROJECT_NAME}
#   ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${TBB_LIBRARIES})
