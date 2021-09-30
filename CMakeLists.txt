cmake_minimum_required(VERSION 3.10)

project(codec)

set(CMAKE_CXX_VERSION 17)

option(ENABLE_LOGGING "Enable logging of syntax elements" OFF)
message("ENABLE_LOGGING = ${ENABLE_LOGGING}")

set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

include_directories(
  googletest/googletest/include
  src
)

add_library(LogLib src/log.cpp)

add_subdirectory(googletest)
add_subdirectory(VVCSoftware_VTM)

add_library(EntropyCoderLib
  src/arith_codec.cpp
  src/bit_stream.cpp
  src/contexts.cpp
)