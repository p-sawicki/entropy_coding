cmake_minimum_required(VERSION 3.10)

project(codec)

set(CMAKE_CXX_VERSION 17)

option(ENABLE_LOGGING "Enable logging of syntax elements" OFF)
option(STANDALONE_ENTROPY_CODEC "Use standalone entropy codec lib instead of the built-in one" OFF)
message("ENABLE_LOGGING = ${ENABLE_LOGGING}")
message("STANDALONE_ENTROPY_CODEC = ${STANDALONE_ENTROPY_CODEC}")

set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

include_directories(
  googletest/googletest/include
  src/binarizer
  src/common
  src/entropy_codec
  src/log
)

add_library(LogLib src/log/log.cpp)
if(ENABLE_LOGGING)
  target_compile_definitions(LogLib PUBLIC ENABLE_LOGGING)
endif()

add_subdirectory(googletest)
add_subdirectory(VVCSoftware_VTM)

add_library(EntropyCodingLib
  src/entropy_codec/arith_codec.cpp
  src/common/bit_stream.cpp
  src/common/contexts.cpp
  src/binarizer/cabac_reader.cpp
  src/binarizer/cabac_writer.cpp
  src/binarizer/coding_structure.cpp
  src/binarizer/context_modelling.cpp
  src/binarizer/mv.cpp
  src/binarizer/rom.cpp
  src/binarizer/slice.cpp
  src/binarizer/unit_partitioner.cpp
  src/binarizer/unit_tools.cpp
  src/binarizer/unit.cpp
)
target_link_libraries(EntropyCodingLib LogLib)