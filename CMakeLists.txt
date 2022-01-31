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

add_subdirectory(VVCSoftware_VTM)

add_library(Common
  src/common/bit_stream.cpp
  src/common/contexts.cpp
  src/common/coding_structure.cpp
  src/common/context_modelling.cpp
  src/common/mv.cpp
  src/common/rom.cpp
  src/common/slice.cpp
  src/common/unit_partitioner.cpp
  src/common/unit_tools.cpp
  src/common/unit.cpp
)

add_library(EntropyCodingLib
  src/entropy_codec/arith_codec.cpp
  src/entropy_codec/cabac_reader.cpp
  src/entropy_codec/cabac_writer.cpp
)

target_link_libraries(EntropyCodingLib LogLib Common)

include(CTest)

function(add_regression_test FILE QP FR CHROMA ENCODE_LOG_HASH DECODE_LOG_HASH ENCODE_OUT_HASH DECODE_OUT_HASH)
  set(ENCODE_TEST ${FILE}_${QP}_ENCODE_TEST)
  set(DECODE_TEST ${FILE}_${QP}_DECODE_TEST)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
  set(TEST_SCRIPT ${TEST_DIR}/run_test.cmake)
  set(ENCODE_ARGS "-i ${TEST_DIR}/${FILE} -c ${TEST_DIR}/encoder_intra_vtm_qp${QP}.cfg -fr ${FR} -f 10 -wdt 176 -hgt 144 --InputChromaFormat=${CHROMA}")
  set(DECODE_ARGS "-b str.bin -o recon.yuv")
  
  add_test(NAME ${ENCODE_TEST}
    COMMAND ${CMAKE_COMMAND}
            -DCMD=$<TARGET_FILE:EncoderApp>
            -DARGS=${ENCODE_ARGS}
            -DLOG_FILE=bin_log.txt
            -DOUT_FILE=str.bin
    -P ${TEST_SCRIPT}
  )
  set_property(TEST ${ENCODE_TEST}
    PROPERTY PASS_REGULAR_EXPRESSION ".*${ENCODE_LOG_HASH}.*bin_log.txt.*${ENCODE_OUT_HASH}.*str.bin.*"
  )

  add_test(NAME ${DECODE_TEST}
    COMMAND ${CMAKE_COMMAND}
            -DCMD=$<TARGET_FILE:DecoderApp>
            -DARGS=${DECODE_ARGS}
            -DLOG_FILE=bin_log.txt
            -DOUT_FILE=recon.yuv
    -P ${TEST_SCRIPT}  
  )
  set_property(TEST ${DECODE_TEST}
    PROPERTY PASS_REGULAR_EXPRESSION ".*${DECODE_LOG_HASH}.*bin_log.txt.*${DECODE_OUT_HASH}.*recon.yuv.*"
  )
endfunction()

# short tests (QP=51)
add_regression_test(akiyo_qcif.y4m           51 60 422 ec240e3870d1573352e77b1e4cbd8456 c0803b5707712acf52976c64b336a5d9 1d05000ad73270b7acf3ac69fdfcf11d a0b9a7b9e415f6b1e60b9c9949d60763)
add_regression_test(mother_daughter_qcif.y4m 51 60 420 87a68e22e0474c9e3fb317df9c4d0e54 55f42d268fc72474557d81a5c8a0c811 8ae78b0364c16c17e9bf8e06c2341be0 4247f7b7f0b9321153b9f80707c44b93)
add_regression_test(foreman_qcif.y4m         51 50 422 b738547ffd20318533742c8129e07d78 a93123cfcfa35f4afac8ef85b6590e3e f6d7aa2946d77d36fa4eba1dca550a8b 49bade579717b9d9a864cbfa46097767)

# long tests (QP=32)
# add_regression_test(akiyo_qcif.y4m           32 60 422 e5d53e83f2b7aadf3022493e5538720c 527d76b47ce224680a528b7b85ac637a 701e32c475d14ced3d87d56155ce6326 d02855846da02a1fad1f49d5ecce3a49)
# add_regression_test(mother_daughter_qcif.y4m 32 60 420 cb07aa6057367e2147c466aa2de8c8a6 8eafc6db6f212cee67fe95eda2ba746d e8bd7c96e6ca65d55ba03614468a7b31 f999dd3c573959b8011a281f18f47170)
# add_regression_test(foreman_qcif.y4m         32 50 422 2030774c80eef75535520d007d6133fe e598fef37767657c46fe704bb98f52c9 4d9b27c3105dd72ca07ca1600ef32072 72d4e806a3ff22948ee29037e964638b)

# very long tests (QP=0)
# add_regression_test(akiyo_qcif.y4m           0  60 422 e966511982d0936a9ca9883a35f517ce c081c2ee883fdc09deb98221fa457c3e 5a61b710458457aabf2e44e3131e088b f4e85576bfe9d5f9add099a634b7b3e6)
# add_regression_test(mother_daughter_qcif.y4m 0  60 420 6c8b606ee72da8d5da2770bb7bcd5e8e 13675b49f48e0b507f5d4ee161e7bc84 08ac3dbfb6567b3ee507c2f915648bb4 cb1fad9a6ca3bce75f57f31e9c73921d)
# add_regression_test(foreman_qcif.y4m         0  50 422 64f62fdbc3beb058a2f67c3e60005f36 4399c57b8943d373171fcbe3237bf2e3 8c11258b22fb565e124edc393b4d786b 7d75ebb1730909beb00fa8fa2237f6c8)
