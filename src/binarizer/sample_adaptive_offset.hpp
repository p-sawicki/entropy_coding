#ifndef ENTROPY_CODEC_SAMPLE_ADAPTIVE_OFFSET
#define ENTROPY_CODEC_SAMPLE_ADAPTIVE_OFFSET

#include "unit.hpp"

namespace EntropyCoding {

#define MAX_SAO_TRUNCATED_BITDEPTH 10

class SampleAdaptiveOffset {
public:
  static int getMaxOffsetQVal(const int channelBitDepth) {
    return (1 << (::std::min<int>(channelBitDepth, MAX_SAO_TRUNCATED_BITDEPTH) -
                  5)) -
           1;
  } // Table 9-32, inclusive
};
} // namespace EntropyCoding

#endif
