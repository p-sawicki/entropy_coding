#ifndef __SAMPLEADAPTIVEOFFSET__
#define __SAMPLEADAPTIVEOFFSET__

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
