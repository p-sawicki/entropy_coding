#ifndef __ALFPARAMETERS__
#define __ALFPARAMETERS__

#include "common_def.hpp"

namespace EntropyCoding {
struct AlfParam {
  bool enabledFlag[MAX_NUM_COMPONENT]; // alf_slice_enable_flag, alf_chroma_idc
  int numAlternativesChroma;           // alf_chroma_num_alts_minus_one + 1
};

struct CcAlfFilterParam {
  bool ccAlfFilterEnabled[2];
  uint8_t ccAlfFilterCount[2];
};
} // namespace EntropyCoding

#endif // end of #ifndef  __ALFPARAMETERS__
