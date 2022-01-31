#ifndef COMMON_ALF_PARAMETERS
#define COMMON_ALF_PARAMETERS

#include "common_def.hpp"

namespace Common {
struct AlfParam {
  std::array<bool, MAX_NUM_COMPONENT> enabledFlag; // alf_slice_enable_flag, alf_chroma_idc
  int numAlternativesChroma;           // alf_chroma_num_alts_minus_one + 1
};

struct CcAlfFilterParam {
  std::array<bool, 2> ccAlfFilterEnabled;
  std::array<uint8_t, 2> ccAlfFilterCount;
};
} // namespace Common

#endif // end of #ifndef  COMMON_ALF_PARAMETERS
