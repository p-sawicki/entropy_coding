#ifndef ENTROPY_CODEC_PICTURE
#define ENTROPY_CODEC_PICTURE

#include <deque>

#include "unit.hpp"

namespace EntropyCoding {

class CodingStructure;

#define M_BUFS(JID, PID) m_bufs[PID]

struct Picture : public UnitArea {
  CodingStructure *cs;
  ::std::array<::std::vector<SAOBlkParam>, 2> m_sao;
  ::std::array<::std::vector<uint8_t>, MAX_NUM_COMPONENT> m_alfCtuEnableFlag;
  ::std::vector<short> m_alfCtbFilterIndex;
  ::std::array<::std::vector<uint8_t>, MAX_NUM_COMPONENT> m_alfCtuAlternative;

  SAOBlkParam *getSAO(int id = 0) { return &m_sao[id][0]; };

  uint8_t *getAlfCtuEnableFlag(int compIdx) {
    return m_alfCtuEnableFlag[compIdx].data();
  }

  short *getAlfCtbFilterIndex() { return m_alfCtbFilterIndex.data(); }

  uint8_t *getAlfCtuAlternativeData(int compIdx) {
    return m_alfCtuAlternative[compIdx].data();
  }
};
} // namespace EntropyCoding

#endif
