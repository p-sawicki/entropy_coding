#ifndef __PICTURE__
#define __PICTURE__

#include <deque>

#include "buffer.hpp"
#include "coding_structure.hpp"
#include "slice.hpp"
#include "unit.hpp"

namespace EntropyCoding {

class SEI;
class AQpLayer;

typedef ::std::list<SEI *> SEIMessages;

#define M_BUFS(JID, PID) m_bufs[PID]

struct Picture : public UnitArea {
  CodingStructure *cs;

  SAOBlkParam *getSAO(int id = 0) { return &m_sao[id][0]; };

  ::std::vector<SAOBlkParam> m_sao[2];

  ::std::vector<uint8_t> m_alfCtuEnableFlag[MAX_NUM_COMPONENT];
  uint8_t *getAlfCtuEnableFlag(int compIdx) {
    return m_alfCtuEnableFlag[compIdx].data();
  }

  ::std::vector<short> m_alfCtbFilterIndex;
  short *getAlfCtbFilterIndex() { return m_alfCtbFilterIndex.data(); }

  ::std::vector<uint8_t> m_alfCtuAlternative[MAX_NUM_COMPONENT];
  uint8_t *getAlfCtuAlternativeData(int compIdx) {
    return m_alfCtuAlternative[compIdx].data();
  }
};
} // namespace EntropyCoding

#endif
