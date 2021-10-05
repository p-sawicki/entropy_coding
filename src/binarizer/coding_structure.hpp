#ifndef __CODINGSTRUCTURE__
#define __CODINGSTRUCTURE__

#include <vector>

#include "buffer.hpp"
#include "slice.hpp"
#include "unit.hpp"
#include "unit_partitioner.hpp"

namespace EntropyCoding {

struct Picture;

enum PictureType {
  PIC_RECONSTRUCTION = 0,
  PIC_ORIGINAL,
  PIC_TRUE_ORIGINAL,
  PIC_FILTERED_ORIGINAL,
  PIC_PREDICTION,
  PIC_RESIDUAL,
  PIC_ORG_RESI,
  PIC_RECON_WRAP,
  PIC_ORIGINAL_INPUT,
  PIC_TRUE_ORIGINAL_INPUT,
  PIC_FILTERED_ORIGINAL_INPUT,
  NUM_PIC_TYPES
};

// ---------------------------------------------------------------------------
// coding structure
// ---------------------------------------------------------------------------

class CodingStructure {
public:
  UnitArea area;

  Picture *picture;
  CodingStructure *parent;
  Slice *slice;

  UnitScale unitScale[MAX_NUM_COMPONENT];
  int chromaQpAdj;
  const SPS *sps;
  const PPS *pps;
  PicHeader *picHeader;
  const PreCalcValues *pcv;

  const CodingUnit *getCU(const Position &pos, const ChannelType _chType) const;
  const PredictionUnit *getPU(const Position &pos,
                              const ChannelType _chType) const;
  const TransformUnit *getTU(const Position &pos, const ChannelType _chType,
                             const int subTuIdx = -1) const;

  CodingUnit *getCU(const Position &pos, const ChannelType _chType);
  CodingUnit *getLumaCU(const Position &pos);
  PredictionUnit *getPU(const Position &pos, const ChannelType _chType);
  TransformUnit *getTU(const Position &pos, const ChannelType _chType,
                       const int subTuIdx = -1);

  const CodingUnit *getCURestricted(const Position &pos, const Position curPos,
                                    const unsigned curSliceIdx,
                                    const unsigned curTileIdx,
                                    const ChannelType _chType) const;
  const CodingUnit *getCURestricted(const Position &pos,
                                    const CodingUnit &curCu,
                                    const ChannelType _chType) const;
  const PredictionUnit *getPURestricted(const Position &pos,
                                        const PredictionUnit &curPu,
                                        const ChannelType _chType) const;

  CodingUnit &addCU(const UnitArea &unit, const ChannelType _chType);
  PredictionUnit &addPU(const UnitArea &unit, const ChannelType _chType);
  TransformUnit &addTU(const UnitArea &unit, const ChannelType _chType);
  void addEmptyTUs(Partitioner &partitioner);

  TreeType treeType; // because partitioner can not go deep to tu and cu coding
                     // (e.g., addCU()), need another variable for indicating
                     // treeType
  ModeType modeType;

  const int signalModeCons(const PartSplit split, Partitioner &partitioner,
                           const ModeType modeTypeParent) const;
  ::std::vector<CodingUnit *> cus;
  ::std::vector<PredictionUnit *> pus;
  ::std::vector<TransformUnit *> tus;

  PLTBuf prevPLT;
  void reorderPrevPLT(PLTBuf &prevPLT, uint8_t curPLTSize[MAX_NUM_CHANNEL_TYPE],
                      Pel curPLT[MAX_NUM_COMPONENT][MAXPLTSIZE],
                      bool reuseflag[MAX_NUM_CHANNEL_TYPE][MAXPLTPREDSIZE],
                      uint32_t compBegin, uint32_t numComp, bool jointPLT);

private:
  // needed for TU encoding
  bool m_isTuEnc;

  unsigned *m_cuIdx[MAX_NUM_CHANNEL_TYPE];
  unsigned *m_puIdx[MAX_NUM_CHANNEL_TYPE];
  unsigned *m_tuIdx[MAX_NUM_CHANNEL_TYPE];

  unsigned m_numCUs;
  unsigned m_numPUs;
  unsigned m_numTUs;

  CUCache &m_cuCache;
  PUCache &m_puCache;
  TUCache &m_tuCache;

  TCoeff *m_coeffs[MAX_NUM_COMPONENT];
  Pel *m_pcmbuf[MAX_NUM_COMPONENT];
  bool *m_runType[MAX_NUM_CHANNEL_TYPE];
  int m_offsets[MAX_NUM_COMPONENT];
};

static inline uint32_t getNumberValidTBlocks(const PreCalcValues &pcv) {
  return (pcv.chrFormat == CHROMA_400)
             ? 1
             : (pcv.multiBlock422 ? MAX_NUM_TBLOCKS : MAX_NUM_COMPONENT);
}
} // namespace EntropyCoding

#endif
