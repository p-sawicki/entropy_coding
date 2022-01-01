#ifndef ENTROPY_CODEC_CODING_STRUCTURE
#define ENTROPY_CODEC_CODING_STRUCTURE

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

  ::std::shared_ptr<Picture> picture;
  ::std::shared_ptr<CodingStructure> parent;
  ::std::shared_ptr<Slice> slice;

  ::std::array<UnitScale, MAX_NUM_COMPONENT> unitScale;
  int chromaQpAdj;
  ::std::shared_ptr<const SPS> sps;
  ::std::shared_ptr<const PPS> pps;
  ::std::shared_ptr<PicHeader> picHeader;
  ::std::shared_ptr<const PreCalcValues> pcv;

  CodingStructure(const UnitArea &_area,
                  const ::std::array<UnitScale, MAX_NUM_COMPONENT> &_unitScale,
                  const int _chromaQpAdj,
                  const TreeType _treeType, const ModeType _modeType,
                  const PLTBuf &_prevPLT, bool isTuEnc,
                  unsigned *const cuIdx[MAX_NUM_CHANNEL_TYPE],
                  unsigned *const puIdx[MAX_NUM_CHANNEL_TYPE],
                  unsigned *const tuIdx[MAX_NUM_CHANNEL_TYPE],
                  const unsigned numCUs, const unsigned numPUs,
                  const unsigned numTUs, TCoeff *const coeffs[MAX_NUM_COMPONENT],
                  Pel *const pcmbuf[MAX_NUM_COMPONENT],
                  bool *const runType[MAX_NUM_CHANNEL_TYPE], const int *offsets)
      : area(_area),
        unitScale(_unitScale), chromaQpAdj(_chromaQpAdj), treeType(_treeType),
        modeType(_modeType), prevPLT(_prevPLT),
        m_isTuEnc(isTuEnc), m_numCUs(numCUs), m_numPUs(numPUs),
        m_numTUs(numTUs) {
    copy_array(cuIdx, m_cuIdx);
    copy_array(puIdx, m_puIdx);
    copy_array(tuIdx, m_tuIdx);

    copy_array(coeffs, m_coeffs);
    copy_array(pcmbuf, m_pcmbuf);
    copy_array(runType, m_runType);
    copy_array(offsets, m_offsets);
  }

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
  ::std::vector<::std::shared_ptr<CodingUnit>> cus;
  ::std::vector<::std::shared_ptr<PredictionUnit>> pus;
  ::std::vector<::std::shared_ptr<TransformUnit>> tus;

  PLTBuf prevPLT;
  void reorderPrevPLT(PLTBuf &prevPLT,
                      ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> &curPLTSize,
                      CurPLT31 &curPLT, ReuseFlag &reuseflag,
                      uint32_t compBegin, uint32_t numComp, bool jointPLT);

  // needed for TU encoding
  bool m_isTuEnc;

  ::std::array<unsigned *, MAX_NUM_CHANNEL_TYPE> m_cuIdx;
  ::std::array<unsigned *, MAX_NUM_CHANNEL_TYPE> m_puIdx;
  ::std::array<unsigned *, MAX_NUM_CHANNEL_TYPE> m_tuIdx;

  unsigned m_numCUs;
  unsigned m_numPUs;
  unsigned m_numTUs;

  ::std::shared_ptr<CUCache> m_cuCache;
  ::std::shared_ptr<PUCache> m_puCache;
  ::std::shared_ptr<TUCache> m_tuCache;

  ::std::array<TCoeff *, MAX_NUM_COMPONENT> m_coeffs;
  ::std::array<Pel *, MAX_NUM_COMPONENT> m_pcmbuf;
  ::std::array<bool *, MAX_NUM_CHANNEL_TYPE> m_runType;
  ::std::array<int, MAX_NUM_COMPONENT> m_offsets;
};

static inline uint32_t getNumberValidTBlocks(const PreCalcValues &pcv) {
  return (pcv.chrFormat == CHROMA_400)
             ? 1
             : (pcv.multiBlock422 ? MAX_NUM_TBLOCKS : MAX_NUM_COMPONENT);
}
} // namespace EntropyCoding

#endif
