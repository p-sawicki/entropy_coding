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

  Picture *picture;
  CodingStructure *parent;
  Slice *slice;

  ::std::array<UnitScale, MAX_NUM_COMPONENT> unitScale;
  int chromaQpAdj;
  const SPS *sps;
  const PPS *pps;
  PicHeader *picHeader;
  const PreCalcValues *pcv;

  CodingStructure(const UnitArea &_area, Picture *_picture,
                  CodingStructure *_parent, Slice *_slice,
                  const ::std::array<UnitScale, MAX_NUM_COMPONENT> &_unitScale,
                  const int _chromaQpAdj, const SPS *_sps, const PPS *_pps,
                  PicHeader *_picHeader, const PreCalcValues *_pcv,
                  const TreeType _treeType, const ModeType _modeType,
                  const ::std::vector<CodingUnit *> &_cus,
                  const ::std::vector<PredictionUnit *> &_pus,
                  const ::std::vector<TransformUnit *> &_tus,
                  const PLTBuf &_prevPLT, bool isTuEnc,
                  unsigned *const cuIdx[MAX_NUM_CHANNEL_TYPE],
                  unsigned *const puIdx[MAX_NUM_CHANNEL_TYPE],
                  unsigned *const tuIdx[MAX_NUM_CHANNEL_TYPE],
                  const unsigned numCUs, const unsigned numPUs,
                  const unsigned numTUs, CUCache *cuCache, PUCache *puCache,
                  TUCache *tuCache, TCoeff *const coeffs[MAX_NUM_COMPONENT],
                  Pel *const pcmbuf[MAX_NUM_COMPONENT],
                  bool *const runType[MAX_NUM_CHANNEL_TYPE], const int *offsets)
      : area(_area), picture(_picture), parent(_parent), slice(_slice),
        unitScale(_unitScale), chromaQpAdj(_chromaQpAdj), sps(_sps), pps(_pps),
        picHeader(_picHeader), pcv(_pcv), treeType(_treeType),
        modeType(_modeType), cus(_cus), pus(_pus), tus(_tus), prevPLT(_prevPLT),
        m_isTuEnc(isTuEnc), m_numCUs(numCUs), m_numPUs(numPUs),
        m_numTUs(numTUs), m_cuCache(cuCache), m_puCache(puCache),
        m_tuCache(tuCache) {
    picture->cs = this;
    slice->setPic(picture);
    slice->setSPS(sps);
    slice->setPPS(pps);
    slice->setPicHeader(picHeader);

    ::std::for_each(cus.begin(), cus.end(), [this](CodingUnit *cu) {
      cu->cs = this;
      cu->slice = slice;
    });

    ::std::for_each(pus.begin(), pus.end(),
                    [this](PredictionUnit *pu) { pu->cs = this; });

    ::std::for_each(tus.begin(), tus.end(),
                    [this](TransformUnit *tu) { tu->cs = this; });

    copy_array(cuIdx, m_cuIdx);
    copy_array(puIdx, m_puIdx);
    copy_array(tuIdx, m_tuIdx);

    copy_array(coeffs, m_coeffs);
    copy_array(pcmbuf, m_pcmbuf);
    copy_array(runType, m_runType);
    copy_array(offsets, m_offsets);
  }

  ~CodingStructure() {
    delete picture;
    delete parent;

    ::std::for_each(cus.begin(), cus.end(), [](CodingUnit *cu) { delete cu; });

    ::std::for_each(pus.begin(), pus.end(),
                    [](PredictionUnit *pu) { delete pu; });

    ::std::for_each(tus.begin(), tus.end(),
                    [](TransformUnit *tu) { delete tu; });

    delete m_cuCache;
    delete m_puCache;
    delete m_tuCache;
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
  ::std::vector<CodingUnit *> cus;
  ::std::vector<PredictionUnit *> pus;
  ::std::vector<TransformUnit *> tus;

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

  CUCache *m_cuCache;
  PUCache *m_puCache;
  TUCache *m_tuCache;

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
