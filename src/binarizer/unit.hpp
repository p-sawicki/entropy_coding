#ifndef __UNIT__
#define __UNIT__

#include <assert.h>

#include "buffer.hpp"
#include "common_def.hpp"
#include "mv.hpp"

namespace EntropyCoding {

// ---------------------------------------------------------------------------
// tools
// ---------------------------------------------------------------------------
struct PLTBuf {
  uint8_t curPLTSize[MAX_NUM_CHANNEL_TYPE];
  Pel curPLT[MAX_NUM_COMPONENT][MAXPLTPREDSIZE];
};
inline Position recalcPosition(const ChromaFormat _cf, const ComponentID srcCId,
                               const ComponentID dstCId, const Position &pos) {
  if (toChannelType(srcCId) == toChannelType(dstCId)) {
    return pos;
  } else if (isLuma(srcCId) && isChroma(dstCId)) {
    return Position(pos.x >> getComponentScaleX(dstCId, _cf),
                    pos.y >> getComponentScaleY(dstCId, _cf));
  } else {
    return Position(pos.x << getComponentScaleX(srcCId, _cf),
                    pos.y << getComponentScaleY(srcCId, _cf));
  }
}

inline Position recalcPosition(const ChromaFormat _cf, const ChannelType srcCHt,
                               const ChannelType dstCHt, const Position &pos) {
  if (srcCHt == dstCHt) {
    return pos;
  } else if (isLuma(srcCHt) && isChroma(dstCHt)) {
    return Position(pos.x >> getChannelTypeScaleX(dstCHt, _cf),
                    pos.y >> getChannelTypeScaleY(dstCHt, _cf));
  } else {
    return Position(pos.x << getChannelTypeScaleX(srcCHt, _cf),
                    pos.y << getChannelTypeScaleY(srcCHt, _cf));
  }
}

inline Size recalcSize(const ChromaFormat _cf, const ComponentID srcCId,
                       const ComponentID dstCId, const Size &size) {
  if (toChannelType(srcCId) == toChannelType(dstCId)) {
    return size;
  } else if (isLuma(srcCId) && isChroma(dstCId)) {
    return Size(size.width >> getComponentScaleX(dstCId, _cf),
                size.height >> getComponentScaleY(dstCId, _cf));
  } else {
    return Size(size.width << getComponentScaleX(srcCId, _cf),
                size.height << getComponentScaleY(srcCId, _cf));
  }
}

inline Size recalcSize(const ChromaFormat _cf, const ChannelType srcCHt,
                       const ChannelType dstCHt, const Size &size) {
  if (srcCHt == dstCHt) {
    return size;
  } else if (isLuma(srcCHt) && isChroma(dstCHt)) {
    return Size(size.width >> getChannelTypeScaleX(dstCHt, _cf),
                size.height >> getChannelTypeScaleY(dstCHt, _cf));
  } else {
    return Size(size.width << getChannelTypeScaleX(srcCHt, _cf),
                size.height << getChannelTypeScaleY(srcCHt, _cf));
  }
}

// ---------------------------------------------------------------------------
// block definition
// ---------------------------------------------------------------------------

struct CompArea : public Area {
  CompArea()
      : Area(), chromaFormat(NUM_CHROMA_FORMAT), compID(MAX_NUM_TBLOCKS) {}
  CompArea(const ComponentID _compID, const ChromaFormat _cf, const Area &_area,
           const bool isLuma = false)
      : Area(_area), chromaFormat(_cf), compID(_compID) {
    if (isLuma)
      xRecalcLumaToChroma();
  }
  CompArea(const ComponentID _compID, const ChromaFormat _cf,
           const Position &_pos, const Size &_size, const bool isLuma = false)
      : Area(_pos, _size), chromaFormat(_cf), compID(_compID) {
    if (isLuma)
      xRecalcLumaToChroma();
  }
  CompArea(const ComponentID _compID, const ChromaFormat _cf, const uint32_t _x,
           const uint32_t _y, const uint32_t _w, const uint32_t _h,
           const bool isLuma = false)
      : Area(_x, _y, _w, _h), chromaFormat(_cf), compID(_compID) {
    if (isLuma)
      xRecalcLumaToChroma();
  }

  ChromaFormat chromaFormat;
  ComponentID compID;

  Position chromaPos() const;
  Position lumaPos() const;

  Size chromaSize() const;
  Size lumaSize() const;

  Position compPos(const ComponentID compID) const;
  Position chanPos(const ChannelType chType) const;

  Position topLeftComp(const ComponentID _compID) const {
    return recalcPosition(chromaFormat, compID, _compID, *this);
  }
  Position topRightComp(const ComponentID _compID) const {
    return recalcPosition(chromaFormat, compID, _compID,
                          {(PosType)(x + width - 1), y});
  }
  Position bottomLeftComp(const ComponentID _compID) const {
    return recalcPosition(chromaFormat, compID, _compID,
                          {x, (PosType)(y + height - 1)});
  }
  Position bottomRightComp(const ComponentID _compID) const {
    return recalcPosition(
        chromaFormat, compID, _compID,
        {(PosType)(x + width - 1), (PosType)(y + height - 1)});
  }

  bool valid() const {
    return chromaFormat < NUM_CHROMA_FORMAT && compID < MAX_NUM_TBLOCKS &&
           width != 0 && height != 0;
  }

  const bool operator==(const CompArea &other) const {
    if (chromaFormat != other.chromaFormat)
      return false;
    if (compID != other.compID)
      return false;

    return Position::operator==(other) && Size::operator==(other);
  }

  const bool operator!=(const CompArea &other) const {
    return !(operator==(other));
  }

#if REUSE_CU_RESULTS_WITH_MULTIPLE_TUS
  void resizeTo(const Size &newSize) { Size::resizeTo(newSize); }
#endif
  void repositionTo(const Position &newPos) { Position::repositionTo(newPos); }
  void positionRelativeTo(const CompArea &origCompArea) {
    Position::relativeTo(origCompArea);
  }

private:
  void xRecalcLumaToChroma();
};

inline CompArea clipArea(const CompArea &compArea, const Area &boundingBox) {
  return CompArea(compArea.compID, compArea.chromaFormat,
                  clipArea((const Area &)compArea, boundingBox));
}

// ---------------------------------------------------------------------------
// unit definition
// ---------------------------------------------------------------------------

typedef static_vector<CompArea, MAX_NUM_TBLOCKS> UnitBlocksType;

struct UnitArea {
  ChromaFormat chromaFormat;
  UnitBlocksType blocks;

  UnitArea() : chromaFormat(NUM_CHROMA_FORMAT) {}
  UnitArea(const ChromaFormat _chromaFormat);
  UnitArea(const ChromaFormat _chromaFormat, const Area &area);
  UnitArea(const ChromaFormat _chromaFormat, const CompArea &blkY);
  UnitArea(const ChromaFormat _chromaFormat, CompArea &&blkY);
  UnitArea(const ChromaFormat _chromaFormat, const CompArea &blkY,
           const CompArea &blkCb, const CompArea &blkCr);
  UnitArea(const ChromaFormat _chromaFormat, CompArea &&blkY, CompArea &&blkCb,
           CompArea &&blkCr);

  CompArea &Y() { return blocks[COMPONENT_Y]; }
  const CompArea &Y() const { return blocks[COMPONENT_Y]; }
  CompArea &Cb() { return blocks[COMPONENT_Cb]; }
  const CompArea &Cb() const { return blocks[COMPONENT_Cb]; }
  CompArea &Cr() { return blocks[COMPONENT_Cr]; }
  const CompArea &Cr() const { return blocks[COMPONENT_Cr]; }

  CompArea &block(const ComponentID comp) { return blocks[comp]; }
  const CompArea &block(const ComponentID comp) const { return blocks[comp]; }

  bool contains(const UnitArea &other) const;
  bool contains(const UnitArea &other, const ChannelType chType) const;

  CompArea &operator[](const int n) { return blocks[n]; }
  const CompArea &operator[](const int n) const { return blocks[n]; }

  const bool operator==(const UnitArea &other) const {
    if (chromaFormat != other.chromaFormat)
      return false;
    if (blocks.size() != other.blocks.size())
      return false;

    for (uint32_t i = 0; i < blocks.size(); i++) {
      if (blocks[i] != other.blocks[i])
        return false;
    }

    return true;
  }

#if REUSE_CU_RESULTS_WITH_MULTIPLE_TUS
  void resizeTo(const UnitArea &unit);
#endif
  void repositionTo(const UnitArea &unit);

  const bool operator!=(const UnitArea &other) const {
    return !(*this == other);
  }

  const Position &lumaPos() const { return Y(); }
  const Size &lumaSize() const { return Y(); }

  const Position &chromaPos() const { return Cb(); }
  const Size &chromaSize() const { return Cb(); }

  const UnitArea singleComp(const ComponentID compID) const;
  const UnitArea singleChan(const ChannelType chType) const;

  const SizeType lwidth() const { return Y().width; }   /*! luma width  */
  const SizeType lheight() const { return Y().height; } /*! luma height */

  const PosType lx() const { return Y().x; } /*! luma x-pos */
  const PosType ly() const { return Y().y; } /*! luma y-pos */

  bool valid() const {
    return chromaFormat != NUM_CHROMA_FORMAT && blocks.size() > 0;
  }
};

inline UnitArea clipArea(const UnitArea &area, const UnitArea &boundingBox) {
  UnitArea ret(area.chromaFormat);

  for (uint32_t i = 0; i < area.blocks.size(); i++) {
    ret.blocks.push_back(clipArea(area.blocks[i], boundingBox.blocks[i]));
  }

  return ret;
}

struct UnitAreaRelative : public UnitArea {
  UnitAreaRelative(const UnitArea &origUnit, const UnitArea &unit) {
    *((UnitArea *)this) = unit;
    for (uint32_t i = 0; i < blocks.size(); i++) {
      blocks[i].positionRelativeTo(origUnit.blocks[i]);
    }
  }
};

class SPS;
class VPS;
class DCI;
class PPS;
class Slice;

// ---------------------------------------------------------------------------
// coding unit
// ---------------------------------------------------------------------------

struct TransformUnit;
struct PredictionUnit;
class CodingStructure;

struct CodingUnit : public UnitArea {
  CodingStructure *cs;
  Slice *slice;
  ChannelType chType;

  PredMode predMode;

  uint8_t depth;   // number of all splits, applied with generalized splits
  uint8_t qtDepth; // number of applied quad-splits, before switching to the
                   // multi-type-tree (mtt)
  // a triple split would increase the mtDepth by 1, but the qtDepth by 2 in the
  // first and last part and by 1 in the middle part (because of the 1-2-1 split
  // proportions)
  uint8_t btDepth; // number of applied binary splits, after switching to the
                   // mtt (or it's equivalent)
  uint8_t mtDepth; // the actual number of splits after switching to mtt (equals
                   // btDepth if only binary splits are allowed)
  int8_t chromaQpAdj;
  int8_t qp;
  SplitSeries splitSeries;
  TreeType treeType;
  ModeType modeType;
  ModeTypeSeries modeTypeSeries;
  bool skip;
  bool mmvdSkip;
  bool affine;
  int affineType;
  bool colorTransform;
  bool geoFlag;
  int bdpcmMode;
  int bdpcmModeChroma;
  uint8_t imv;
  bool rootCbf;
  uint8_t sbtInfo;
  uint32_t tileIdx;
  uint32_t lfnstIdx;
  uint8_t BcwIdx;
  bool mipFlag;

  // // needed for fast imv mode decisions
  uint8_t smvdMode;
  uint8_t ispMode;
  bool useEscape[MAX_NUM_CHANNEL_TYPE];
  bool useRotation[MAX_NUM_CHANNEL_TYPE];
  bool reuseflag[MAX_NUM_CHANNEL_TYPE][MAXPLTPREDSIZE];
  uint8_t lastPLTSize[MAX_NUM_CHANNEL_TYPE];
  uint8_t reusePLTSize[MAX_NUM_CHANNEL_TYPE];
  uint8_t curPLTSize[MAX_NUM_CHANNEL_TYPE];
  Pel curPLT[MAX_NUM_COMPONENT][MAXPLTSIZE];

  void initData();

  unsigned idx;
  CodingUnit *next;

  PredictionUnit *firstPU;
  PredictionUnit *lastPU;

  TransformUnit *firstTU;
  TransformUnit *lastTU;

  const uint8_t getSbtIdx() const {
    assert(((sbtInfo >> 0) & 0xf) < NUMBER_SBT_IDX);
    return (sbtInfo >> 0) & 0xf;
  }
  const uint8_t getSbtPos() const { return (sbtInfo >> 4) & 0x3; }
  void setSbtIdx(uint8_t idx) {
    CHECK(idx >= NUMBER_SBT_IDX, "sbt_idx wrong");
    sbtInfo = (idx << 0) + (sbtInfo & 0xf0);
  }
  void setSbtPos(uint8_t pos) {
    CHECK(pos >= 4, "sbt_pos wrong");
    sbtInfo = (pos << 4) + (sbtInfo & 0xcf);
  }
  uint8_t getSbtTuSplit() const;
  const uint8_t checkAllowedSbt() const;
  const bool checkCCLMAllowed() const;
  const bool isSepTree() const;
  const bool isLocalSepTree() const;
  const bool isConsInter() const { return modeType == MODE_TYPE_INTER; }
  const bool isConsIntra() const { return modeType == MODE_TYPE_INTRA; }
};

// ---------------------------------------------------------------------------
// prediction unit
// ---------------------------------------------------------------------------

struct IntraPredictionData {
  uint32_t intraDir[MAX_NUM_CHANNEL_TYPE];
  bool mipTransposedFlag;
  int multiRefIdx;
};

struct InterPredictionData {
  bool mergeFlag;
  bool regularMergeFlag;
  uint8_t mergeIdx;
  uint8_t geoSplitDir;
  uint8_t geoMergeIdx0;
  uint8_t geoMergeIdx1;
  bool mmvdMergeFlag;
  uint32_t mmvdMergeIdx;
  uint8_t interDir;
  uint8_t mvpIdx[NUM_REF_PIC_LIST_01];
  Mv mvd[NUM_REF_PIC_LIST_01];
  Mv mv[NUM_REF_PIC_LIST_01];
  int16_t refIdx[NUM_REF_PIC_LIST_01];
  MergeType mergeType;
  Mv mvdAffi[NUM_REF_PIC_LIST_01][3];
  bool ciipFlag;
};

struct PredictionUnit : public UnitArea,
                        public IntraPredictionData,
                        public InterPredictionData {
  CodingUnit *cu;
  CodingStructure *cs;
  ChannelType chType;

  void initData();
  unsigned idx;
  PredictionUnit *next;
};

// ---------------------------------------------------------------------------
// transform unit
// ---------------------------------------------------------------------------

struct TransformUnit : public UnitArea {
  CodingUnit *cu;
  CodingStructure *cs;
  ChannelType chType;

  uint8_t depth;
  uint8_t mtsIdx[MAX_NUM_TBLOCKS];
  bool noResidual;
  uint8_t jointCbCr;
  uint8_t cbf[MAX_NUM_TBLOCKS];

  void initData();

  unsigned idx;
  TransformUnit *next;
  TransformUnit *prev;
  void init(TCoeff **coeffs, Pel **pcmbuf, bool **runType);

  void checkTuNoResidual(unsigned idx);
  int getTbAreaAfterCoefZeroOut(ComponentID compID) const;

  CoeffBuf getCoeffs(const ComponentID id);
  const CCoeffBuf getCoeffs(const ComponentID id) const;
  PelBuf getPcmbuf(const ComponentID id);
  const CPelBuf getPcmbuf(const ComponentID id) const;
  PelBuf getcurPLTIdx(const ComponentID id);
  const CPelBuf getcurPLTIdx(const ComponentID id) const;
  PLTtypeBuf getrunType(const ComponentID id);
  const CPLTtypeBuf getrunType(const ComponentID id) const;
  PLTescapeBuf getescapeValue(const ComponentID id);
  const CPLTescapeBuf getescapeValue(const ComponentID id) const;

private:
  TCoeff *m_coeffs[MAX_NUM_TBLOCKS];
  Pel *m_pcmbuf[MAX_NUM_TBLOCKS];
  bool *m_runType[MAX_NUM_TBLOCKS - 1];
};

// ---------------------------------------------------------------------------
// Utility class for easy for-each like unit traversing
// ---------------------------------------------------------------------------

#include <iterator>

template <typename T>
class UnitIterator : public ::std::iterator<::std::forward_iterator_tag, T> {
private:
  T *m_punit;

public:
  UnitIterator() : m_punit(nullptr) {}
  UnitIterator(T *_punit) : m_punit(_punit) {}

  typedef T &reference;
  typedef T const &const_reference;
  typedef T *pointer;
  typedef T const *const_pointer;

  reference operator*() { return *m_punit; }
  const_reference operator*() const { return *m_punit; }
  pointer operator->() { return m_punit; }
  const_pointer operator->() const { return m_punit; }

  UnitIterator<T> &operator++() {
    m_punit = m_punit->next;
    return *this;
  }
  UnitIterator<T> operator++(int) {
    auto x = *this;
    ++(*this);
    return x;
  }
  bool operator!=(const UnitIterator<T> &other) const {
    return m_punit != other.m_punit;
  }
  bool operator==(const UnitIterator<T> &other) const {
    return m_punit == other.m_punit;
  }
};

template <typename T> class UnitTraverser {
private:
  T *m_begin;
  T *m_end;

public:
  UnitTraverser() : m_begin(nullptr), m_end(nullptr) {}
  UnitTraverser(T *_begin, T *_end) : m_begin(_begin), m_end(_end) {}

  typedef T value_type;
  typedef size_t size_type;
  typedef T &reference;
  typedef T const &const_reference;
  typedef T *pointer;
  typedef T const *const_pointer;
  typedef UnitIterator<T> iterator;
  typedef UnitIterator<const T> const_iterator;

  iterator begin() { return UnitIterator<T>(m_begin); }
  const_iterator begin() const { return UnitIterator<T>(m_begin); }
  const_iterator cbegin() const { return UnitIterator<T>(m_begin); }
  iterator end() { return UnitIterator<T>(m_end); }
  const_iterator end() const { return UnitIterator<T>(m_end); }
  const_iterator cend() const { return UnitIterator<T>(m_end); }
};

typedef UnitTraverser<CodingUnit> CUTraverser;
typedef UnitTraverser<PredictionUnit> PUTraverser;
typedef UnitTraverser<TransformUnit> TUTraverser;

typedef UnitTraverser<const CodingUnit> cCUTraverser;
typedef UnitTraverser<const PredictionUnit> cPUTraverser;
typedef UnitTraverser<const TransformUnit> cTUTraverser;
} // namespace EntropyCoding

#endif
