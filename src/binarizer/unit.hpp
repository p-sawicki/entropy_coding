#ifndef ENTROPY_CODEC_UNIT
#define ENTROPY_CODEC_UNIT

#include <assert.h>

#include "buffer.hpp"
#include "common_def.hpp"
#include "mv.hpp"

namespace EntropyCoding {
typedef ::std::array<::std::array<Pel, MAXPLTSIZE>, MAX_NUM_COMPONENT> CurPLT31;
typedef ::std::array<::std::array<Pel, MAXPLTPREDSIZE>, MAX_NUM_COMPONENT>
    CurPLT63;

// ---------------------------------------------------------------------------
// tools
// ---------------------------------------------------------------------------
struct PLTBuf {
  ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> curPLTSize;
  CurPLT63 curPLT;
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
  // void resizeTo(const Size &newSize) { Size::resizeTo(newSize); }
  void resizeTo(const Size &newSize) {}
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
  UnitArea(const ChromaFormat _chromaFormat, const UnitBlocksType &_blocks);
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

typedef ::std::array<::std::array<bool, MAXPLTPREDSIZE>, MAX_NUM_CHANNEL_TYPE>
    ReuseFlag;

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
  ::std::array<bool, MAX_NUM_CHANNEL_TYPE> useEscape;
  ::std::array<bool, MAX_NUM_CHANNEL_TYPE> useRotation;
  ReuseFlag reuseflag;
  ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> lastPLTSize;
  ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> reusePLTSize;
  ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> curPLTSize;
  CurPLT31 curPLT;

  void initData();

  unsigned idx;
  CodingUnit *next;

  PredictionUnit *firstPU;
  PredictionUnit *lastPU;

  TransformUnit *firstTU;
  TransformUnit *lastTU;

  CodingUnit() {}

  CodingUnit(const UnitArea &unitArea, const ChannelType _chType,
             const PredMode _predMode, const uint8_t _depth,
             const uint8_t _qtDepth, const uint8_t _btDepth,
             const uint8_t _mtDepth, const uint8_t _chromaQpAdj,
             const uint8_t _qp, const SplitSeries _splitSeries,
             const TreeType _treeType, const ModeType _modeType,
             const ModeTypeSeries _modeTypeSeries, const bool _skip,
             const bool _mmvdSkip, const bool _affine, const int _affineType,
             const bool _colorTransform, const bool _geoFlag,
             const int _bdpcmMode, const int _bdpcmModeChroma,
             const uint8_t _imv, const bool _rootCbf, const uint8_t _sbtInfo,
             const uint32_t _tileIdx, const uint32_t _lfnstIdx,
             const uint8_t _BcwIdx, const bool _mipFlag,
             const uint8_t _smvdMode, const uint8_t _ispMode,
             const bool *_useEscape, const bool *_useRotation,
             const bool _reuseflag[][MAXPLTPREDSIZE],
             const uint8_t *_lastPLTSize, const uint8_t *_reusePLTSize,
             const uint8_t *_curPLTSize, const Pel _curPLT[][MAXPLTSIZE],
             const unsigned _idx)
      : UnitArea(unitArea), chType(_chType), predMode(_predMode), depth(_depth),
        qtDepth(_qtDepth), btDepth(_btDepth), mtDepth(_mtDepth),
        chromaQpAdj(_chromaQpAdj), qp(_qp), splitSeries(_splitSeries),
        treeType(_treeType), modeType(_modeType),
        modeTypeSeries(_modeTypeSeries), skip(_skip), mmvdSkip(_mmvdSkip),
        affine(_affine), affineType(_affineType),
        colorTransform(_colorTransform), geoFlag(_geoFlag),
        bdpcmMode(_bdpcmMode), bdpcmModeChroma(_bdpcmModeChroma), imv(_imv),
        rootCbf(_rootCbf), sbtInfo(_sbtInfo), tileIdx(_tileIdx),
        lfnstIdx(_lfnstIdx), BcwIdx(_BcwIdx), mipFlag(_mipFlag),
        smvdMode(_smvdMode), ispMode(_ispMode), idx(_idx) {
    copy_array(_useEscape, useEscape);
    copy_array(_useRotation, useRotation);
    for (int i = 0; i < reuseflag.size(); ++i) {
      copy_array(_reuseflag[i], reuseflag[i]);
    }
    copy_array(_lastPLTSize, lastPLTSize);
    copy_array(_reusePLTSize, reusePLTSize);
    copy_array(_curPLTSize, curPLTSize);
    for (int i = 0; i < curPLT.size(); ++i) {
      copy_array(_curPLT[i], curPLT[i]);
    }
  }

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
  ::std::array<uint32_t, MAX_NUM_CHANNEL_TYPE> intraDir;
  bool mipTransposedFlag;
  int multiRefIdx;

  IntraPredictionData() {}
  IntraPredictionData(const uint32_t *_intraDir, const bool _mipTransposedFlag,
                      const int _multiRefIdx)
      : mipTransposedFlag(_mipTransposedFlag), multiRefIdx(_multiRefIdx) {
    copy_array(_intraDir, intraDir);
  }
};

typedef ::std::array<::std::array<Mv, 3>, NUM_REF_PIC_LIST_01> MvdAffi;

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
  ::std::array<uint8_t, NUM_REF_PIC_LIST_01> mvpIdx;
  ::std::array<Mv, NUM_REF_PIC_LIST_01> mvd;
  ::std::array<Mv, NUM_REF_PIC_LIST_01> mv;
  ::std::array<int16_t, NUM_REF_PIC_LIST_01> refIdx;
  MergeType mergeType;
  MvdAffi mvdAffi;
  bool ciipFlag;

  InterPredictionData() {}
  InterPredictionData(const bool _mergeFlag, const bool _regularMergeFlag,
                      const uint8_t _mergeIdx, const uint8_t _geoSplitDir,
                      const uint8_t _geoMergeIdx0, const uint8_t _geoMergeIdx1,
                      const bool _mmvdMergeFlag, const uint32_t _mmvdMergeIdx,
                      const uint8_t _interDir, const uint8_t *_mvpIdx,
                      const ::std::array<Mv, NUM_REF_PIC_LIST_01> &_mvd,
                      const ::std::array<Mv, NUM_REF_PIC_LIST_01> &_mv,
                      const int16_t *_refIdx, const MergeType _mergeType,
                      const MvdAffi &_mvdAffi, const bool _ciipFlag)
      : mergeFlag(_mergeFlag), regularMergeFlag(_regularMergeFlag),
        mergeIdx(_mergeIdx), geoSplitDir(_geoSplitDir),
        geoMergeIdx0(_geoMergeIdx0), geoMergeIdx1(_geoMergeIdx1),
        mmvdMergeFlag(_mmvdMergeFlag), mmvdMergeIdx(_mmvdMergeIdx),
        interDir(_interDir), mvd(_mvd), mv(_mv), mergeType(_mergeType),
        mvdAffi(_mvdAffi), ciipFlag(_ciipFlag) {
    copy_array(_mvpIdx, mvpIdx);
    copy_array(_refIdx, refIdx);
  }
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

  PredictionUnit() {}
  PredictionUnit(const UnitArea &_unitArea, const IntraPredictionData &_intraPd,
                 const InterPredictionData &_interPd, const ChannelType _chType,
                 const unsigned _idx)
      : UnitArea(_unitArea), IntraPredictionData(_intraPd),
        InterPredictionData(_interPd), chType(_chType), idx(_idx) {}
};

// ---------------------------------------------------------------------------
// transform unit
// ---------------------------------------------------------------------------

struct TransformUnit : public UnitArea {
  CodingUnit *cu;
  CodingStructure *cs;
  ChannelType chType;

  uint8_t depth;
  ::std::array<uint8_t, MAX_NUM_TBLOCKS> mtsIdx;
  bool noResidual;
  uint8_t jointCbCr;
  ::std::array<uint8_t, MAX_NUM_TBLOCKS> cbf;

  TransformUnit() {}
  TransformUnit(const UnitArea &_unitArea, const ChannelType _chType,
                const uint8_t _depth, const uint8_t *_mtsIdx,
                const bool _noResidual, const uint8_t _jointCbCr,
                const uint8_t *_cbf, const unsigned _idx,
                TCoeff *const coeffs[MAX_NUM_TBLOCKS],
                Pel *const pcmbuf[MAX_NUM_TBLOCKS],
                bool *const runType[MAX_NUM_TBLOCKS - 1])
      : UnitArea(_unitArea), chType(_chType), depth(_depth),
        noResidual(_noResidual), jointCbCr(_jointCbCr), idx(_idx) {
    copy_array(_mtsIdx, mtsIdx);
    copy_array(_cbf, cbf);
    copy_array(coeffs, m_coeffs);
    copy_array(pcmbuf, m_pcmbuf);
    copy_array(runType, m_runType);
  }

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

  const ::std::array<TCoeff *, MAX_NUM_TBLOCKS> getCoeffs() const {
    return m_coeffs;
  }
  const ::std::array<Pel *, MAX_NUM_TBLOCKS> getPcmBuf() const {
    return m_pcmbuf;
  }
  const ::std::array<bool *, MAX_NUM_TBLOCKS - 1> getRunType() const {
    return m_runType;
  }

private:
  ::std::array<TCoeff *, MAX_NUM_TBLOCKS> m_coeffs;
  ::std::array<Pel *, MAX_NUM_TBLOCKS> m_pcmbuf;
  ::std::array<bool *, MAX_NUM_TBLOCKS - 1> m_runType;
};
} // namespace EntropyCoding

// ---------------------------------------------------------------------------
// Utility class for easy for-each like unit traversing
// ---------------------------------------------------------------------------

#include <iterator>

namespace EntropyCoding {
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
