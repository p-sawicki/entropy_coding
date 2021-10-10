#include <cstring>

#include "slice.hpp"
#include "unit.hpp"
#include "unit_partitioner.hpp"
#include "unit_tools.hpp"
#include "coding_structure.hpp"

namespace EntropyCoding {

Position CompArea::chromaPos() const {
  if (isLuma(compID)) {
    uint32_t scaleX = getComponentScaleX(compID, chromaFormat);
    uint32_t scaleY = getComponentScaleY(compID, chromaFormat);

    return Position(x >> scaleX, y >> scaleY);
  } else {
    return *this;
  }
}

Size CompArea::lumaSize() const {
  if (isChroma(compID)) {
    uint32_t scaleX = getComponentScaleX(compID, chromaFormat);
    uint32_t scaleY = getComponentScaleY(compID, chromaFormat);

    return Size(width << scaleX, height << scaleY);
  } else {
    return *this;
  }
}

Size CompArea::chromaSize() const {
  if (isLuma(compID)) {
    uint32_t scaleX = getComponentScaleX(compID, chromaFormat);
    uint32_t scaleY = getComponentScaleY(compID, chromaFormat);

    return Size(width >> scaleX, height >> scaleY);
  } else {
    return *this;
  }
}

Position CompArea::lumaPos() const {
  if (isChroma(compID)) {
    uint32_t scaleX = getComponentScaleX(compID, chromaFormat);
    uint32_t scaleY = getComponentScaleY(compID, chromaFormat);

    return Position(x << scaleX, y << scaleY);
  } else {
    return *this;
  }
}

Position CompArea::compPos(const ComponentID compID) const {
  return isLuma(compID) ? lumaPos() : chromaPos();
}

Position CompArea::chanPos(const ChannelType chType) const {
  return isLuma(chType) ? lumaPos() : chromaPos();
}

void CompArea::xRecalcLumaToChroma() {
  const uint32_t csx = getComponentScaleX(compID, chromaFormat);
  const uint32_t csy = getComponentScaleY(compID, chromaFormat);

  x >>= csx;
  y >>= csy;
  width >>= csx;
  height >>= csy;
}

UnitArea::UnitArea(const ChromaFormat _chromaFormat)
    : chromaFormat(_chromaFormat) {}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, const UnitBlocksType &_blocks)
    : chromaFormat(_chromaFormat), blocks(_blocks) {}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, const Area &_area)
    : chromaFormat(_chromaFormat),
      blocks(getNumberValidComponents(_chromaFormat)) {
  const uint32_t numCh = getNumberValidComponents(chromaFormat);

  for (uint32_t i = 0; i < numCh; i++) {
    blocks[i] = CompArea(ComponentID(i), chromaFormat, _area, true);
  }
}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, const CompArea &blkY)
    : chromaFormat(_chromaFormat), blocks{blkY} {}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, CompArea &&blkY)
    : chromaFormat(_chromaFormat), blocks{::std::forward<CompArea>(blkY)} {}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, const CompArea &blkY,
                   const CompArea &blkCb, const CompArea &blkCr)
    : chromaFormat(_chromaFormat), blocks{blkY, blkCb, blkCr} {}

UnitArea::UnitArea(const ChromaFormat _chromaFormat, CompArea &&blkY,
                   CompArea &&blkCb, CompArea &&blkCr)
    : chromaFormat(_chromaFormat), blocks{::std::forward<CompArea>(blkY),
                                          ::std::forward<CompArea>(blkCb),
                                          ::std::forward<CompArea>(blkCr)} {}

bool UnitArea::contains(const UnitArea &other) const {
  bool ret = true;
  bool any = false;

  for (const auto &blk : other.blocks) {
    if (blk.valid() && blocks[blk.compID].valid()) {
      ret &= blocks[blk.compID].contains(blk);
      any = true;
    }
  }

  return any && ret;
}

bool UnitArea::contains(const UnitArea &other, const ChannelType chType) const {
  bool ret = true;
  bool any = false;

  for (const auto &blk : other.blocks) {
    if (toChannelType(blk.compID) == chType && blk.valid() &&
        blocks[blk.compID].valid()) {
      ret &= blocks[blk.compID].contains(blk);
      any = true;
    }
  }

  return any && ret;
}

void UnitArea::repositionTo(const UnitArea &unitArea) {
  for (uint32_t i = 0; i < blocks.size(); i++) {
    blocks[i].repositionTo(unitArea.blocks[i]);
  }
}

const UnitArea UnitArea::singleComp(const ComponentID compID) const {
  UnitArea ret(chromaFormat);

  for (const auto &blk : blocks) {
    if (blk.compID == compID) {
      ret.blocks.push_back(blk);
    } else {
      ret.blocks.push_back(CompArea());
    }
  }

  return ret;
}

const UnitArea UnitArea::singleChan(const ChannelType chType) const {
  UnitArea ret(chromaFormat);

  for (const auto &blk : blocks) {
    if (toChannelType(blk.compID) == chType) {
      ret.blocks.push_back(blk);
    } else {
      ret.blocks.push_back(CompArea());
    }
  }

  return ret;
}

void CodingUnit::initData() {
  predMode = NUMBER_OF_PREDICTION_MODES;
  qtDepth = 0;
  depth = 0;
  btDepth = 0;
  mtDepth = 0;
  splitSeries = 0;
  skip = false;
  mmvdSkip = false;
  affine = false;
  affineType = 0;
  colorTransform = false;
  geoFlag = false;
  bdpcmMode = 0;
  bdpcmModeChroma = 0;
  qp = 0;
  chromaQpAdj = 0;
  rootCbf = true;
  sbtInfo = 0;
  // mtsFlag           = 0;
  lfnstIdx = 0;
  tileIdx = 0;
  imv = 0;
  // imvNumCand        = 0;
  BcwIdx = BCW_DEFAULT;
  for (int i = 0; i < 2; i++)
    // refIdxBi[i] = -1;
    smvdMode = 0;
  ispMode = 0;
  mipFlag = false;

  for (int idx = 0; idx < MAX_NUM_CHANNEL_TYPE; idx++) {
    curPLTSize[idx] = 0;
    reusePLTSize[idx] = 0;
    lastPLTSize[idx] = 0;
    useEscape[idx] = false;
    useRotation[idx] = false;
    memset(reuseflag[idx].data(), false, MAXPLTPREDSIZE * sizeof(bool));
  }

  for (int idx = 0; idx < MAX_NUM_COMPONENT; idx++) {
    memset(curPLT[idx].data(), 0, MAXPLTSIZE * sizeof(Pel));
  }

  treeType = TREE_D;
  modeType = MODE_TYPE_ALL;
  modeTypeSeries = 0;
}

uint8_t CodingUnit::getSbtTuSplit() const {
  uint8_t sbtTuSplitType = 0;

  switch (getSbtIdx()) {
  case SBT_VER_HALF:
    sbtTuSplitType =
        (getSbtPos() == SBT_POS0 ? 0 : 1) + SBT_VER_HALF_POS0_SPLIT;
    break;
  case SBT_HOR_HALF:
    sbtTuSplitType =
        (getSbtPos() == SBT_POS0 ? 0 : 1) + SBT_HOR_HALF_POS0_SPLIT;
    break;
  case SBT_VER_QUAD:
    sbtTuSplitType =
        (getSbtPos() == SBT_POS0 ? 0 : 1) + SBT_VER_QUAD_POS0_SPLIT;
    break;
  case SBT_HOR_QUAD:
    sbtTuSplitType =
        (getSbtPos() == SBT_POS0 ? 0 : 1) + SBT_HOR_QUAD_POS0_SPLIT;
    break;
  default:
    assert(0);
    break;
  }

  assert(sbtTuSplitType <= SBT_HOR_QUAD_POS1_SPLIT &&
         sbtTuSplitType >= SBT_VER_HALF_POS0_SPLIT);
  return sbtTuSplitType;
}

const uint8_t CodingUnit::checkAllowedSbt() const {
  if (!slice->getSPS()->getUseSBT()) {
    return 0;
  }

  // check on prediction mode
  if (predMode == MODE_INTRA || predMode == MODE_IBC ||
      predMode == MODE_PLT) // intra, palette or IBC
  {
    return 0;
  }
  if (firstPU->ciipFlag) {
    return 0;
  }

  uint8_t sbtAllowed = 0;
  int cuWidth = lwidth();
  int cuHeight = lheight();
  bool allow_type[NUMBER_SBT_IDX];
  memset(allow_type, false, NUMBER_SBT_IDX * sizeof(bool));

  // parameter
  int maxSbtCUSize = cs->sps->getMaxTbSize();
  int minSbtCUSize = 1 << (MIN_CU_LOG2 + 1);

  // check on size
  if (cuWidth > maxSbtCUSize || cuHeight > maxSbtCUSize) {
    return 0;
  }

  allow_type[SBT_VER_HALF] = cuWidth >= minSbtCUSize;
  allow_type[SBT_HOR_HALF] = cuHeight >= minSbtCUSize;
  allow_type[SBT_VER_QUAD] = cuWidth >= (minSbtCUSize << 1);
  allow_type[SBT_HOR_QUAD] = cuHeight >= (minSbtCUSize << 1);

  for (int i = 0; i < NUMBER_SBT_IDX; i++) {
    sbtAllowed += (uint8_t)allow_type[i] << i;
  }

  return sbtAllowed;
}

const bool CodingUnit::checkCCLMAllowed() const {
  bool allowCCLM = false;

  if (!CS::isDualITree(
          *cs)) // single tree I slice or non-I slice (Note: judging chType is
                // no longer equivalent to checking dual-tree I slice since the
                // local dual-tree is introduced)
  {
    allowCCLM = true;
  } else if (slice->getSPS()->getCTUSize() <= 32) // dual tree, CTUsize < 64
  {
    allowCCLM = true;
  } else // dual tree, CTU size 64 or 128
  {
    int depthFor64x64Node = slice->getSPS()->getCTUSize() == 128 ? 1 : 0;
    const PartSplit cuSplitTypeDepth1 =
        CU::getSplitAtDepth(*this, depthFor64x64Node);
    const PartSplit cuSplitTypeDepth2 =
        CU::getSplitAtDepth(*this, depthFor64x64Node + 1);

    // allow CCLM if 64x64 chroma tree node uses QT split or HBT+VBT split
    // combination
    if (cuSplitTypeDepth1 == CU_QUAD_SPLIT ||
        (cuSplitTypeDepth1 == CU_HORZ_SPLIT &&
         cuSplitTypeDepth2 == CU_VERT_SPLIT)) {
      if (chromaFormat == CHROMA_420) {
        CHECK(!(blocks[COMPONENT_Cb].width <= 16 &&
                blocks[COMPONENT_Cb].height <= 16),
              "chroma cu size shall be <= 16x16 for YUV420 format");
      }
      allowCCLM = true;
    }
    // allow CCLM if 64x64 chroma tree node uses NS (No Split) and becomes a
    // chroma CU containing 32x32 chroma blocks
    else if (cuSplitTypeDepth1 == CU_DONT_SPLIT) {
      if (chromaFormat == CHROMA_420) {
        CHECK(!(blocks[COMPONENT_Cb].width == 32 &&
                blocks[COMPONENT_Cb].height == 32),
              "chroma cu size shall be 32x32 for YUV420 format");
      }
      allowCCLM = true;
    }
    // allow CCLM if 64x32 chroma tree node uses NS and becomes a chroma CU
    // containing 32x16 chroma blocks
    else if (cuSplitTypeDepth1 == CU_HORZ_SPLIT &&
             cuSplitTypeDepth2 == CU_DONT_SPLIT) {
      if (chromaFormat == CHROMA_420) {
        CHECK(!(blocks[COMPONENT_Cb].width == 32 &&
                blocks[COMPONENT_Cb].height == 16),
              "chroma cu size shall be 32x16 for YUV420 format");
      }
      allowCCLM = true;
    }

    // further check luma conditions
    if (allowCCLM) {
      // disallow CCLM if luma 64x64 block uses BT or TT or NS with ISP
      const Position lumaRefPos(
          chromaPos().x << getComponentScaleX(COMPONENT_Cb, chromaFormat),
          chromaPos().y << getComponentScaleY(COMPONENT_Cb, chromaFormat));
      const CodingUnit *colLumaCu =
          cs->picture->cs->getCU(lumaRefPos, CHANNEL_TYPE_LUMA);

      if (colLumaCu->lwidth() < 64 ||
          colLumaCu->lheight() < 64) // further split at 64x64 luma node
      {
        const PartSplit cuSplitTypeDepth1Luma =
            CU::getSplitAtDepth(*colLumaCu, depthFor64x64Node);
        CHECK(!(cuSplitTypeDepth1Luma >= CU_QUAD_SPLIT &&
                cuSplitTypeDepth1Luma <= CU_TRIV_SPLIT),
              "split mode shall be BT, TT or QT");
        if (cuSplitTypeDepth1Luma != CU_QUAD_SPLIT) {
          allowCCLM = false;
        }
      } else if (colLumaCu->lwidth() == 64 && colLumaCu->lheight() == 64 &&
                 colLumaCu
                     ->ispMode) // not split at 64x64 luma node and use ISP mode
      {
        allowCCLM = false;
      }
    }
  }

  return allowCCLM;
}

const bool CodingUnit::isSepTree() const {
  return treeType != TREE_D || CS::isDualITree(*cs);
}

const bool CodingUnit::isLocalSepTree() const {
  return treeType != TREE_D && !CS::isDualITree(*cs);
}

void PredictionUnit::initData() {
  // intra data - need this default initialization for PCM
  intraDir[0] = DC_IDX;
  intraDir[1] = PLANAR_IDX;
  mipTransposedFlag = false;
  multiRefIdx = 0;

  // inter data
  mergeFlag = false;
  regularMergeFlag = false;
  mergeIdx = MAX_UCHAR;
  geoSplitDir = MAX_UCHAR;
  geoMergeIdx0 = MAX_UCHAR;
  geoMergeIdx1 = MAX_UCHAR;
  mmvdMergeFlag = false;
  mmvdMergeIdx = MAX_UINT;
  interDir = MAX_UCHAR;
  mergeType = MRG_TYPE_DEFAULT_N;
  // bv.setZero();
  // bvd.setZero();
  // mvRefine = false;
  for (uint32_t i = 0; i < MAX_NUM_SUBCU_DMVR; i++) {
    // mvdL0SubPu[i].setZero();
  }
  for (uint32_t i = 0; i < NUM_REF_PIC_LIST_01; i++) {
    mvpIdx[i] = MAX_UCHAR;
    // mvpNum[i] = MAX_UCHAR;
    refIdx[i] = -1;
    mv[i].setZero();
    mvd[i].setZero();
    for (uint32_t j = 0; j < 3; j++) {
      mvdAffi[i][j].setZero();
    }
    for (uint32_t j = 0; j < 3; j++) {
      // mvAffi[i][j].setZero();
// #if GDR_ENABLED
//       mvAffiSolid[i][j] = true;
//       mvAffiValid[i][j] = true;
// #endif
    }
  }
  ciipFlag = false;
  // mmvdEncOptMode = 0;
}

void TransformUnit::initData() {
  for (unsigned i = 0; i < MAX_NUM_TBLOCKS; i++) {
    cbf[i] = 0;
    mtsIdx[i] = MTS_DCT2_DCT2;
  }
  depth = 0;
  noResidual = false;
  jointCbCr = 0;
  // m_chromaResScaleInv = 0;
}

void TransformUnit::init(TCoeff **coeffs, Pel **pcmbuf, bool **runType) {
  uint32_t numBlocks = getNumberValidTBlocks(*cs->pcv);

  for (uint32_t i = 0; i < numBlocks; i++) {
    m_coeffs[i] = coeffs[i];
    m_pcmbuf[i] = pcmbuf[i];
  }

  // numBlocks is either 1 for 4:0:0, or 3 otherwise. It would perhaps be better
  // to loop over getNumberValidChannels(*cs->pcv.chrFormat) for m_runType.
  for (uint32_t i = 0; i < ::std::max<uint32_t>(2, numBlocks) - 1; i++) {
    m_runType[i] = runType[i];
  }
}

void TransformUnit::checkTuNoResidual(unsigned idx) {
  if (CU::getSbtIdx(cu->sbtInfo) == SBT_OFF_DCT) {
    return;
  }

  if ((CU::getSbtPos(cu->sbtInfo) == SBT_POS0 && idx == 1) ||
      (CU::getSbtPos(cu->sbtInfo) == SBT_POS1 && idx == 0)) {
    noResidual = true;
  }
}

int TransformUnit::getTbAreaAfterCoefZeroOut(ComponentID compID) const {
  int tbArea = blocks[compID].width * blocks[compID].height;
  int tbZeroOutWidth = blocks[compID].width;
  int tbZeroOutHeight = blocks[compID].height;

  if (cs->sps->getUseMTS() && cu->sbtInfo != 0 && blocks[compID].width <= 32 &&
      blocks[compID].height <= 32 && compID == COMPONENT_Y) {
    tbZeroOutWidth = (blocks[compID].width == 32) ? 16 : tbZeroOutWidth;
    tbZeroOutHeight = (blocks[compID].height == 32) ? 16 : tbZeroOutHeight;
  }
  tbZeroOutWidth = ::std::min<int>(JVET_C0024_ZERO_OUT_TH, tbZeroOutWidth);
  tbZeroOutHeight = ::std::min<int>(JVET_C0024_ZERO_OUT_TH, tbZeroOutHeight);
  tbArea = tbZeroOutWidth * tbZeroOutHeight;
  return tbArea;
}

CoeffBuf TransformUnit::getCoeffs(const ComponentID id) {
  return CoeffBuf(m_coeffs[id], blocks[id]);
}
const CCoeffBuf TransformUnit::getCoeffs(const ComponentID id) const {
  return CCoeffBuf(m_coeffs[id], blocks[id]);
}

PelBuf TransformUnit::getPcmbuf(const ComponentID id) {
  return PelBuf(m_pcmbuf[id], blocks[id]);
}
const CPelBuf TransformUnit::getPcmbuf(const ComponentID id) const {
  return CPelBuf(m_pcmbuf[id], blocks[id]);
}

PelBuf TransformUnit::getcurPLTIdx(const ComponentID id) {
  return PelBuf(m_pcmbuf[id], blocks[id]);
}
const CPelBuf TransformUnit::getcurPLTIdx(const ComponentID id) const {
  return CPelBuf(m_pcmbuf[id], blocks[id]);
}

PLTtypeBuf TransformUnit::getrunType(const ComponentID id) {
  return PLTtypeBuf(m_runType[id], blocks[id]);
}
const CPLTtypeBuf TransformUnit::getrunType(const ComponentID id) const {
  return CPLTtypeBuf(m_runType[id], blocks[id]);
}

PLTescapeBuf TransformUnit::getescapeValue(const ComponentID id) {
  return PLTescapeBuf(m_coeffs[id], blocks[id]);
}
const CPLTescapeBuf TransformUnit::getescapeValue(const ComponentID id) const {
  return CPLTescapeBuf(m_coeffs[id], blocks[id]);
}
} // namespace EntropyCoding
