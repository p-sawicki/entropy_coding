#include "unit_tools.hpp"
#include "coding_structure.hpp"

namespace EntropyCoding {

uint32_t getCtuAddr(const Position &pos, const PreCalcValues &pcv) {
  return (pos.x >> pcv.maxCUWidthLog2) +
         (pos.y >> pcv.maxCUHeightLog2) * pcv.widthInCtus;
}

UnitArea CS::getArea(const CodingStructure &cs, const UnitArea &area,
                     const ChannelType chType) {
  return isDualITree(cs) || cs.treeType != TREE_D ? area.singleChan(chType)
                                                  : area;
}

bool CS::isDualITree(const CodingStructure &cs) {
  return cs.slice->isIntra() && !cs.pcv->ISingleTree;
}

bool CU::isIntra(const CodingUnit &cu) { return cu.predMode == MODE_INTRA; }

bool CU::isInter(const CodingUnit &cu) { return cu.predMode == MODE_INTER; }

bool CU::isIBC(const CodingUnit &cu) { return cu.predMode == MODE_IBC; }

bool CU::isPLT(const CodingUnit &cu) { return cu.predMode == MODE_PLT; }

bool CU::isSameCtu(const CodingUnit &cu, const CodingUnit &cu2) {
  uint32_t ctuSizeBit = floorLog2(cu.cs->sps->getMaxCUWidth());

  Position pos1Ctu(cu.lumaPos().x >> ctuSizeBit, cu.lumaPos().y >> ctuSizeBit);
  Position pos2Ctu(cu2.lumaPos().x >> ctuSizeBit,
                   cu2.lumaPos().y >> ctuSizeBit);

  return pos1Ctu.x == pos2Ctu.x && pos1Ctu.y == pos2Ctu.y;
}

bool CU::isSameSliceAndTile(const CodingUnit &cu, const CodingUnit &cu2) {
  return (cu.slice->getIndependentSliceIdx() ==
          cu2.slice->getIndependentSliceIdx()) &&
         (cu.tileIdx == cu2.tileIdx);
}

bool CU::isLastSubCUOfCtu(const CodingUnit &cu) {
  const Area cuAreaY =
      cu.isSepTree()
          ? Area(recalcPosition(cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA,
                                cu.blocks[cu.chType].pos()),
                 recalcSize(cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA,
                            cu.blocks[cu.chType].size()))
          : (const Area &)cu.Y();

  return (
      (((cuAreaY.x + cuAreaY.width) & cu.cs->pcv->maxCUWidthMask) == 0 ||
       cuAreaY.x + cuAreaY.width == cu.cs->pps->getPicWidthInLumaSamples()) &&
      (((cuAreaY.y + cuAreaY.height) & cu.cs->pcv->maxCUHeightMask) == 0 ||
       cuAreaY.y + cuAreaY.height == cu.cs->pps->getPicHeightInLumaSamples()));
}

uint32_t CU::getCtuAddr(const CodingUnit &cu) {
  return getCtuAddr(cu.blocks[cu.chType].lumaPos(), *cu.cs->pcv);
}

int CU::predictQP(const CodingUnit &cu, const int prevQP) {
  const CodingStructure &cs = *cu.cs;

  uint32_t ctuRsAddr = getCtuAddr(cu);
  uint32_t ctuXPosInCtus = ctuRsAddr % cs.pcv->widthInCtus;
  uint32_t tileColIdx = cu.slice->getPPS()->ctuToTileCol(ctuXPosInCtus);
  uint32_t tileXPosInCtus = cu.slice->getPPS()->getTileColumnBd(tileColIdx);
  if (ctuXPosInCtus == tileXPosInCtus &&
      !(cu.blocks[cu.chType].x &
        (cs.pcv->maxCUWidthMask >>
         getChannelTypeScaleX(cu.chType, cu.chromaFormat))) &&
      !(cu.blocks[cu.chType].y &
        (cs.pcv->maxCUHeightMask >>
         getChannelTypeScaleY(cu.chType, cu.chromaFormat))) &&
      (cs.getCU(cu.blocks[cu.chType].pos().offset(0, -1), cu.chType) != NULL) &&
      CU::isSameSliceAndTile(
          *cs.getCU(cu.blocks[cu.chType].pos().offset(0, -1), cu.chType), cu)) {
    return (
        (cs.getCU(cu.blocks[cu.chType].pos().offset(0, -1), cu.chType))->qp);
  } else {
    const int a =
        (cu.blocks[cu.chType].y &
         (cs.pcv->maxCUHeightMask >>
          getChannelTypeScaleY(cu.chType, cu.chromaFormat)))
            ? (cs.getCU(cu.blocks[cu.chType].pos().offset(0, -1), cu.chType))
                  ->qp
            : prevQP;
    const int b =
        (cu.blocks[cu.chType].x &
         (cs.pcv->maxCUWidthMask >>
          getChannelTypeScaleX(cu.chType, cu.chromaFormat)))
            ? (cs.getCU(cu.blocks[cu.chType].pos().offset(-1, 0), cu.chType))
                  ->qp
            : prevQP;

    return (a + b + 1) >> 1;
  }
}

uint32_t CU::getNumPUs(const CodingUnit &cu) {
  uint32_t cnt = 0;
  PredictionUnit *pu = cu.firstPU;

  do {
    cnt++;
  } while ((pu != cu.lastPU) && (pu = pu->next));

  return cnt;
}

PartSplit CU::getSplitAtDepth(const CodingUnit &cu, const unsigned depth) {
  if (depth >= cu.depth) {
    return CU_DONT_SPLIT;
  }

  const PartSplit cuSplitType =
      PartSplit((cu.splitSeries >> (depth * SPLIT_DMULT)) & SPLIT_MASK);

  if (cuSplitType == CU_QUAD_SPLIT) {
    return CU_QUAD_SPLIT;
  }

  else if (cuSplitType == CU_HORZ_SPLIT) {
    return CU_HORZ_SPLIT;
  }

  else if (cuSplitType == CU_VERT_SPLIT) {
    return CU_VERT_SPLIT;
  }

  else if (cuSplitType == CU_TRIH_SPLIT) {
    return CU_TRIH_SPLIT;
  } else if (cuSplitType == CU_TRIV_SPLIT) {
    return CU_TRIV_SPLIT;
  } else {
    THROW("Unknown split mode");
    return CU_QUAD_SPLIT;
  }
}

ModeType CU::getModeTypeAtDepth(const CodingUnit &cu, const unsigned depth) {
  ModeType modeType = ModeType((cu.modeTypeSeries >> (depth * 3)) & 0x07);
  CHECK(depth > cu.depth, " depth is wrong");
  return modeType;
}

bool CU::isBcwIdxCoded(const CodingUnit &cu) {
  if (cu.cs->sps->getUseBcw() == false) {
    CHECK(cu.BcwIdx != BCW_DEFAULT, "Error: cu.BcwIdx != BCW_DEFAULT");
    return false;
  }

  if (cu.predMode == MODE_IBC) {
    return false;
  }

  if (cu.predMode == MODE_INTRA || cu.cs->slice->isInterP()) {
    return false;
  }

  if (cu.lwidth() * cu.lheight() < BCW_SIZE_CONSTRAINT) {
    return false;
  }

  if (!cu.firstPU->mergeFlag) {
    if (cu.firstPU->interDir == 3) {
      const int refIdx0 = cu.firstPU->refIdx[REF_PIC_LIST_0];
      const int refIdx1 = cu.firstPU->refIdx[REF_PIC_LIST_1];

      const WPScalingParam *wp0 =
          cu.cs->slice->getWpScaling(REF_PIC_LIST_0, refIdx0);
      const WPScalingParam *wp1 =
          cu.cs->slice->getWpScaling(REF_PIC_LIST_1, refIdx1);

      return !(WPScalingParam::isWeighted(wp0) ||
               WPScalingParam::isWeighted(wp1));
    }
  }

  return false;
}

uint8_t CU::getValidBcwIdx(const CodingUnit &cu) {
  if (cu.firstPU->interDir == 3 && !cu.firstPU->mergeFlag) {
    return cu.BcwIdx;
  } else if (cu.firstPU->interDir == 3 && cu.firstPU->mergeFlag &&
             cu.firstPU->mergeType == MRG_TYPE_DEFAULT_N) {
    // This is intended to do nothing here.
  } else if (cu.firstPU->mergeFlag &&
             cu.firstPU->mergeType == MRG_TYPE_SUBPU_ATMVP) {
    CHECK(cu.BcwIdx != BCW_DEFAULT, " cu.BcwIdx != BCW_DEFAULT ");
  } else {
    CHECK(cu.BcwIdx != BCW_DEFAULT, " cu.BcwIdx != BCW_DEFAULT ");
  }

  return BCW_DEFAULT;
}

void CU::setBcwIdx(CodingUnit &cu, uint8_t uh) {
  int8_t uhCnt = 0;

  if (cu.firstPU->interDir == 3 && !cu.firstPU->mergeFlag) {
    cu.BcwIdx = uh;
    ++uhCnt;
  } else if (cu.firstPU->interDir == 3 && cu.firstPU->mergeFlag &&
             cu.firstPU->mergeType == MRG_TYPE_DEFAULT_N) {
    // This is intended to do nothing here.
  } else if (cu.firstPU->mergeFlag &&
             cu.firstPU->mergeType == MRG_TYPE_SUBPU_ATMVP) {
    cu.BcwIdx = BCW_DEFAULT;
  } else {
    cu.BcwIdx = BCW_DEFAULT;
  }

  CHECK(uhCnt <= 0, " uhCnt <= 0 ");
}

uint8_t CU::getSbtIdx(const uint8_t sbtInfo) { return (sbtInfo >> 0) & 0xf; }

uint8_t CU::getSbtPos(const uint8_t sbtInfo) { return (sbtInfo >> 4) & 0x3; }

bool CU::bdpcmAllowed(const CodingUnit &cu, const ComponentID compID) {
  SizeType transformSkipMaxSize =
      1 << cu.cs->sps->getLog2MaxTransformSkipBlockSize();

  bool bdpcmAllowed = cu.cs->sps->getBDPCMEnabledFlag();
  bdpcmAllowed &= CU::isIntra(cu);
  if (isLuma(compID)) {
    bdpcmAllowed &= (cu.lwidth() <= transformSkipMaxSize &&
                     cu.lheight() <= transformSkipMaxSize);
  } else {
    bdpcmAllowed &= (cu.chromaSize().width <= transformSkipMaxSize &&
                     cu.chromaSize().height <= transformSkipMaxSize) &&
                    !cu.colorTransform;
  }
  return bdpcmAllowed;
}

bool CU::isMTSAllowed(const CodingUnit &cu, const ComponentID compID) {
  SizeType tsMaxSize = 1 << cu.cs->sps->getLog2MaxTransformSkipBlockSize();
  const int maxSize =
      CU::isIntra(cu) ? MTS_INTRA_MAX_CU_SIZE : MTS_INTER_MAX_CU_SIZE;
  const int cuWidth = cu.blocks[0].lumaSize().width;
  const int cuHeight = cu.blocks[0].lumaSize().height;
  bool mtsAllowed = cu.chType == CHANNEL_TYPE_LUMA && compID == COMPONENT_Y;

  mtsAllowed &= CU::isIntra(cu)
                    ? cu.cs->sps->getUseIntraMTS()
                    : cu.cs->sps->getUseInterMTS() && CU::isInter(cu);
  mtsAllowed &= cuWidth <= maxSize && cuHeight <= maxSize;
  mtsAllowed &= !cu.ispMode;
  mtsAllowed &= !cu.sbtInfo;
  mtsAllowed &=
      !(cu.bdpcmMode && cuWidth <= tsMaxSize && cuHeight <= tsMaxSize);
  return mtsAllowed;
}

bool CU::divideTuInRows(const CodingUnit &cu) {
  CHECK(cu.ispMode != HOR_INTRA_SUBPARTITIONS &&
            cu.ispMode != VER_INTRA_SUBPARTITIONS,
        "Intra Subpartitions type not recognized!");
  return cu.ispMode == HOR_INTRA_SUBPARTITIONS ? true : false;
}

PartSplit CU::getISPType(const CodingUnit &cu, const ComponentID compID) {
  if (cu.ispMode && isLuma(compID)) {
    const bool tuIsDividedInRows = CU::divideTuInRows(cu);

    return tuIsDividedInRows ? TU_1D_HORZ_SPLIT : TU_1D_VERT_SPLIT;
  }
  return TU_NO_ISP;
}

bool CU::isISPFirst(const CodingUnit &cu, const CompArea &tuArea,
                    const ComponentID compID) {
  return tuArea == cu.firstTU->blocks[compID];
}

bool CU::canUseISP(const CodingUnit &cu, const ComponentID compID) {
  const int width = cu.blocks[compID].width;
  const int height = cu.blocks[compID].height;
  const int maxTrSize = cu.cs->sps->getMaxTbSize();
  return CU::canUseISP(width, height, maxTrSize);
}

bool CU::canUseISP(const int width, const int height, const int maxTrSize) {
  bool notEnoughSamplesToSplit =
      (floorLog2(width) + floorLog2(height) <= (floorLog2(MIN_TB_SIZEY) << 1));
  bool cuSizeLargerThanMaxTrSize = width > maxTrSize || height > maxTrSize;
  if (notEnoughSamplesToSplit || cuSizeLargerThanMaxTrSize) {
    return false;
  }
  return true;
}

bool CU::canUseLfnstWithISP(const CompArea &cuArea,
                            const ISPType ispSplitType) {
  if (ispSplitType == NOT_INTRA_SUBPARTITIONS) {
    return false;
  }
  Size tuSize =
      (ispSplitType == HOR_INTRA_SUBPARTITIONS)
          ? Size(cuArea.width, CU::getISPSplitDim(cuArea.width, cuArea.height,
                                                  TU_1D_HORZ_SPLIT))
          : Size(CU::getISPSplitDim(cuArea.width, cuArea.height,
                                    TU_1D_VERT_SPLIT),
                 cuArea.height);

  if (!(tuSize.width >= MIN_TB_SIZEY && tuSize.height >= MIN_TB_SIZEY)) {
    return false;
  }
  return true;
}

bool CU::canUseLfnstWithISP(const CodingUnit &cu, const ChannelType chType) {
  CHECK(!isLuma(chType), "Wrong ISP mode!");
  return CU::canUseLfnstWithISP(cu.blocks[chType == CHANNEL_TYPE_LUMA ? 0 : 1],
                                (ISPType)cu.ispMode);
}

uint32_t CU::getISPSplitDim(const int width, const int height,
                            const PartSplit ispType) {
  bool divideTuInRows = ispType == TU_1D_HORZ_SPLIT;
  uint32_t splitDimensionSize, nonSplitDimensionSize, partitionSize,
      divShift = 2;

  if (divideTuInRows) {
    splitDimensionSize = height;
    nonSplitDimensionSize = width;
  } else {
    splitDimensionSize = width;
    nonSplitDimensionSize = height;
  }

  const int minNumberOfSamplesPerCu = 1 << ((floorLog2(MIN_TB_SIZEY) << 1));
  const int factorToMinSamples =
      nonSplitDimensionSize < minNumberOfSamplesPerCu
          ? minNumberOfSamplesPerCu >> floorLog2(nonSplitDimensionSize)
          : 1;
  partitionSize = (splitDimensionSize >> divShift) < factorToMinSamples
                      ? factorToMinSamples
                      : (splitDimensionSize >> divShift);

  CHECK(floorLog2(partitionSize) + floorLog2(nonSplitDimensionSize) <
            floorLog2(minNumberOfSamplesPerCu),
        "A partition has less than the minimum amount of samples!");
  return partitionSize;
}

PUTraverser CU::traversePUs(CodingUnit &cu) {
  return PUTraverser(cu.firstPU, cu.lastPU->next);
}

TUTraverser CU::traverseTUs(CodingUnit &cu) {
  return TUTraverser(cu.firstTU, cu.lastTU->next);
}

cPUTraverser CU::traversePUs(const CodingUnit &cu) {
  return cPUTraverser(cu.firstPU, cu.lastPU->next);
}

cTUTraverser CU::traverseTUs(const CodingUnit &cu) {
  return cTUTraverser(cu.firstTU, cu.lastTU->next);
}

bool CU::hasSubCUNonZeroMVd(const CodingUnit &cu) {
  bool bNonZeroMvd = false;

  for (const auto &pu : CU::traversePUs(cu)) {
    if ((!pu.mergeFlag) && (!cu.skip)) {
      if (pu.interDir != 2 /* PRED_L1 */) {
        bNonZeroMvd |= pu.mvd[REF_PIC_LIST_0].getHor() != 0;
        bNonZeroMvd |= pu.mvd[REF_PIC_LIST_0].getVer() != 0;
      }
      if (pu.interDir != 1 /* PRED_L0 */) {
        if (!pu.cu->cs->picHeader->getMvdL1ZeroFlag() ||
            pu.interDir != 3 /* PRED_BI */) {
          bNonZeroMvd |= pu.mvd[REF_PIC_LIST_1].getHor() != 0;
          bNonZeroMvd |= pu.mvd[REF_PIC_LIST_1].getVer() != 0;
        }
      }
    }
  }

  return bNonZeroMvd;
}

bool CU::hasSubCUNonZeroAffineMVd(const CodingUnit &cu) {
  bool nonZeroAffineMvd = false;

  if (!cu.affine || cu.firstPU->mergeFlag) {
    return false;
  }

  for (const auto &pu : EntropyCoding::CU::traversePUs(cu)) {
    if ((!pu.mergeFlag) && (!cu.skip)) {
      if (pu.interDir != 2 /* PRED_L1 */) {
        for (int i = 0; i < (cu.affineType == AFFINEMODEL_6PARAM ? 3 : 2);
             i++) {
          nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_0][i].getHor() != 0;
          nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_0][i].getVer() != 0;
        }
      }

      if (pu.interDir != 1 /* PRED_L0 */) {
        if (!pu.cu->cs->picHeader->getMvdL1ZeroFlag() ||
            pu.interDir != 3 /* PRED_BI */) {
          for (int i = 0; i < (cu.affineType == AFFINEMODEL_6PARAM ? 3 : 2);
               i++) {
            nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_1][i].getHor() != 0;
            nonZeroAffineMvd |= pu.mvdAffi[REF_PIC_LIST_1][i].getVer() != 0;
          }
        }
      }
    }
  }

  return nonZeroAffineMvd;
}

uint8_t CU::targetSbtAllowed(uint8_t sbtIdx, uint8_t sbtAllowed) {
  uint8_t val = 0;
  switch (sbtIdx) {
  case SBT_VER_HALF:
    val = ((sbtAllowed >> SBT_VER_HALF) & 0x1);
    break;
  case SBT_HOR_HALF:
    val = ((sbtAllowed >> SBT_HOR_HALF) & 0x1);
    break;
  case SBT_VER_QUAD:
    val = ((sbtAllowed >> SBT_VER_QUAD) & 0x1);
    break;
  case SBT_HOR_QUAD:
    val = ((sbtAllowed >> SBT_HOR_QUAD) & 0x1);
    break;
  default:
    THROW("unknown SBT type");
  }
  return val;
}

int PU::getLMSymbolList(const PredictionUnit &pu, int *modeList) {
  int idx = 0;

  modeList[idx++] = LM_CHROMA_IDX;
  modeList[idx++] = MDLM_L_IDX;
  modeList[idx++] = MDLM_T_IDX;
  return idx;
}

uint32_t PU::getCoLocatedIntraLumaMode(const PredictionUnit &pu)
{
  return PU::getIntraDirLuma(PU::getCoLocatedLumaPU(pu));
}

void PU::getIntraChromaCandModes(const PredictionUnit &pu,
                                 unsigned modeList[NUM_CHROMA_MODE]) {
  modeList[0] = PLANAR_IDX;
  modeList[1] = VER_IDX;
  modeList[2] = HOR_IDX;
  modeList[3] = DC_IDX;
  modeList[4] = LM_CHROMA_IDX;
  modeList[5] = MDLM_L_IDX;
  modeList[6] = MDLM_T_IDX;
  modeList[7] = DM_CHROMA_IDX;

  // If Direct Mode is MIP, mode cannot be already in the list.
  if (isDMChromaMIP(pu)) {
    return;
  }

  const uint32_t lumaMode = getCoLocatedIntraLumaMode(pu);
  for (int i = 0; i < 4; i++) {
    if (lumaMode == modeList[i]) {
      modeList[i] = VDIA_IDX;
      break;
    }
  }
}

int PU::getIntraMPMs(const PredictionUnit &pu, unsigned *mpm,
                     const ChannelType &channelType /*= CHANNEL_TYPE_LUMA*/) {
  const int numMPMs = NUM_MOST_PROBABLE_MODES;
  {
    CHECK(channelType != CHANNEL_TYPE_LUMA, "Not harmonized yet");
    int numCand = -1;
    int leftIntraDir = PLANAR_IDX, aboveIntraDir = PLANAR_IDX;

    const CompArea &area = pu.block(getFirstComponentOfChannel(channelType));
    const Position posRT = area.topRight();
    const Position posLB = area.bottomLeft();

    // Get intra direction of left PU
    const PredictionUnit *puLeft =
        pu.cs->getPURestricted(posLB.offset(-1, 0), pu, channelType);
    if (puLeft && CU::isIntra(*puLeft->cu)) {
      leftIntraDir = PU::getIntraDirLuma(*puLeft);
    }

    // Get intra direction of above PU
    const PredictionUnit *puAbove =
        pu.cs->getPURestricted(posRT.offset(0, -1), pu, channelType);
    if (puAbove && CU::isIntra(*puAbove->cu) &&
        CU::isSameCtu(*pu.cu, *puAbove->cu)) {
      aboveIntraDir = PU::getIntraDirLuma(*puAbove);
    }

    CHECK(2 >= numMPMs, "Invalid number of most probable modes");

    const int offset = (int)NUM_LUMA_MODE - 6;
    const int mod = offset + 3;

    mpm[0] = PLANAR_IDX;
    mpm[1] = DC_IDX;
    mpm[2] = VER_IDX;
    mpm[3] = HOR_IDX;
    mpm[4] = VER_IDX - 4;
    mpm[5] = VER_IDX + 4;

    if (leftIntraDir == aboveIntraDir) {
      numCand = 1;
      if (leftIntraDir > DC_IDX) {
        mpm[0] = PLANAR_IDX;
        mpm[1] = leftIntraDir;
        mpm[2] = ((leftIntraDir + offset) % mod) + 2;
        mpm[3] = ((leftIntraDir - 1) % mod) + 2;
        mpm[4] = ((leftIntraDir + offset - 1) % mod) + 2;
        mpm[5] = (leftIntraDir % mod) + 2;
      }
    } else // L!=A
    {
      numCand = 2;
      int maxCandModeIdx = mpm[0] > mpm[1] ? 0 : 1;

      if ((leftIntraDir > DC_IDX) && (aboveIntraDir > DC_IDX)) {
        mpm[0] = PLANAR_IDX;
        mpm[1] = leftIntraDir;
        mpm[2] = aboveIntraDir;
        maxCandModeIdx = mpm[1] > mpm[2] ? 1 : 2;
        int minCandModeIdx = mpm[1] > mpm[2] ? 2 : 1;
        if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] == 1) {
          mpm[3] = ((mpm[minCandModeIdx] + offset) % mod) + 2;
          mpm[4] = ((mpm[maxCandModeIdx] - 1) % mod) + 2;
          mpm[5] = ((mpm[minCandModeIdx] + offset - 1) % mod) + 2;
        } else if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] >= 62) {
          mpm[3] = ((mpm[minCandModeIdx] - 1) % mod) + 2;
          mpm[4] = ((mpm[maxCandModeIdx] + offset) % mod) + 2;
          mpm[5] = (mpm[minCandModeIdx] % mod) + 2;
        } else if (mpm[maxCandModeIdx] - mpm[minCandModeIdx] == 2) {
          mpm[3] = ((mpm[minCandModeIdx] - 1) % mod) + 2;
          mpm[4] = ((mpm[minCandModeIdx] + offset) % mod) + 2;
          mpm[5] = ((mpm[maxCandModeIdx] - 1) % mod) + 2;
        } else {
          mpm[3] = ((mpm[minCandModeIdx] + offset) % mod) + 2;
          mpm[4] = ((mpm[minCandModeIdx] - 1) % mod) + 2;
          mpm[5] = ((mpm[maxCandModeIdx] + offset) % mod) + 2;
        }
      } else if (leftIntraDir + aboveIntraDir >= 2) {
        mpm[0] = PLANAR_IDX;
        mpm[1] = (leftIntraDir < aboveIntraDir) ? aboveIntraDir : leftIntraDir;
        maxCandModeIdx = 1;
        mpm[2] = ((mpm[maxCandModeIdx] + offset) % mod) + 2;
        mpm[3] = ((mpm[maxCandModeIdx] - 1) % mod) + 2;
        mpm[4] = ((mpm[maxCandModeIdx] + offset - 1) % mod) + 2;
        mpm[5] = (mpm[maxCandModeIdx] % mod) + 2;
      }
    }
    for (int i = 0; i < numMPMs; i++) {
      CHECK(mpm[i] >= NUM_LUMA_MODE, "Invalid MPM");
    }
    CHECK(numCand == 0, "No candidates found");
    return numCand;
  }
}

bool PU::isMIP(const PredictionUnit &pu, const ChannelType &chType) {
  if (chType == CHANNEL_TYPE_LUMA) {
    // Default case if chType is omitted.
    return pu.cu->mipFlag;
  } else {
    return isDMChromaMIP(pu) &&
           (pu.intraDir[CHANNEL_TYPE_CHROMA] == DM_CHROMA_IDX);
  }
}

bool PU::isDMChromaMIP(const PredictionUnit &pu) {
  return !pu.cu->isSepTree() && (pu.chromaFormat == CHROMA_444) &&
         getCoLocatedLumaPU(pu).cu->mipFlag;
}

uint32_t PU::getIntraDirLuma(const PredictionUnit &pu) {
  if (isMIP(pu)) {
    return PLANAR_IDX;
  } else {
    return pu.intraDir[CHANNEL_TYPE_LUMA];
  }
}

const PredictionUnit &PU::getCoLocatedLumaPU(const PredictionUnit &pu) {
  Position topLeftPos = pu.blocks[pu.chType].lumaPos();
  Position refPos =
      topLeftPos.offset(pu.blocks[pu.chType].lumaSize().width >> 1,
                        pu.blocks[pu.chType].lumaSize().height >> 1);
  const PredictionUnit &lumaPU =
      pu.cu->isSepTree() ? *pu.cs->picture->cs->getPU(refPos, CHANNEL_TYPE_LUMA)
                         : *pu.cs->getPU(topLeftPos, CHANNEL_TYPE_LUMA);

  return lumaPU;
}

bool PU::isBipredRestriction(const PredictionUnit &pu) {
  if (pu.cu->lumaSize().width == 4 && pu.cu->lumaSize().height == 4) {
    return true;
  }
  /* disable bi-prediction for 4x8/8x4 */
  if (pu.cu->lumaSize().width + pu.cu->lumaSize().height == 12) {
    return true;
  }
  return false;
}

bool PU::isLMCMode(unsigned mode) {
  return (mode >= LM_CHROMA_IDX && mode <= MDLM_T_IDX);
}

bool TU::getCbf(const TransformUnit &tu, const ComponentID &compID) {
  return getCbfAtDepth(tu, compID, tu.depth);
}

bool TU::getCbfAtDepth(const TransformUnit &tu, const ComponentID &compID,
                       const unsigned &depth) {
  if (!tu.blocks[compID].valid()) {
    CHECK(tu.cbf[compID] != 0,
          "cbf must be 0 if the component is not available");
  }
  return ((tu.cbf[compID] >> depth) & 1) == 1;
}

void TU::setCbfAtDepth(TransformUnit &tu, const ComponentID &compID,
                       const unsigned &depth, const bool &cbf) {
  // first clear the CBF at the depth
  tu.cbf[compID] &= ~(1 << depth);
  // then set the CBF
  tu.cbf[compID] |= ((cbf ? 1 : 0) << depth);
}

bool TU::isTSAllowed(const TransformUnit &tu, const ComponentID compID) {
  const int maxSize = tu.cs->sps->getLog2MaxTransformSkipBlockSize();

  bool tsAllowed = tu.cs->sps->getTransformSkipEnabledFlag();
  tsAllowed &= (!tu.cu->ispMode || !isLuma(compID));
  SizeType transformSkipMaxSize = 1 << maxSize;
  tsAllowed &= !(tu.cu->bdpcmMode && isLuma(compID));
  tsAllowed &= !(tu.cu->bdpcmModeChroma && isChroma(compID));
  tsAllowed &= tu.blocks[compID].width <= transformSkipMaxSize &&
               tu.blocks[compID].height <= transformSkipMaxSize;
  tsAllowed &= !tu.cu->sbtInfo;

  return tsAllowed;
}

TransformUnit *TU::getPrevTU(const TransformUnit &tu,
                             const ComponentID compID) {
  TransformUnit *prevTU = tu.prev;

  if (prevTU != nullptr &&
      (prevTU->cu != tu.cu || !prevTU->blocks[compID].valid())) {
    prevTU = nullptr;
  }

  return prevTU;
}

bool TU::getPrevTuCbfAtDepth(const TransformUnit &currentTu,
                             const ComponentID compID, const int trDepth) {
  const TransformUnit *prevTU = getPrevTU(currentTu, compID);
  return (prevTU != nullptr) ? TU::getCbfAtDepth(*prevTU, compID, trDepth)
                             : false;
}

bool allowLfnstWithMip(const Size &block) {
  if (block.width >= 16 && block.height >= 16) {
    return true;
  }
  return false;
}

int getNumModesMip(const Size &block) {
  switch (getMipSizeId(block)) {
  case 0:
    return 16;
  case 1:
    return 8;
  case 2:
    return 6;
  default:
    THROW("Invalid mipSizeId");
  }
}

int getMipSizeId(const Size &block) {
  if (block.width == 4 && block.height == 4) {
    return 0;
  } else if (block.width == 4 || block.height == 4 ||
             (block.width == 8 && block.height == 8)) {
    return 1;
  } else {
    return 2;
  }
}
} // namespace EntropyCoding