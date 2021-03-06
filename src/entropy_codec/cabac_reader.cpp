#include <algorithm>

#include "cabac_reader.hpp"
#include "sample_adaptive_offset.hpp"
#include "unit_tools.hpp"
#include "coding_structure.hpp"
#include "log.hpp"

using namespace EntropyCoding;
using namespace Common;

#if RExt__DECODER_DEBUG_BIT_STATISTICS
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(x)                       \
  const CodingStatisticsClassType CSCT(x);                                     \
  m_BinDecoder.set(CSCT)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET2(x, y)                   \
  const CodingStatisticsClassType CSCT(x, y);                                  \
  m_BinDecoder.set(CSCT)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE(x, s)               \
  const CodingStatisticsClassType CSCT(x, s.width, s.height);                  \
  m_BinDecoder.set(CSCT)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(x, s, z)           \
  const CodingStatisticsClassType CSCT(x, s.width, s.height, z);               \
  m_BinDecoder.set(CSCT)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_SET(x) m_BinDecoder.set(x);
#else
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(x)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET2(x, y)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE(x, s)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(x, s, z)
#define RExt__DECODER_DEBUG_BIT_STATISTICS_SET(x)
#endif

void CABACReader::initCtxModels(Slice &slice) {
  SliceType sliceType = slice.getSliceType();
  int qp = slice.getSliceQp();
  if (slice.getPPS()->getCabacInitPresentFlag() && slice.getCabacInitFlag()) {
    switch (sliceType) {
    case P_SLICE: // change initialization table to B_SLICE initialization
      sliceType = B_SLICE;
      break;
    case B_SLICE: // change initialization table to P_SLICE initialization
      sliceType = P_SLICE;
      break;
    default: // should not occur
      THROW("Invalid slice type");
      break;
    }
  }
  m_BinDecoder.reset(qp, (int)sliceType);
  // m_BinDecoder.setBaseLevel(slice.getRiceBaseLevel());
#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
  m_BinDecoder.riceStatReset(slice.getSPS()->getBitDepth(CHANNEL_TYPE_LUMA),
                             slice.getSPS()
                                 ->getSpsRangeExtension()
                                 .getPersistentRiceAdaptationEnabledFlag());
#else
  m_BinDecoder.riceStatReset(slice.getSPS()->getBitDepth(CHANNEL_TYPE_LUMA));
#endif
}

//================================================================================
//  clause 7.3.8.1
//--------------------------------------------------------------------------------
//    bool  terminating_bit()
//    void  remaining_bytes( noTrailingBytesExpected )
//================================================================================

bool CABACReader::terminating_bit() {
  if (m_BinDecoder.decodeBinTrm()) {
    m_BinDecoder.finish();
#if RExt__DECODER_DEBUG_BIT_STATISTICS
    CodingStatistics::IncrementStatisticEP(
        STATS__TRAILING_BITS, m_Bitstream->readOutTrailingBits(), 0);
#else
    m_Bitstream->readOutTrailingBits();
#endif
    return true;
  }
  return false;
}

void CABACReader::remaining_bytes(bool noTrailingBytesExpected) {
  if (noTrailingBytesExpected) {
    CHECK(0 != m_Bitstream->getNumBitsLeft(), "Bits left when not supposed");
  } else {
    while (m_Bitstream->getNumBitsLeft()) {
      unsigned trailingNullByte = m_Bitstream->readByte();
      if (trailingNullByte != 0) {
        THROW("Trailing byte should be '0', but has a value of "
              << std::hex << trailingNullByte << std::dec << "\n");
      }
    }
  }
}

//================================================================================
//  clause 7.3.8.2
//--------------------------------------------------------------------------------
//    void  coding_tree_unit( cs, area, qpL, qpC, ctuRsAddr )
//================================================================================

void CABACReader::coding_tree_unit(CodingStructure &cs, const UnitArea &area,
                                   int (&qps)[2], unsigned ctuRsAddr) {
  CUCtx cuCtx(qps[CH_L]);
  QTBTPartitioner partitioner;

  partitioner.initCtu(area, CH_L, *cs.slice);
  cs.treeType = partitioner.treeType = TREE_D;
  cs.modeType = partitioner.modeType = MODE_TYPE_ALL;

  sao(cs, ctuRsAddr);
  if (cs.sps->getALFEnabledFlag() &&
      (cs.slice->getAlfEnabledFlag(COMPONENT_Y))) {
    const PreCalcValues &pcv = *cs.pcv;
    int frame_width_in_ctus = pcv.widthInCtus;
    int ry = ctuRsAddr / frame_width_in_ctus;
    int rx = ctuRsAddr - ry * frame_width_in_ctus;
    const Position pos(rx * cs.pcv->maxCUWidth, ry * cs.pcv->maxCUHeight);
    const uint32_t curSliceIdx = cs.slice->getIndependentSliceIdx();
    const uint32_t curTileIdx = cs.pps->getTileIdx(pos);
    bool leftAvail = cs.getCURestricted(pos.offset(-(int)pcv.maxCUWidth, 0),
                                        pos, curSliceIdx, curTileIdx, CH_L)
                         ? true
                         : false;
    bool aboveAvail = cs.getCURestricted(pos.offset(0, -(int)pcv.maxCUHeight),
                                         pos, curSliceIdx, curTileIdx, CH_L)
                          ? true
                          : false;

    int leftCTUAddr = leftAvail ? ctuRsAddr - 1 : -1;
    int aboveCTUAddr = aboveAvail ? ctuRsAddr - frame_width_in_ctus : -1;

    for (int compIdx = 0; compIdx < MAX_NUM_COMPONENT; compIdx++) {
      if (cs.slice->getAlfEnabledFlag((ComponentID)compIdx)) {
        uint8_t *ctbAlfFlag = cs.slice->getPic()->getAlfCtuEnableFlag(compIdx);
        int ctx = 0;
        ctx += leftCTUAddr > -1 ? (ctbAlfFlag[leftCTUAddr] ? 1 : 0) : 0;
        ctx += aboveCTUAddr > -1 ? (ctbAlfFlag[aboveCTUAddr] ? 1 : 0) : 0;

        RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__ALF);
        ctbAlfFlag[ctuRsAddr] =
            m_BinDecoder.decodeBin(Ctx::ctbAlfFlag(compIdx * 3 + ctx));
        binLogger.LogElements(SyntaxElement::alf_ctb_flag, ctbAlfFlag[ctuRsAddr]);

        if (isLuma((ComponentID)compIdx) && ctbAlfFlag[ctuRsAddr]) {
          readAlfCtuFilterIndex(cs, ctuRsAddr);
        }
        if (isChroma((ComponentID)compIdx)) {
          int apsIdx = cs.slice->getAlfApsIdChroma();
          CHECK(cs.slice->getAlfAPSs()[apsIdx] == nullptr,
                "APS not initialized");
          const AlfParam &alfParam =
              cs.slice->getAlfAPSs()[apsIdx]->getAlfAPSParam();
          const int numAlts = alfParam.numAlternativesChroma;
          uint8_t *ctbAlfAlternative =
              cs.slice->getPic()->getAlfCtuAlternativeData(compIdx);
          ctbAlfAlternative[ctuRsAddr] = 0;
          if (ctbAlfFlag[ctuRsAddr]) {
            uint8_t decoded = 0;
            while (decoded < numAlts - 1 &&
                   m_BinDecoder.decodeBin(Ctx::ctbAlfAlternative(compIdx - 1))) {
              ++decoded;
              binLogger.LogElement(SyntaxElement::alf_ctb_filter_alt_idx);
                   }
            ctbAlfAlternative[ctuRsAddr] = decoded;
          }
        }
      }
    }
  }
  if (cs.sps->getCCALFEnabledFlag()) {
    for (int compIdx = 1; compIdx < getNumberValidComponents(cs.pcv->chrFormat);
         compIdx++) {
      if (cs.slice->m_ccAlfFilterParam.ccAlfFilterEnabled[compIdx - 1]) {
        const int filterCount =
            cs.slice->m_ccAlfFilterParam.ccAlfFilterCount[compIdx - 1];

        const int ry = ctuRsAddr / cs.pcv->widthInCtus;
        const int rx = ctuRsAddr % cs.pcv->widthInCtus;
        const Position lumaPos(rx * cs.pcv->maxCUWidth,
                               ry * cs.pcv->maxCUHeight);

        ccAlfFilterControlIdc(cs, ComponentID(compIdx), ctuRsAddr,
                              cs.slice->m_ccAlfFilterControl[compIdx - 1],
                              lumaPos, filterCount);
      }
    }
  }

  if (CS::isDualITree(cs) && cs.pcv->chrFormat != CHROMA_400 &&
      cs.pcv->maxCUWidth > 64) {
    QTBTPartitioner chromaPartitioner;
    chromaPartitioner.initCtu(area, CH_C, *cs.slice);
    CUCtx cuCtxChroma(qps[CH_C]);
    coding_tree(cs, partitioner, cuCtx, &chromaPartitioner, &cuCtxChroma);
    qps[CH_L] = cuCtx.qp;
    qps[CH_C] = cuCtxChroma.qp;
  } else {
    coding_tree(cs, partitioner, cuCtx);
    qps[CH_L] = cuCtx.qp;
    if (CS::isDualITree(cs) && cs.pcv->chrFormat != CHROMA_400) {
      CUCtx cuCtxChroma(qps[CH_C]);
      partitioner.initCtu(area, CH_C, *cs.slice);
      coding_tree(cs, partitioner, cuCtxChroma);
      qps[CH_C] = cuCtxChroma.qp;
    }
  }
}

void CABACReader::readAlfCtuFilterIndex(CodingStructure &cs,
                                        unsigned ctuRsAddr) {
  short *alfCtbFilterSetIndex = cs.slice->getPic()->getAlfCtbFilterIndex();
  unsigned numAps = cs.slice->getNumAlfApsIdsLuma();
  unsigned numAvailableFiltSets = numAps + NUM_FIXED_FILTER_SETS;
  uint32_t filtIndex = 0;
  if (numAvailableFiltSets > NUM_FIXED_FILTER_SETS) {
    unsigned usePrevFilt = m_BinDecoder.decodeBin(Ctx::AlfUseTemporalFilt());
    binLogger.LogElements(SyntaxElement::alf_use_aps_flag, usePrevFilt);
    if (usePrevFilt) {
      if (numAps > 1) {
        xReadTruncBinCode(filtIndex,
                          numAvailableFiltSets - NUM_FIXED_FILTER_SETS);
        binLogger.LogElements(SyntaxElement::alf_luma_fixed_filter_idx, filtIndex);
      }
      filtIndex += (unsigned)(NUM_FIXED_FILTER_SETS);
    } else {
      xReadTruncBinCode(filtIndex, NUM_FIXED_FILTER_SETS);
      binLogger.LogElements(SyntaxElement::alf_luma_fixed_filter_idx, filtIndex);
    }
  } else {
    xReadTruncBinCode(filtIndex, NUM_FIXED_FILTER_SETS);
    binLogger.LogElements(SyntaxElement::alf_luma_fixed_filter_idx, filtIndex);
  }
  alfCtbFilterSetIndex[ctuRsAddr] = filtIndex;
}
void CABACReader::ccAlfFilterControlIdc(CodingStructure &cs,
                                        const ComponentID compID,
                                        const int curIdx,
                                        uint8_t *filterControlIdc,
                                        Position lumaPos, int filterCount) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
      STATS__CABAC_BITS__CROSS_COMPONENT_ALF_BLOCK_LEVEL_IDC);

  Position leftLumaPos = lumaPos.offset(-(int)cs.pcv->maxCUWidth, 0);
  Position aboveLumaPos = lumaPos.offset(0, -(int)cs.pcv->maxCUWidth);
  const uint32_t curSliceIdx = cs.slice->getIndependentSliceIdx();
  const uint32_t curTileIdx = cs.pps->getTileIdx(lumaPos);
  bool leftAvail =
      cs.getCURestricted(leftLumaPos, lumaPos, curSliceIdx, curTileIdx, CH_L)
          ? true
          : false;
  bool aboveAvail =
      cs.getCURestricted(aboveLumaPos, lumaPos, curSliceIdx, curTileIdx, CH_L)
          ? true
          : false;
  int ctxt = 0;

  if (leftAvail) {
    ctxt += (filterControlIdc[curIdx - 1]) ? 1 : 0;
  }
  if (aboveAvail) {
    ctxt += (filterControlIdc[curIdx - cs.pcv->widthInCtus]) ? 1 : 0;
  }
  ctxt += (compID == COMPONENT_Cr) ? 3 : 0;

  int idcVal = m_BinDecoder.decodeBin(Ctx::CcAlfFilterControlFlag(ctxt));
  binLogger.LogElements(SyntaxElement::alf_ctb_filter_alt_idx, idcVal);
  if (idcVal) {
    while ((idcVal != filterCount) && m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::alf_ctb_filter_alt_idx);
      idcVal++;
    }
  }
  filterControlIdc[curIdx] = idcVal;
}

//================================================================================
//  clause 7.3.8.3
//--------------------------------------------------------------------------------
//    void  sao( slice, ctuRsAddr )
//================================================================================

void CABACReader::sao(CodingStructure &cs, unsigned ctuRsAddr) {
  const SPS &sps = *cs.sps;

  if (!sps.getSAOEnabledFlag()) {
    return;
  }

  const Slice &slice = *cs.slice;
  SAOBlkParam &sao_ctu_pars = cs.picture->getSAO()[ctuRsAddr];
  bool slice_sao_luma_flag = (slice.getSaoEnabledFlag(CHANNEL_TYPE_LUMA));
  bool slice_sao_chroma_flag = (slice.getSaoEnabledFlag(CHANNEL_TYPE_CHROMA) &&
                                sps.getChromaFormatIdc() != CHROMA_400);
  sao_ctu_pars[COMPONENT_Y].modeIdc = SAO_MODE_OFF;
  sao_ctu_pars[COMPONENT_Cb].modeIdc = SAO_MODE_OFF;
  sao_ctu_pars[COMPONENT_Cr].modeIdc = SAO_MODE_OFF;
  if (!slice_sao_luma_flag && !slice_sao_chroma_flag) {
    return;
  }

  // merge
  int frame_width_in_ctus = cs.pcv->widthInCtus;
  int ry = ctuRsAddr / frame_width_in_ctus;
  int rx = ctuRsAddr - ry * frame_width_in_ctus;
  int sao_merge_type = -1;
  const Position pos(rx * cs.pcv->maxCUWidth, ry * cs.pcv->maxCUHeight);
  const unsigned curSliceIdx = cs.slice->getIndependentSliceIdx();

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__SAO);

  const unsigned curTileIdx = cs.pps->getTileIdx(pos);
  if (cs.getCURestricted(pos.offset(-(int)cs.pcv->maxCUWidth, 0), pos,
                         curSliceIdx, curTileIdx, CH_L)) {
    // sao_merge_left_flag
    sao_merge_type += int(m_BinDecoder.decodeBin(Ctx::SaoMergeFlag()));
    binLogger.LogElements(SyntaxElement::sao_merge_left_flag, sao_merge_type);
  }

  if (sao_merge_type < 0 &&
      cs.getCURestricted(pos.offset(0, -(int)cs.pcv->maxCUHeight), pos,
                         curSliceIdx, curTileIdx, CH_L)) {
    // sao_merge_above_flag
    sao_merge_type += int(m_BinDecoder.decodeBin(Ctx::SaoMergeFlag())) << 1;
    binLogger.LogElements(SyntaxElement::sao_merge_up_flag, sao_merge_type);
  }
  if (sao_merge_type >= 0) {
    if (slice_sao_luma_flag || slice_sao_chroma_flag) {
      sao_ctu_pars[COMPONENT_Y].modeIdc = SAO_MODE_MERGE;
      sao_ctu_pars[COMPONENT_Y].typeIdc = sao_merge_type;
    }
    if (slice_sao_chroma_flag) {
      sao_ctu_pars[COMPONENT_Cb].modeIdc = SAO_MODE_MERGE;
      sao_ctu_pars[COMPONENT_Cr].modeIdc = SAO_MODE_MERGE;
      sao_ctu_pars[COMPONENT_Cb].typeIdc = sao_merge_type;
      sao_ctu_pars[COMPONENT_Cr].typeIdc = sao_merge_type;
    }
    return;
  }

  // explicit parameters
  ComponentID firstComp = (slice_sao_luma_flag ? COMPONENT_Y : COMPONENT_Cb);
  ComponentID lastComp = (slice_sao_chroma_flag ? COMPONENT_Cr : COMPONENT_Y);
  for (ComponentID compID = firstComp; compID <= lastComp;
       compID = ComponentID(compID + 1)) {
    SAOOffset &sao_pars = sao_ctu_pars[compID];

    // sao_type_idx_luma / sao_type_idx_chroma
    if (compID != COMPONENT_Cr) {
      if (m_BinDecoder.decodeBin(Ctx::SaoTypeIdx())) {
        binLogger.LogElement(SyntaxElement::sao_type_idx_luma);
        if (m_BinDecoder.decodeBinEP()) {
          // edge offset
          binLogger.LogElements(SyntaxElement::sao_type_idx_luma, 1);
          sao_pars.modeIdc = SAO_MODE_NEW;
          sao_pars.typeIdc = SAO_TYPE_START_EO;
        } else {
          // band offset
          binLogger.LogElements(SyntaxElement::sao_type_idx_luma, 0);
          sao_pars.modeIdc = SAO_MODE_NEW;
          sao_pars.typeIdc = SAO_TYPE_START_BO;
        }
      }
    } else // Cr, follow Cb SAO type
    {
      binLogger.LogElement(SyntaxElement::sao_type_idx_chroma);
      sao_pars.modeIdc = sao_ctu_pars[COMPONENT_Cb].modeIdc;
      sao_pars.typeIdc = sao_ctu_pars[COMPONENT_Cb].typeIdc;
    }
    if (sao_pars.modeIdc == SAO_MODE_OFF) {
      continue;
    }

    // sao_offset_abs
    int offset[4];
    const int maxOffsetQVal = SampleAdaptiveOffset::getMaxOffsetQVal(
        sps.getBitDepth(toChannelType(compID)));
    offset[0] = (int)unary_max_eqprob(maxOffsetQVal);
    offset[1] = (int)unary_max_eqprob(maxOffsetQVal);
    offset[2] = (int)unary_max_eqprob(maxOffsetQVal);
    offset[3] = (int)unary_max_eqprob(maxOffsetQVal);
    binLogger.LogElements(SyntaxElement::sao_offset_abs, offset[0], offset[1], offset[2], offset[3]);

    // band offset mode
    if (sao_pars.typeIdc == SAO_TYPE_START_BO) {
      // sao_offset_sign
      for (int k = 0; k < 4; k++) {
        if (offset[k] && m_BinDecoder.decodeBinEP()) {
          binLogger.LogElement(SyntaxElement::sao_offset_sign_flag);
          offset[k] = -offset[k];
        }
      }
      // sao_band_position
      sao_pars.typeAuxInfo = m_BinDecoder.decodeBinsEP(NUM_SAO_BO_CLASSES_LOG2);
      binLogger.LogElements(SyntaxElement::sao_band_position, sao_pars.typeAuxInfo);
      for (int k = 0; k < 4; k++) {
        sao_pars.offset[(sao_pars.typeAuxInfo + k) % MAX_NUM_SAO_CLASSES] =
            offset[k];
      }
      continue;
    }

    // edge offset mode
    sao_pars.typeAuxInfo = 0;
    if (compID != COMPONENT_Cr) {
      // sao_eo_class_luma / sao_eo_class_chroma
      sao_pars.typeIdc += m_BinDecoder.decodeBinsEP(NUM_SAO_EO_TYPES_LOG2);
      binLogger.LogElement(SyntaxElement::sao_type_idx_luma);
    } else {
      sao_pars.typeIdc = sao_ctu_pars[COMPONENT_Cb].typeIdc;
    }
    sao_pars.offset[SAO_CLASS_EO_FULL_VALLEY] = offset[0];
    sao_pars.offset[SAO_CLASS_EO_HALF_VALLEY] = offset[1];
    sao_pars.offset[SAO_CLASS_EO_PLAIN] = 0;
    sao_pars.offset[SAO_CLASS_EO_HALF_PEAK] = -offset[2];
    sao_pars.offset[SAO_CLASS_EO_FULL_PEAK] = -offset[3];
  }
}

//================================================================================
//  clause 7.3.8.4
//--------------------------------------------------------------------------------
//    void  coding_tree       ( cs, partitioner, cuCtx )
//    bool  split_cu_flag     ( cs, partitioner )
//    split split_cu_mode_mt  ( cs, partitioner )
//================================================================================

void CABACReader::coding_tree(CodingStructure &cs, Partitioner &partitioner,
                              CUCtx &cuCtx, Partitioner *pPartitionerChroma,
                              CUCtx *pCuCtxChroma) {
  const PPS &pps = *cs.pps;
  const UnitArea &currArea = partitioner.currArea();

  // Reset delta QP coding flag and ChromaQPAdjustemt coding flag
  // Note: do not reset qg at chroma CU
  if (pps.getUseDQP() && partitioner.currQgEnable() &&
      !isChroma(partitioner.chType)) {
    cuCtx.qgStart = true;
    cuCtx.isDQPCoded = false;
  }
  if (cs.slice->getUseChromaQpAdj() && partitioner.currQgChromaEnable()) {
    cuCtx.isChromaQpAdjCoded = false;
    cs.chromaQpAdj = 0;
  }

  // Reset delta QP coding flag and ChromaQPAdjustemt coding flag
  if (CS::isDualITree(cs) && pPartitionerChroma != nullptr) {
    if (pps.getUseDQP() && pPartitionerChroma->currQgEnable()) {
      pCuCtxChroma->qgStart = true;
      pCuCtxChroma->isDQPCoded = false;
    }
    if (cs.slice->getUseChromaQpAdj() &&
        pPartitionerChroma->currQgChromaEnable()) {
      pCuCtxChroma->isChromaQpAdjCoded = false;
      cs.chromaQpAdj = 0;
    }
  }

  const PartSplit splitMode = split_cu_mode(cs, partitioner);

  CHECK(!partitioner.canSplit(splitMode, cs), "Got an invalid split!");

  if (splitMode != CU_DONT_SPLIT) {
    if (CS::isDualITree(cs) && pPartitionerChroma != nullptr &&
        (partitioner.currArea().lwidth() >= 64 ||
         partitioner.currArea().lheight() >= 64)) {
      partitioner.splitCurrArea(CU_QUAD_SPLIT, cs);
      pPartitionerChroma->splitCurrArea(CU_QUAD_SPLIT, cs);
      bool beContinue = true;
      bool lumaContinue = true;
      bool chromaContinue = true;

      while (beContinue) {
        if (partitioner.currArea().lwidth() > 64 ||
            partitioner.currArea().lheight() > 64) {
          if (cs.area.blocks[partitioner.chType].contains(
                  partitioner.currArea().blocks[partitioner.chType].pos())) {
            coding_tree(cs, partitioner, cuCtx, pPartitionerChroma,
                        pCuCtxChroma);
          }
          lumaContinue = partitioner.nextPart(cs);
          chromaContinue = pPartitionerChroma->nextPart(cs);
          CHECK(lumaContinue != chromaContinue,
                "luma chroma partition should be matched");
          beContinue = lumaContinue;
        } else {
          // dual tree coding under 64x64 block
          if (cs.area.blocks[partitioner.chType].contains(
                  partitioner.currArea().blocks[partitioner.chType].pos())) {
            coding_tree(cs, partitioner, cuCtx);
          }
          lumaContinue = partitioner.nextPart(cs);
          if (cs.area.blocks[pPartitionerChroma->chType].contains(
                  pPartitionerChroma->currArea()
                      .blocks[pPartitionerChroma->chType]
                      .pos())) {
            coding_tree(cs, *pPartitionerChroma, *pCuCtxChroma);
          }
          chromaContinue = pPartitionerChroma->nextPart(cs);
          CHECK(lumaContinue != chromaContinue,
                "luma chroma partition should be matched");
          beContinue = lumaContinue;
        }
      }
      partitioner.exitCurrSplit();
      pPartitionerChroma->exitCurrSplit();

      // cat the chroma CUs together
      CodingUnit *currentCu =
          cs.getCU(partitioner.currArea().lumaPos(), CHANNEL_TYPE_LUMA);
      CodingUnit *nextCu = nullptr;
      CodingUnit *tempLastLumaCu = nullptr;
      CodingUnit *tempLastChromaCu = nullptr;
      ChannelType currentChType = currentCu->chType;
      while (currentCu->next != nullptr) {
        nextCu = currentCu->next;
        if (currentChType != nextCu->chType &&
            currentChType == CHANNEL_TYPE_LUMA) {
          tempLastLumaCu = currentCu;
          if (tempLastChromaCu != nullptr) // swap
          {
            tempLastChromaCu->next = nextCu;
          }
        } else if (currentChType != nextCu->chType &&
                   currentChType == CHANNEL_TYPE_CHROMA) {
          tempLastChromaCu = currentCu;
          if (tempLastLumaCu != nullptr) // swap
          {
            tempLastLumaCu->next = nextCu;
          }
        }
        currentCu = nextCu;
        currentChType = currentCu->chType;
      }

      CodingUnit *chromaFirstCu = cs.getCU(
          pPartitionerChroma->currArea().chromaPos(), CHANNEL_TYPE_CHROMA);
      tempLastLumaCu->next = chromaFirstCu;
    } else {
      const ModeType modeTypeParent = partitioner.modeType;
      cs.modeType = partitioner.modeType =
          mode_constraint(cs, partitioner, splitMode); // change for child nodes
      // decide chroma split or not
      bool chromaNotSplit = modeTypeParent == MODE_TYPE_ALL &&
                            partitioner.modeType == MODE_TYPE_INTRA;
      CHECK(chromaNotSplit && partitioner.chType != CHANNEL_TYPE_LUMA,
            "chType must be luma");
      if (partitioner.treeType == TREE_D) {
        cs.treeType = partitioner.treeType = chromaNotSplit ? TREE_L : TREE_D;
      }
      partitioner.splitCurrArea(splitMode, cs);
      do {
        if (cs.area.blocks[partitioner.chType].contains(
                partitioner.currArea().blocks[partitioner.chType].pos())) {
          coding_tree(cs, partitioner, cuCtx);
        }
      } while (partitioner.nextPart(cs));

      partitioner.exitCurrSplit();
      if (chromaNotSplit) {
        CHECK(partitioner.chType != CHANNEL_TYPE_LUMA, "must be luma status");
        partitioner.chType = CHANNEL_TYPE_CHROMA;
        cs.treeType = partitioner.treeType = TREE_C;

        if (cs.picture->blocks[partitioner.chType].contains(
                partitioner.currArea().blocks[partitioner.chType].pos())) {
          coding_tree(cs, partitioner, cuCtx);
        }

        // recover treeType
        partitioner.chType = CHANNEL_TYPE_LUMA;
        cs.treeType = partitioner.treeType = TREE_D;
      }

      // recover ModeType
      cs.modeType = partitioner.modeType = modeTypeParent;
    }
    return;
  }

  CodingUnit &cu = cs.addCU(CS::getArea(cs, currArea, partitioner.chType),
                            partitioner.chType);

  partitioner.setCUData(cu);
  cu.slice = cs.slice.get();
  cu.tileIdx = cs.pps->getTileIdx(currArea.lumaPos());
  CHECK(cu.cs->treeType != partitioner.treeType, "treeType mismatch");
  int lumaQPinLocalDualTree = -1;

  // Predict QP on start of quantization group
  if (cuCtx.qgStart) {
    cuCtx.qgStart = false;
    cuCtx.qp = CU::predictQP(cu, cuCtx.qp);
  }

  if (pps.getUseDQP() && partitioner.isSepTree(cs) && isChroma(cu.chType)) {
    const Position chromaCentral(cu.chromaPos().offset(
        cu.chromaSize().width >> 1, cu.chromaSize().height >> 1));
    const Position lumaRefPos(
        chromaCentral.x << getComponentScaleX(COMPONENT_Cb, cu.chromaFormat),
        chromaCentral.y << getComponentScaleY(COMPONENT_Cb, cu.chromaFormat));
    // derive chroma qp, but the chroma qp is saved in cuCtx.qp which is used
    // for luma qp therefore, after decoding the chroma CU, the cuCtx.qp shall
    // be recovered to luma qp in order to decode next luma cu qp
    const CodingUnit *colLumaCu = cs.getLumaCU(lumaRefPos);
    CHECK(colLumaCu == nullptr, "colLumaCU shall exist");
    lumaQPinLocalDualTree = cuCtx.qp;

    if (colLumaCu) {
      cuCtx.qp = colLumaCu->qp;
    }
  }

  cu.qp =
      cuCtx.qp; // NOTE: CU QP can be changed by deltaQP signaling at TU level
  cu.chromaQpAdj =
      cs.chromaQpAdj; // NOTE: CU chroma QP adjustment can be
                      // changed by adjustment signaling at TU level

  // coding unit

  coding_unit(cu, partitioner, cuCtx);
  // recover cuCtx.qp to luma qp after decoding the chroma CU
  if (pps.getUseDQP() && partitioner.isSepTree(cs) && isChroma(cu.chType)) {
    cuCtx.qp = lumaQPinLocalDualTree;
  }

  uint32_t compBegin;
  uint32_t numComp;
  bool jointPLT = false;
  if (cu.isSepTree()) {
    if (cu.isLocalSepTree()) {
      compBegin = COMPONENT_Y;
      numComp = (cu.chromaFormat != CHROMA_400) ? 3 : 1;
      jointPLT = true;
    } else {
      if (isLuma(partitioner.chType)) {
        compBegin = COMPONENT_Y;
        numComp = 1;
      } else {
        compBegin = COMPONENT_Cb;
        numComp = 2;
      }
    }
  } else {
    compBegin = COMPONENT_Y;
    numComp = (cu.chromaFormat != CHROMA_400) ? 3 : 1;
    jointPLT = true;
  }
  if (CU::isPLT(cu)) {
    cs.reorderPrevPLT(cs.prevPLT, cu.curPLTSize, cu.curPLT, cu.reuseflag,
                      compBegin, numComp, jointPLT);
  }
}

ModeType CABACReader::mode_constraint(CodingStructure &cs,
                                      Partitioner &partitioner,
                                      PartSplit splitMode) {
  int val = cs.signalModeCons(splitMode, partitioner, partitioner.modeType);
  if (val == LDT_MODE_TYPE_SIGNAL) {
    int ctxIdx = DeriveCtx::CtxModeConsFlag(cs, partitioner);
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
        STATS__CABAC_BITS__MODE_CONSTRAINT_FLAG,
        partitioner.currArea().blocks[partitioner.chType].size(),
        partitioner.chType);
    bool flag = m_BinDecoder.decodeBin(Ctx::ModeConsFlag(ctxIdx));
    binLogger.LogElements(SyntaxElement::non_inter_flag, flag);
    return flag ? MODE_TYPE_INTRA : MODE_TYPE_INTER;
  } else if (val == LDT_MODE_TYPE_INFER) {
    return MODE_TYPE_INTRA;
  } else {
    return partitioner.modeType;
  }
}

PartSplit CABACReader::split_cu_mode(CodingStructure &cs,
                                     Partitioner &partitioner) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__SPLIT_FLAG,
      partitioner.currArea().blocks[partitioner.chType].size(),
      partitioner.chType);

  PartSplit mode = CU_DONT_SPLIT;

  bool canNo, canQt, canBh, canBv, canTh, canTv;
  partitioner.canSplit(cs, canNo, canQt, canBh, canBv, canTh, canTv);

  bool canSpl[6] = {canNo, canQt, canBh, canBv, canTh, canTv};

  unsigned ctxSplit = 0, ctxQtSplit = 0, ctxBttHV = 0, ctxBttH12 = 0, ctxBttV12;
  DeriveCtx::CtxSplit(cs, partitioner, ctxSplit, ctxQtSplit, ctxBttHV,
                      ctxBttH12, ctxBttV12, canSpl);

  bool isSplit = canBh || canBv || canTh || canTv || canQt;

  if (canNo && isSplit) {
    isSplit = m_BinDecoder.decodeBin(Ctx::SplitFlag(ctxSplit));
    binLogger.LogElements(SyntaxElement::split_cu_flag, isSplit);
  }

  if (!isSplit) {
    return CU_DONT_SPLIT;
  }

  const bool canBtt = canBh || canBv || canTh || canTv;
  bool isQt = canQt;

  if (isQt && canBtt) {
    isQt = m_BinDecoder.decodeBin(Ctx::SplitQtFlag(ctxQtSplit));
    binLogger.LogElements(SyntaxElement::split_qt_flag, isQt);
  }

  if (isQt) {
    return CU_QUAD_SPLIT;
  }

  const bool canHor = canBh || canTh;
  bool isVer = canBv || canTv;

  if (isVer && canHor) {
    isVer = m_BinDecoder.decodeBin(Ctx::SplitHvFlag(ctxBttHV));
    binLogger.LogElements(SyntaxElement::mtt_split_cu_vertical_flag, isVer);
  }

  const bool can14 = isVer ? canTv : canTh;
  bool is12 = isVer ? canBv : canBh;

  if (is12 && can14) {
    is12 =
        m_BinDecoder.decodeBin(Ctx::Split12Flag(isVer ? ctxBttV12 : ctxBttH12));
    binLogger.LogElements(SyntaxElement::mtt_split_cu_binary_flag, is12);
  }

  if (isVer && is12) {
    mode = CU_VERT_SPLIT;
  } else if (isVer && !is12) {
    mode = CU_TRIV_SPLIT;
  } else if (!isVer && is12) {
    mode = CU_HORZ_SPLIT;
  } else {
    mode = CU_TRIH_SPLIT;
  }

  return mode;
}

//================================================================================
//  clause 7.3.8.5
//--------------------------------------------------------------------------------
//    void  coding_unit               ( cu, partitioner, cuCtx )
//    void  cu_skip_flag              ( cu )
//    void  pred_mode                 ( cu )
//    void  part_mode                 ( cu )
//    void  cu_pred_data              ( pus )
//    void  cu_lic_flag               ( cu )
//    void  intra_luma_pred_modes     ( pus )
//    void  intra_chroma_pred_mode    ( pu )
//    void  cu_residual               ( cu, partitioner, cuCtx )
//    void  rqt_root_cbf              ( cu )
//    void  end_of_ctu                ( cu, cuCtx )
//================================================================================

void CABACReader::coding_unit(CodingUnit &cu, Partitioner &partitioner,
                              CUCtx &cuCtx) {
  CodingStructure &cs = *cu.cs;
  CHECK(cu.treeType != partitioner.treeType ||
            cu.modeType != partitioner.modeType,
        "treeType or modeType mismatch");
  PredictionUnit &pu = cs.addPU(cu, partitioner.chType);
  // skip flag
  if ((!cs.slice->isIntra() || cs.slice->getSPS()->getIBCFlag()) &&
      cu.Y().valid()) {
    cu_skip_flag(cu);
  }

  // skip data
  if (cu.skip) {
    cu.colorTransform = false;
    cs.addEmptyTUs(partitioner);
    MergeCtx mrgCtx;
    prediction_unit(pu, mrgCtx);
    end_of_ctu(cu, cuCtx);
    return;
  }

  // prediction mode and partitioning data
  pred_mode(cu);
  if (CU::isIntra(cu)) {
    adaptive_color_transform(cu);
  }
  if (CU::isPLT(cu)) {
    cu.colorTransform = false;
    cs.addTU(cu, partitioner.chType);
    if (cu.isSepTree()) {
      if (isLuma(partitioner.chType)) {
        cu_palette_info(cu, COMPONENT_Y, 1, cuCtx);
      }
      if (cu.chromaFormat != CHROMA_400 &&
          (partitioner.chType == CHANNEL_TYPE_CHROMA)) {
        cu_palette_info(cu, COMPONENT_Cb, 2, cuCtx);
      }
    } else {
      if (cu.chromaFormat != CHROMA_400) {
        cu_palette_info(cu, COMPONENT_Y, 3, cuCtx);
      } else {
        cu_palette_info(cu, COMPONENT_Y, 1, cuCtx);
      }
    }
    end_of_ctu(cu, cuCtx);
    return;
  }

  // --> create PUs

  // prediction data ( intra prediction modes / reference indexes + motion
  // vectors )
  cu_pred_data(cu);

  // residual data ( coded block flags + transform coefficient levels )
  cu_residual(cu, partitioner, cuCtx);

  // check end of cu
  end_of_ctu(cu, cuCtx);
}

void CABACReader::cu_skip_flag(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__SKIP_FLAG);

  if ((cu.slice->isIntra() || cu.isConsIntra()) &&
      cu.cs->slice->getSPS()->getIBCFlag()) {
    cu.skip = false;
    cu.rootCbf = false;
    cu.predMode = MODE_INTRA;
    cu.mmvdSkip = false;
    if (cu.lwidth() < 128 &&
        cu.lheight() < 128) // disable IBC mode larger than 64x64
    {
      unsigned ctxId = DeriveCtx::CtxSkipFlag(cu);
      unsigned skip = m_BinDecoder.decodeBin(Ctx::SkipFlag(ctxId));
      binLogger.LogElements(SyntaxElement::cu_skip_flag, skip);
      if (skip) {
        cu.skip = true;
        cu.rootCbf = false;
        cu.predMode = MODE_IBC;
        cu.mmvdSkip = false;
      }
    }
    return;
  }
  if (!cu.cs->slice->getSPS()->getIBCFlag() && cu.lwidth() == 4 &&
      cu.lheight() == 4) {
    return;
  }
  if (!cu.cs->slice->getSPS()->getIBCFlag() && cu.isConsIntra()) {
    return;
  }
  unsigned ctxId = DeriveCtx::CtxSkipFlag(cu);
  unsigned skip = m_BinDecoder.decodeBin(Ctx::SkipFlag(ctxId));
  binLogger.LogElements(SyntaxElement::cu_skip_flag, skip);

  if (skip && cu.cs->slice->getSPS()->getIBCFlag()) {
    if (cu.lwidth() < 128 && cu.lheight() < 128 &&
        !cu.isConsInter()) // disable IBC mode larger than 64x64 and disable IBC
                           // when only allowing inter mode
    {
      if (cu.lwidth() == 4 && cu.lheight() == 4) {
        cu.skip = true;
        cu.rootCbf = false;
        cu.predMode = MODE_IBC;
        cu.mmvdSkip = false;
        return;
      }
      unsigned ctxidx = DeriveCtx::CtxIBCFlag(cu);
      if (m_BinDecoder.decodeBin(Ctx::IBCFlag(ctxidx))) {
        binLogger.LogElement(SyntaxElement::pred_mode_ibc_flag);
        cu.skip = true;
        cu.rootCbf = false;
        cu.predMode = MODE_IBC;
        cu.mmvdSkip = false;
        cu.firstPU->regularMergeFlag = false;
      } else {
        cu.predMode = MODE_INTER;
      }
    } else {
      cu.predMode = MODE_INTER;
    }
  }
  if ((skip && CU::isInter(cu) && cu.cs->slice->getSPS()->getIBCFlag()) ||
      (skip && !cu.cs->slice->getSPS()->getIBCFlag())) {
    cu.skip = true;
    cu.rootCbf = false;
    cu.predMode = MODE_INTER;
  }
}

void CABACReader::imv_mode(CodingUnit &cu, MergeCtx &mrgCtx) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__OTHER);

  if (!cu.cs->sps->getAMVREnabledFlag()) {
    return;
  }

  bool bNonZeroMvd = CU::hasSubCUNonZeroMVd(cu);
  if (!bNonZeroMvd) {
    return;
  }

  if (cu.affine) {
    return;
  }

  const SPS *sps = cu.cs->sps.get();

  unsigned value = 0;
  if (CU::isIBC(cu)) {
    value = 1;
  } else {
    value = m_BinDecoder.decodeBin(Ctx::ImvFlag(0));
    binLogger.LogElements(SyntaxElement::amvr_flag, value);
  }

  cu.imv = value;
  if (sps->getAMVREnabledFlag() && value) {
    if (!CU::isIBC(cu)) {
      value = m_BinDecoder.decodeBin(Ctx::ImvFlag(4));
      binLogger.LogElements(SyntaxElement::amvr_precision_idx, value);
      cu.imv = value ? 1 : IMV_HPEL;
    }
    if (value) {
      value = m_BinDecoder.decodeBin(Ctx::ImvFlag(1));
      binLogger.LogElements(SyntaxElement::amvr_precision_idx, value);
      value++;
      cu.imv = value;
    }
  }
}

void CABACReader::affine_amvr_mode(CodingUnit &cu, MergeCtx &mrgCtx) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__OTHER);

  const SPS *sps = cu.slice->getSPS();

  if (!sps->getAffineAmvrEnabledFlag() || !cu.affine) {
    return;
  }

  if (!CU::hasSubCUNonZeroAffineMVd(cu)) {
    return;
  }

  unsigned value = 0;
  value = m_BinDecoder.decodeBin(Ctx::ImvFlag(2));
  binLogger.LogElements(SyntaxElement::amvr_flag, value);

  if (value) {
    value = m_BinDecoder.decodeBin(Ctx::ImvFlag(3));
    binLogger.LogElements(SyntaxElement::amvr_precision_idx, value);
    value++;
  }

  cu.imv = value;
}

void CABACReader::pred_mode(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__PRED_MODE);
  if (cu.cs->slice->getSPS()->getIBCFlag() &&
      cu.chType != CHANNEL_TYPE_CHROMA) {
    if (cu.isConsInter()) {
      cu.predMode = MODE_INTER;
      return;
    }

    if (cu.cs->slice->isIntra() || (cu.lwidth() == 4 && cu.lheight() == 4) ||
        cu.isConsIntra()) {
      cu.predMode = MODE_INTRA;
      if (cu.lwidth() < 128 &&
          cu.lheight() < 128) // disable IBC mode larger than 64x64
      {
        unsigned ctxidx = DeriveCtx::CtxIBCFlag(cu);
        if (m_BinDecoder.decodeBin(Ctx::IBCFlag(ctxidx))) {
          binLogger.LogElement(SyntaxElement::pred_mode_ibc_flag);
          cu.predMode = MODE_IBC;
        }
      }
      if (!CU::isIBC(cu) && cu.cs->slice->getSPS()->getPLTMode() &&
          cu.lwidth() <= 64 && cu.lheight() <= 64 &&
          (cu.lumaSize().width * cu.lumaSize().height > 16)) {
        if (m_BinDecoder.decodeBin(Ctx::PLTFlag(0))) {
          binLogger.LogElement(SyntaxElement::pred_mode_plt_flag);
          cu.predMode = MODE_PLT;
        }
      }
    } else {
      if (m_BinDecoder.decodeBin(
              Ctx::PredMode(DeriveCtx::CtxPredModeFlag(cu)))) {
        binLogger.LogElement(SyntaxElement::pred_mode_flag);
        cu.predMode = MODE_INTRA;
        if (cu.cs->slice->getSPS()->getPLTMode() && cu.lwidth() <= 64 &&
            cu.lheight() <= 64 &&
            (cu.lumaSize().width * cu.lumaSize().height > 16)) {
          if (m_BinDecoder.decodeBin(Ctx::PLTFlag(0))) {
            binLogger.LogElement(SyntaxElement::pred_mode_plt_flag);
            cu.predMode = MODE_PLT;
          }
        }
      } else {
        cu.predMode = MODE_INTER;
        if (cu.lwidth() < 128 &&
            cu.lheight() < 128) // disable IBC mode larger than 64x64
        {
          unsigned ctxidx = DeriveCtx::CtxIBCFlag(cu);
          if (m_BinDecoder.decodeBin(Ctx::IBCFlag(ctxidx))) {
            binLogger.LogElement(SyntaxElement::pred_mode_ibc_flag);
            cu.predMode = MODE_IBC;
          }
        }
      }
    }
  } else {
    if (cu.isConsInter()) {
      cu.predMode = MODE_INTER;
      return;
    }

    if (cu.cs->slice->isIntra() || (cu.lwidth() == 4 && cu.lheight() == 4) ||
        cu.isConsIntra()) {
      cu.predMode = MODE_INTRA;
      if (cu.cs->slice->getSPS()->getPLTMode() && cu.lwidth() <= 64 &&
          cu.lheight() <= 64 &&
          (((!isLuma(cu.chType)) &&
            (cu.chromaSize().width * cu.chromaSize().height > 16)) ||
           ((isLuma(cu.chType)) &&
            ((cu.lumaSize().width * cu.lumaSize().height) > 16))) &&
          (!cu.isLocalSepTree() || isLuma(cu.chType))) {
        if (m_BinDecoder.decodeBin(Ctx::PLTFlag(0))) {
          binLogger.LogElement(SyntaxElement::pred_mode_plt_flag);
          cu.predMode = MODE_PLT;
        }
      }
    } else {
      cu.predMode =
          m_BinDecoder.decodeBin(Ctx::PredMode(DeriveCtx::CtxPredModeFlag(cu)))
              ? MODE_INTRA
              : MODE_INTER;
      binLogger.LogElements(SyntaxElement::pred_mode_flag, cu.predMode);
      if (CU::isIntra(cu) && cu.cs->slice->getSPS()->getPLTMode() &&
          cu.lwidth() <= 64 && cu.lheight() <= 64 &&
          (((!isLuma(cu.chType)) &&
            (cu.chromaSize().width * cu.chromaSize().height > 16)) ||
           ((isLuma(cu.chType)) &&
            ((cu.lumaSize().width * cu.lumaSize().height) > 16))) &&
          (!cu.isLocalSepTree() || isLuma(cu.chType))) {
        if (m_BinDecoder.decodeBin(Ctx::PLTFlag(0))) {
          binLogger.LogElement(SyntaxElement::pred_mode_plt_flag);
          cu.predMode = MODE_PLT;
        }
      }
    }
  }
}

void CABACReader::bdpcm_mode(CodingUnit &cu, const ComponentID compID) {
  if (!CU::bdpcmAllowed(cu, compID)) {
    if (isLuma(compID)) {
      cu.bdpcmMode = 0;
      if (!CS::isDualITree(*cu.cs)) {
        cu.bdpcmModeChroma = 0;
      }
    } else {
      cu.bdpcmModeChroma = 0;
    }
    return;
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__BDPCM_MODE, cu.block(compID).lumaSize(), compID);

  int bdpcmMode;
  unsigned ctxId = isLuma(compID) ? 0 : 2;
  bdpcmMode = m_BinDecoder.decodeBin(Ctx::BDPCMMode(ctxId));
  binLogger.LogElements(isLuma(compID) ? SyntaxElement::intra_bdpcm_luma_flag
                                       : SyntaxElement::intra_bdpcm_chroma_flag, bdpcmMode);
  if (bdpcmMode) {
    bdpcmMode += m_BinDecoder.decodeBin(Ctx::BDPCMMode(ctxId + 1));
    binLogger.LogElement(isLuma(compID)
                              ? SyntaxElement::intra_bdpcm_luma_dir_flag
                              : SyntaxElement::intra_bdpcm_chroma_dir_flag);
  }
  if (isLuma(compID)) {
    cu.bdpcmMode = bdpcmMode;
  } else {
    cu.bdpcmModeChroma = bdpcmMode;
  }
}

void CABACReader::cu_pred_data(CodingUnit &cu) {
  if (CU::isIntra(cu)) {
    if (cu.Y().valid()) {
      bdpcm_mode(cu, COMPONENT_Y);
    }
    intra_luma_pred_modes(cu);
    if ((!cu.Y().valid() || (!cu.isSepTree() && cu.Y().valid())) &&
        isChromaEnabled(cu.chromaFormat)) {
      bdpcm_mode(cu, ComponentID(CHANNEL_TYPE_CHROMA));
    }
    intra_chroma_pred_modes(cu);
    return;
  }
  if (!cu.Y().valid()) // dual tree chroma CU
  {
    cu.predMode = MODE_IBC;
    return;
  }
  MergeCtx mrgCtx;

  for (auto &pu : CU::traversePUs(cu)) {
    prediction_unit(pu, mrgCtx);
  }

  imv_mode(cu, mrgCtx);
  affine_amvr_mode(cu, mrgCtx);
  cu_bcw_flag(cu);
}

void CABACReader::cu_bcw_flag(CodingUnit &cu) {
  if (!CU::isBcwIdxCoded(cu)) {
    return;
  }

  CHECK(!(BCW_NUM > 1 && (BCW_NUM == 2 || (BCW_NUM & 0x01) == 1)),
        " !( BCW_NUM > 1 && ( BCW_NUM == 2 || ( BCW_NUM & 0x01 ) == 1 ) ) ");

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__BCW_IDX);

  uint32_t idx = 0;

  uint32_t symbol = m_BinDecoder.decodeBin(Ctx::BcwIdx(0));
  binLogger.LogElements(SyntaxElement::bcw_idx, symbol);

  int32_t numBcw = (cu.slice->getCheckLDC()) ? 5 : 3;
  if (symbol == 1) {
    uint32_t prefixNumBits = numBcw - 2;
    uint32_t step = 1;

    idx = 1;

    for (int ui = 0; ui < prefixNumBits; ++ui) {
      symbol = m_BinDecoder.decodeBinEP();
      binLogger.LogElements(SyntaxElement::bcw_idx, symbol);
      if (symbol == 0) {
        break;
      }
      idx += step;
    }
  }

  uint8_t bcwIdx = (uint8_t)g_BcwParsingOrder[idx];
  CU::setBcwIdx(cu, bcwIdx);
}

void CABACReader::xReadTruncBinCode(uint32_t &symbol, uint32_t maxSymbol) {
  int thresh;
  if (maxSymbol > 256) {
    int threshVal = 1 << 8;
    thresh = 8;
    while (threshVal <= maxSymbol) {
      thresh++;
      threshVal <<= 1;
    }
    thresh--;
  } else {
    thresh = g_tbMax[maxSymbol];
  }

  int val = 1 << thresh;
  int b = maxSymbol - val;
  symbol = m_BinDecoder.decodeBinsEP(thresh);
  if (symbol >= val - b) {
    uint32_t altSymbol;
    altSymbol = m_BinDecoder.decodeBinEP();
    symbol <<= 1;
    symbol += altSymbol;
    symbol -= (val - b);
  }
}

void CABACReader::extend_ref_line(CodingUnit &cu) {
  if (!cu.Y().valid() || cu.predMode != MODE_INTRA || !isLuma(cu.chType) ||
      cu.bdpcmMode) {
    cu.firstPU->multiRefIdx = 0;
    return;
  }
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
      STATS__CABAC_BITS__MULTI_REF_LINE);

  const int numBlocks = CU::getNumPUs(cu);
  PredictionUnit *pu = cu.firstPU;

  for (int k = 0; k < numBlocks; k++) {
    if (!cu.cs->sps->getUseMRL()) {
      pu->multiRefIdx = 0;
      pu = pu->next;
      continue;
    }
    bool isFirstLineOfCtu = (((cu.block(COMPONENT_Y).y) &
                              ((cu.cs->sps)->getMaxCUWidth() - 1)) == 0);
    if (isFirstLineOfCtu) {
      pu->multiRefIdx = 0;
      continue;
    }
    int multiRefIdx = 0;

    if (MRL_NUM_REF_LINES > 1) {
      multiRefIdx = m_BinDecoder.decodeBin(Ctx::MultiRefLineIdx(0)) == 1
                        ? MULTI_REF_LINE_IDX[1]
                        : MULTI_REF_LINE_IDX[0];
      binLogger.LogElements(SyntaxElement::ref_idx_l0, multiRefIdx);
      if (MRL_NUM_REF_LINES > 2 && multiRefIdx != MULTI_REF_LINE_IDX[0]) {
        multiRefIdx = m_BinDecoder.decodeBin(Ctx::MultiRefLineIdx(1)) == 1
                          ? MULTI_REF_LINE_IDX[2]
                          : MULTI_REF_LINE_IDX[1];
        binLogger.LogElements(SyntaxElement::ref_idx_l1, multiRefIdx);
      }
    }
    pu->multiRefIdx = multiRefIdx;
    pu = pu->next;
  }
}

void CABACReader::intra_luma_pred_modes(CodingUnit &cu) {
  if (!cu.Y().valid()) {
    return;
  }

  if (cu.bdpcmMode) {
    cu.firstPU->intraDir[0] = cu.bdpcmMode == 2 ? VER_IDX : HOR_IDX;
    return;
  }

  mip_flag(cu);
  if (cu.mipFlag) {
    mip_pred_modes(cu);
    return;
  }
  extend_ref_line(cu);
  isp_mode(cu);

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__INTRA_DIR_ANG, cu.lumaSize(), CHANNEL_TYPE_LUMA);

  // prev_intra_luma_pred_flag
  int numBlocks = CU::getNumPUs(cu);
  int mpmFlag[4];
  for (int k = 0; k < numBlocks; k++) {
    CHECK(numBlocks != 1, "not supported yet");
    if (cu.firstPU->multiRefIdx) {
      mpmFlag[0] = true;
    } else {
      mpmFlag[k] = m_BinDecoder.decodeBin(Ctx::IntraLumaMpmFlag());
      binLogger.LogElements(SyntaxElement::intra_luma_mpm_flag, mpmFlag[k]);
    }
  }

  PredictionUnit *pu = cu.firstPU;

  unsigned
      mpm_pred[NUM_MOST_PROBABLE_MODES]; // mpm_idx / rem_intra_luma_pred_mode
  for (int k = 0; k < numBlocks; k++) {
    PU::getIntraMPMs(*pu, mpm_pred);

    if (mpmFlag[k]) {
      uint32_t ipred_idx = 0;
      {
        unsigned ctx = (pu->cu->ispMode == NOT_INTRA_SUBPARTITIONS ? 1 : 0);
        if (pu->multiRefIdx == 0) {
          ipred_idx = m_BinDecoder.decodeBin(Ctx::IntraLumaPlanarFlag(ctx));
          binLogger.LogElements(SyntaxElement::intra_luma_not_planar_flag, ipred_idx);
        } else {
          ipred_idx = 1;
        }
        if (ipred_idx) {
          ipred_idx += m_BinDecoder.decodeBinEP();
          binLogger.LogElement(SyntaxElement::intra_luma_mpm_idx);
        }
        if (ipred_idx > 1) {
          ipred_idx += m_BinDecoder.decodeBinEP();
          binLogger.LogElement(SyntaxElement::intra_luma_mpm_idx);
        }
        if (ipred_idx > 2) {
          ipred_idx += m_BinDecoder.decodeBinEP();
          binLogger.LogElement(SyntaxElement::intra_luma_mpm_idx);
        }
        if (ipred_idx > 3) {
          ipred_idx += m_BinDecoder.decodeBinEP();
          binLogger.LogElement(SyntaxElement::intra_luma_mpm_idx);
        }
      }
      pu->intraDir[0] = mpm_pred[ipred_idx];
    } else {
      unsigned ipred_mode = 0;

      xReadTruncBinCode(ipred_mode, NUM_LUMA_MODE - NUM_MOST_PROBABLE_MODES);
      binLogger.LogElements(SyntaxElement::intra_luma_mpm_remainder, ipred_mode);
      // postponed sorting of MPMs (only in remaining branch)
      std::sort(mpm_pred, mpm_pred + NUM_MOST_PROBABLE_MODES);

      for (uint32_t i = 0; i < NUM_MOST_PROBABLE_MODES; i++) {
        ipred_mode += (ipred_mode >= mpm_pred[i]);
      }

      pu->intraDir[0] = ipred_mode;
    }
    pu = pu->next;
  }
}

void CABACReader::intra_chroma_pred_modes(CodingUnit &cu) {
  if (cu.chromaFormat == CHROMA_400 ||
      (cu.isSepTree() && cu.chType == CHANNEL_TYPE_LUMA)) {
    return;
  }

  if (cu.bdpcmModeChroma) {
    cu.firstPU->intraDir[1] = cu.bdpcmModeChroma == 2 ? VER_IDX : HOR_IDX;
    return;
  }
  PredictionUnit *pu = cu.firstPU;

  CHECK(pu->cu != &cu, "Inkonsistent PU-CU mapping");
  intra_chroma_pred_mode(*pu);
}

bool CABACReader::intra_chroma_lmc_mode(PredictionUnit &pu) {
  int lmModeList[10];
  PU::getLMSymbolList(pu, lmModeList);

  int symbol = m_BinDecoder.decodeBin(Ctx::CclmModeIdx(0));
  binLogger.LogElements(SyntaxElement::cclm_mode_idx, symbol);

  if (symbol == 0) {
    pu.intraDir[1] = lmModeList[symbol];
    CHECK(pu.intraDir[1] != LM_CHROMA_IDX, "should be LM_CHROMA");
  } else {
    symbol += m_BinDecoder.decodeBinEP();
    binLogger.LogElement(SyntaxElement::cclm_mode_idx);
    pu.intraDir[1] = lmModeList[symbol];
  }
  return true; // it will only enter this function for LMC modes, so always
               // return true ;
}

void CABACReader::intra_chroma_pred_mode(PredictionUnit &pu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__INTRA_DIR_ANG, pu.cu->blocks[pu.chType].lumaSize(),
      CHANNEL_TYPE_CHROMA);
  if (pu.cu->colorTransform) {
    pu.intraDir[CHANNEL_TYPE_CHROMA] = DM_CHROMA_IDX;
    return;
  }

  // LM chroma mode

  if (pu.cs->sps->getUseLMChroma() && pu.cu->checkCCLMAllowed()) {
    bool isLMCMode =
        m_BinDecoder.decodeBin(Ctx::CclmModeFlag(0)) ? true : false;
    binLogger.LogElements(SyntaxElement::cclm_mode_flag, isLMCMode);
    if (isLMCMode) {
      intra_chroma_lmc_mode(pu);
      return;
    }
  }

  if (m_BinDecoder.decodeBin(Ctx::IntraChromaPredMode(0)) == 0) {
    binLogger.LogElement(SyntaxElement::intra_chroma_pred_mode);
    pu.intraDir[1] = DM_CHROMA_IDX;
    return;
  }

  unsigned candId = m_BinDecoder.decodeBinsEP(2);
  binLogger.LogElements(SyntaxElement::intra_chroma_pred_mode, candId);

  unsigned chromaCandModes[NUM_CHROMA_MODE];
  PU::getIntraChromaCandModes(pu, chromaCandModes);

  CHECK(candId >= NUM_CHROMA_MODE,
        "Chroma prediction mode index out of bounds");
  CHECK(PU::isLMCMode(chromaCandModes[candId]),
        "The intra dir cannot be LM_CHROMA for this path");
  CHECK(chromaCandModes[candId] == DM_CHROMA_IDX,
        "The intra dir cannot be DM_CHROMA for this path");

  pu.intraDir[1] = chromaCandModes[candId];
}

void CABACReader::cu_residual(CodingUnit &cu, Partitioner &partitioner,
                              CUCtx &cuCtx) {
  if (!CU::isIntra(cu)) {
    PredictionUnit &pu = *cu.firstPU;
    if (!pu.mergeFlag) {
      rqt_root_cbf(cu);
    } else {
      cu.rootCbf = true;
    }
    if (cu.rootCbf) {
      sbt_mode(cu);
    }
    if (!cu.rootCbf) {
      cu.colorTransform = false;
      cu.cs->addEmptyTUs(partitioner);
      return;
    }
  }

  if (CU::isInter(cu) || CU::isIBC(cu)) {
    adaptive_color_transform(cu);
  }

  cuCtx.violatesLfnstConstrained[CHANNEL_TYPE_LUMA] = false;
  cuCtx.violatesLfnstConstrained[CHANNEL_TYPE_CHROMA] = false;
  cuCtx.lfnstLastScanPos = false;
  cuCtx.violatesMtsCoeffConstraint = false;
  cuCtx.mtsLastScanPos = false;

  ChromaCbfs chromaCbfs;
  if (cu.ispMode && isLuma(partitioner.chType)) {
    TUIntraSubPartitioner subTuPartitioner(partitioner);
    transform_tree(
        *cu.cs, subTuPartitioner, cuCtx,
        CU::getISPType(cu, getFirstComponentOfChannel(partitioner.chType)), 0);
  } else {
    transform_tree(*cu.cs, partitioner, cuCtx);
  }

  residual_lfnst_mode(cu, cuCtx);
  mts_idx(cu, cuCtx);
}

void CABACReader::rqt_root_cbf(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__QT_ROOT_CBF);

  cu.rootCbf = (m_BinDecoder.decodeBin(Ctx::QtRootCbf()));
  binLogger.LogElements(SyntaxElement::cu_coded_flag, cu.rootCbf);
}

void CABACReader::adaptive_color_transform(CodingUnit &cu) {
  if (!cu.slice->getSPS()->getUseColorTrans()) {
    return;
  }

  if (cu.isSepTree()) {
    return;
  }

  if (CU::isInter(cu) || CU::isIBC(cu) || CU::isIntra(cu)) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__ACT);
    cu.colorTransform = (m_BinDecoder.decodeBin(Ctx::ACTFlag()));
    binLogger.LogElements(SyntaxElement::cu_act_enabled_flag, cu.colorTransform);
  }
}

void CABACReader::sbt_mode(CodingUnit &cu) {
  const uint8_t sbtAllowed = cu.checkAllowedSbt();
  if (!sbtAllowed) {
    return;
  }

  SizeType cuWidth = cu.lwidth();
  SizeType cuHeight = cu.lheight();

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__SBT_MODE);
  // bin - flag
  uint8_t ctxIdx = (cuWidth * cuHeight <= 256) ? 1 : 0;
  bool sbtFlag = m_BinDecoder.decodeBin(Ctx::SbtFlag(ctxIdx));
  binLogger.LogElements(SyntaxElement::cu_sbt_flag, sbtFlag);
  if (!sbtFlag) {
    return;
  }

  uint8_t sbtVerHalfAllow = CU::targetSbtAllowed(SBT_VER_HALF, sbtAllowed);
  uint8_t sbtHorHalfAllow = CU::targetSbtAllowed(SBT_HOR_HALF, sbtAllowed);
  uint8_t sbtVerQuadAllow = CU::targetSbtAllowed(SBT_VER_QUAD, sbtAllowed);
  uint8_t sbtHorQuadAllow = CU::targetSbtAllowed(SBT_HOR_QUAD, sbtAllowed);

  // bin - type
  bool sbtQuadFlag = false;
  if ((sbtHorHalfAllow || sbtVerHalfAllow) &&
      (sbtHorQuadAllow || sbtVerQuadAllow)) {
    sbtQuadFlag = m_BinDecoder.decodeBin(Ctx::SbtQuadFlag(0));
    binLogger.LogElements(SyntaxElement::cu_sbt_quad_flag, sbtQuadFlag);
  } else {
    sbtQuadFlag = 0;
  }

  // bin - dir
  bool sbtHorFlag = false;
  if ((sbtQuadFlag && sbtVerQuadAllow && sbtHorQuadAllow) ||
      (!sbtQuadFlag && sbtVerHalfAllow &&
       sbtHorHalfAllow)) // both direction allowed
  {
    uint8_t ctxIdx = (cuWidth == cuHeight) ? 0 : (cuWidth < cuHeight ? 1 : 2);
    sbtHorFlag = m_BinDecoder.decodeBin(Ctx::SbtHorFlag(ctxIdx));
    binLogger.LogElements(SyntaxElement::cu_sbt_horizontal_flag, sbtHorFlag);
  } else {
    sbtHorFlag =
        (sbtQuadFlag && sbtHorQuadAllow) || (!sbtQuadFlag && sbtHorHalfAllow);
  }
  cu.setSbtIdx(sbtHorFlag ? (sbtQuadFlag ? SBT_HOR_QUAD : SBT_HOR_HALF)
                          : (sbtQuadFlag ? SBT_VER_QUAD : SBT_VER_HALF));

  // bin - pos
  bool sbtPosFlag = m_BinDecoder.decodeBin(Ctx::SbtPosFlag(0));
  binLogger.LogElements(SyntaxElement::cu_sbt_pos_flag, sbtPosFlag);
  cu.setSbtPos(sbtPosFlag ? SBT_POS1 : SBT_POS0);
}

void CABACReader::end_of_ctu(CodingUnit &cu, CUCtx &cuCtx) {
  const Position rbPos =
      recalcPosition(cu.chromaFormat, cu.chType, CHANNEL_TYPE_LUMA,
                     cu.blocks[cu.chType].bottomRight().offset(1, 1));

  if (((rbPos.x & cu.cs->pcv->maxCUWidthMask) == 0 ||
       rbPos.x == cu.cs->pps->getPicWidthInLumaSamples()) &&
      ((rbPos.y & cu.cs->pcv->maxCUHeightMask) == 0 ||
       rbPos.y == cu.cs->pps->getPicHeightInLumaSamples()) &&
      (!cu.isSepTree() || cu.chromaFormat == CHROMA_400 ||
       isChroma(cu.chType))) {
    cuCtx.isDQPCoded = (cu.cs->pps->getUseDQP() && !cuCtx.isDQPCoded);
  }
}

void CABACReader::cu_palette_info(CodingUnit &cu, ComponentID compBegin,
                                  uint32_t numComp, CUCtx &cuCtx) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__PLT_MODE);

  const SPS &sps = *(cu.cs->sps);
  TransformUnit &tu = *cu.firstTU;
  int curPLTidx = 0;

  if (cu.isLocalSepTree())
    cu.cs->prevPLT.curPLTSize[compBegin] =
        cu.cs->prevPLT.curPLTSize[COMPONENT_Y];
  cu.lastPLTSize[compBegin] = cu.cs->prevPLT.curPLTSize[compBegin];

  int maxPltSize = cu.isSepTree() ? MAXPLTSIZE_DUALTREE : MAXPLTSIZE;

  if (cu.lastPLTSize[compBegin]) {
    xDecodePLTPredIndicator(cu, maxPltSize, compBegin);
  }

  for (int idx = 0; idx < cu.lastPLTSize[compBegin]; idx++) {
    if (cu.reuseflag[compBegin][idx]) {
      if (cu.isLocalSepTree()) {
        for (int comp = COMPONENT_Y; comp < MAX_NUM_COMPONENT; comp++) {
          cu.curPLT[comp][curPLTidx] = cu.cs->prevPLT.curPLT[comp][idx];
        }
      } else {
        for (int comp = compBegin; comp < (compBegin + numComp); comp++) {
          cu.curPLT[comp][curPLTidx] = cu.cs->prevPLT.curPLT[comp][idx];
        }
      }
      curPLTidx++;
    }
  }

  uint32_t recievedPLTnum = 0;
  if (curPLTidx < maxPltSize) {
    recievedPLTnum = exp_golomb_eqprob(0);
    binLogger.LogElements(SyntaxElement::new_palette_entries, recievedPLTnum);
  }

  cu.curPLTSize[compBegin] = curPLTidx + recievedPLTnum;
  if (cu.isLocalSepTree())
    cu.curPLTSize[COMPONENT_Y] = cu.curPLTSize[compBegin];
  for (int comp = compBegin; comp < (compBegin + numComp); comp++) {
    for (int idx = curPLTidx; idx < cu.curPLTSize[compBegin]; idx++) {
      ComponentID compID = (ComponentID)comp;
      const int channelBitDepth = sps.getBitDepth(toChannelType(compID));
      cu.curPLT[compID][idx] = m_BinDecoder.decodeBinsEP(channelBitDepth);
      binLogger.LogElements(SyntaxElement::palette_idx_idc, cu.curPLT[compID][idx]);
      if (cu.isLocalSepTree()) {
        if (isLuma(cu.chType)) {
          cu.curPLT[COMPONENT_Cb][idx] =
              1 << (cu.cs->sps->getBitDepth(CHANNEL_TYPE_CHROMA) - 1);
          cu.curPLT[COMPONENT_Cr][idx] =
              1 << (cu.cs->sps->getBitDepth(CHANNEL_TYPE_CHROMA) - 1);
        } else {
          cu.curPLT[COMPONENT_Y][idx] =
              1 << (cu.cs->sps->getBitDepth(CHANNEL_TYPE_LUMA) - 1);
        }
      }
    }
  }
  cu.useEscape[compBegin] = true;
  if (cu.curPLTSize[compBegin] > 0) {
    uint32_t escCode = 0;
    escCode = m_BinDecoder.decodeBinEP();
    binLogger.LogElements(SyntaxElement::palette_escape_val_present_flag, escCode);
    cu.useEscape[compBegin] = (escCode != 0);
  }
  uint32_t indexMaxSize = cu.useEscape[compBegin]
                              ? (cu.curPLTSize[compBegin] + 1)
                              : cu.curPLTSize[compBegin];
  // encode index map
  uint32_t height = cu.block(compBegin).height;
  uint32_t width = cu.block(compBegin).width;

  uint32_t total = height * width;
  if (indexMaxSize > 1) {
    parseScanRotationModeFlag(cu, compBegin);
  } else {
    cu.useRotation[compBegin] = false;
  }

  if (cu.useEscape[compBegin] && cu.cs->pps->getUseDQP() && !cuCtx.isDQPCoded) {
    if (!cu.isSepTree() || isLuma(tu.chType)) {
      cu_qp_delta(cu, cuCtx.qp, cu.qp);
      cuCtx.qp = cu.qp;
      cuCtx.isDQPCoded = true;
    }
  }
  if (cu.useEscape[compBegin] && cu.cs->slice->getUseChromaQpAdj() &&
      !cuCtx.isChromaQpAdjCoded) {
    if (!cu.isSepTree() || isChroma(tu.chType)) {
      cu_chroma_qp_offset(cu);
      cuCtx.isChromaQpAdjCoded = true;
    }
  }

  m_scanOrder =
      g_scanOrder[SCAN_UNGROUPED]
                 [(cu.useRotation[compBegin]) ? SCAN_TRAV_VER : SCAN_TRAV_HOR]
                 [gp_sizeIdxInfo->idxFrom(width)]
                 [gp_sizeIdxInfo->idxFrom(height)];
  uint32_t prevRunPos = 0;
  unsigned prevRunType = 0;
  for (int subSetId = 0; subSetId <= (total - 1) >> LOG2_PALETTE_CG_SIZE;
       subSetId++) {
    cuPaletteSubblockInfo(cu, compBegin, numComp, subSetId, prevRunPos,
                          prevRunType);
  }
  CHECK(cu.curPLTSize[compBegin] > maxPltSize,
        " Current palette size is larger than maximum palette size");
}

void CABACReader::cuPaletteSubblockInfo(CodingUnit &cu, ComponentID compBegin,
                                        uint32_t numComp, int subSetId,
                                        uint32_t &prevRunPos,
                                        unsigned &prevRunType) {
  const SPS &sps = *(cu.cs->sps);
  TransformUnit &tu = *cu.firstTU;
  PLTtypeBuf runType = tu.getrunType(compBegin);
  PelBuf curPLTIdx = tu.getcurPLTIdx(compBegin);
  uint32_t indexMaxSize = cu.useEscape[compBegin]
                              ? (cu.curPLTSize[compBegin] + 1)
                              : cu.curPLTSize[compBegin];
  uint32_t totalPel = cu.block(compBegin).height * cu.block(compBegin).width;

  int minSubPos = subSetId << LOG2_PALETTE_CG_SIZE;
  int maxSubPos = minSubPos + (1 << LOG2_PALETTE_CG_SIZE);
  maxSubPos = (maxSubPos > totalPel)
                  ? totalPel
                  : maxSubPos; // if last position is out of the current CU size

  unsigned runCopyFlag[(1 << LOG2_PALETTE_CG_SIZE)];
  for (int i = 0; i < (1 << LOG2_PALETTE_CG_SIZE); i++) {
    runCopyFlag[i] = MAX_INT;
  }
  if (minSubPos == 0) {
    runCopyFlag[0] = 0;
  }

  // PLT runCopy flag and runType - context coded
  int curPos = minSubPos;
  for (; curPos < maxSubPos && indexMaxSize > 1; curPos++) {
    uint32_t posy = m_scanOrder[curPos].y;
    uint32_t posx = m_scanOrder[curPos].x;
    uint32_t posyprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].y;
    uint32_t posxprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].x;
    unsigned identityFlag = 1;

    const CtxSet &ctxSet =
        (prevRunType == PLT_RUN_INDEX) ? Ctx::IdxRunModel : Ctx::CopyRunModel;
    if (curPos > 0) {
      int dist = curPos - prevRunPos - 1;
      const unsigned ctxId = DeriveCtx::CtxPltCopyFlag(prevRunType, dist);
      identityFlag = m_BinDecoder.decodeBin(ctxSet(ctxId));
      binLogger.LogElements(SyntaxElement::run_copy_flag, identityFlag);
      runCopyFlag[curPos - minSubPos] = identityFlag;
    }

    if (identityFlag == 0 || curPos == 0) {
      if (((posy == 0) && !cu.useRotation[compBegin]) ||
          ((posx == 0) && cu.useRotation[compBegin])) {
        runType.at(posx, posy) = PLT_RUN_INDEX;
      } else if (curPos != 0 &&
                 runType.at(posxprev, posyprev) == PLT_RUN_COPY) {
        runType.at(posx, posy) = PLT_RUN_INDEX;
      } else {
        runType.at(posx, posy) = (m_BinDecoder.decodeBin(Ctx::RunTypeFlag()));
        binLogger.LogElements(SyntaxElement::copy_above_palette_indices_flag, runType.at(posx, posy));
      }
      prevRunType = runType.at(posx, posy);
      prevRunPos = curPos;
    } else // assign run information
    {
      runType.at(posx, posy) = runType.at(posxprev, posyprev);
    }
  }

  // PLT index values - bypass coded
  uint32_t adjust;
  uint32_t symbol = 0;
  curPos = minSubPos;
  if (indexMaxSize > 1) {
    for (; curPos < maxSubPos; curPos++) {
      if (curPos > 0) {
        adjust = 1;
      } else {
        adjust = 0;
      }

      uint32_t posy = m_scanOrder[curPos].y;
      uint32_t posx = m_scanOrder[curPos].x;
      uint32_t posyprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].y;
      uint32_t posxprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].x;
      if (runCopyFlag[curPos - minSubPos] == 0 &&
          runType.at(posx, posy) == PLT_RUN_INDEX) {
        xReadTruncBinCode(symbol, indexMaxSize - adjust);
        binLogger.LogElements(SyntaxElement::pred_mode_plt_flag, symbol);
        xAdjustPLTIndex(cu, symbol, curPos, curPLTIdx, runType, indexMaxSize,
                        compBegin);
      } else if (runType.at(posx, posy) == PLT_RUN_INDEX) {
        curPLTIdx.at(posx, posy) = curPLTIdx.at(posxprev, posyprev);
      } else {
        curPLTIdx.at(posx, posy) = (cu.useRotation[compBegin])
                                       ? curPLTIdx.at(posx - 1, posy)
                                       : curPLTIdx.at(posx, posy - 1);
      }
    }
  } else {
    for (; curPos < maxSubPos; curPos++) {
      uint32_t posy = m_scanOrder[curPos].y;
      uint32_t posx = m_scanOrder[curPos].x;
      uint32_t posyprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].y;
      uint32_t posxprev = (curPos == 0) ? 0 : m_scanOrder[curPos - 1].x;
      runType.at(posx, posy) = PLT_RUN_INDEX;
      if (runCopyFlag[curPos - minSubPos] == 0 &&
          runType.at(posx, posy) == PLT_RUN_INDEX) {
        curPLTIdx.at(posx, posy) = 0;
      } else {
        curPLTIdx.at(posx, posy) = curPLTIdx.at(posxprev, posyprev);
      }
    }
  }

  // Quantized escape colors - bypass coded
  uint32_t scaleX = getComponentScaleX(COMPONENT_Cb, sps.getChromaFormatIdc());
  uint32_t scaleY = getComponentScaleY(COMPONENT_Cb, sps.getChromaFormatIdc());
  for (int comp = compBegin; comp < (compBegin + numComp); comp++) {
    ComponentID compID = (ComponentID)comp;
    for (curPos = minSubPos; curPos < maxSubPos; curPos++) {
      uint32_t posy = m_scanOrder[curPos].y;
      uint32_t posx = m_scanOrder[curPos].x;
      if (curPLTIdx.at(posx, posy) == cu.curPLTSize[compBegin]) {
        PLTescapeBuf escapeValue = tu.getescapeValue((ComponentID)comp);
        if (compID == COMPONENT_Y || compBegin != COMPONENT_Y) {
          escapeValue.at(posx, posy) = exp_golomb_eqprob(5);
          binLogger.LogElements(SyntaxElement::palette_escape_val, escapeValue.at(posx, posy));
          assert(escapeValue.at(posx, posy) <
                 (TCoeff(1) << (cu.cs->sps->getBitDepth(
                                    toChannelType((ComponentID)comp)) +
                                1)));
        }
        if (compBegin == COMPONENT_Y && compID != COMPONENT_Y &&
            posy % (1 << scaleY) == 0 && posx % (1 << scaleX) == 0) {
          uint32_t posxC = posx >> scaleX;
          uint32_t posyC = posy >> scaleY;
          escapeValue.at(posxC, posyC) = exp_golomb_eqprob(5);
          binLogger.LogElements(SyntaxElement::palette_escape_val, escapeValue.at(posxC, posyC));
          assert(escapeValue.at(posxC, posyC) <
                 (TCoeff(1)
                  << (cu.cs->sps->getBitDepth(toChannelType(compID)) + 1)));
        }
      }
    }
  }
}

void CABACReader::parseScanRotationModeFlag(CodingUnit &cu,
                                            ComponentID compBegin) {
  cu.useRotation[compBegin] = m_BinDecoder.decodeBin(Ctx::RotationFlag());
  binLogger.LogElements(SyntaxElement::palette_transpose_flag, cu.useRotation[compBegin]);
}

void CABACReader::xDecodePLTPredIndicator(CodingUnit &cu, uint32_t maxPLTSize,
                                          ComponentID compBegin) {
  uint32_t symbol, numPltPredicted = 0, idx = 0;

  symbol = exp_golomb_eqprob(0);
  binLogger.LogElements(SyntaxElement::palette_predictor_run, symbol);

  if (symbol != 1) {
    while (idx < cu.lastPLTSize[compBegin] && numPltPredicted < maxPLTSize) {
      if (idx > 0) {
        symbol = exp_golomb_eqprob(0);
        binLogger.LogElements(SyntaxElement::palette_predictor_run, symbol);
      }
      if (symbol == 1) {
        break;
      }

      if (symbol) {
        idx += symbol - 1;
      }
      cu.reuseflag[compBegin][idx] = 1;
      if (cu.isLocalSepTree()) {
        cu.reuseflag[COMPONENT_Y][idx] = 1;
      }
      numPltPredicted++;
      idx++;
    }
  }
}
void CABACReader::xAdjustPLTIndex(CodingUnit &cu, Pel curLevel, uint32_t idx,
                                  PelBuf &paletteIdx,
                                  PLTtypeBuf &paletteRunType, int maxSymbol,
                                  ComponentID compBegin) {
  uint32_t symbol;
  int refLevel = MAX_INT;
  uint32_t posy = m_scanOrder[idx].y;
  uint32_t posx = m_scanOrder[idx].x;
  if (idx) {
    uint32_t prevposy = m_scanOrder[idx - 1].y;
    uint32_t prevposx = m_scanOrder[idx - 1].x;
    if (paletteRunType.at(prevposx, prevposy) == PLT_RUN_INDEX) {
      refLevel = paletteIdx.at(prevposx, prevposy);
      if (paletteIdx.at(prevposx, prevposy) ==
          cu.curPLTSize[compBegin]) // escape
      {
        refLevel = maxSymbol - 1;
      }
    } else {
      if (cu.useRotation[compBegin]) {
        assert(prevposx > 0);
        refLevel = paletteIdx.at(posx - 1, posy);
        if (paletteIdx.at(posx - 1, posy) ==
            cu.curPLTSize[compBegin]) // escape mode
        {
          refLevel = maxSymbol - 1;
        }
      } else {
        assert(prevposy > 0);
        refLevel = paletteIdx.at(posx, posy - 1);
        if (paletteIdx.at(posx, posy - 1) ==
            cu.curPLTSize[compBegin]) // escape mode
        {
          refLevel = maxSymbol - 1;
        }
      }
    }
    maxSymbol--;
  }
  symbol = curLevel;
  if (curLevel >= refLevel) // include escape mode
  {
    symbol++;
  }
  paletteIdx.at(posx, posy) = symbol;
}

//================================================================================
//  clause 7.3.8.6
//--------------------------------------------------------------------------------
//    void  prediction_unit ( pu, mrgCtx );
//    void  merge_flag      ( pu );
//    void  merge_data      ( pu, mrgCtx );
//    void  merge_idx       ( pu );
//    void  inter_pred_idc  ( pu );
//    void  ref_idx         ( pu, refList );
//    void  mvp_flag        ( pu, refList );
//================================================================================

void CABACReader::prediction_unit(PredictionUnit &pu, MergeCtx &mrgCtx) {
  if (pu.cu->skip) {
    pu.mergeFlag = true;
  } else {
    merge_flag(pu);
  }
  if (pu.mergeFlag) {
    merge_data(pu);
  } else if (CU::isIBC(*pu.cu)) {
    pu.interDir = 1;
    pu.cu->affine = false;
    pu.refIdx[REF_PIC_LIST_0] = MAX_NUM_REF;
    mvd_coding(pu.mvd[REF_PIC_LIST_0]);
    if (pu.cs->sps->getMaxNumIBCMergeCand() == 1) {
      pu.mvpIdx[REF_PIC_LIST_0] = 0;
    } else {
      mvp_flag(pu, REF_PIC_LIST_0);
    }
  } else {
    inter_pred_idc(pu);
    affine_flag(*pu.cu);
    smvd_mode(pu);

    if (pu.interDir != 2 /* PRED_L1 */) {
      ref_idx(pu, REF_PIC_LIST_0);
      if (pu.cu->affine) {
        mvd_coding(pu.mvdAffi[REF_PIC_LIST_0][0]);
        mvd_coding(pu.mvdAffi[REF_PIC_LIST_0][1]);
        if (pu.cu->affineType == AFFINEMODEL_6PARAM) {
          mvd_coding(pu.mvdAffi[REF_PIC_LIST_0][2]);
        }
      } else {
        mvd_coding(pu.mvd[REF_PIC_LIST_0]);
      }
      mvp_flag(pu, REF_PIC_LIST_0);
    }

    if (pu.interDir != 1 /* PRED_L0 */) {
      if (pu.cu->smvdMode != 1) {
        ref_idx(pu, REF_PIC_LIST_1);
        if (pu.cu->cs->picHeader->getMvdL1ZeroFlag() &&
            pu.interDir == 3 /* PRED_BI */) {
          pu.mvd[REF_PIC_LIST_1] = Mv();
          pu.mvdAffi[REF_PIC_LIST_1][0] = Mv();
          pu.mvdAffi[REF_PIC_LIST_1][1] = Mv();
          pu.mvdAffi[REF_PIC_LIST_1][2] = Mv();
        } else if (pu.cu->affine) {
          mvd_coding(pu.mvdAffi[REF_PIC_LIST_1][0]);
          mvd_coding(pu.mvdAffi[REF_PIC_LIST_1][1]);
          if (pu.cu->affineType == AFFINEMODEL_6PARAM) {
            mvd_coding(pu.mvdAffi[REF_PIC_LIST_1][2]);
          }
        } else {
          mvd_coding(pu.mvd[REF_PIC_LIST_1]);
        }
      }
      mvp_flag(pu, REF_PIC_LIST_1);
    }
  }
  if (pu.interDir == 3 /* PRED_BI */ && PU::isBipredRestriction(pu)) {
    pu.mv[REF_PIC_LIST_1] = Mv(0, 0);
    pu.refIdx[REF_PIC_LIST_1] = -1;
    pu.interDir = 1;
    pu.cu->BcwIdx = BCW_DEFAULT;
  }

  if (pu.cu->smvdMode) {
    RefPicList eCurRefList = (RefPicList)(pu.cu->smvdMode - 1);
    pu.mvd[1 - eCurRefList].set(-pu.mvd[eCurRefList].hor,
                                -pu.mvd[eCurRefList].ver);
    CHECK(!((pu.mvd[1 - eCurRefList].getHor() >= MVD_MIN) &&
            (pu.mvd[1 - eCurRefList].getHor() <= MVD_MAX)) ||
              !((pu.mvd[1 - eCurRefList].getVer() >= MVD_MIN) &&
                (pu.mvd[1 - eCurRefList].getVer() <= MVD_MAX)),
          "Illegal MVD value");
    pu.refIdx[1 - eCurRefList] = pu.cs->slice->getSymRefIdx(1 - eCurRefList);
  }
}

void CABACReader::smvd_mode(PredictionUnit &pu) {
  pu.cu->smvdMode = 0;
  if (pu.interDir != 3 || pu.cu->affine) {
    return;
  }

  if (pu.cs->slice->getBiDirPred() == false) {
    return;
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__SYMMVD_FLAG);

  pu.cu->smvdMode = m_BinDecoder.decodeBin(Ctx::SmvdFlag()) ? 1 : 0;
  binLogger.LogElements(SyntaxElement::sym_mvd_flag, pu.cu->smvdMode);
}

void CABACReader::subblock_merge_flag(CodingUnit &cu) {
  cu.affine = false;

  if (!cu.cs->slice->isIntra() &&
      (cu.slice->getPicHeader()->getMaxNumAffineMergeCand() > 0) &&
      cu.lumaSize().width >= 8 && cu.lumaSize().height >= 8) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
        STATS__CABAC_BITS__AFFINE_FLAG);

    unsigned ctxId = DeriveCtx::CtxAffineFlag(cu);
    cu.affine = m_BinDecoder.decodeBin(Ctx::SubblockMergeFlag(ctxId));
    binLogger.LogElements(SyntaxElement::merge_subblock_flag, cu.affine);
  }
}

void CABACReader::affine_flag(CodingUnit &cu) {
  if (!cu.cs->slice->isIntra() && cu.cs->sps->getUseAffine() &&
      cu.lumaSize().width > 8 && cu.lumaSize().height > 8) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
        STATS__CABAC_BITS__AFFINE_FLAG);

    unsigned ctxId = DeriveCtx::CtxAffineFlag(cu);
    cu.affine = m_BinDecoder.decodeBin(Ctx::AffineFlag(ctxId));
    binLogger.LogElements(SyntaxElement::inter_affine_flag, cu.affine);

    if (cu.affine && cu.cs->sps->getUseAffineType()) {
      ctxId = 0;
      cu.affineType = m_BinDecoder.decodeBin(Ctx::AffineType(ctxId));
      binLogger.LogElements(SyntaxElement::cu_affine_type_flag, cu.affineType);
    } else {
      cu.affineType = AFFINEMODEL_4PARAM;
    }
  }
}

void CABACReader::merge_flag(PredictionUnit &pu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__MERGE_FLAG);

  pu.mergeFlag = (m_BinDecoder.decodeBin(Ctx::MergeFlag()));
  binLogger.LogElements(SyntaxElement::general_merge_flag, pu.mergeFlag);

  if (pu.mergeFlag && CU::isIBC(*pu.cu)) {
    pu.mmvdMergeFlag = false;
    pu.regularMergeFlag = false;
    return;
  }
}

void CABACReader::merge_data(PredictionUnit &pu) {
  if (CU::isIBC(*pu.cu)) {
    merge_idx(pu);
    return;
  } else {
    CodingUnit cu = *pu.cu;
    subblock_merge_flag(*pu.cu);
    if (pu.cu->affine) {
      merge_idx(pu);
      cu.firstPU->regularMergeFlag = false;
      return;
    }

    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
        STATS__CABAC_BITS__MERGE_FLAG);

    const bool ciipAvailable = pu.cs->sps->getUseCiip() && !pu.cu->skip &&
                               pu.cu->lwidth() < MAX_CU_SIZE &&
                               pu.cu->lheight() < MAX_CU_SIZE &&
                               pu.cu->lwidth() * pu.cu->lheight() >= 64;
    const bool geoAvailable = pu.cu->cs->slice->getSPS()->getUseGeo() &&
                              pu.cu->cs->slice->isInterB() &&
                              pu.cs->sps->getMaxNumGeoCand() > 1 &&
                              pu.cu->lwidth() >= GEO_MIN_CU_SIZE &&
                              pu.cu->lheight() >= GEO_MIN_CU_SIZE &&
                              pu.cu->lwidth() <= GEO_MAX_CU_SIZE &&
                              pu.cu->lheight() <= GEO_MAX_CU_SIZE &&
                              pu.cu->lwidth() < 8 * pu.cu->lheight() &&
                              pu.cu->lheight() < 8 * pu.cu->lwidth();
    if (geoAvailable || ciipAvailable) {
      cu.firstPU->regularMergeFlag =
          m_BinDecoder.decodeBin(Ctx::RegularMergeFlag(cu.skip ? 0 : 1));
      binLogger.LogElements(SyntaxElement::regular_merge_flag, cu.firstPU->regularMergeFlag);
    } else {
      cu.firstPU->regularMergeFlag = true;
    }
    if (cu.firstPU->regularMergeFlag) {
      if (cu.cs->slice->getSPS()->getUseMMVD()) {
        cu.firstPU->mmvdMergeFlag = m_BinDecoder.decodeBin(Ctx::MmvdFlag(0));
        binLogger.LogElements(SyntaxElement::mmvd_merge_flag, cu.firstPU->mmvdMergeFlag);
      } else {
        cu.firstPU->mmvdMergeFlag = false;
      }
      if (cu.skip) {
        cu.mmvdSkip = cu.firstPU->mmvdMergeFlag;
      }
    } else {
      pu.mmvdMergeFlag = false;
      pu.cu->mmvdSkip = false;
      if (geoAvailable && ciipAvailable) {
        Ciip_flag(pu);
      } else if (ciipAvailable) {
        pu.ciipFlag = true;
      } else {
        pu.ciipFlag = false;
      }
      if (pu.ciipFlag) {
        pu.intraDir[0] = PLANAR_IDX;
        pu.intraDir[1] = DM_CHROMA_IDX;
      } else {
        pu.cu->geoFlag = true;
      }
    }
  }
  if (pu.mmvdMergeFlag || pu.cu->mmvdSkip) {
    mmvd_merge_idx(pu);
  } else {
    merge_idx(pu);
  }
}

void CABACReader::merge_idx(PredictionUnit &pu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__MERGE_INDEX);

  if (pu.cu->affine) {
    int numCandminus1 = int(pu.cs->picHeader->getMaxNumAffineMergeCand()) - 1;
    pu.mergeIdx = 0;
    if (numCandminus1 > 0) {
      if (m_BinDecoder.decodeBin(Ctx::AffMergeIdx())) {
        binLogger.LogElement(SyntaxElement::merge_idx);
        pu.mergeIdx++;
        for (; pu.mergeIdx < numCandminus1; pu.mergeIdx++) {
          if (!m_BinDecoder.decodeBinEP()) {
            break;
          }
          binLogger.LogElement(SyntaxElement::merge_idx);
        }
      }
    }
  } else {
    int numCandminus1 = int(pu.cs->sps->getMaxNumMergeCand()) - 1;
    pu.mergeIdx = 0;

    if (pu.cu->geoFlag) {
      RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
          STATS__CABAC_BITS__GEO_INDEX);
      uint32_t splitDir = 0;
      xReadTruncBinCode(splitDir, GEO_NUM_PARTITION_MODE);
      binLogger.LogElements(SyntaxElement::merge_idx, splitDir);
      pu.geoSplitDir = splitDir;
      const int maxNumGeoCand = pu.cs->sps->getMaxNumGeoCand();
      CHECK(maxNumGeoCand < 2, "Incorrect max number of geo candidates");
      CHECK(pu.cu->lheight() > 64 || pu.cu->lwidth() > 64,
            "Incorrect block size of geo flag");
      int numCandminus2 = maxNumGeoCand - 2;
      pu.mergeIdx = 0;
      int mergeCand0 = 0;
      int mergeCand1 = 0;
      if (m_BinDecoder.decodeBin(Ctx::MergeIdx())) {
        binLogger.LogElement(SyntaxElement::merge_idx);
        mergeCand0 += unary_max_eqprob(numCandminus2) + 1;
        binLogger.LogElements(SyntaxElement::amvr_precision_idx, mergeCand0);
      }
      if (numCandminus2 > 0) {
        if (m_BinDecoder.decodeBin(Ctx::MergeIdx())) {
          binLogger.LogElement(SyntaxElement::merge_idx);
          mergeCand1 += unary_max_eqprob(numCandminus2 - 1) + 1;
          binLogger.LogElements(SyntaxElement::amvr_precision_idx, mergeCand1);
        }
      }
      mergeCand1 += mergeCand1 >= mergeCand0 ? 1 : 0;
      pu.geoMergeIdx0 = mergeCand0;
      pu.geoMergeIdx1 = mergeCand1;
      return;
    }

    if (pu.cu->predMode == MODE_IBC) {
      numCandminus1 = int(pu.cs->sps->getMaxNumIBCMergeCand()) - 1;
    }
    if (numCandminus1 > 0) {
      if (m_BinDecoder.decodeBin(Ctx::MergeIdx())) {
        binLogger.LogElement(SyntaxElement::merge_idx);
        pu.mergeIdx++;
        for (; pu.mergeIdx < numCandminus1; pu.mergeIdx++) {
          if (!m_BinDecoder.decodeBinEP()) {
            break;
          }
          binLogger.LogElement(SyntaxElement::merge_idx);
        }
      }
    }
  }
}

void CABACReader::mmvd_merge_idx(PredictionUnit &pu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__MERGE_INDEX);

  int var0 = 0;
  if (pu.cs->sps->getMaxNumMergeCand() > 1) {
    static_assert(MMVD_BASE_MV_NUM == 2, "");
    var0 = m_BinDecoder.decodeBin(Ctx::MmvdMergeIdx());
    binLogger.LogElements(SyntaxElement::mmvd_merge_flag, var0);
  }
  int numCandminus1_step = MMVD_REFINE_STEP - 1;
  int var1 = 0;
  if (m_BinDecoder.decodeBin(Ctx::MmvdStepMvpIdx())) {
    binLogger.LogElement(SyntaxElement::mmvd_distance_idx);
    var1++;
    for (; var1 < numCandminus1_step; var1++) {
      if (!m_BinDecoder.decodeBinEP()) {
        break;
      }
      binLogger.LogElement(SyntaxElement::mmvd_distance_idx);
    }
  }
  int var2 = 0;
  if (m_BinDecoder.decodeBinEP()) {
    binLogger.LogElement(SyntaxElement::mmvd_distance_idx);
    var2 += 2;
    if (m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::mmvd_distance_idx);
      var2 += 1;
    }
  } else {
    var2 += 0;
    if (m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::mmvd_distance_idx);
      var2 += 1;
    }
  }
  int mvpIdx = (var0 * MMVD_MAX_REFINE_NUM + var1 * 4 + var2);
  pu.mmvdMergeIdx = mvpIdx;
}

void CABACReader::inter_pred_idc(PredictionUnit &pu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__INTER_DIR);

  if (pu.cs->slice->isInterP()) {
    pu.interDir = 1;
    return;
  }
  if (!(PU::isBipredRestriction(pu))) {
    unsigned ctxId = DeriveCtx::CtxInterDir(pu);
    if (m_BinDecoder.decodeBin(Ctx::InterDir(ctxId))) {
      binLogger.LogElement(SyntaxElement::inter_pred_idc);
      pu.interDir = 3;
      return;
    }
  }
  if (m_BinDecoder.decodeBin(Ctx::InterDir(5))) {
    binLogger.LogElement(SyntaxElement::inter_pred_idc);
    pu.interDir = 2;
    return;
  }
  pu.interDir = 1;
  return;
}

void CABACReader::ref_idx(PredictionUnit &pu, RefPicList eRefList) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__REF_FRM_IDX);

  if (pu.cu->smvdMode) {
    pu.refIdx[eRefList] = pu.cs->slice->getSymRefIdx(eRefList);
    return;
  }

  int numRef = pu.cs->slice->getNumRefIdx(eRefList);

  if (numRef <= 1 || !m_BinDecoder.decodeBin(Ctx::RefPic())) {
    binLogger.LogElement(SyntaxElement::ref_idx_l0);
    pu.refIdx[eRefList] = 0;
    return;
  }
  if (numRef <= 2 || !m_BinDecoder.decodeBin(Ctx::RefPic(1))) {
    binLogger.LogElement(SyntaxElement::ref_idx_l1);
    pu.refIdx[eRefList] = 1;
    return;
  }
  for (int idx = 3;; idx++) {
    if (numRef <= idx || !m_BinDecoder.decodeBinEP()) {
      pu.refIdx[eRefList] = (signed char)(idx - 1);
      return;
    }
    binLogger.LogElement(SyntaxElement::ref_idx_l0);
  }
}

void CABACReader::mvp_flag(PredictionUnit &pu, RefPicList eRefList) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__MVP_IDX);

  unsigned mvp_idx = m_BinDecoder.decodeBin(Ctx::MVPIdx());
  binLogger.LogElements(SyntaxElement::mvp_l0_flag, mvp_idx);
  pu.mvpIdx[eRefList] = mvp_idx;
}

void CABACReader::Ciip_flag(PredictionUnit &pu) {
  if (!pu.cs->sps->getUseCiip()) {
    pu.ciipFlag = false;
    return;
  }
  if (pu.cu->skip) {
    pu.ciipFlag = false;
    return;
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
      STATS__CABAC_BITS__MH_INTRA_FLAG);

  pu.ciipFlag = (m_BinDecoder.decodeBin(Ctx::CiipFlag()));
  binLogger.LogElements(SyntaxElement::ciip_flag, pu.ciipFlag);
}

//================================================================================
//  clause 7.3.8.8
//--------------------------------------------------------------------------------
//    void  transform_tree      ( cs, area, cuCtx, chromaCbfs )
//    bool  split_transform_flag( depth )
//    bool  cbf_comp            ( area, depth )
//================================================================================

void CABACReader::transform_tree(CodingStructure &cs, Partitioner &partitioner,
                                 CUCtx &cuCtx, const PartSplit ispType,
                                 const int subTuIdx) {
  const UnitArea &area = partitioner.currArea();
  CodingUnit &cu =
      *cs.getCU(area.blocks[partitioner.chType], partitioner.chType);
  int subTuCounter = subTuIdx;

  // split_transform_flag
  bool split = partitioner.canSplit(TU_MAX_TR_SPLIT, cs);
  const unsigned trDepth = partitioner.currTrDepth;

  if (cu.sbtInfo && partitioner.canSplit(PartSplit(cu.getSbtTuSplit()), cs)) {
    split = true;
  }

  if (!split && cu.ispMode) {
    split = partitioner.canSplit(ispType, cs);
  }

  if (split) {
    if (partitioner.canSplit(TU_MAX_TR_SPLIT, cs)) {
      partitioner.splitCurrArea(TU_MAX_TR_SPLIT, cs);
    } else if (cu.ispMode) {
      partitioner.splitCurrArea(ispType, cs);
    } else if (cu.sbtInfo &&
               partitioner.canSplit(PartSplit(cu.getSbtTuSplit()), cs)) {
      partitioner.splitCurrArea(PartSplit(cu.getSbtTuSplit()), cs);
    } else {
      THROW("Implicit TU split not available!");
    }

    do {
      transform_tree(cs, partitioner, cuCtx, ispType, subTuCounter);
      subTuCounter += subTuCounter != -1 ? 1 : 0;
    } while (partitioner.nextPart(cs));

    partitioner.exitCurrSplit();
  } else {
    TransformUnit &tu =
        cs.addTU(CS::getArea(cs, area, partitioner.chType), partitioner.chType);
    unsigned numBlocks = getNumberValidTBlocks(*cs.pcv);
    tu.checkTuNoResidual(partitioner.currPartIdx());

    for (unsigned compID = COMPONENT_Y; compID < numBlocks; compID++) {
      if (tu.blocks[compID].valid()) {
        tu.getCoeffs(ComponentID(compID)).fill(0);
        tu.getPcmbuf(ComponentID(compID)).fill(0);
      }
    }
    tu.depth = trDepth;

    transform_unit(tu, cuCtx, partitioner, subTuCounter);
  }
}

bool CABACReader::cbf_comp(CodingStructure &cs, const CompArea &area,
                           unsigned depth, const bool prevCbf,
                           const bool useISP) {
  unsigned ctxId =
      DeriveCtx::CtxQtCbf(area.compID, prevCbf, useISP && isLuma(area.compID));
  const CtxSet &ctxSet = Ctx::QtCbf[area.compID];

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(STATS__CABAC_BITS__QT_CBF,
                                                      area.size(), area.compID);

  unsigned cbf = 0;
  if ((area.compID == COMPONENT_Y &&
       cs.getCU(area.pos(), toChannelType(area.compID))->bdpcmMode) ||
      (area.compID != COMPONENT_Y &&
       cs.getCU(area.pos(), toChannelType(area.compID))->bdpcmModeChroma)) {
    if (area.compID == COMPONENT_Y) {
      ctxId = 1;
    } else if (area.compID == COMPONENT_Cb) {
      ctxId = 1;
    } else {
      ctxId = 2;
    }
    cbf = m_BinDecoder.decodeBin(ctxSet(ctxId));
    binLogger.LogElements(area.compID == COMPONENT_Y
                              ? SyntaxElement::intra_bdpcm_luma_flag
                              : SyntaxElement::intra_bdpcm_chroma_flag,
                          cbf);
  } else {
    cbf = m_BinDecoder.decodeBin(ctxSet(ctxId));
    binLogger.LogElements(area.compID == COMPONENT_Y
                              ? SyntaxElement::intra_bdpcm_luma_flag
                              : SyntaxElement::intra_bdpcm_chroma_flag,
                          cbf);
  }

  return cbf;
}

//================================================================================
//  clause 7.3.8.9
//--------------------------------------------------------------------------------
//    void  mvd_coding( pu, refList )
//================================================================================

void CABACReader::mvd_coding(Mv &rMvd) {
#if RExt__DECODER_DEBUG_BIT_STATISTICS
  CodingStatisticsClassType ctype_mvd(STATS__CABAC_BITS__MVD);
  CodingStatisticsClassType ctype_mvd_ep(STATS__CABAC_BITS__MVD_EP);
#endif

  RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_mvd);

  // abs_mvd_greater0_flag[ 0 | 1 ]
  int horAbs = (int)m_BinDecoder.decodeBin(Ctx::Mvd());
  int verAbs = (int)m_BinDecoder.decodeBin(Ctx::Mvd());
  binLogger.LogElements(SyntaxElement::abs_mvd_greater0_flag, horAbs, verAbs);

  // abs_mvd_greater1_flag[ 0 | 1 ]
  if (horAbs) {
    horAbs += (int)m_BinDecoder.decodeBin(Ctx::Mvd(1));
    binLogger.LogElement(SyntaxElement::abs_mvd_greater1_flag);
  }
  if (verAbs) {
    verAbs += (int)m_BinDecoder.decodeBin(Ctx::Mvd(1));
    binLogger.LogElement(SyntaxElement::abs_mvd_greater1_flag);
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_mvd_ep);

  // abs_mvd_minus2[ 0 | 1 ] and mvd_sign_flag[ 0 | 1 ]
  if (horAbs) {
    if (horAbs > 1) {
      horAbs += m_BinDecoder.decodeRemAbsEP(1, 0, MV_BITS - 1);
      binLogger.LogElement(SyntaxElement::abs_mvd_minus2);
    }
    if (m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::mvd_sign_flag);
      horAbs = -horAbs;
    }
  }
  if (verAbs) {
    if (verAbs > 1) {
      verAbs += m_BinDecoder.decodeRemAbsEP(1, 0, MV_BITS - 1);
      binLogger.LogElement(SyntaxElement::abs_mvd_minus2);
    }
    if (m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::mvd_sign_flag);
      verAbs = -verAbs;
    }
  }
  rMvd = Mv(horAbs, verAbs);
  CHECK(!((horAbs >= MVD_MIN) && (horAbs <= MVD_MAX)) ||
            !((verAbs >= MVD_MIN) && (verAbs <= MVD_MAX)),
        "Illegal MVD value");
}

//================================================================================
//  clause 7.3.8.10
//--------------------------------------------------------------------------------
//    void  transform_unit      ( tu, cuCtx, chromaCbfs )
//    void  cu_qp_delta         ( cu )
//    void  cu_chroma_qp_offset ( cu )
//================================================================================
void CABACReader::transform_unit(TransformUnit &tu, CUCtx &cuCtx,
                                 Partitioner &partitioner,
                                 const int subTuCounter) {
  const UnitArea &area = partitioner.currArea();
  const unsigned trDepth = partitioner.currTrDepth;

  CodingStructure &cs = *tu.cs;
  CodingUnit &cu = *tu.cu;
  ChromaCbfs chromaCbfs;
  chromaCbfs.Cb = chromaCbfs.Cr = false;

  const bool chromaCbfISP = area.chromaFormat != CHROMA_400 &&
                            area.blocks[COMPONENT_Cb].valid() && cu.ispMode;

  // cbf_cb & cbf_cr
  if (area.chromaFormat != CHROMA_400 && area.blocks[COMPONENT_Cb].valid() &&
      (!cu.isSepTree() || partitioner.chType == CHANNEL_TYPE_CHROMA) &&
      (!cu.ispMode || chromaCbfISP)) {
    const int cbfDepth = chromaCbfISP ? trDepth - 1 : trDepth;
    if (!(cu.sbtInfo && tu.noResidual)) {
      chromaCbfs.Cb = cbf_comp(cs, area.blocks[COMPONENT_Cb], cbfDepth);
    }

    if (!(cu.sbtInfo && tu.noResidual)) {
      chromaCbfs.Cr =
          cbf_comp(cs, area.blocks[COMPONENT_Cr], cbfDepth, chromaCbfs.Cb);
    }
  } else if (cu.isSepTree()) {
    chromaCbfs = ChromaCbfs(false);
  }

  if (!isChroma(partitioner.chType)) {
    if (!CU::isIntra(cu) && trDepth == 0 &&
        !chromaCbfs.sigChroma(area.chromaFormat)) {
      TU::setCbfAtDepth(tu, COMPONENT_Y, trDepth, 1);
    } else if (cu.sbtInfo && tu.noResidual) {
      TU::setCbfAtDepth(tu, COMPONENT_Y, trDepth, 0);
    } else if (cu.sbtInfo && !chromaCbfs.sigChroma(area.chromaFormat)) {
      assert(!tu.noResidual);
      TU::setCbfAtDepth(tu, COMPONENT_Y, trDepth, 1);
    } else {
      bool lumaCbfIsInferredACT =
          (cu.colorTransform && cu.predMode == MODE_INTRA && trDepth == 0 &&
           !chromaCbfs.sigChroma(area.chromaFormat));
      bool lastCbfIsInferred =
          lumaCbfIsInferredACT; // ISP and ACT are mutually exclusive
      bool previousCbf = false;
      bool rootCbfSoFar = false;
      if (cu.ispMode) {
        uint32_t nTus = cu.ispMode == HOR_INTRA_SUBPARTITIONS
                            ? cu.lheight() >> floorLog2(tu.lheight())
                            : cu.lwidth() >> floorLog2(tu.lwidth());
        if (subTuCounter == nTus - 1) {
          TransformUnit *tuPointer = cu.firstTU;
          for (int tuIdx = 0; tuIdx < nTus - 1; tuIdx++) {
            rootCbfSoFar |= TU::getCbfAtDepth(*tuPointer, COMPONENT_Y, trDepth);
            tuPointer = tuPointer->next;
          }
          if (!rootCbfSoFar) {
            lastCbfIsInferred = true;
          }
        }
        if (!lastCbfIsInferred) {
          previousCbf = TU::getPrevTuCbfAtDepth(tu, COMPONENT_Y, trDepth);
        }
      }
      bool cbfY = lastCbfIsInferred
                      ? true
                      : cbf_comp(cs, tu.Y(), trDepth, previousCbf, cu.ispMode);
      TU::setCbfAtDepth(tu, COMPONENT_Y, trDepth, (cbfY ? 1 : 0));
    }
  }
  if (area.chromaFormat != CHROMA_400 && (!cu.ispMode || chromaCbfISP)) {
    TU::setCbfAtDepth(tu, COMPONENT_Cb, trDepth, (chromaCbfs.Cb ? 1 : 0));
    TU::setCbfAtDepth(tu, COMPONENT_Cr, trDepth, (chromaCbfs.Cr ? 1 : 0));
  }
  bool lumaOnly =
      (cu.chromaFormat == CHROMA_400 || !tu.blocks[COMPONENT_Cb].valid());
  bool cbfLuma = (tu.cbf[COMPONENT_Y] != 0);
  bool cbfChroma = (lumaOnly ? false : (chromaCbfs.Cb || chromaCbfs.Cr));

  if ((cu.lwidth() > 64 || cu.lheight() > 64 || cbfLuma || cbfChroma) &&
      (!tu.cu->isSepTree() || isLuma(tu.chType))) {
    if (cu.cs->pps->getUseDQP() && !cuCtx.isDQPCoded) {
      cu_qp_delta(cu, cuCtx.qp, cu.qp);
      cuCtx.qp = cu.qp;
      cuCtx.isDQPCoded = true;
    }
  }
  if (!cu.isSepTree() || isChroma(tu.chType)) // !DUAL_TREE_LUMA
  {
    SizeType channelWidth =
        !cu.isSepTree() ? cu.lwidth() : cu.chromaSize().width;
    SizeType channelHeight =
        !cu.isSepTree() ? cu.lheight() : cu.chromaSize().height;

    if (cu.cs->slice->getUseChromaQpAdj() &&
        (channelWidth > 64 || channelHeight > 64 || cbfChroma) &&
        !cuCtx.isChromaQpAdjCoded) {
      cu_chroma_qp_offset(cu);
      cuCtx.isChromaQpAdjCoded = true;
    }
  }

  if (!lumaOnly) {
    joint_cb_cr(tu, (tu.cbf[COMPONENT_Cb] ? 2 : 0) +
                        (tu.cbf[COMPONENT_Cr] ? 1 : 0));
  }

  if (cbfLuma) {
    residual_coding(tu, COMPONENT_Y, cuCtx);
  }
  if (!lumaOnly) {
    for (ComponentID compID = COMPONENT_Cb; compID <= COMPONENT_Cr;
         compID = ComponentID(compID + 1)) {
      if (tu.cbf[compID]) {
        residual_coding(tu, compID, cuCtx);
      }
    }
  }
}

void CABACReader::cu_qp_delta(CodingUnit &cu, int predQP, int8_t &qp) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__DELTA_QP_EP);

  CHECK(predQP == std::numeric_limits<int>::max(), "Invalid predicted QP");
  int qpY = predQP;
  int DQp = unary_max_symbol(Ctx::DeltaQP(), Ctx::DeltaQP(1), CU_DQP_TU_CMAX);
  binLogger.LogElements(SyntaxElement::cu_qp_delta_abs, DQp);
  if (DQp >= CU_DQP_TU_CMAX) {
    DQp += exp_golomb_eqprob(CU_DQP_EG_k);
    binLogger.LogElement(SyntaxElement::cu_qp_delta_abs);
  }
  if (DQp > 0) {
    if (m_BinDecoder.decodeBinEP()) {
      binLogger.LogElement(SyntaxElement::cu_qp_delta_sign_flag);
      DQp = -DQp;
    }
    int qpBdOffsetY = cu.cs->sps->getQpBDOffset(CHANNEL_TYPE_LUMA);
    qpY = ((predQP + DQp + (MAX_QP + 1) + 2 * qpBdOffsetY) %
           ((MAX_QP + 1) + qpBdOffsetY)) -
          qpBdOffsetY;
  }
  qp = (int8_t)qpY;
}

void CABACReader::cu_chroma_qp_offset(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__CHROMA_QP_ADJUSTMENT, cu.blocks[cu.chType].lumaSize(),
      CHANNEL_TYPE_CHROMA);

  // cu_chroma_qp_offset_flag
  int length = cu.cs->pps->getChromaQpOffsetListLen();
  unsigned qpAdj = m_BinDecoder.decodeBin(Ctx::ChromaQpAdjFlag());
  binLogger.LogElements(SyntaxElement::cu_chroma_qp_offset_flag, qpAdj);
  if (qpAdj && length > 1) {
    // cu_chroma_qp_offset_idx
    qpAdj += unary_max_symbol(Ctx::ChromaQpAdjIdc(), Ctx::ChromaQpAdjIdc(),
                              length - 1);
    binLogger.LogElement(SyntaxElement::cu_chroma_qp_offset_idx);
  }
  /* NB, symbol = 0 if outer flag is not set,
   *              1 if outer flag is set and there is no inner flag
   *              1+ otherwise */
  cu.chromaQpAdj = cu.cs->chromaQpAdj = qpAdj;
}

//================================================================================
//  clause 7.3.8.11
//--------------------------------------------------------------------------------
//    void        residual_coding         ( tu, compID )
//    bool        transform_skip_flag     ( tu, compID )
//    int         last_sig_coeff          ( coeffCtx )
//    void        residual_coding_subblock( coeffCtx )
//================================================================================

void CABACReader::joint_cb_cr(TransformUnit &tu, const int cbfMask) {
  if (!tu.cu->slice->getSPS()->getJointCbCrEnabledFlag()) {
    return;
  }

  if ((CU::isIntra(*tu.cu) && cbfMask) || (cbfMask == 3)) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
        STATS__CABAC_BITS__JOINT_CB_CR, tu.blocks[COMPONENT_Cr].lumaSize(),
        CHANNEL_TYPE_CHROMA);
    tu.jointCbCr =
        (m_BinDecoder.decodeBin(Ctx::JointCbCrFlag(cbfMask - 1)) ? cbfMask : 0);
    binLogger.LogElements(SyntaxElement::tu_joint_cbcr_residual_flag, tu.jointCbCr);
  }
}

void CABACReader::residual_coding(TransformUnit &tu, ComponentID compID,
                                  CUCtx &cuCtx) {
  const CodingUnit &cu = *tu.cu;

  if (compID == COMPONENT_Cr && tu.jointCbCr == 3) {
    return;
  }

  ts_flag(tu, compID);

  if (tu.mtsIdx[compID] == MTS_SKIP &&
      !tu.cs->slice->getTSResidualCodingDisabledFlag()) {
    residual_codingTS(tu, compID);
    return;
  }

  // determine sign hiding
  bool signHiding = cu.cs->slice->getSignDataHidingEnabledFlag();

  // init coeff coding context
  CoeffCodingContext cctx(tu, compID, signHiding);
  TCoeff *coeff = tu.getCoeffs(compID).buf;

  // parse last coeff position
  cctx.setScanPosLast(last_sig_coeff(cctx, tu, compID));
  if (tu.mtsIdx[compID] != MTS_SKIP && tu.blocks[compID].height >= 4 &&
      tu.blocks[compID].width >= 4) {
    const int maxLfnstPos =
        ((tu.blocks[compID].height == 4 && tu.blocks[compID].width == 4) ||
         (tu.blocks[compID].height == 8 && tu.blocks[compID].width == 8))
            ? 7
            : 15;
    cuCtx.violatesLfnstConstrained[toChannelType(compID)] |=
        cctx.scanPosLast() > maxLfnstPos;
  }
  if (tu.mtsIdx[compID] != MTS_SKIP && tu.blocks[compID].height >= 4 &&
      tu.blocks[compID].width >= 4) {
    const int lfnstLastScanPosTh =
        isLuma(compID) ? LFNST_LAST_SIG_LUMA : LFNST_LAST_SIG_CHROMA;
    cuCtx.lfnstLastScanPos |= cctx.scanPosLast() >= lfnstLastScanPosTh;
  }
  if (isLuma(compID) && tu.mtsIdx[compID] != MTS_SKIP) {
    cuCtx.mtsLastScanPos |= cctx.scanPosLast() >= 1;
  }

  // parse subblocks
  const int stateTransTab =
      (tu.cs->slice->getDepQuantEnabledFlag() ? 32040 : 0);
  int state = 0;

  int ctxBinSampleRatio = (compID == COMPONENT_Y)
                              ? MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA
                              : MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA;
  cctx.regBinLimit =
      (tu.getTbAreaAfterCoefZeroOut(compID) * ctxBinSampleRatio) >> 4;

  // int baseLevel = m_BinDecoder.getCtx().getBaseLevel();
  // cctx.setBaseLevel(baseLevel);
  if (tu.cs->slice->getSPS()
          ->getSpsRangeExtension()
          .getPersistentRiceAdaptationEnabledFlag()) {
    cctx.setUpdateHist(1);
    unsigned riceStats =
        m_BinDecoder.getCtx().getGRAdaptStats((unsigned)compID);
    TCoeff historyValue = (TCoeff)1 << riceStats;
    cctx.setHistValue(historyValue);
  }
  for (int subSetId = (cctx.scanPosLast() >> cctx.log2CGSize()); subSetId >= 0;
       subSetId--) {
    cctx.initSubblock(subSetId);

    if (tu.cs->sps->getUseMTS() && tu.cu->sbtInfo != 0 &&
        tu.blocks[compID].height <= 32 && tu.blocks[compID].width <= 32 &&
        compID == COMPONENT_Y) {
      if ((tu.blocks[compID].height == 32 &&
           cctx.cgPosY() >= (16 >> cctx.log2CGHeight())) ||
          (tu.blocks[compID].width == 32 &&
           cctx.cgPosX() >= (16 >> cctx.log2CGWidth()))) {
        continue;
      }
    }
    residual_coding_subblock(cctx, coeff, stateTransTab, state);

    if (isLuma(compID) && cctx.isSigGroup() &&
        (cctx.cgPosY() > 3 || cctx.cgPosX() > 3)) {
      cuCtx.violatesMtsCoeffConstraint = true;
    }
  }
}

void CABACReader::ts_flag(TransformUnit &tu, ComponentID compID) {
  int tsFlag = ((tu.cu->bdpcmMode && isLuma(compID)) ||
                (tu.cu->bdpcmModeChroma && isChroma(compID)))
                   ? 1
                   : tu.mtsIdx[compID] == MTS_SKIP ? 1 : 0;
  int ctxIdx = isLuma(compID) ? 0 : 1;

  if (TU::isTSAllowed(tu, compID)) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
        STATS__CABAC_BITS__MTS_FLAGS, tu.blocks[compID], compID);
    tsFlag = m_BinDecoder.decodeBin(Ctx::TransformSkipFlag(ctxIdx));
    binLogger.LogElements(SyntaxElement::transform_skip_flag, tsFlag);
  }

  tu.mtsIdx[compID] = tsFlag ? MTS_SKIP : MTS_DCT2_DCT2;
}

void CABACReader::mts_idx(CodingUnit &cu, CUCtx &cuCtx) {
  TransformUnit &tu = *cu.firstTU;
  int mtsIdx =
      tu.mtsIdx[COMPONENT_Y]; // Transform skip flag has already been decoded

  if (CU::isMTSAllowed(cu, COMPONENT_Y) && !cuCtx.violatesMtsCoeffConstraint &&
      cuCtx.mtsLastScanPos && cu.lfnstIdx == 0 && mtsIdx != MTS_SKIP) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
        STATS__CABAC_BITS__MTS_FLAGS, tu.blocks[COMPONENT_Y], COMPONENT_Y);
    int ctxIdx = 0;
    int symbol = m_BinDecoder.decodeBin(Ctx::MTSIdx(ctxIdx));
    binLogger.LogElements(SyntaxElement::mts_idx, symbol);

    if (symbol) {
      ctxIdx = 1;
      mtsIdx = MTS_DST7_DST7; // mtsIdx = 2 -- 4
      for (int i = 0; i < 3; i++, ctxIdx++) {
        symbol = m_BinDecoder.decodeBin(Ctx::MTSIdx(ctxIdx));
        binLogger.LogElements(SyntaxElement::mts_idx, symbol);
        mtsIdx += symbol;

        if (!symbol) {
          break;
        }
      }
    }
  }

  tu.mtsIdx[COMPONENT_Y] = mtsIdx;
}

void CABACReader::isp_mode(CodingUnit &cu) {
  if (!CU::isIntra(cu) || !isLuma(cu.chType) || cu.firstPU->multiRefIdx ||
      !cu.cs->sps->getUseISP() || cu.bdpcmMode ||
      !CU::canUseISP(cu, getFirstComponentOfChannel(cu.chType)) ||
      cu.colorTransform) {
    cu.ispMode = NOT_INTRA_SUBPARTITIONS;
    return;
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
      STATS__CABAC_BITS__ISP_MODE_FLAG);

  int symbol = m_BinDecoder.decodeBin(Ctx::ISPMode(0));
  binLogger.LogElements(SyntaxElement::intra_subpartitions_mode_flag, symbol);

  if (symbol) {
    RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(
        STATS__CABAC_BITS__ISP_SPLIT_FLAG);
    cu.ispMode = 1 + m_BinDecoder.decodeBin(Ctx::ISPMode(1));
    binLogger.LogElements(SyntaxElement::intra_subpartitions_mode_flag, cu.ispMode);
  }
}

void CABACReader::residual_lfnst_mode(CodingUnit &cu, CUCtx &cuCtx) {
  int chIdx = cu.isSepTree() && cu.chType == CHANNEL_TYPE_CHROMA ? 1 : 0;
  if ((cu.ispMode && !CU::canUseLfnstWithISP(cu, cu.chType)) ||
      (cu.cs->sps->getUseLFNST() && CU::isIntra(cu) && cu.mipFlag &&
       !allowLfnstWithMip(cu.firstPU->lumaSize())) ||
      (cu.isSepTree() && cu.chType == CHANNEL_TYPE_CHROMA &&
       std::min(cu.blocks[1].width, cu.blocks[1].height) < 4) ||
      (cu.blocks[chIdx].lumaSize().width > cu.cs->sps->getMaxTbSize() ||
       cu.blocks[chIdx].lumaSize().height > cu.cs->sps->getMaxTbSize())) {
    return;
  }

  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__LFNST);

  if (cu.cs->sps->getUseLFNST() && CU::isIntra(cu)) {
    const bool lumaFlag =
        cu.isSepTree() ? (isLuma(cu.chType) ? true : false) : true;
    const bool chromaFlag =
        cu.isSepTree() ? (isChroma(cu.chType) ? true : false) : true;
    bool nonZeroCoeffNonTsCorner8x8 =
        (lumaFlag && cuCtx.violatesLfnstConstrained[CHANNEL_TYPE_LUMA]) ||
        (chromaFlag && cuCtx.violatesLfnstConstrained[CHANNEL_TYPE_CHROMA]);
    bool isTrSkip = false;
    for (auto &currTU : CU::traverseTUs(cu)) {
      const uint32_t numValidComp = getNumberValidComponents(cu.chromaFormat);
      for (uint32_t compID = COMPONENT_Y; compID < numValidComp; compID++) {
        if (currTU.blocks[compID].valid() &&
            TU::getCbf(currTU, (ComponentID)compID) &&
            currTU.mtsIdx[compID] == MTS_SKIP) {
          isTrSkip = true;
          break;
        }
      }
    }
    if ((!cuCtx.lfnstLastScanPos && !cu.ispMode) ||
        nonZeroCoeffNonTsCorner8x8 || isTrSkip) {
      cu.lfnstIdx = 0;
      return;
    }
  } else {
    cu.lfnstIdx = 0;
    return;
  }

  unsigned cctx = 0;
  if (cu.isSepTree())
    cctx++;

  uint32_t idxLFNST = m_BinDecoder.decodeBin(Ctx::LFNSTIdx(cctx));
  binLogger.LogElements(SyntaxElement::lfnst_idx, idxLFNST);
  if (idxLFNST) {
    idxLFNST += m_BinDecoder.decodeBin(Ctx::LFNSTIdx(2));
    binLogger.LogElement(SyntaxElement::lfnst_idx);
  }
  cu.lfnstIdx = idxLFNST;
}

int CABACReader::last_sig_coeff(CoeffCodingContext &cctx, TransformUnit &tu,
                                ComponentID compID) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__LAST_SIG_X_Y, Size(cctx.width(), cctx.height()),
      cctx.compID());

  unsigned PosLastX = 0, PosLastY = 0;
  unsigned maxLastPosX = cctx.maxLastPosX();
  unsigned maxLastPosY = cctx.maxLastPosY();
#if JVET_W0046_RLSCP
  unsigned zoTbWdith = std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, cctx.width());
  unsigned zoTbHeight =
      std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, cctx.height());
#endif

  if (tu.cs->sps->getUseMTS() && tu.cu->sbtInfo != 0 &&
      tu.blocks[compID].width <= 32 && tu.blocks[compID].height <= 32 &&
      compID == COMPONENT_Y) {
    maxLastPosX =
        (tu.blocks[compID].width == 32) ? g_groupIdx[15] : maxLastPosX;
    maxLastPosY =
        (tu.blocks[compID].height == 32) ? g_groupIdx[15] : maxLastPosY;
#if JVET_W0046_RLSCP
    zoTbWdith = (tu.blocks[compID].width == 32) ? 16 : zoTbWdith;
    zoTbHeight = (tu.blocks[compID].height == 32) ? 16 : zoTbHeight;
#endif
  }

  for (; PosLastX < maxLastPosX; PosLastX++) {
    if (!m_BinDecoder.decodeBin(cctx.lastXCtxId(PosLastX))) {
      break;
    }
    binLogger.LogElement(SyntaxElement::last_sig_coeff_x_prefix);
  }
  for (; PosLastY < maxLastPosY; PosLastY++) {
    if (!m_BinDecoder.decodeBin(cctx.lastYCtxId(PosLastY))) {
      break;
    }
    binLogger.LogElement(SyntaxElement::last_sig_coeff_y_prefix);
  }
  if (PosLastX > 3) {
    uint32_t temp = 0;
    uint32_t uiCount = (PosLastX - 2) >> 1;
    for (int i = uiCount - 1; i >= 0; i--) {
      temp += m_BinDecoder.decodeBinEP() << i;
      binLogger.LogElement(SyntaxElement::last_sig_coeff_x_suffix);
    }
    PosLastX = g_minInGroup[PosLastX] + temp;
  }
  if (PosLastY > 3) {
    uint32_t temp = 0;
    uint32_t uiCount = (PosLastY - 2) >> 1;
    for (int i = uiCount - 1; i >= 0; i--) {
      temp += m_BinDecoder.decodeBinEP() << i;
      binLogger.LogElement(SyntaxElement::last_sig_coeff_y_suffix);
    }
    PosLastY = g_minInGroup[PosLastY] + temp;
  }

// #if JVET_W0046_RLSCP
//   if (tu.cu->slice->getReverseLastSigCoeffFlag()) {
//     PosLastX = zoTbWdith - 1 - PosLastX;
//     PosLastY = zoTbHeight - 1 - PosLastY;
//   }
// #endif
  int blkPos;
  blkPos = PosLastX + (PosLastY * cctx.width());

  int scanPos = 0;
  for (; scanPos < cctx.maxNumCoeff() - 1; scanPos++) {
    if (blkPos == cctx.blockPos(scanPos)) {
      break;
    }
  }
  return scanPos;
}

static void check_coeff_conformance(const CoeffCodingContext &cctx,
                                    const TCoeff coeff) {
  CHECK(coeff < cctx.minCoeff() || coeff > cctx.maxCoeff(),
        "TransCoeffLevel outside allowable range");
}

void CABACReader::residual_coding_subblock(CoeffCodingContext &cctx,
                                           TCoeff *coeff,
                                           const int stateTransTable,
                                           int &state) {
  // NOTE: All coefficients of the subblock must be set to zero before calling
  // this function
#if RExt__DECODER_DEBUG_BIT_STATISTICS
  CodingStatisticsClassType ctype_group(STATS__CABAC_BITS__SIG_COEFF_GROUP_FLAG,
                                        cctx.width(), cctx.height(),
                                        cctx.compID());
  CodingStatisticsClassType ctype_map(STATS__CABAC_BITS__SIG_COEFF_MAP_FLAG,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_par(STATS__CABAC_BITS__PAR_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_gt1(STATS__CABAC_BITS__GT1_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_gt2(STATS__CABAC_BITS__GT2_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_escs(STATS__CABAC_BITS__ESCAPE_BITS,
                                       cctx.width(), cctx.height(),
                                       cctx.compID());
#endif

  //===== init =====
  const int minSubPos = cctx.minSubPos();
  const bool isLast = cctx.isLast();
  int firstSigPos = (isLast ? cctx.scanPosLast() : cctx.maxSubPos());
  int nextSigPos = firstSigPos;
  int baseLevel = cctx.getBaseLevel();
  bool updateHistory = cctx.getUpdateHist();

  //===== decode significant_coeffgroup_flag =====
  RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_group);
  bool sigGroup = (isLast || !minSubPos);
  if (!sigGroup) {
    sigGroup = m_BinDecoder.decodeBin(cctx.sigGroupCtxId());
    binLogger.LogElements(SyntaxElement::sig_coeff_flag, sigGroup);
  }
  if (sigGroup) {
    cctx.setSigGroup();
  } else {
    return;
  }

  uint8_t ctxOffset[16];

  //===== decode absolute values =====
  const int inferSigPos = nextSigPos != cctx.scanPosLast()
                              ? (cctx.isNotFirst() ? minSubPos : -1)
                              : nextSigPos;
  int firstNZPos = nextSigPos;
  int lastNZPos = -1;
  int numNonZero = 0;
  int remRegBins = cctx.regBinLimit;
  int firstPosMode2 = minSubPos - 1;
  int sigBlkPos[1 << MLS_CG_SIZE];

  for (; nextSigPos >= minSubPos && remRegBins >= 4; nextSigPos--) {
    int blkPos = cctx.blockPos(nextSigPos);
    unsigned sigFlag = (!numNonZero && nextSigPos == inferSigPos);
    if (!sigFlag) {
      RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_map);
      const unsigned sigCtxId = cctx.sigCtxIdAbs(nextSigPos, coeff, state);
      sigFlag = m_BinDecoder.decodeBin(sigCtxId);
      binLogger.LogElements(SyntaxElement::sig_coeff_flag, sigFlag);
      remRegBins--;
    } else if (nextSigPos != cctx.scanPosLast()) {
      cctx.sigCtxIdAbs(nextSigPos, coeff,
                       state); // required for setting variables that are needed
                               // for gtx/par context selection
    }

    if (sigFlag) {
      uint8_t &ctxOff = ctxOffset[nextSigPos - minSubPos];
      ctxOff = cctx.ctxOffsetAbs();
      sigBlkPos[numNonZero++] = blkPos;
      firstNZPos = nextSigPos;
      lastNZPos = std::max<int>(lastNZPos, nextSigPos);

      RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_gt1);
      unsigned gt1Flag = m_BinDecoder.decodeBin(cctx.greater1CtxIdAbs(ctxOff));
      binLogger.LogElements(SyntaxElement::abs_mvd_greater0_flag, gt1Flag);
      remRegBins--;

      unsigned parFlag = 0;
      unsigned gt2Flag = 0;
      if (gt1Flag) {
        RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_par);
        parFlag = m_BinDecoder.decodeBin(cctx.parityCtxIdAbs(ctxOff));
        binLogger.LogElements(SyntaxElement::par_level_flag, parFlag);

        remRegBins--;
        RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_gt2);
        gt2Flag = m_BinDecoder.decodeBin(cctx.greater2CtxIdAbs(ctxOff));
        binLogger.LogElements(SyntaxElement::abs_mvd_greater1_flag, gt2Flag);
        remRegBins--;
      }
      coeff[blkPos] += 1 + parFlag + gt1Flag + (gt2Flag << 1);
    }

    state =
        (stateTransTable >> ((state << 2) + ((coeff[blkPos] & 1) << 1))) & 3;
  }
  firstPosMode2 = nextSigPos;
  cctx.regBinLimit = remRegBins;

  //===== 2nd PASS: Go-rice codes =====
  unsigned ricePar = 0;
  for (int scanPos = firstSigPos; scanPos > firstPosMode2; scanPos--) {
    ricePar = (cctx.*(cctx.deriveRiceRRC))(scanPos, coeff, baseLevel);

    TCoeff &tcoeff = coeff[cctx.blockPos(scanPos)];
    if (tcoeff >= 4) {
      RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_escs);
      int rem = m_BinDecoder.decodeRemAbsEP(ricePar, COEF_REMAIN_BIN_REDUCTION,
                                            cctx.maxLog2TrDRange());
      binLogger.LogElements(SyntaxElement::abs_remainder, rem);
      tcoeff += (rem << 1);
      if ((updateHistory) && (rem > 0)) {
        unsigned &riceStats =
            m_BinDecoder.getCtx().getGRAdaptStats((unsigned)(cctx.compID()));
        cctx.updateRiceStat(riceStats, rem, 1);
        cctx.setUpdateHist(0);
        updateHistory = 0;
      }
    }
  }

  //===== coeff bypass ====
  for (int scanPos = firstPosMode2; scanPos >= minSubPos; scanPos--) {
    int rice = (cctx.*(cctx.deriveRiceRRC))(scanPos, coeff, 0);
    int pos0 = g_goRicePosCoeff0(state, rice);
    RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_escs);
    int rem = m_BinDecoder.decodeRemAbsEP(rice, COEF_REMAIN_BIN_REDUCTION,
                                          cctx.maxLog2TrDRange());
    binLogger.LogElements(SyntaxElement::abs_remainder, rem);
    TCoeff tcoeff = (rem == pos0 ? 0 : rem < pos0 ? rem + 1 : rem);
    state = (stateTransTable >> ((state << 2) + ((tcoeff & 1) << 1))) & 3;
    if ((updateHistory) && (rem > 0)) {
      unsigned &riceStats =
          m_BinDecoder.getCtx().getGRAdaptStats((unsigned)(cctx.compID()));
      cctx.updateRiceStat(riceStats, rem, 0);
      cctx.setUpdateHist(0);
      updateHistory = 0;
    }
    if (tcoeff) {
      int blkPos = cctx.blockPos(scanPos);
      sigBlkPos[numNonZero++] = blkPos;
      firstNZPos = scanPos;
      lastNZPos = std::max<int>(lastNZPos, scanPos);
      coeff[blkPos] = tcoeff;
    }
  }

  //===== decode sign's =====
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
      STATS__CABAC_BITS__SIGN_BIT, Size(cctx.width(), cctx.height()),
      cctx.compID());
  const unsigned numSigns =
      (cctx.hideSign(firstNZPos, lastNZPos) ? numNonZero - 1 : numNonZero);
  unsigned signPattern = m_BinDecoder.decodeBinsEP(numSigns) << (32 - numSigns);
  binLogger.LogElements(SyntaxElement::num_signalled_palette_entries, signPattern);

  //===== set final coefficents =====
  TCoeff sumAbs = 0;
  for (unsigned k = 0; k < numSigns; k++) {
    TCoeff AbsCoeff = coeff[sigBlkPos[k]];
    sumAbs += AbsCoeff;
    coeff[sigBlkPos[k]] = (signPattern & (1u << 31) ? -AbsCoeff : AbsCoeff);
    signPattern <<= 1;
    check_coeff_conformance(cctx, coeff[sigBlkPos[k]]);
  }
  if (numNonZero > numSigns) {
    int k = numSigns;
    TCoeff AbsCoeff = coeff[sigBlkPos[k]];
    sumAbs += AbsCoeff;
    coeff[sigBlkPos[k]] = (sumAbs & 1 ? -AbsCoeff : AbsCoeff);
    check_coeff_conformance(cctx, coeff[sigBlkPos[k]]);
  }
}

void CABACReader::residual_codingTS(TransformUnit &tu, ComponentID compID) {
  // init coeff coding context
  CoeffCodingContext cctx(tu, compID, false,
                          isLuma(compID) ? tu.cu->bdpcmMode
                                         : tu.cu->bdpcmModeChroma);
  TCoeff *coeff = tu.getCoeffs(compID).buf;
  int maxCtxBins = (cctx.maxNumCoeff() * 7) >> 2;
  cctx.setNumCtxBins(maxCtxBins);

  for (int subSetId = 0;
       subSetId <= (cctx.maxNumCoeff() - 1) >> cctx.log2CGSize(); subSetId++) {
    cctx.initSubblock(subSetId);
    int goRiceParam = 1;
    if (tu.cu->slice->getSPS()
            ->getSpsRangeExtension()
            .getTSRCRicePresentFlag() &&
        tu.mtsIdx[compID] == MTS_SKIP) {
      goRiceParam = goRiceParam + tu.cu->slice->get_tsrc_index();
    }
    residual_coding_subblockTS(cctx, coeff, goRiceParam);
  }
}

void CABACReader::residual_coding_subblockTS(CoeffCodingContext &cctx,
                                             TCoeff *coeff, int riceParam) {
  // NOTE: All coefficients of the subblock must be set to zero before calling
  // this function
#if RExt__DECODER_DEBUG_BIT_STATISTICS
  CodingStatisticsClassType ctype_group(STATS__CABAC_BITS__SIG_COEFF_GROUP_FLAG,
                                        cctx.width(), cctx.height(),
                                        cctx.compID());
#if TR_ONLY_COEFF_STATS
  CodingStatisticsClassType ctype_map(STATS__CABAC_BITS__SIG_COEFF_MAP_FLAG_TS,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_par(STATS__CABAC_BITS__PAR_FLAG_TS,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_gt1(STATS__CABAC_BITS__GT1_FLAG_TS,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_gt2(STATS__CABAC_BITS__GT2_FLAG_TS,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_escs(STATS__CABAC_BITS__ESCAPE_BITS_TS,
                                       cctx.width(), cctx.height(),
                                       cctx.compID());
#else
  CodingStatisticsClassType ctype_map(STATS__CABAC_BITS__SIG_COEFF_MAP_FLAG,
                                      cctx.width(), cctx.height(),
                                      cctx.compID());
  CodingStatisticsClassType ctype_par(STATS__CABAC_BITS__PAR_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_gt1(STATS__CABAC_BITS__GT1_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_gt2(STATS__CABAC_BITS__GT2_FLAG, cctx.width(),
                                      cctx.height(), cctx.compID());
  CodingStatisticsClassType ctype_escs(STATS__CABAC_BITS__ESCAPE_BITS,
                                       cctx.width(), cctx.height(),
                                       cctx.compID());
#endif

#endif

  //===== init =====
  const int minSubPos = cctx.maxSubPos();
  int firstSigPos = cctx.minSubPos();
  int nextSigPos = firstSigPos;
  unsigned signPattern = 0;

  //===== decode significant_coeffgroup_flag =====
  RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_group);
  bool sigGroup = cctx.isLastSubSet() && cctx.noneSigGroup();
  if (!sigGroup) {
    sigGroup = m_BinDecoder.decodeBin(cctx.sigGroupCtxId(true));
    binLogger.LogElements(SyntaxElement::sig_coeff_flag, sigGroup);
  }
  if (sigGroup) {
    cctx.setSigGroup();
  } else {
    return;
  }

  //===== decode absolute values =====
  const int inferSigPos = minSubPos;
  int numNonZero = 0;
  int sigBlkPos[1 << MLS_CG_SIZE];

  int lastScanPosPass1 = -1;
  int lastScanPosPass2 = -1;
  for (; nextSigPos <= minSubPos && cctx.numCtxBins() >= 4; nextSigPos++) {
    int blkPos = cctx.blockPos(nextSigPos);
    unsigned sigFlag = (!numNonZero && nextSigPos == inferSigPos);
    if (!sigFlag) {
      RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_map);
      const unsigned sigCtxId = cctx.sigCtxIdAbsTS(nextSigPos, coeff);
      sigFlag = m_BinDecoder.decodeBin(sigCtxId);
      binLogger.LogElements(SyntaxElement::sig_coeff_flag, sigFlag);
      cctx.decimateNumCtxBins(1);
    }

    if (sigFlag) {
      //===== decode sign's =====
#if TR_ONLY_COEFF_STATS
      RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
          STATS__CABAC_BITS__SIGN_BIT_TS, Size(cctx.width(), cctx.height()),
          cctx.compID());
#else
      RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET_SIZE2(
          STATS__CABAC_BITS__SIGN_BIT, Size(cctx.width(), cctx.height()),
          cctx.compID());
#endif
      int sign;
      const unsigned signCtxId =
          cctx.signCtxIdAbsTS(nextSigPos, coeff, cctx.bdpcm());
      sign = m_BinDecoder.decodeBin(signCtxId);
      binLogger.LogElements(SyntaxElement::sig_coeff_flag, sign);
      cctx.decimateNumCtxBins(1);

      signPattern += (sign << numNonZero);

      sigBlkPos[numNonZero++] = blkPos;

      RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_gt1);
      unsigned gt1Flag;
      const unsigned gt1CtxId =
          cctx.lrg1CtxIdAbsTS(nextSigPos, coeff, cctx.bdpcm());
      gt1Flag = m_BinDecoder.decodeBin(gt1CtxId);
      binLogger.LogElements(SyntaxElement::abs_mvd_greater0_flag, gt1Flag);
      cctx.decimateNumCtxBins(1);

      unsigned parFlag = 0;
      if (gt1Flag) {
        RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_par);
        parFlag = m_BinDecoder.decodeBin(cctx.parityCtxIdAbsTS());
        binLogger.LogElements(SyntaxElement::par_level_flag, parFlag);
        cctx.decimateNumCtxBins(1);
      }
      coeff[blkPos] = (sign ? -1 : 1) * (TCoeff)(1 + parFlag + gt1Flag);
    }
    lastScanPosPass1 = nextSigPos;
  }

  int cutoffVal = 2;
  const int numGtBins = 4;

  //===== 2nd PASS: gt2 =====
  for (int scanPos = firstSigPos;
       scanPos <= minSubPos && cctx.numCtxBins() >= 4; scanPos++) {
    TCoeff &tcoeff = coeff[cctx.blockPos(scanPos)];
    cutoffVal = 2;
    for (int i = 0; i < numGtBins; i++) {
      if (tcoeff < 0) {
        tcoeff = -tcoeff;
      }
      if (tcoeff >= cutoffVal) {
        RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_gt2);
        unsigned gt2Flag;
        gt2Flag =
            m_BinDecoder.decodeBin(cctx.greaterXCtxIdAbsTS(cutoffVal >> 1));
        binLogger.LogElements(SyntaxElement::abs_mvd_greater1_flag, gt2Flag);
        tcoeff += (gt2Flag << 1);
        cctx.decimateNumCtxBins(1);
      }
      cutoffVal += 2;
    }
    lastScanPosPass2 = scanPos;
  }
  //===== 3rd PASS: Go-rice codes =====
  for (int scanPos = firstSigPos; scanPos <= minSubPos; scanPos++) {
    TCoeff &tcoeff = coeff[cctx.blockPos(scanPos)];
    RExt__DECODER_DEBUG_BIT_STATISTICS_SET(ctype_escs);

    cutoffVal =
        (scanPos <= lastScanPosPass2 ? 10
                                     : (scanPos <= lastScanPosPass1 ? 2 : 0));
    if (tcoeff < 0) {
      tcoeff = -tcoeff;
    }
    if (tcoeff >= cutoffVal) {
      int rice = riceParam;
      int rem = m_BinDecoder.decodeRemAbsEP(rice, COEF_REMAIN_BIN_REDUCTION,
                                            cctx.maxLog2TrDRange());
      binLogger.LogElements(SyntaxElement::abs_remainder, rem);
      tcoeff += (scanPos <= lastScanPosPass1) ? (rem << 1) : rem;
      if (tcoeff && scanPos > lastScanPosPass1) {
        int blkPos = cctx.blockPos(scanPos);
        int sign = m_BinDecoder.decodeBinEP();
        binLogger.LogElements(SyntaxElement::coeff_sign_flag, sign);
        signPattern += (sign << numNonZero);
        sigBlkPos[numNonZero++] = blkPos;
      }
    }
    if (!cctx.bdpcm() && cutoffVal) {
      if (tcoeff > 0) {
        int rightPixel, belowPixel;
        cctx.neighTS(rightPixel, belowPixel, scanPos, coeff);
        tcoeff = cctx.decDeriveModCoeff(rightPixel, belowPixel, tcoeff);
      }
    }
  }

  //===== set final coefficents =====
  for (unsigned k = 0; k < numNonZero; k++) {
    TCoeff AbsCoeff = coeff[sigBlkPos[k]];
    coeff[sigBlkPos[k]] = (signPattern & 1 ? -AbsCoeff : AbsCoeff);
    signPattern >>= 1;
    check_coeff_conformance(cctx, coeff[sigBlkPos[k]]);
  }
}

//================================================================================
//  helper functions
//--------------------------------------------------------------------------------
//    unsigned  unary_max_symbol ( ctxId0, ctxId1, maxSymbol )
//    unsigned  unary_max_eqprob (                 maxSymbol )
//    unsigned  exp_golomb_eqprob( count )
//================================================================================

unsigned CABACReader::unary_max_symbol(unsigned ctxId0, unsigned ctxIdN,
                                       unsigned maxSymbol) {
  unsigned onesRead = 0;
  while (onesRead < maxSymbol &&
         m_BinDecoder.decodeBin(onesRead == 0 ? ctxId0 : ctxIdN) == 1) {
    ++onesRead;
  }
  return onesRead;
}

unsigned CABACReader::unary_max_eqprob(unsigned maxSymbol) {
  for (unsigned k = 0; k < maxSymbol; k++) {
    if (!m_BinDecoder.decodeBinEP()) {
      return k;
    }
  }
  return maxSymbol;
}

unsigned CABACReader::exp_golomb_eqprob(unsigned count) {
  unsigned symbol = 0;
  unsigned bit = 1;
  while (bit) {
    bit = m_BinDecoder.decodeBinEP();
    symbol += bit << count++;
  }
  if (--count) {
    symbol += m_BinDecoder.decodeBinsEP(count);
  }
  return symbol;
}

void CABACReader::mip_flag(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__OTHER);

  if (!cu.Y().valid()) {
    return;
  }
  if (!cu.cs->sps->getUseMIP()) {
    cu.mipFlag = false;
    return;
  }

  unsigned ctxId = DeriveCtx::CtxMipFlag(cu);
  cu.mipFlag = m_BinDecoder.decodeBin(Ctx::MipFlag(ctxId));
  binLogger.LogElements(SyntaxElement::intra_mip_flag, cu.mipFlag);
}

void CABACReader::mip_pred_modes(CodingUnit &cu) {
  RExt__DECODER_DEBUG_BIT_STATISTICS_CREATE_SET(STATS__CABAC_BITS__OTHER);

  if (!cu.Y().valid()) {
    return;
  }
  for (auto &pu : CU::traversePUs(cu)) {
    mip_pred_mode(pu);
  }
}

void CABACReader::mip_pred_mode(PredictionUnit &pu) {
  pu.mipTransposedFlag = bool(m_BinDecoder.decodeBinEP());
  binLogger.LogElements(SyntaxElement::intra_mip_transposed_flag, pu.mipTransposedFlag);

  uint32_t mipMode;
  const int numModes = getNumModesMip(pu.Y());
  xReadTruncBinCode(mipMode, numModes);
  binLogger.LogElements(SyntaxElement::intra_mip_mode, mipMode);
  pu.intraDir[CHANNEL_TYPE_LUMA] = mipMode;
  CHECKD(pu.intraDir[CHANNEL_TYPE_LUMA] < 0 ||
             pu.intraDir[CHANNEL_TYPE_LUMA] >= numModes,
         "Invalid MIP mode");
}
