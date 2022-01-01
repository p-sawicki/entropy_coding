#include "coding_structure.hpp"
#include "unit_tools.hpp"

namespace EntropyCoding {

CodingUnit *CodingStructure::getCU(const Position &pos,
                                   const ChannelType effChType) {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos) ||
      (treeType == TREE_C && effChType == CHANNEL_TYPE_LUMA)) {
    // keep this check, which is helpful to identify bugs
    if (treeType == TREE_C && effChType == CHANNEL_TYPE_LUMA) {
      CHECK(parent == nullptr,
            "parent shall be valid; consider using function getLumaCU()");
      CHECK(parent->treeType != TREE_D, "wrong parent treeType ");
    }
    if (parent) {
      return parent->getCU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_cuIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];

    if (idx != 0) {
      return cus[idx - 1].get();
    } else {
      return nullptr;
    }
  }
}

const CodingUnit *CodingStructure::getCU(const Position &pos,
                                         const ChannelType effChType) const {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos) ||
      (treeType == TREE_C && effChType == CHANNEL_TYPE_LUMA)) {
    if (treeType == TREE_C && effChType == CHANNEL_TYPE_LUMA) {
      CHECK(parent == nullptr,
            "parent shall be valid; consider using function getLumaCU()");
      CHECK(parent->treeType != TREE_D, "wrong parent treeType");
    }
    if (parent) {
      return parent->getCU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_cuIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];

    if (idx != 0) {
      return cus[idx - 1].get();
    } else {
      return nullptr;
    }
  }
}

PredictionUnit *CodingStructure::getPU(const Position &pos,
                                       const ChannelType effChType) {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos)) {
    if (parent) {
      return parent->getPU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_puIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];

    if (idx != 0) {
      return pus[idx - 1].get();
    } else {
      return nullptr;
    }
  }
}

const PredictionUnit *
CodingStructure::getPU(const Position &pos, const ChannelType effChType) const {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos)) {
    if (parent) {
      return parent->getPU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_puIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];

    if (idx != 0) {
      return pus[idx - 1].get();
    } else {
      return nullptr;
    }
  }
}

TransformUnit *CodingStructure::getTU(const Position &pos,
                                      const ChannelType effChType,
                                      const int subTuIdx) {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos)) {
    if (parent) {
      return parent->getTU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_tuIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];

    if (idx != 0) {
      unsigned extraIdx = 0;
      if (isLuma(effChType)) {
        const TransformUnit &tu = *tus[idx - 1];

        if (tu.cu->ispMode) // Intra SubPartitions mode
        {
          // we obtain the offset to index the corresponding sub-partition
          if (subTuIdx != -1) {
            extraIdx = subTuIdx;
          } else {
            while (!tus[idx - 1 + extraIdx]
                        ->blocks[getFirstComponentOfChannel(effChType)]
                        .contains(pos)) {
              extraIdx++;
              CHECK(tus[idx - 1 + extraIdx]->cu->treeType == TREE_C,
                    "tu searched by position points to a chroma tree CU");
              CHECK(extraIdx > 3, "extraIdx > 3");
            }
          }
        }
      }
      return tus[idx - 1 + extraIdx].get();
    } else if (m_isTuEnc) {
      return parent->getTU(pos, effChType);
    } else {
      return nullptr;
    }
  }
}

const TransformUnit *CodingStructure::getTU(const Position &pos,
                                            const ChannelType effChType,
                                            const int subTuIdx) const {
  const CompArea &_blk = area.blocks[effChType];

  if (!_blk.contains(pos)) {
    if (parent) {
      return parent->getTU(pos, effChType);
    } else {
      return nullptr;
    }
  } else {
    const unsigned idx = m_tuIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                   unitScale[effChType])];
    if (idx != 0) {
      unsigned extraIdx = 0;
      if (isLuma(effChType)) {
        const TransformUnit &tu = *tus[idx - 1];
        if (tu.cu->ispMode) // Intra SubPartitions mode
        {
          // we obtain the offset to index the corresponding sub-partition
          if (subTuIdx != -1) {
            extraIdx = subTuIdx;
          } else {
            while (!tus[idx - 1 + extraIdx]
                        ->blocks[getFirstComponentOfChannel(effChType)]
                        .contains(pos)) {
              extraIdx++;
              CHECK(tus[idx - 1 + extraIdx]->cu->treeType == TREE_C,
                    "tu searched by position points to a chroma tree CU");
              CHECK(extraIdx > 3, "extraIdx > 3");
            }
          }
        }
      }
      return tus[idx - 1 + extraIdx].get();
    } else if (m_isTuEnc) {
      return parent->getTU(pos, effChType);
    } else {
      return nullptr;
    }
  }
}

CodingUnit *CodingStructure::getLumaCU(const Position &pos) {
  const ChannelType effChType = CHANNEL_TYPE_LUMA;
  const CompArea &_blk = area.blocks[effChType];
  CHECK(!_blk.contains(pos), "must contain the pos");

  const unsigned idx = m_cuIdx[effChType][rsAddr(pos, _blk.pos(), _blk.width,
                                                 unitScale[effChType])];

  if (idx != 0) {
    return cus[idx - 1].get();
  } else {
    return nullptr;
  }
}

const CodingUnit *
CodingStructure::getCURestricted(const Position &pos, const CodingUnit &curCu,
                                 const ChannelType _chType) const {
  const CodingUnit *cu = getCU(pos, _chType);
  // exists       same slice and tile                  cu precedes curCu in
  // encoding order
  //                                                  (thus, is either from
  //                                                  parent CS in RD-search or
  //                                                  its index is lower)
  const bool wavefrontsEnabled =
      curCu.slice->getSPS()->getEntropyCodingSyncEnabledFlag();
  int ctuSizeBit = floorLog2(curCu.cs->sps->getMaxCUWidth());
  int xNbY = pos.x << getChannelTypeScaleX(_chType, curCu.chromaFormat);
  int xCurr = curCu.blocks[_chType].x
              << getChannelTypeScaleX(_chType, curCu.chromaFormat);
  bool addCheck =
      (wavefrontsEnabled && (xNbY >> ctuSizeBit) >= (xCurr >> ctuSizeBit) + 1)
          ? false
          : true;
  if (cu && CU::isSameSliceAndTile(*cu, curCu) &&
      (cu->cs != curCu.cs || cu->idx <= curCu.idx) && addCheck) {
    return cu;
  } else {
    return nullptr;
  }
}

const CodingUnit *CodingStructure::getCURestricted(
    const Position &pos, const Position curPos, const unsigned curSliceIdx,
    const unsigned curTileIdx, const ChannelType _chType) const {
  const CodingUnit *cu = getCU(pos, _chType);
  const bool wavefrontsEnabled =
      this->slice->getSPS()->getEntropyCodingSyncEnabledFlag();
  int ctuSizeBit = floorLog2(this->sps->getMaxCUWidth());
  int xNbY = pos.x << getChannelTypeScaleX(_chType, this->area.chromaFormat);
  int xCurr =
      curPos.x << getChannelTypeScaleX(_chType, this->area.chromaFormat);
  bool addCheck =
      (wavefrontsEnabled && (xNbY >> ctuSizeBit) >= (xCurr >> ctuSizeBit) + 1)
          ? false
          : true;
  return (cu && cu->slice->getIndependentSliceIdx() == curSliceIdx &&
          cu->tileIdx == curTileIdx && addCheck)
             ? cu
             : nullptr;
}

CodingUnit &CodingStructure::addCU(const UnitArea &unit,
                                   const ChannelType chType) {
  CodingUnit *cu = m_cuCache->get();

  cu->UnitArea::operator=(unit);
  cu->initData();
  cu->cs = this;
  cu->slice = nullptr;
  cu->next = nullptr;
  cu->firstPU = nullptr;
  cu->lastPU = nullptr;
  cu->firstTU = nullptr;
  cu->lastTU = nullptr;
  cu->chType = chType;
  cu->treeType = treeType;
  cu->modeType = modeType;

  CodingUnit *prevCU = m_numCUs > 0 ? cus.back().get() : nullptr;

  if (prevCU) {
    prevCU->next = cu;
  }

  cus.push_back(::std::shared_ptr<CodingUnit>(cu));

  uint32_t idx = ++m_numCUs;
  cu->idx = idx;

  uint32_t numCh = getNumberValidChannels(area.chromaFormat);

  for (uint32_t i = 0; i < numCh; i++) {
    if (!cu->blocks[i].valid()) {
      continue;
    }

    const CompArea &_selfBlk = area.blocks[i];
    const CompArea &_blk = cu->blocks[i];

    const UnitScale &scale = unitScale[_blk.compID];
    const Area scaledSelf = scale.scale(_selfBlk);
    const Area scaledBlk = scale.scale(_blk);
    unsigned *idxPtr = m_cuIdx[i] + rsAddr(scaledBlk.pos(), scaledSelf.pos(),
                                           scaledSelf.width);
    CHECK(*idxPtr, "Overwriting a pre-existing value, should be '0'!");
    AreaBuf<uint32_t>(idxPtr, scaledSelf.width, scaledBlk.size()).fill(idx);
  }

  return *cu;
}

const PredictionUnit *
CodingStructure::getPURestricted(const Position &pos,
                                 const PredictionUnit &curPu,
                                 const ChannelType _chType) const {
  const PredictionUnit *pu = getPU(pos, _chType);
  // exists       same slice and tile                  pu precedes curPu in
  // encoding order
  //                                                  (thus, is either from
  //                                                  parent CS in RD-search or
  //                                                  its index is lower)
  const bool wavefrontsEnabled =
      curPu.cu->slice->getSPS()->getEntropyCodingSyncEnabledFlag();
  int ctuSizeBit = floorLog2(curPu.cs->sps->getMaxCUWidth());
  int xNbY = pos.x << getChannelTypeScaleX(_chType, curPu.chromaFormat);
  int xCurr = curPu.blocks[_chType].x
              << getChannelTypeScaleX(_chType, curPu.chromaFormat);
  bool addCheck =
      (wavefrontsEnabled && (xNbY >> ctuSizeBit) >= (xCurr >> ctuSizeBit) + 1)
          ? false
          : true;
  if (pu && CU::isSameSliceAndTile(*pu->cu, *curPu.cu) &&
      (pu->cs != curPu.cs || pu->idx <= curPu.idx) && addCheck) {
    return pu;
  } else {
    return nullptr;
  }
}

PredictionUnit &CodingStructure::addPU(const UnitArea &unit,
                                       const ChannelType chType) {
  PredictionUnit *pu = m_puCache->get();

  pu->UnitArea::operator=(unit);
  pu->initData();
  pu->next = nullptr;
  pu->cs = this;
  pu->cu = m_isTuEnc ? cus[0].get() : getCU(unit.blocks[chType].pos(), chType);
  pu->chType = chType;

  PredictionUnit *prevPU = m_numPUs > 0 ? pus.back().get() : nullptr;

  if (prevPU && prevPU->cu == pu->cu) {
    prevPU->next = pu;
  }

  pus.push_back(::std::shared_ptr<PredictionUnit>(pu));

  if (pu->cu->firstPU == nullptr) {
    pu->cu->firstPU = pu;
  }
  pu->cu->lastPU = pu;

  uint32_t idx = ++m_numPUs;
  pu->idx = idx;

  uint32_t numCh = getNumberValidChannels(area.chromaFormat);
  for (uint32_t i = 0; i < numCh; i++) {
    if (!pu->blocks[i].valid()) {
      continue;
    }

    const CompArea &_selfBlk = area.blocks[i];
    const CompArea &_blk = pu->blocks[i];

    const UnitScale &scale = unitScale[_blk.compID];
    const Area scaledSelf = scale.scale(_selfBlk);
    const Area scaledBlk = scale.scale(_blk);
    unsigned *idxPtr = m_puIdx[i] + rsAddr(scaledBlk.pos(), scaledSelf.pos(),
                                           scaledSelf.width);
    CHECK(*idxPtr, "Overwriting a pre-existing value, should be '0'!");
    AreaBuf<uint32_t>(idxPtr, scaledSelf.width, scaledBlk.size()).fill(idx);
  }

  return *pu;
}

TransformUnit &CodingStructure::addTU(const UnitArea &unit,
                                      const ChannelType chType) {
  TransformUnit *tu = m_tuCache->get();

  tu->UnitArea::operator=(unit);
  tu->initData();
  tu->next = nullptr;
  tu->prev = nullptr;
  tu->cs = this;
  tu->cu = m_isTuEnc ? cus[0].get() : getCU(unit.blocks[chType].pos(), chType);
  tu->chType = chType;

  TransformUnit *prevTU = m_numTUs > 0 ? tus.back().get() : nullptr;

  if (prevTU && prevTU->cu == tu->cu) {
    prevTU->next = tu;
    tu->prev = prevTU;
  }

  tus.push_back(::std::shared_ptr<TransformUnit>(tu));

  if (tu->cu) {
    if (tu->cu->firstTU == nullptr) {
      tu->cu->firstTU = tu;
    }
    tu->cu->lastTU = tu;
  }

  uint32_t idx = ++m_numTUs;
  tu->idx = idx;

  TCoeff *coeffs[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  Pel *pcmbuf[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  bool *runType[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};

  uint32_t numCh = getNumberValidComponents(area.chromaFormat);

  for (uint32_t i = 0; i < numCh; i++) {
    if (!tu->blocks[i].valid()) {
      continue;
    }

    if (i < getNumberValidChannels(area.chromaFormat)) {
      const CompArea &_selfBlk = area.blocks[i];
      const CompArea &_blk = tu->blocks[i];

      bool isIspTu =
          tu->cu != nullptr && tu->cu->ispMode && isLuma(_blk.compID);

      bool isFirstIspTu = false;
      if (isIspTu) {
        isFirstIspTu = CU::isISPFirst(
            *tu->cu, _blk, getFirstComponentOfChannel(ChannelType(i)));
      }
      if (!isIspTu || isFirstIspTu) {
        const UnitScale &scale = unitScale[_blk.compID];

        const Area scaledSelf = scale.scale(_selfBlk);
        const Area scaledBlk =
            isIspTu ? scale.scale(tu->cu->blocks[i]) : scale.scale(_blk);
        unsigned *idxPtr =
            m_tuIdx[i] +
            rsAddr(scaledBlk.pos(), scaledSelf.pos(), scaledSelf.width);
        CHECK(*idxPtr, "Overwriting a pre-existing value, should be '0'!");
        AreaBuf<uint32_t>(idxPtr, scaledSelf.width, scaledBlk.size()).fill(idx);
      }
    }

    coeffs[i] = m_coeffs[i] + m_offsets[i];
    pcmbuf[i] = m_pcmbuf[i] + m_offsets[i];

    if (i < MAX_NUM_CHANNEL_TYPE) {
      if (m_runType[i] != nullptr) {
        runType[i] = m_runType[i] + m_offsets[i];
      }
    }

    unsigned areaSize = tu->blocks[i].area();
    m_offsets[i] += areaSize;
  }
  tu->init(coeffs, pcmbuf, runType);

  return *tu;
}

void CodingStructure::addEmptyTUs(Partitioner &partitioner) {
  const UnitArea &area = partitioner.currArea();
  bool split = partitioner.canSplit(TU_MAX_TR_SPLIT, *this);
  const unsigned trDepth = partitioner.currTrDepth;

  if (split) {
    partitioner.splitCurrArea(TU_MAX_TR_SPLIT, *this);
    do {
      addEmptyTUs(partitioner);
    } while (partitioner.nextPart(*this));

    partitioner.exitCurrSplit();
  } else {
    TransformUnit &tu = this->addTU(
        CS::getArea(*this, area, partitioner.chType), partitioner.chType);
    unsigned numBlocks = getNumberValidTBlocks(*this->pcv);
    for (unsigned compID = COMPONENT_Y; compID < numBlocks; compID++) {
      if (tu.blocks[compID].valid()) {
        tu.getCoeffs(ComponentID(compID)).fill(0);
        tu.getPcmbuf(ComponentID(compID)).fill(0);
      }
    }
    tu.depth = trDepth;
  }
}

const int CodingStructure::signalModeCons(const PartSplit split,
                                          Partitioner &partitioner,
                                          const ModeType modeTypeParent) const {
  if (CS::isDualITree(*this) || modeTypeParent != MODE_TYPE_ALL ||
      partitioner.currArea().chromaFormat == CHROMA_444 ||
      partitioner.currArea().chromaFormat == CHROMA_400) {
    return LDT_MODE_TYPE_INHERIT;
  }
  int minLumaArea = partitioner.currArea().lumaSize().area();
  if (split == CU_QUAD_SPLIT || split == CU_TRIH_SPLIT ||
      split == CU_TRIV_SPLIT) // the area is split into 3 or 4 parts
  {
    minLumaArea = minLumaArea >> 2;
  } else if (split == CU_VERT_SPLIT ||
             split == CU_HORZ_SPLIT) // the area is split into 2 parts
  {
    minLumaArea = minLumaArea >> 1;
  }
  int minChromaBlock =
      minLumaArea >>
      (getChannelTypeScaleX(CHANNEL_TYPE_CHROMA,
                            partitioner.currArea().chromaFormat) +
       getChannelTypeScaleY(CHANNEL_TYPE_CHROMA,
                            partitioner.currArea().chromaFormat));
  bool is2xNChroma = (partitioner.currArea().chromaSize().width == 4 &&
                      split == CU_VERT_SPLIT) ||
                     (partitioner.currArea().chromaSize().width == 8 &&
                      split == CU_TRIV_SPLIT);
  return minChromaBlock >= 16 && !is2xNChroma
             ? LDT_MODE_TYPE_INHERIT
             : ((minLumaArea < 32) || slice->isIntra()) ? LDT_MODE_TYPE_INFER
                                                        : LDT_MODE_TYPE_SIGNAL;
}

void CodingStructure::reorderPrevPLT(
    PLTBuf &prevPLT, ::std::array<uint8_t, MAX_NUM_CHANNEL_TYPE> &curPLTSize,
    CurPLT31 &curPLT, ReuseFlag &reuseflag, uint32_t compBegin,
    uint32_t numComp, bool jointPLT) {
  Pel stuffedPLT[MAX_NUM_COMPONENT][MAXPLTPREDSIZE];
  uint8_t tempCurPLTsize[MAX_NUM_CHANNEL_TYPE];
  uint8_t stuffPLTsize[MAX_NUM_COMPONENT];

  uint32_t maxPredPltSize = jointPLT ? MAXPLTPREDSIZE : MAXPLTPREDSIZE_DUALTREE;

  for (int i = compBegin; i < (compBegin + numComp); i++) {
    ComponentID comID = jointPLT ? (ComponentID)compBegin
                                 : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
    tempCurPLTsize[comID] = curPLTSize[comID];
    stuffPLTsize[i] = 0;
    memcpy(stuffedPLT[i], curPLT[i].data(), curPLTSize[comID] * sizeof(Pel));
  }

  for (int ch = compBegin; ch < (compBegin + numComp); ch++) {
    ComponentID comID = jointPLT ? (ComponentID)compBegin
                                 : ((ch > 0) ? COMPONENT_Cb : COMPONENT_Y);
    if (ch > 1)
      break;
    for (int i = 0; i < prevPLT.curPLTSize[comID]; i++) {
      if (tempCurPLTsize[comID] + stuffPLTsize[ch] >= maxPredPltSize) {
        break;
      }

      if (!reuseflag[comID][i]) {
        if (ch == COMPONENT_Y) {
          stuffedPLT[0][tempCurPLTsize[comID] + stuffPLTsize[ch]] =
              prevPLT.curPLT[0][i];
        } else {
          stuffedPLT[1][tempCurPLTsize[comID] + stuffPLTsize[ch]] =
              prevPLT.curPLT[1][i];
          stuffedPLT[2][tempCurPLTsize[comID] + stuffPLTsize[ch]] =
              prevPLT.curPLT[2][i];
        }
        stuffPLTsize[ch]++;
      }
    }
  }

  for (int i = compBegin; i < (compBegin + numComp); i++) {
    ComponentID comID = jointPLT ? (ComponentID)compBegin
                                 : ((i > 0) ? COMPONENT_Cb : COMPONENT_Y);
    prevPLT.curPLTSize[comID] = curPLTSize[comID] + stuffPLTsize[comID];
    memcpy(prevPLT.curPLT[i].data(), stuffedPLT[i],
           prevPLT.curPLTSize[comID] * sizeof(Pel));
    CHECK(prevPLT.curPLTSize[comID] > maxPredPltSize,
          " Maximum palette predictor size exceed limit");
  }
}
} // namespace EntropyCoding