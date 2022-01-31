#include "slice.hpp"

using namespace Common;

uint32_t PreCalcValues::getValIdx(const Slice &slice,
                                  const ChannelType chType) const {
  return slice.isIntra() ? (ISingleTree ? 0 : (chType << 1)) : 1;
}

uint32_t PreCalcValues::getMaxBtDepth(const Slice &slice,
                                      const ChannelType chType) const {
  if (slice.getPicHeader()->getSplitConsOverrideFlag()) {
    return slice.getPicHeader()->getMaxMTTHierarchyDepth(
        slice.getSliceType(), ISingleTree ? CHANNEL_TYPE_LUMA : chType);
  } else {
    return maxBtDepth[getValIdx(slice, chType)];
  }
}

uint32_t PreCalcValues::getMinBtSize(const Slice &slice,
                                     const ChannelType chType) const {
  return minBtSize[getValIdx(slice, chType)];
}

uint32_t PreCalcValues::getMaxBtSize(const Slice &slice,
                                     const ChannelType chType) const {
  if (slice.getPicHeader()->getSplitConsOverrideFlag()) {
    return slice.getPicHeader()->getMaxBTSize(
        slice.getSliceType(), ISingleTree ? CHANNEL_TYPE_LUMA : chType);
  } else {
    return maxBtSize[getValIdx(slice, chType)];
  }
}

uint32_t PreCalcValues::getMinTtSize(const Slice &slice,
                                     const ChannelType chType) const {
  return minTtSize[getValIdx(slice, chType)];
}

uint32_t PreCalcValues::getMaxTtSize(const Slice &slice,
                                     const ChannelType chType) const {
  if (slice.getPicHeader()->getSplitConsOverrideFlag()) {
    return slice.getPicHeader()->getMaxTTSize(
        slice.getSliceType(), ISingleTree ? CHANNEL_TYPE_LUMA : chType);
  } else {
    return maxTtSize[getValIdx(slice, chType)];
  }
}
uint32_t PreCalcValues::getMinQtSize(const Slice &slice,
                                     const ChannelType chType) const {
  if (slice.getPicHeader()->getSplitConsOverrideFlag()) {
    return slice.getPicHeader()->getMinQTSize(
        slice.getSliceType(), ISingleTree ? CHANNEL_TYPE_LUMA : chType);
  } else {
    return minQtSize[getValIdx(slice, chType)];
  }
}

//! get tables for weighted prediction
const WPScalingParam *Slice::getWpScaling(const RefPicList refPicList,
                                          const int refIdx) const {
  CHECK(refPicList >= NUM_REF_PIC_LIST_01, "Invalid picture reference list");
  if (refIdx < 0) {
    return nullptr;
  } else {
    return m_weightPredTable[refPicList][refIdx].data();
  }
}

WPScalingParam *Slice::getWpScaling(const RefPicList refPicList,
                                    const int refIdx) {
  CHECK(refPicList >= NUM_REF_PIC_LIST_01, "Invalid picture reference list");
  if (refIdx < 0) {
    return nullptr;
  } else {
    return m_weightPredTable[refPicList][refIdx].data();
  }
}