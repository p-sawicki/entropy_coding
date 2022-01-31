#ifndef COMMON_SLICE
#define COMMON_SLICE

#include <cstring>
#include <list>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "alf_parameters.hpp"
#include "picture.hpp"

namespace Common {

struct MotionInfo;

struct Picture;
class Pic;
class TrQuant;
// ====================================================================================================================
// Constants
// ====================================================================================================================
class PreCalcValues;
static const uint32_t REF_PIC_LIST_NUM_IDX = 32;

typedef std::list<Picture *> PicList;

// ====================================================================================================================
// Class definition
// ====================================================================================================================

/// SPS RExt class
class SPSRExt // Names aligned to text specification
{
private:
  bool m_extendedPrecisionProcessingFlag;
  bool m_tsrcRicePresentFlag;
  bool m_persistentRiceAdaptationEnabledFlag;
  bool m_rrcRiceExtensionEnableFlag;

public:
  SPSRExt(const bool extendedPrecisionProcessingFlag,
          const bool tsrcRicePresentFlag,
          const bool persistentRiceAdaptationEnabledFlag,
          const bool rrcRiceExtensionEnableFlag)
      : m_extendedPrecisionProcessingFlag(extendedPrecisionProcessingFlag),
        m_tsrcRicePresentFlag(tsrcRicePresentFlag),
        m_persistentRiceAdaptationEnabledFlag(
            persistentRiceAdaptationEnabledFlag),
        m_rrcRiceExtensionEnableFlag(rrcRiceExtensionEnableFlag) {}

  bool getExtendedPrecisionProcessingFlag() const {
    return m_extendedPrecisionProcessingFlag;
  }

  bool getTSRCRicePresentFlag() const { return m_tsrcRicePresentFlag; }

  bool getPersistentRiceAdaptationEnabledFlag() const {
    return m_persistentRiceAdaptationEnabledFlag;
  }

  bool getRrcRiceExtensionEnableFlag() const {
    return m_rrcRiceExtensionEnableFlag;
  }
};

/// SPS class
class SPS {
private:
  bool m_affineAmvrEnabledFlag;
  bool m_MMVD;
  bool m_SBT;
  bool m_ISP;
  ChromaFormat m_chromaFormatIdc;

  int m_log2MinCodingBlockSize;
  unsigned m_CTUSize;
  uint32_t m_uiMaxCUWidth;

  // Tool list

  bool m_transformSkipEnabledFlag;
  int m_log2MaxTransformSkipBlockSize;
  bool m_BDPCMEnabledFlag;
  bool m_JointCbCrEnabledFlag;
  // Parameter
  BitDepths m_bitDepths;
  bool m_EntropyCodingSyncEnabledFlag; //!< Flag for enabling WPP
  std::array<int, MAX_NUM_CHANNEL_TYPE> m_qpBDOffset;
  uint32_t m_log2MaxTbSize;

  bool m_saoEnabledFlag;
  SPSRExt m_spsRangeExtension;
  bool m_alfEnabledFlag;
  bool m_ccalfEnabledFlag;
  unsigned m_IBCFlag;
  bool m_useColorTrans;
  unsigned m_PLTMode;

  bool m_AMVREnabledFlag;
  bool m_LMChroma;
  bool m_MTS;
  bool m_IntraMTS; // 18
  bool m_InterMTS; // 19
  bool m_LFNST;
  bool m_Affine;
  bool m_AffineType;
  bool m_bcw; //
  bool m_ciip;
  bool m_Geo;
  bool m_MRL;
  bool m_MIP;
  uint32_t m_maxNumMergeCand;
  uint32_t m_maxNumIBCMergeCand;
  uint32_t m_maxNumGeoCand;

public:
  SPS(const bool affineAmvrEnabledFlag, const bool MMVD, const bool SBT,
      const bool ISP, const ChromaFormat chromaFormatIdc,
      const int log2MinCodingBlockSize, const unsigned CTUSize,
      const uint32_t uiMaxCUWidth, const bool transformSkipEnabledFlag,
      const int log2MaxTransformSkipBlockSize, const bool BDPCMEnabledFlag,
      const bool JointCbCrEnabledFlag, const BitDepths &bitDepths,
      const bool entropyCodingSyncEnabledFlag, const int *qpBDOffset,
      const uint32_t log2MaxTbSize, const bool saoEnabledFlag,
      const SPSRExt &spsRangeExtension, const bool alfEnabledFlag,
      const bool ccalfEnabledFlag, const unsigned IBCFlag,
      const bool useColorTrans, const unsigned PLTMode,
      const bool AMVREnabledFlag, const bool LMChroma, const bool MTS,
      const bool IntraMTS, const bool InterMTS, const bool LFNST,
      const bool Affine, const bool AffineType, const bool bcw, const bool ciip,
      const bool Geo, const bool MRL, const bool MIP,
      const uint32_t maxNumMergeCand, const uint32_t maxNumIBCMergeCand,
      const uint32_t maxNumGeoCand)
      : m_affineAmvrEnabledFlag(affineAmvrEnabledFlag), m_MMVD(MMVD),
        m_SBT(SBT), m_ISP(ISP), m_chromaFormatIdc(chromaFormatIdc),
        m_log2MinCodingBlockSize(log2MinCodingBlockSize), m_CTUSize(CTUSize),
        m_uiMaxCUWidth(uiMaxCUWidth),
        m_transformSkipEnabledFlag(transformSkipEnabledFlag),
        m_log2MaxTransformSkipBlockSize(log2MaxTransformSkipBlockSize),
        m_BDPCMEnabledFlag(BDPCMEnabledFlag),
        m_JointCbCrEnabledFlag(JointCbCrEnabledFlag), m_bitDepths(bitDepths),
        m_EntropyCodingSyncEnabledFlag(entropyCodingSyncEnabledFlag),
        m_log2MaxTbSize(log2MaxTbSize), m_saoEnabledFlag(saoEnabledFlag),
        m_spsRangeExtension(spsRangeExtension),
        m_alfEnabledFlag(alfEnabledFlag), m_ccalfEnabledFlag(ccalfEnabledFlag),
        m_IBCFlag(IBCFlag), m_useColorTrans(useColorTrans), m_PLTMode(PLTMode),
        m_AMVREnabledFlag(AMVREnabledFlag), m_LMChroma(LMChroma), m_MTS(MTS),
        m_IntraMTS(IntraMTS), m_InterMTS(InterMTS), m_LFNST(LFNST),
        m_Affine(Affine), m_AffineType(AffineType), m_bcw(bcw), m_ciip(ciip),
        m_Geo(Geo), m_MRL(MRL), m_MIP(MIP), m_maxNumMergeCand(maxNumMergeCand),
        m_maxNumIBCMergeCand(maxNumIBCMergeCand),
        m_maxNumGeoCand(maxNumGeoCand) {
    std::copy(qpBDOffset, qpBDOffset + MAX_NUM_CHANNEL_TYPE,
                m_qpBDOffset.begin());
  }

  ChromaFormat getChromaFormatIdc() const { return m_chromaFormatIdc; }
  int getLog2MinCodingBlockSize() const { return m_log2MinCodingBlockSize; }
  unsigned getCTUSize() const { return m_CTUSize; }
  uint32_t getMaxCUWidth() const { return m_uiMaxCUWidth; }
  bool getTransformSkipEnabledFlag() const {
    return m_transformSkipEnabledFlag;
  }
  uint32_t getLog2MaxTransformSkipBlockSize() const {
    return m_log2MaxTransformSkipBlockSize;
  }
  bool getBDPCMEnabledFlag() const { return m_BDPCMEnabledFlag; }
  uint32_t getMaxTbSize() const { return 1 << m_log2MaxTbSize; }
  uint32_t getLog2MaxTbSize() const {return m_log2MaxTbSize; }
  // Bit-depth
  int getBitDepth(ChannelType type) const { return m_bitDepths.recon[type]; }
  const BitDepths &getBitDepths() const { return m_bitDepths; }

  bool getEntropyCodingSyncEnabledFlag() const {
    return m_EntropyCodingSyncEnabledFlag;
  }
#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
  int getMaxLog2TrDynamicRange(ChannelType channelType) const {
    return getSpsRangeExtension().getExtendedPrecisionProcessingFlag()
               ? std::min<int>(20, int(m_bitDepths.recon[channelType] + 6))
               : 15;
  }
#else
  int getMaxLog2TrDynamicRange(ChannelType channelType) const {
    return getSpsRangeExtension().getExtendedPrecisionProcessingFlag() &&
                   int(m_bitDepths.recon[channelType]) > 10
               ? std::min<int>(20, int(m_bitDepths.recon[channelType] + 6))
               : 15;
  }
#endif
  int getQpBDOffset(ChannelType type) const { return m_qpBDOffset[type]; }
  bool getSAOEnabledFlag() const { return m_saoEnabledFlag; }

  bool getALFEnabledFlag() const { return m_alfEnabledFlag; }
  bool getCCALFEnabledFlag() const { return m_ccalfEnabledFlag; }
  bool getJointCbCrEnabledFlag() const { return m_JointCbCrEnabledFlag; }

  bool getUseMMVD() const { return m_MMVD; }
  uint32_t getMaxNumMergeCand() const { return m_maxNumMergeCand; }
  uint32_t getMaxNumIBCMergeCand() const { return m_maxNumIBCMergeCand; }
  uint32_t getMaxNumGeoCand() const { return m_maxNumGeoCand; }
  bool getAffineAmvrEnabledFlag() const { return m_affineAmvrEnabledFlag; }

  const SPSRExt &getSpsRangeExtension() const { return m_spsRangeExtension; }
  SPSRExt &getSpsRangeExtension() { return m_spsRangeExtension; }

  unsigned getIBCFlag() const { return m_IBCFlag; }
  bool getUseColorTrans() const { return m_useColorTrans; }
  unsigned getPLTMode() const { return m_PLTMode; }
  bool getUseSBT() const { return m_SBT; }
  bool getUseISP() const { return m_ISP; }

  bool getAMVREnabledFlag() const { return m_AMVREnabledFlag; }
  bool getUseAffine() const { return m_Affine; }
  bool getUseAffineType() const { return m_AffineType; }
  bool getUseLMChroma() const { return m_LMChroma; }
  bool getUseMTS() const { return m_MTS; }
  bool getUseIntraMTS() const { return m_IntraMTS; }
  bool getUseInterMTS() const { return m_InterMTS; }
  bool getUseLFNST() const { return m_LFNST; }
  bool getUseBcw() const { return m_bcw; }
  bool getUseCiip() const { return m_ciip; }
  bool getUseGeo() const { return m_Geo; }
  bool getUseMRL() const { return m_MRL; }
  bool getUseMIP() const { return m_MIP; }
};

/// PPS class
class PPS {
private:
  bool m_useDQP;

  // Chroma QP Adjustments
  int m_chromaQpOffsetListLen; // size (excludes the null entry used in the
                               // following array).
  uint8_t m_ctuSize;           //!< CTU size
  uint32_t m_numTileCols;      //!< number of tile columns
  std::vector<uint32_t>
      m_tileColBd; //!< tile column left-boundaries in units of CTUs
  std::vector<uint32_t> m_ctuToTileCol; //!< mapping between CTU horizontal
                                          //!< address and tile column index
  std::vector<uint32_t> m_ctuToTileRow; //!< mapping between CTU vertical
                                          //!< address and tile row index

  bool m_cabacInitPresentFlag;
  uint32_t m_picWidthInLumaSamples;
  uint32_t m_picHeightInLumaSamples;

public:
  PPS(const bool useDQP, const int chromaQpOffsetListLen, const uint8_t ctuSize,
      const uint32_t numTileCols, const std::vector<uint32_t> &tileColBd,
      const std::vector<uint32_t> &ctuToTileCol,
      const std::vector<uint32_t> &ctuToTileRow,
      const bool cabacInitPresentFlag, const uint32_t picWidthInLumaSamples,
      const uint32_t picHeightInLumaSamples)
      : m_useDQP(useDQP), m_chromaQpOffsetListLen(chromaQpOffsetListLen),
        m_ctuSize(ctuSize), m_numTileCols(numTileCols), m_tileColBd(tileColBd),
        m_ctuToTileCol(ctuToTileCol), m_ctuToTileRow(ctuToTileRow),
        m_cabacInitPresentFlag(cabacInitPresentFlag),
        m_picWidthInLumaSamples(picWidthInLumaSamples),
        m_picHeightInLumaSamples(picHeightInLumaSamples) {}

  bool getUseDQP() const { return m_useDQP; }
  int getChromaQpOffsetListLen() const { return m_chromaQpOffsetListLen; }
  uint32_t getNumTileColumns() const { return m_numTileCols; }
  uint32_t getTileColumnBd(int idx) const {
    CHECK(idx >= m_tileColBd.size(), "Tile column index exceeds valid range");
    return m_tileColBd[idx];
  }
  uint32_t ctuToTileCol(int ctuX) const {
    CHECK(ctuX >= m_ctuToTileCol.size(),
          "CTU address index exceeds valid range");
    return m_ctuToTileCol[ctuX];
  }
  uint32_t ctuToTileRow(int ctuY) const {
    CHECK(ctuY >= m_ctuToTileRow.size(),
          "CTU address index exceeds valid range");
    return m_ctuToTileRow[ctuY];
  }
  uint32_t getTileIdx(uint32_t ctuX, uint32_t ctuY) const {
    return (ctuToTileRow(ctuY) * getNumTileColumns()) + ctuToTileCol(ctuX);
  }
  uint32_t getTileIdx(const Position &pos) const {
    return getTileIdx(pos.x / m_ctuSize, pos.y / m_ctuSize);
  }

  bool getCabacInitPresentFlag() const { return m_cabacInitPresentFlag; }
  uint32_t getPicWidthInLumaSamples() const { return m_picWidthInLumaSamples; }
  uint32_t getPicHeightInLumaSamples() const {
    return m_picHeightInLumaSamples;
  }
};

class APS {
private:
  AlfParam m_alfAPSParam;

public:
  APS() {}
  APS(const AlfParam &alfAPSParam) : m_alfAPSParam(alfAPSParam) {}

  AlfParam &getAlfAPSParam() { return m_alfAPSParam; }
  const AlfParam &getAlfAPSParam() const {return m_alfAPSParam; }
};

struct WPScalingParam {
  // Explicit weighted prediction parameters parsed in slice header,
  // or Implicit weighted prediction parameters (8 bits depth values).
  bool presentFlag;
  uint32_t log2WeightDenom;
  int codedWeight;
  int codedOffset;

  // Weighted prediction scaling values built from above parameters (bitdepth
  // scaled):
  int w;
  int o;
  int offset;
  int shift;
  int round;

  static bool isWeighted(const WPScalingParam *wp);
};

inline bool WPScalingParam::isWeighted(const WPScalingParam *wp) {
  return wp != nullptr &&
         (wp[COMPONENT_Y].presentFlag || wp[COMPONENT_Cb].presentFlag ||
          wp[COMPONENT_Cr].presentFlag);
}

// picture header class
class PicHeader {
private:
  bool m_splitConsOverrideFlag;    //!< partitioning constraint override flag
  uint32_t m_cuQpDeltaSubdivIntra; //!< CU QP delta maximum subdivision for
                                   //!< intra slices
  uint32_t m_cuQpDeltaSubdivInter; //!< CU QP delta maximum subdivision for
                                   //!< inter slices
  uint32_t m_cuChromaQpOffsetSubdivIntra; //!< CU chroma QP offset maximum
                                          //!< subdivision for intra slices
  uint32_t m_cuChromaQpOffsetSubdivInter; //!< CU chroma QP offset maximum
                                          //!< subdivision for inter slices
  bool m_mvdL1ZeroFlag;                   //!< L1 MVD set to zero flag
  uint32_t
      m_maxNumAffineMergeCand; //!< max number of sub-block merge candidates
  std::array<unsigned, 3> m_minQT; //!< minimum quad-tree size  0: I slice luma; 1: P/B
                       //!< slice; 2: I slice chroma
  std::array<unsigned, 3> m_maxMTTHierarchyDepth; //!< maximum MTT depth
  std::array<unsigned, 3> m_maxBTSize;            //!< maximum BT size
  std::array<unsigned, 3> m_maxTTSize;            //!< maximum TT size

public:
  PicHeader(const bool splitConsOverrideFlag,
            const uint32_t cuQpDeltaSubdivIntra,
            const uint32_t cuQpDeltaSubdivInter,
            const uint32_t cuChromaQpOffsetSubdivIntra,
            const uint32_t cuChromaQpOffsetSubdivInter,
            const bool mvdL1ZeroFlag, const uint32_t maxNumAffineMergeCand,
            const unsigned *minQT, const unsigned *maxMTTHierarchyDepth,
            const unsigned *maxBTSize, const unsigned *maxTTSize)
      : m_splitConsOverrideFlag(splitConsOverrideFlag),
        m_cuQpDeltaSubdivIntra(cuQpDeltaSubdivIntra),
        m_cuQpDeltaSubdivInter(cuQpDeltaSubdivInter),
        m_cuChromaQpOffsetSubdivIntra(cuChromaQpOffsetSubdivIntra),
        m_cuChromaQpOffsetSubdivInter(cuChromaQpOffsetSubdivInter),
        m_mvdL1ZeroFlag(mvdL1ZeroFlag),
        m_maxNumAffineMergeCand(maxNumAffineMergeCand) {
    std::copy(minQT, minQT + 3, m_minQT.begin());
    std::copy(maxMTTHierarchyDepth, maxMTTHierarchyDepth + 3,
                m_maxMTTHierarchyDepth.begin());
    std::copy(maxBTSize, maxBTSize + 3, m_maxBTSize.begin());
    std::copy(maxTTSize, maxTTSize + 3, m_maxTTSize.begin());
  }

  bool getSplitConsOverrideFlag() const { return m_splitConsOverrideFlag; }
  uint32_t getCuQpDeltaSubdivIntra() const { return m_cuQpDeltaSubdivIntra; }
  uint32_t getCuQpDeltaSubdivInter() const { return m_cuQpDeltaSubdivInter; }
  uint32_t getCuChromaQpOffsetSubdivIntra() const {
    return m_cuChromaQpOffsetSubdivIntra;
  }
  uint32_t getCuChromaQpOffsetSubdivInter() const {
    return m_cuChromaQpOffsetSubdivInter;
  }
  bool getMvdL1ZeroFlag() const { return m_mvdL1ZeroFlag; }
  uint32_t getMaxNumAffineMergeCand() const { return m_maxNumAffineMergeCand; }

  unsigned getMinQTSize(SliceType slicetype,
                        ChannelType chType = CHANNEL_TYPE_LUMA) const {
    return slicetype == I_SLICE
               ? (chType == CHANNEL_TYPE_LUMA ? m_minQT[0] : m_minQT[2])
               : m_minQT[1];
  }
  const std::array<unsigned, 3> &getMinQt() const {return m_minQT; }
  unsigned
  getMaxMTTHierarchyDepth(SliceType slicetype,
                          ChannelType chType = CHANNEL_TYPE_LUMA) const {
    return slicetype == I_SLICE
               ? (chType == CHANNEL_TYPE_LUMA ? m_maxMTTHierarchyDepth[0]
                                              : m_maxMTTHierarchyDepth[2])
               : m_maxMTTHierarchyDepth[1];
  }
  const std::array<unsigned, 3> &getMaxMTTHierarchyDepth() const {return m_maxMTTHierarchyDepth;}
  unsigned getMaxBTSize(SliceType slicetype,
                        ChannelType chType = CHANNEL_TYPE_LUMA) const {
    return slicetype == I_SLICE
               ? (chType == CHANNEL_TYPE_LUMA ? m_maxBTSize[0] : m_maxBTSize[2])
               : m_maxBTSize[1];
  }
  const std::array<unsigned, 3> &getMaxBTSize() const {return m_maxBTSize;}
  unsigned getMaxTTSize(SliceType slicetype,
                        ChannelType chType = CHANNEL_TYPE_LUMA) const {
    return slicetype == I_SLICE
               ? (chType == CHANNEL_TYPE_LUMA ? m_maxTTSize[0] : m_maxTTSize[2])
               : m_maxTTSize[1];
  }
  const std::array<unsigned, 3> &getMaxTTSize() const {return m_maxTTSize;}
};

typedef std::array<
    std::array<std::array<WPScalingParam, MAX_NUM_COMPONENT>, MAX_NUM_REF>,
    NUM_REF_PIC_LIST_01>
    WeightPredTable;

/// slice header class
class Slice {

private:
  //  Bitstream writing
  std::array<bool, MAX_NUM_CHANNEL_TYPE>  m_saoEnabledFlag;
  SliceType m_eSliceType;
  int m_iSliceQp;
  bool m_ChromaQpAdjEnabled;
  bool m_depQuantEnabledFlag;       //!< dependent quantization enabled flag
  int m_riceBaseLevelValue;         //< baseLevel value for abs_remainder
  bool m_signDataHidingEnabledFlag; //!< sign data hiding enabled flag
  bool m_tsResidualCodingDisabledFlag;
  std::array<int, NUM_REF_PIC_LIST_01> m_aiNumRefIdx; //  for multiple reference of current
                                          //  slice

  bool m_bCheckLDC;
  bool m_biDirPred;
  std::array<int, 2> m_symRefIdx;

  // access channel
  const SPS *m_pcSPS;
  const PPS *m_pcPPS;
  Picture *m_pcPic;
  const PicHeader *m_pcPicHeader; //!< pointer to picture header structure
  uint32_t m_independentSliceIdx;
  WeightPredTable m_weightPredTable; // [REF_PIC_LIST_0 or
                                     // REF_PIC_LIST_1][refIdx][0:Y,
                                     // 1:U, 2:V]

  bool m_cabacInitFlag;

  SliceType
      m_encCABACTableIdx; // Used to transmit table selection across slices.
  std::array<APS*, ALF_CTB_MAX_NUM_APS> m_alfApss;
  std::array<bool, MAX_NUM_COMPONENT> m_alfEnabledFlag;
  int m_numAlfApsIdsLuma;
  int m_alfApsIdChroma;
  int m_tsrc_index;
  std::array<unsigned, 8> m_riceBit;

public:
  Slice(const bool *saoEnabledFlag, const SliceType eSliceType,
        const int iSliceQp, const bool ChromaQpAdjEnabled,
        const bool depQuantEnabledFlag, const int riceBaseLevelValue,
        const bool signDataHidingEnabledFlag,
        const bool tsResidualCodingDisabledFlag, const int *aiNumRefIdx,
        const bool bCheckLDC, const bool biDirPred, const int *symRefIdx,
        const uint32_t independentSliceIdx,
        const WeightPredTable &weightPredTable, const bool cabacInitFlag,
        const SliceType encCABACTableIdx,
        const std::array<APS *, ALF_CTB_MAX_NUM_APS> &alfApss,
        const bool *alfEnabledFlag, const int numAlfApsIdsLuma,
        const int alfApsIdChroma, const int tsrc_index, const unsigned *riceBit,
        const CcAlfFilterParam &ccAlfFilterParam, uint8_t* const *ccAlfFilterControl)
      : m_eSliceType(eSliceType), m_iSliceQp(iSliceQp),
        m_ChromaQpAdjEnabled(ChromaQpAdjEnabled),
        m_depQuantEnabledFlag(depQuantEnabledFlag),
        m_riceBaseLevelValue(riceBaseLevelValue),
        m_signDataHidingEnabledFlag(signDataHidingEnabledFlag),
        m_tsResidualCodingDisabledFlag(tsResidualCodingDisabledFlag),
        m_bCheckLDC(bCheckLDC), m_biDirPred(biDirPred),
        m_independentSliceIdx(independentSliceIdx),
        m_weightPredTable(weightPredTable), m_cabacInitFlag(cabacInitFlag),
        m_encCABACTableIdx(encCABACTableIdx), m_alfApss(alfApss),
        m_numAlfApsIdsLuma(numAlfApsIdsLuma), m_alfApsIdChroma(alfApsIdChroma),
        m_tsrc_index(tsrc_index), m_ccAlfFilterParam(ccAlfFilterParam) {
    std::copy(saoEnabledFlag, saoEnabledFlag + m_saoEnabledFlag.size(),
                m_saoEnabledFlag.begin());
    std::copy(aiNumRefIdx, aiNumRefIdx + m_aiNumRefIdx.size(),
                m_aiNumRefIdx.begin());
    std::copy(symRefIdx, symRefIdx + m_symRefIdx.size(), m_symRefIdx.begin());
    std::copy(alfEnabledFlag, alfEnabledFlag + m_alfEnabledFlag.size(),
                m_alfEnabledFlag.begin());
    std::copy(riceBit, riceBit + m_riceBit.size(), m_riceBit.begin());
    copy_array(ccAlfFilterControl, m_ccAlfFilterControl);
  }

  ~Slice() {
    std::for_each(m_alfApss.begin(), m_alfApss.end(), [](APS *aps) {
      delete aps;
    });
  }

  const PicHeader *getPicHeader() const { return m_pcPicHeader; }
  void setPicHeader(const PicHeader *picHeader) { m_pcPicHeader = picHeader; }
  const SPS *getSPS() const { return m_pcSPS; }
  void setSPS(const SPS *sps) { m_pcSPS = sps; }
  const PPS *getPPS() const { return m_pcPPS; }
  void setPPS(const PPS *pps) { m_pcPPS = pps; }

  APS **getAlfAPSs() { return m_alfApss.data(); }
  const std::array<APS *, ALF_CTB_MAX_NUM_APS> &getAlfAPSs() const {
    return m_alfApss;
  }
  bool getSaoEnabledFlag(ChannelType chType) const {
    return m_saoEnabledFlag[chType];
  }
  SliceType getSliceType() const { return m_eSliceType; }
  int getSliceQp() const { return m_iSliceQp; }
  bool getUseChromaQpAdj() const { return m_ChromaQpAdjEnabled; }

  int getNumRefIdx(RefPicList e) const { return m_aiNumRefIdx[e]; }
  Picture *getPic() { return m_pcPic; }
  const Picture *getPic() const { return m_pcPic; }
  void setPic(Picture *pic) { m_pcPic = pic; }
  bool getCheckLDC() const { return m_bCheckLDC; }
  bool getDepQuantEnabledFlag() const { return m_depQuantEnabledFlag; }
  int getRiceBaseLevel() const { return m_riceBaseLevelValue; }
  bool getSignDataHidingEnabledFlag() const {
    return m_signDataHidingEnabledFlag;
  }
  bool getTSResidualCodingDisabledFlag() const {
    return m_tsResidualCodingDisabledFlag;
  }
  bool getBiDirPred() const { return m_biDirPred; }
  int getSymRefIdx(int refList) const { return m_symRefIdx[refList]; }

  bool isIntra() const { return m_eSliceType == I_SLICE; }
  bool isInterB() const { return m_eSliceType == B_SLICE; }
  bool isInterP() const { return m_eSliceType == P_SLICE; }

  uint32_t getCuQpDeltaSubdiv() const {
    return this->isIntra() ? m_pcPicHeader->getCuQpDeltaSubdivIntra()
                           : m_pcPicHeader->getCuQpDeltaSubdivInter();
  }
  uint32_t getCuChromaQpOffsetSubdiv() const {
    return this->isIntra() ? m_pcPicHeader->getCuChromaQpOffsetSubdivIntra()
                           : m_pcPicHeader->getCuChromaQpOffsetSubdivInter();
  }
  uint32_t getIndependentSliceIdx() const { return m_independentSliceIdx; }
  const WeightPredTable &getWeightPredTable() const {return m_weightPredTable; }
  WPScalingParam *getWpScaling(const RefPicList refPicList, const int refIdx);
  const WPScalingParam *getWpScaling(const RefPicList refPicList,
                                     const int refIdx) const;
  bool getCabacInitFlag() const {
    return m_cabacInitFlag;
  } //!< get CABAC initial flag

  SliceType getEncCABACTableIdx() const { return m_encCABACTableIdx; }

  bool getAlfEnabledFlag(ComponentID compId) const {
    return m_alfEnabledFlag[compId];
  }
  int getNumAlfApsIdsLuma() const { return m_numAlfApsIdsLuma; }
  int getAlfApsIdChroma() const { return m_alfApsIdChroma; }

  CcAlfFilterParam m_ccAlfFilterParam;
  std::array<uint8_t *, 2> m_ccAlfFilterControl;
  int get_tsrc_index() const { return m_tsrc_index; }
  void setRiceBit(int idx, int i) { m_riceBit[idx] = i; }
  unsigned getRiceBit(int idx) const { return m_riceBit[idx]; }
}; // END CLASS DEFINITION Slice

class PreCalcValues {
public:
  const ChromaFormat chrFormat;
  const bool multiBlock422;
  const unsigned maxCUWidth;
  const unsigned maxCUHeight;
  // to get CTU position, use (x & maxCUWidthMask) rather than (x % maxCUWidth)
  const unsigned maxCUWidthMask;
  const unsigned maxCUHeightMask;
  const unsigned maxCUWidthLog2;
  const unsigned maxCUHeightLog2;
  const unsigned widthInCtus;
  const unsigned sizeInCtus;
  const bool noChroma2x2;
  const bool ISingleTree;

private:
  std::array<unsigned, 3> maxBtDepth;
  std::array<unsigned, 3> minBtSize;
  std::array<unsigned, 3> maxBtSize;
  std::array<unsigned, 3> minTtSize;
  std::array<unsigned, 3> maxTtSize;
  std::array<unsigned, 3> minQtSize;

  unsigned getValIdx(const Slice &slice, const ChannelType chType) const;

public:
  PreCalcValues(const ChromaFormat _chrFormat, const bool _multiBlock422,
                const unsigned _maxCUWidth, const unsigned _maxCUHeight,
                const unsigned _maxCUWidthMask, const unsigned _maxCUHeightMask,
                const unsigned _maxCUWidthLog2, const unsigned _maxCUHeightLog2,
                const unsigned _widthInCtus, const unsigned _sizeInCtus,
                const bool _noChroma2x2, const bool _ISingleTree,
                const unsigned *_maxBtDepth, const unsigned *_minBtSize,
                const unsigned *_maxBtSize, const unsigned *_minTtSize,
                const unsigned *_maxTtSize, const unsigned *_minQtSize)
      : chrFormat(_chrFormat), multiBlock422(_multiBlock422),
        maxCUWidth(_maxCUWidth), maxCUHeight(_maxCUHeight),
        maxCUWidthMask(_maxCUWidthMask), maxCUHeightMask(_maxCUHeightMask),
        maxCUWidthLog2(_maxCUWidthLog2), maxCUHeightLog2(_maxCUHeightLog2),
        widthInCtus(_widthInCtus), sizeInCtus(_sizeInCtus),
        noChroma2x2(_noChroma2x2), ISingleTree(_ISingleTree) {
    copy_array(_maxBtDepth, maxBtDepth);
    copy_array(_minBtSize, minBtSize);
    copy_array(_maxBtSize, maxBtSize);
    copy_array(_minTtSize, minTtSize);
    copy_array(_maxTtSize, maxTtSize);
    copy_array(_minQtSize, minQtSize);
  }

  unsigned getMaxBtDepth(const Slice &slice, const ChannelType chType) const;
  unsigned getMinBtSize(const Slice &slice, const ChannelType chType) const;
  unsigned getMaxBtSize(const Slice &slice, const ChannelType chType) const;
  unsigned getMinTtSize(const Slice &slice, const ChannelType chType) const;
  unsigned getMaxTtSize(const Slice &slice, const ChannelType chType) const;
  unsigned getMinQtSize(const Slice &slice, const ChannelType chType) const;
};
} // namespace Common

#endif // COMMON_SLICE
