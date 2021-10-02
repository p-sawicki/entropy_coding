#include "common_def.hpp"
#include "mv.hpp"

/// class for motion vector with reference index
struct MvField {
  Mv mv;
  int16_t refIdx;

  MvField() : refIdx(NOT_VALID) {}
  MvField(Mv const &cMv, const int iRefIdx) : mv(cMv), refIdx(iRefIdx) {}

  void setMvField(Mv const &cMv, const int iRefIdx) {
    CHECK(iRefIdx == -1 && cMv != Mv(0, 0), "Must not happen.");
    mv = cMv;
    refIdx = iRefIdx;
  }

  bool operator==(const MvField &other) const {
    CHECK(refIdx == -1 && mv != Mv(0, 0), "Error in operator== of MvField.");
    CHECK(other.refIdx == -1 && other.mv != Mv(0, 0),
          "Error in operator== of MvField.");
    return refIdx == other.refIdx && mv == other.mv;
  }
  bool operator!=(const MvField &other) const {
    CHECK(refIdx == -1 && mv != Mv(0, 0), "Error in operator!= of MvField.");
    CHECK(other.refIdx == -1 && other.mv != Mv(0, 0),
          "Error in operator!= of MvField.");
    return refIdx != other.refIdx || mv != other.mv;
  }
};

struct MotionInfo {
  bool isInter;
  bool isIBCmot;
  char interDir;
  bool useAltHpelIf;
  uint16_t sliceIdx;
  Mv mv[NUM_REF_PIC_LIST_01];
  int16_t refIdx[NUM_REF_PIC_LIST_01];
  uint8_t BcwIdx;
  Mv bv;
#if GDR_ENABLED
  bool sourceClean;   // source Position is clean/dirty
  Position sourcePos; // source Position of Mv
#endif

  MotionInfo()
      : isInter(false), isIBCmot(false), interDir(0), useAltHpelIf(false),
        sliceIdx(0), refIdx{NOT_VALID, NOT_VALID}, BcwIdx(0) {}
  // ensure that MotionInfo(0) produces '\x000....' bit pattern - needed to work
  // with AreaBuf - don't use this constructor for anything else
  MotionInfo(int i)
      : isInter(i != 0), isIBCmot(false), interDir(0), useAltHpelIf(false),
        sliceIdx(0), refIdx{0, 0}, BcwIdx(0) {
    CHECKD(i != 0, "The argument for this constructor has to be '0'");
  }

  bool operator==(const MotionInfo &mi) const {
    if (isInter != mi.isInter)
      return false;
    if (isIBCmot != mi.isIBCmot)
      return false;
    if (isInter) {
      if (sliceIdx != mi.sliceIdx)
        return false;
      if (interDir != mi.interDir)
        return false;

      if (interDir != 2) {
        if (refIdx[0] != mi.refIdx[0])
          return false;
        if (mv[0] != mi.mv[0])
          return false;
      }

      if (interDir != 1) {
        if (refIdx[1] != mi.refIdx[1])
          return false;
        if (mv[1] != mi.mv[1])
          return false;
      }
    }

    return true;
  }

  bool operator!=(const MotionInfo &mi) const { return !(*this == mi); }
};

struct LutMotionCand {
  static_vector<MotionInfo, MAX_NUM_HMVP_CANDS> lut;
  static_vector<MotionInfo, MAX_NUM_HMVP_CANDS> lutIbc;
};

/// parameters for AMVP
struct AMVPInfo {
  Mv mvCand[AMVP_MAX_NUM_CANDS_MEM]; ///< array of motion vector predictor
                                     ///< candidates
  unsigned numCand; ///< number of motion vector predictor candidates
#if GDR_ENABLED
  bool allCandSolidInAbove;
  bool mvSolid[AMVP_MAX_NUM_CANDS_MEM];
  bool mvValid[AMVP_MAX_NUM_CANDS_MEM];

  Position mvPos[AMVP_MAX_NUM_CANDS_MEM];
  MvpType mvType[AMVP_MAX_NUM_CANDS_MEM];
#endif
};

struct AffineAMVPInfo {
  Mv mvCandLT[AMVP_MAX_NUM_CANDS_MEM]; ///< array of affine motion vector
                                       ///< predictor candidates for left-top
                                       ///< corner
  Mv mvCandRT[AMVP_MAX_NUM_CANDS_MEM]; ///< array of affine motion vector
                                       ///< predictor candidates for right-top
                                       ///< corner
  Mv mvCandLB[AMVP_MAX_NUM_CANDS_MEM]; ///< array of affine motion vector
                                       ///< predictor candidates for left-bottom
                                       ///< corner
  unsigned numCand; ///< number of motion vector predictor candidates
#if GDR_ENABLED
  bool allCandSolidInAbove;

  bool mvSolidLT[AMVP_MAX_NUM_CANDS_MEM];
  bool mvSolidRT[AMVP_MAX_NUM_CANDS_MEM];
  bool mvSolidLB[AMVP_MAX_NUM_CANDS_MEM];

  bool mvValidLT[AMVP_MAX_NUM_CANDS_MEM];
  bool mvValidRT[AMVP_MAX_NUM_CANDS_MEM];
  bool mvValidLB[AMVP_MAX_NUM_CANDS_MEM];

  MvpType mvTypeLT[AMVP_MAX_NUM_CANDS_MEM];
  MvpType mvTypeRT[AMVP_MAX_NUM_CANDS_MEM];
  MvpType mvTypeLB[AMVP_MAX_NUM_CANDS_MEM];

  Position mvPosLT[AMVP_MAX_NUM_CANDS_MEM];
  Position mvPosRT[AMVP_MAX_NUM_CANDS_MEM];
  Position mvPosLB[AMVP_MAX_NUM_CANDS_MEM];
#endif
};