#ifndef ENTROPY_CODEC_CONTEXTS
#define ENTROPY_CODEC_CONTEXTS

#include <cstdint>
#include <initializer_list>
#include <vector>

#include "common_def.hpp"

namespace EntropyCoding {

static constexpr int PROB_BITS =
    15; // Nominal number of bits to represent probabilities
static constexpr int PROB_BITS_0 =
    10; // Number of bits to represent 1st estimate
static constexpr int PROB_BITS_1 =
    14; // Number of bits to represent 2nd estimate
static constexpr int MASK_0 = ~(~0u << PROB_BITS_0)
                              << (PROB_BITS - PROB_BITS_0);
static constexpr int MASK_1 = ~(~0u << PROB_BITS_1)
                              << (PROB_BITS - PROB_BITS_1);
static constexpr uint8_t DWS = 8; // 0x47 Default window sizes

enum BPMType { BPM_Undefined = 0, BPM_Std, BPM_NUM };

class CtxSet {
public:
  CtxSet(uint16_t offset, uint16_t size);
  CtxSet(const CtxSet &ctxSet);
  CtxSet(::std::initializer_list<CtxSet> ctxSets);

public:
  uint16_t operator()() const;
  uint16_t operator()(uint16_t inc) const;
  bool operator==(const CtxSet &ctxSet) const;
  bool operator!=(const CtxSet &ctxSet) const;

public:
  uint16_t Offset;
  uint16_t Size;
};

class ContextSetCfg {
public:
  // context sets: specify offset and size
  static const CtxSet SplitFlag;
  static const CtxSet SplitQtFlag;
  static const CtxSet SplitHvFlag;
  static const CtxSet Split12Flag;
  static const CtxSet ModeConsFlag;
  static const CtxSet SkipFlag;
  static const CtxSet MergeFlag;
  static const CtxSet RegularMergeFlag;
  static const CtxSet MergeIdx;
  static const CtxSet PredMode;
  static const CtxSet MultiRefLineIdx;
  static const CtxSet IntraLumaMpmFlag;
  static const CtxSet IntraLumaPlanarFlag;
  static const CtxSet CclmModeFlag;
  static const CtxSet CclmModeIdx;
  static const CtxSet IntraChromaPredMode;
  static const CtxSet MipFlag;
  static const CtxSet DeltaQP;
  static const CtxSet InterDir;
  static const CtxSet RefPic;
  static const CtxSet MmvdFlag;
  static const CtxSet MmvdMergeIdx;
  static const CtxSet MmvdStepMvpIdx;
  static const CtxSet SubblockMergeFlag;
  static const CtxSet AffineFlag;
  static const CtxSet AffineType;
  static const CtxSet AffMergeIdx;
  static const CtxSet Mvd;
  static const CtxSet BDPCMMode;
  static const CtxSet QtRootCbf;
  static const CtxSet ACTFlag;
  static const CtxSet QtCbf[3];         // [ channel ]
  static const CtxSet SigCoeffGroup[2]; // [ ChannelType ]
  static const CtxSet LastX[2];         // [ ChannelType ]
  static const CtxSet LastY[2];         // [ ChannelType ]
  static const CtxSet SigFlag[6];       // [ ChannelType + State ]
  static const CtxSet ParFlag[2];       // [ ChannelType ]
  static const CtxSet GtxFlag[4];       // [ ChannelType + x ]
  static const CtxSet TsSigCoeffGroup;
  static const CtxSet TsSigFlag;
  static const CtxSet TsParFlag;
  static const CtxSet TsGtxFlag;
  static const CtxSet TsLrg1Flag;
  static const CtxSet TsResidualSign;
  static const CtxSet MVPIdx;
  static const CtxSet SaoMergeFlag;
  static const CtxSet SaoTypeIdx;
  static const CtxSet TransformSkipFlag;
  static const CtxSet MTSIdx;
  static const CtxSet LFNSTIdx;
  static const CtxSet PLTFlag;
  static const CtxSet RotationFlag;
  static const CtxSet RunTypeFlag;
  static const CtxSet IdxRunModel;
  static const CtxSet CopyRunModel;
  static const CtxSet SbtFlag;
  static const CtxSet SbtQuadFlag;
  static const CtxSet SbtHorFlag;
  static const CtxSet SbtPosFlag;
  static const CtxSet ChromaQpAdjFlag;
  static const CtxSet ChromaQpAdjIdc;
  static const CtxSet ImvFlag;
  static const CtxSet BcwIdx;
  static const CtxSet ctbAlfFlag;
  static const CtxSet ctbAlfAlternative;
  static const CtxSet AlfUseTemporalFilt;
  static const CtxSet CcAlfFilterControlFlag;
  static const CtxSet CiipFlag;
  static const CtxSet SmvdFlag;
  static const CtxSet IBCFlag;
  static const CtxSet ISPMode;
  static const CtxSet JointCbCrFlag;
  static const unsigned NumberOfContexts;

  // combined sets for less complex copying
  // NOTE: The contained CtxSet's should directly follow each other in the
  // initalization list;
  //       otherwise, you will copy more elements than you want !!!
  static const CtxSet Sao;
  static const CtxSet Alf;
  static const CtxSet Palette;

public:
  static const ::std::vector<uint8_t> &getInitTable(unsigned initId);

private:
  static ::std::vector<::std::vector<uint8_t>> sm_InitTables;
  static CtxSet addCtxSet(
      ::std::initializer_list<::std::initializer_list<uint8_t>> initSet2d);
};

struct BinFracBits {
  uint32_t intBits[2];
};

class ProbModelTables {
protected:
  static const BinFracBits m_binFracBits[256];
  static const uint8_t m_RenormTable_32[32]; // Std         MP   MPI
};

class BinProbModelBase : public ProbModelTables {
public:
  BinProbModelBase() = default;
  ~BinProbModelBase() = default;
  static uint32_t estFracBitsEP();
  static uint32_t estFracBitsEP(unsigned numBins);
};

class BinProbModel_Std : public BinProbModelBase {
public:
  BinProbModel_Std();
  BinProbModel_Std(const uint16_t state0, const uint16_t state1,
                   const uint8_t rate)
      : m_state{state0, state1}, m_rate(rate) {}
  ~BinProbModel_Std() = default;

public:
  void init(int qp, int initId);
  void update(unsigned bin);
  void setLog2WindowSize(uint8_t log2WindowSize);
  void estFracBitsUpdate(unsigned bin, uint64_t &b);
  uint32_t estFracBits(unsigned bin) const;
  static uint32_t estFracBitsTrm(unsigned bin);
  BinFracBits getFracBitsArray() const;

public:
  uint8_t state() const;
  uint8_t mps() const;
  uint8_t getLPS(unsigned range) const;
  static uint8_t getRenormBitsLPS(unsigned LPS);
  static uint8_t getRenormBitsRange(unsigned);
  uint16_t getState() const;
  void setState(uint16_t pState);

public:
  uint64_t estFracExcessBits(const BinProbModel_Std &r) const;

  uint16_t getState0() const { return m_state[0]; }
  uint16_t getState1() const { return m_state[1]; }
  uint8_t getRate() const { return m_rate; }

private:
  uint16_t m_state[2];
  uint8_t m_rate;
};

class FracBitsAccess {
public:
  virtual BinFracBits getFracBitsArray(unsigned ctxId) const = 0;
};

template <class BinProbModel> class CtxStore : public FracBitsAccess {
public:
  CtxStore();
  CtxStore(bool dummy);
  CtxStore(const CtxStore<BinProbModel> &ctxStore);
  CtxStore(const ::std::vector<BinProbModel> &ctxBuffer)
      : m_CtxBuffer(ctxBuffer), m_Ctx(m_CtxBuffer.data()) {}

public:
  void copyFrom(const CtxStore<BinProbModel> &src);
  void copyFrom(const CtxStore<BinProbModel> &src, const CtxSet &ctxSet);
  void init(int qp, int initId);
  void setWinSizes(const ::std::vector<uint8_t> &log2WindowSizes);
  void loadPStates(const ::std::vector<uint16_t> &probStates);
  void savePStates(::std::vector<uint16_t> &probStates) const;

  const BinProbModel &operator[](unsigned ctxId) const;
  BinProbModel &operator[](unsigned ctxId);
  uint32_t estFracBits(unsigned bin, unsigned ctxId) const;

  BinFracBits getFracBitsArray(unsigned ctxId) const;
  const ::std::vector<BinProbModel> &getCtxBuffer() const {
    return m_CtxBuffer;
  }

private:
  void checkInit();

private:
  ::std::vector<BinProbModel> m_CtxBuffer;
  BinProbModel *m_Ctx;
};

class Ctx;
class SubCtx {
  friend class Ctx;

public:
  SubCtx(const CtxSet &ctxSet, const Ctx &ctx);
  SubCtx(const SubCtx &subCtx);
  const SubCtx &operator=(const SubCtx &) = delete;

private:
  const CtxSet m_CtxSet;
  const Ctx &m_Ctx;
};

class Ctx : public ContextSetCfg {
public:
  Ctx();
  Ctx(const BinProbModel_Std *dummy);
  Ctx(const Ctx &ctx);
  Ctx(const BPMType bpmType, const CtxStore<BinProbModel_Std> &ctxStore)
      : m_BPMType(bpmType), m_CtxStore_Std(ctxStore) {}

public:
  const Ctx &operator=(const Ctx &ctx);

  SubCtx operator=(SubCtx &&subCtx);

  void init(int qp, int initId);
  void init(const Ctx& ctx, const CtxSet &ctxSet);

#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
  void riceStatReset(int bitDepth, bool persistentRiceAdaptationEnabledFlag);
#else
  void riceStatReset(int bitDepth);
#endif

  void loadPStates(const ::std::vector<uint16_t> &probStates);

  void savePStates(::std::vector<uint16_t> &probStates) const;

  void initCtxAndWinSize(unsigned ctxId, const Ctx &ctx, const uint8_t winSize);

  const unsigned &getGRAdaptStats(unsigned id) const;
  unsigned &getGRAdaptStats(unsigned id);

  const unsigned getBaseLevel() const;
  void setBaseLevel(int value);

public:
  unsigned getBPMType() const;
  const Ctx &getCtx() const;
  Ctx &getCtx();

  explicit operator const CtxStore<BinProbModel_Std> &() const;
  explicit operator CtxStore<BinProbModel_Std> &();

  const FracBitsAccess &getFracBitsAcess() const;

private:
  BPMType m_BPMType;
  CtxStore<BinProbModel_Std> m_CtxStore_Std;

protected:
  unsigned m_GRAdaptStats[RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS];
#if JVET_V0106_RRC_RICE
  int m_baseLevel;
#endif
};
} // namespace EntropyCoding

#endif // ENTROPY_CODEC_CONTEXTS