#ifndef ENTROPY_CODEC_ARITH_CODEC
#define ENTROPY_CODEC_ARITH_CODEC

#include "bit_stream.hpp"
#include "contexts.hpp"

namespace EntropyCoding {

class BinStore {
public:
  BinStore();
  ~BinStore() {}

  void reset();
  void addBin(unsigned bin, unsigned ctxId);

  void setUse(bool useStore);
  bool inUse() const;
  const std::vector<bool> &getBinVector(unsigned ctxId) const;

private:
  void xCheckAlloc();

private:
  static const std::size_t m_maxNumBins = 100000;
  bool m_inUse;
  bool m_allocated;
  std::vector<std::vector<bool>> m_binBuffer;
};

class BinEncIf : public Common::Ctx {
protected:
  template <class BinProbModel>
  BinEncIf(const BinProbModel *dummy) : Common::Ctx(dummy) {}

public:
  virtual ~BinEncIf() = default;

public:
  virtual void init(Common::OutputBitstream *bitstream) = 0;
  virtual void uninit() = 0;
  virtual void start() = 0;
  virtual void finish() = 0;
  virtual void restart() = 0;
  virtual void reset(int qp, int initId) = 0;

public:
  virtual void resetBits() = 0;
  virtual uint64_t getEstFracBits() const = 0;
  virtual unsigned getNumBins(unsigned ctxId) const = 0;

public:
  virtual void encodeBin(unsigned bin, unsigned ctxId) = 0;
  virtual void encodeBinEP(unsigned bin) = 0;
  virtual void encodeBinsEP(unsigned bins, unsigned numBins) = 0;
  virtual void encodeRemAbsEP(unsigned bins, unsigned goRicePar,
                              unsigned cutoff, int maxLog2TrDynamicRange) = 0;
  virtual void encodeBinTrm(unsigned bin) = 0;
  virtual void align() = 0;

public:
  virtual uint32_t getNumBins() = 0;
  virtual bool isEncoding() = 0;
  virtual unsigned getNumWrittenBits() = 0;

public:
  virtual void setBinStorage(bool b) = 0;
  virtual const BinStore *getBinStore() const = 0;
  virtual BinEncIf *getTestBinEncoder() const = 0;
};

class BinCounter {
public:
  BinCounter();
  ~BinCounter() = default;

public:
  void reset();
  void addCtx(unsigned ctxId);
  void addEP(unsigned num);
  void addEP();
  void addTrm();
  uint32_t getAll() const;
  uint32_t getCtx(unsigned ctxId) const;
  uint32_t getEP() const;
  uint32_t getTrm() const;

private:
  std::vector<uint32_t> m_CtxBinsCodedBuffer;
  uint32_t *m_NumBinsCtx;
  uint32_t m_NumBinsEP;
  uint32_t m_NumBinsTrm;
};

class BinEncoderBase : public BinEncIf, public BinCounter {
protected:
  template <class BinProbModel> BinEncoderBase(const BinProbModel *dummy);

public:
  ~BinEncoderBase() {}

public:
  void init(Common::OutputBitstream *bitstream);
  void uninit();
  void start();
  void finish();
  void restart();
  void reset(int qp, int initId);
#if JVET_V0106_RRC_RICE
  void riceStatReset(int bitDepth);
#endif
public:
  void resetBits();
  uint64_t getEstFracBits() const;
  unsigned getNumBins(unsigned ctxId) const;

public:
  void encodeBinEP(unsigned bin);
  void encodeBinsEP(unsigned bins, unsigned numBins);
  void encodeRemAbsEP(unsigned bins, unsigned goRicePar, unsigned cutoff,
                      int maxLog2TrDynamicRange);
  void encodeBinTrm(unsigned bin);
  void align();
  unsigned getNumWrittenBits();

public:
  uint32_t getNumBins();
  bool isEncoding();

protected:
  void encodeAlignedBinsEP(unsigned bins, unsigned numBins);
  void writeOut();

protected:
  Common::OutputBitstream *m_Bitstream;
  uint32_t m_Low;
  uint32_t m_Range;
  uint32_t m_bufferedByte;
  int32_t m_numBufferedBytes;
  int32_t m_bitsLeft;
  BinStore m_BinStore;
};

template <class BinProbModel> class TBinEncoder : public BinEncoderBase {
public:
  TBinEncoder();
  ~TBinEncoder() = default;
  void encodeBin(unsigned bin, unsigned ctxId);

public:
  void setBinStorage(bool b);
  const BinStore *getBinStore() const;
  BinEncIf *getTestBinEncoder() const;

private:
  Common::CtxStore<BinProbModel> &m_Ctx;
};

class BitEstimatorBase : public BinEncIf {
protected:
  template <class BinProbModel> BitEstimatorBase(const BinProbModel *dummy);

public:
  ~BitEstimatorBase() = default;

public:
  void init(Common::OutputBitstream *bitstream);
  void uninit();
  void start();
  void finish();
  void restart();
  void reset(int qp, int initId);

public:
  void resetBits();

  uint64_t getEstFracBits() const;
  unsigned getNumBins(unsigned ctxId) const;

public:
  void encodeBinEP(unsigned bin);
  void encodeBinsEP(unsigned bins, unsigned numBins);
  void encodeRemAbsEP(unsigned bins, unsigned goRicePar, unsigned cutoff,
                      int maxLog2TrDynamicRange);
  void align();

public:
  uint32_t getNumBins();
  bool isEncoding();
  unsigned getNumWrittenBits();

protected:
  uint64_t m_EstFracBits;
};

template <class BinProbModel> class TBitEstimator : public BitEstimatorBase {
public:
  TBitEstimator();
  ~TBitEstimator() = default;

  void encodeBin(unsigned bin, unsigned ctxId);
  void encodeBinTrm(unsigned bin);
  void setBinStorage(bool b);
  const BinStore *getBinStore() const;
  BinEncIf *getTestBinEncoder() const;

private:
  Common::CtxStore<BinProbModel> &m_Ctx;
};

typedef TBinEncoder<Common::BinProbModel_Std> BinEncoder_Std;

typedef TBitEstimator<Common::BinProbModel_Std> BitEstimator_Std;

class BinDecoderBase : public Common::Ctx {
protected:
  template <class BinProbModel> BinDecoderBase(const BinProbModel *dummy);

public:
  ~BinDecoderBase() = default;

public:
  void init(Common::InputBitstream *bitstream);
  void uninit();
  void start();
  void finish();
  void reset(int qp, int initId);
#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
  void riceStatReset(int bitDepth, bool persistentRiceAdaptationEnabledFlag);
#else
  void riceStatReset(int bitDepth);
#endif
#if RExt__DECODER_DEBUG_BIT_STATISTICS
  void set(const CodingStatisticsClassType &type);
#endif

public:
  virtual unsigned decodeBin(unsigned ctxId) = 0;

public:
  unsigned decodeBinEP();
  unsigned decodeBinsEP(unsigned numBins);
  unsigned decodeRemAbsEP(unsigned goRicePar, unsigned cutoff,
                          int maxLog2TrDynamicRange);
  unsigned decodeBinTrm();
  void align();
  unsigned getNumBitsRead();

private:
  unsigned decodeAlignedBinsEP(unsigned numBins);

protected:
  Common::InputBitstream *m_Bitstream;
  uint32_t m_Range;
  uint32_t m_Value;
  int32_t m_bitsNeeded;
#if RExt__DECODER_DEBUG_BIT_STATISTICS
  const CodingStatisticsClassType *ptype;
#endif
};

template <class BinProbModel> class TBinDecoder : public BinDecoderBase {
public:
  TBinDecoder();
  ~TBinDecoder() = default;
  unsigned decodeBin(unsigned ctxId);

private:
  Common::CtxStore<BinProbModel> &m_Ctx;
};

typedef TBinDecoder<Common::BinProbModel_Std> BinDecoder_Std;
} // namespace EntropyCoding

#endif // ENTROPY_CODEC_ARITH_CODEC