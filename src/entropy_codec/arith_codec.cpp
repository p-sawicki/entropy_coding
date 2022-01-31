#include "arith_codec.hpp"

using namespace EntropyCoding;
using namespace Common;

#define CNT_OFFSET 0

BinStore::BinStore() : m_inUse(false), m_allocated(false) {}

void BinStore::reset() {
  if (m_inUse) {
    for (unsigned n = 0; n < Ctx::NumberOfContexts; n++) {
      m_binBuffer[n].clear();
    }
  }
}

void BinStore::addBin(unsigned bin, unsigned ctxId) {
  if (m_inUse) {
    std::vector<bool> &binBuffer = m_binBuffer[ctxId];
    if (binBuffer.size() < m_maxNumBins) {
      binBuffer.push_back(bin == 1);
    }
  }
}

void BinStore::setUse(bool useStore) {
  m_inUse = useStore;
  if (m_inUse) {
    xCheckAlloc();
  }
}

bool BinStore::inUse() const { return m_inUse; }

const std::vector<bool> &BinStore::getBinVector(unsigned ctxId) const {
  return m_binBuffer[ctxId];
}

void BinStore::xCheckAlloc() {
  if (!m_allocated) {
    m_binBuffer.resize(Ctx::NumberOfContexts);
    for (unsigned n = 0; n < Ctx::NumberOfContexts; n++) {
      m_binBuffer[n].reserve(m_maxNumBins);
    }
    m_allocated = true;
  }
}

template <class BinProbModel>
BinDecoderBase::BinDecoderBase(const BinProbModel *dummy)
    : Ctx(dummy), m_Bitstream(0), m_Range(0), m_Value(0), m_bitsNeeded(0) {}

void BinDecoderBase::init(InputBitstream *bitstream) {
  m_Bitstream = bitstream;
}

void BinDecoderBase::uninit() { m_Bitstream = 0; }

void BinDecoderBase::start() {
  CHECK(m_Bitstream->getNumBitsUntilByteAligned(),
        "Bitstream is not byte aligned.");
  m_Range = 510;
  m_Value = (m_Bitstream->readByte() << 8) + m_Bitstream->readByte();
  m_bitsNeeded = -8;
}

void BinDecoderBase::finish() {
  unsigned lastByte;
  m_Bitstream->peekPreviousByte(lastByte);
  CHECK(((lastByte << (8 + m_bitsNeeded)) & 0xff) != 0x80,
        "No proper stop/alignment pattern at end of CABAC stream.");
}

void BinDecoderBase::reset(int qp, int initId) {
  Ctx::init(qp, initId);
  start();
}

#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
void BinDecoderBase::riceStatReset(int bitDepth,
                                   bool persistentRiceAdaptationEnabledFlag)
#else
void BinDecoderBase::riceStatReset(int bitDepth)
#endif
{
#if JVET_W0178_CONSTRAINTS_ON_REXT_TOOLS
  Ctx::riceStatReset(bitDepth, persistentRiceAdaptationEnabledFlag);
#else
  Ctx::riceStatReset(bitDepth);
#endif
}

#if RExt__DECODER_DEBUG_BIT_STATISTICS
void BinDecoderBase::set(const CodingStatisticsClassType &type) {
  ptype = &type;
}
#endif

unsigned BinDecoderBase::decodeBinEP() {
  m_Value += m_Value;
  if (++m_bitsNeeded >= 0) {
    m_Value += m_Bitstream->readByte();
    m_bitsNeeded = -8;
  }

  unsigned bin = 0;
  unsigned SR = m_Range << 7;
  if (m_Value >= SR) {
    m_Value -= SR;
    bin = 1;
  }
  return bin;
}

unsigned BinDecoderBase::decodeBinsEP(unsigned numBins) {
  if (m_Range == 256) {
    return decodeAlignedBinsEP(numBins);
  }
  unsigned remBins = numBins;
  unsigned bins = 0;
  while (remBins > 8) {
    m_Value = (m_Value << 8) + (m_Bitstream->readByte() << (8 + m_bitsNeeded));
    unsigned SR = m_Range << 15;
    for (int i = 0; i < 8; i++) {
      bins += bins;
      SR >>= 1;
      if (m_Value >= SR) {
        bins++;
        m_Value -= SR;
      }
    }
    remBins -= 8;
  }
  m_bitsNeeded += remBins;
  m_Value <<= remBins;
  if (m_bitsNeeded >= 0) {
    m_Value += m_Bitstream->readByte() << m_bitsNeeded;
    m_bitsNeeded -= 8;
  }
  unsigned SR = m_Range << (remBins + 7);
  for (unsigned int i = 0; i < remBins; i++) {
    bins += bins;
    SR >>= 1;
    if (m_Value >= SR) {
      bins++;
      m_Value -= SR;
    }
  }
  return bins;
}

unsigned BinDecoderBase::decodeRemAbsEP(unsigned goRicePar, unsigned cutoff,
                                        int maxLog2TrDynamicRange) {
  unsigned prefix = 0;
  {
    const unsigned maxPrefix = 32 - maxLog2TrDynamicRange;
    unsigned codeWord = 0;
    do {
      prefix++;
      codeWord = decodeBinEP();
    } while (codeWord && prefix < maxPrefix);
    prefix -= 1 - codeWord;
  }

  unsigned length = goRicePar, offset;
  if (prefix < cutoff) {
    offset = prefix << goRicePar;
  } else {
    offset = (((1 << (prefix - cutoff)) + cutoff - 1) << goRicePar);
    {
      length +=
          (prefix == (32u - static_cast<unsigned int>(maxLog2TrDynamicRange))
               ? maxLog2TrDynamicRange - goRicePar
               : prefix - cutoff);
    }
  }
  return offset + decodeBinsEP(length);
}

unsigned BinDecoderBase::decodeBinTrm() {
  m_Range -= 2;
  unsigned SR = m_Range << 7;
  if (m_Value >= SR) {
    return 1;
  } else {
    if (m_Range < 256) {
      m_Range += m_Range;
      m_Value += m_Value;
      if (++m_bitsNeeded == 0) {
        m_Value += m_Bitstream->readByte();
        m_bitsNeeded = -8;
      }
    }
    return 0;
  }
}

void BinDecoderBase::align() { m_Range = 256; }

unsigned BinDecoderBase::getNumBitsRead() {
  return m_Bitstream->getNumBitsRead() + m_bitsNeeded;
}

unsigned BinDecoderBase::decodeAlignedBinsEP(unsigned numBins) {
  unsigned remBins = numBins;
  unsigned bins = 0;
  while (remBins > 0) {
    // The MSB of m_Value is known to be 0 because range is 256. Therefore:
    //   > The comparison against the symbol range of 128 is simply a test on
    //   the next-most-significant bit > "Subtracting" the symbol range if the
    //   decoded bin is 1 simply involves clearing that bit.
    //  As a result, the required bins are simply the <binsToRead>
    //  next-most-significant bits of m_Value (m_Value is stored MSB-aligned in
    //  a 16-bit buffer - hence the shift of 15)
    //
    //    m_Value = |0|V|V|V|V|V|V|V|V|B|B|B|B|B|B|B|        (V = usable bit, B
    //    = potential buffered bit (buffer refills when m_bitsNeeded >= 0))
    //
    unsigned binsToRead =
        std::min<unsigned>(remBins, 8); // read bytes if able to take advantage
                                        // of the system's byte-read function
    unsigned binMask = (1 << binsToRead) - 1;
    unsigned newBins = (m_Value >> (15 - binsToRead)) & binMask;
    bins = (bins << binsToRead) | newBins;
    m_Value = (m_Value << binsToRead) & 0x7FFF;
    remBins -= binsToRead;
    m_bitsNeeded += binsToRead;
    if (m_bitsNeeded >= 0) {
      m_Value |= m_Bitstream->readByte() << m_bitsNeeded;
      m_bitsNeeded -= 8;
    }
  }
  return bins;
}

template <class BinProbModel>
TBinDecoder<BinProbModel>::TBinDecoder()
    : BinDecoderBase(static_cast<const BinProbModel *>(nullptr)),
      m_Ctx(static_cast<CtxStore<BinProbModel> &>(*this)) {}

template <class BinProbModel>
unsigned TBinDecoder<BinProbModel>::decodeBin(unsigned ctxId) {
  BinProbModel &rcProbModel = m_Ctx[ctxId];
  unsigned bin = rcProbModel.mps();
  uint32_t LPS = rcProbModel.getLPS(m_Range);

  m_Range -= LPS;
  uint32_t SR = m_Range << 7;
  if (m_Value < SR) {
    // MPS path
    if (m_Range < 256) {
      int numBits = rcProbModel.getRenormBitsRange(m_Range);
      m_Range <<= numBits;
      m_Value <<= numBits;
      m_bitsNeeded += numBits;
      if (m_bitsNeeded >= 0) {
        m_Value += m_Bitstream->readByte() << m_bitsNeeded;
        m_bitsNeeded -= 8;
      }
    }
  } else {
    bin = 1 - bin;
    // LPS path
    int numBits = rcProbModel.getRenormBitsLPS(LPS);
    m_Value -= SR;
    m_Value = m_Value << numBits;
    m_Range = LPS << numBits;
    m_bitsNeeded += numBits;
    if (m_bitsNeeded >= 0) {
      m_Value += m_Bitstream->readByte() << m_bitsNeeded;
      m_bitsNeeded -= 8;
    }
  }
  rcProbModel.update(bin);
  return bin;
}

template class TBinDecoder<BinProbModel_Std>;

BinCounter::BinCounter()
    : m_CtxBinsCodedBuffer(Ctx::NumberOfContexts),
      m_NumBinsCtx(m_CtxBinsCodedBuffer.data()), m_NumBinsEP(0),
      m_NumBinsTrm(0) {}

void BinCounter::reset() {
  for (std::size_t k = 0; k < m_CtxBinsCodedBuffer.size(); k++) {
    m_NumBinsCtx[k] = 0;
  }
  m_NumBinsEP = 0;
  m_NumBinsTrm = 0;
}

void BinCounter::addCtx(unsigned ctxId) { m_NumBinsCtx[ctxId]++; }

void BinCounter::addEP(unsigned num) { m_NumBinsEP += num; }

void BinCounter::addEP() { m_NumBinsEP++; }

void BinCounter::addTrm() { m_NumBinsTrm++; }

uint32_t BinCounter::getAll() const {
  uint32_t count = m_NumBinsEP + m_NumBinsTrm;
  for (std::size_t k = 0; k < m_CtxBinsCodedBuffer.size(); k++) {
    count += m_NumBinsCtx[k];
  }
  return count;
}

uint32_t BinCounter::getCtx(unsigned ctxId) const {
  return m_NumBinsCtx[ctxId];
}

uint32_t BinCounter::getEP() const { return m_NumBinsEP; }

uint32_t BinCounter::getTrm() const { return m_NumBinsTrm; }

template <class BinProbModel>
BinEncoderBase::BinEncoderBase(const BinProbModel *dummy)
    : BinEncIf(dummy), m_Bitstream(0), m_Low(0), m_Range(0), m_bufferedByte(0),
      m_numBufferedBytes(0), m_bitsLeft(0) {}

void BinEncoderBase::init(OutputBitstream *bitstream) {
  m_Bitstream = bitstream;
}

void BinEncoderBase::uninit() { m_Bitstream = 0; }

void BinEncoderBase::start() {
  m_Low = 0;
  m_Range = 510;
  m_bufferedByte = 0xff;
  m_numBufferedBytes = 0;
  m_bitsLeft = 23;
  BinCounter::reset();
  m_BinStore.reset();
}

void BinEncoderBase::finish() {
  if (m_Low >> (32 - m_bitsLeft)) {
    m_Bitstream->write(m_bufferedByte + 1, 8);
    while (m_numBufferedBytes > 1) {
      m_Bitstream->write(0x00, 8);
      m_numBufferedBytes--;
    }
    m_Low -= 1 << (32 - m_bitsLeft);
  } else {
    if (m_numBufferedBytes > 0) {
      m_Bitstream->write(m_bufferedByte, 8);
    }
    while (m_numBufferedBytes > 1) {
      m_Bitstream->write(0xff, 8);
      m_numBufferedBytes--;
    }
  }
  m_Bitstream->write(m_Low >> 8, 24 - m_bitsLeft);
}

void BinEncoderBase::restart() {
  m_Low = 0;
  m_Range = 510;
  m_bufferedByte = 0xff;
  m_numBufferedBytes = 0;
  m_bitsLeft = 23;
}

void BinEncoderBase::reset(int qp, int initId) {
  Ctx::init(qp, initId);
  start();
}

void BinEncoderBase::resetBits() {
  m_Low = 0;
  m_bufferedByte = 0xff;
  m_numBufferedBytes = 0;
  m_bitsLeft = 23;
  BinCounter::reset();
}

uint64_t BinEncoderBase::getEstFracBits() const {
  THROW("not supported");
  return 0;
}

unsigned BinEncoderBase::getNumBins(unsigned ctxId) const {
  return BinCounter::getCtx(ctxId);
}

void BinEncoderBase::encodeBinEP(unsigned bin) {
  BinCounter::addEP();
  m_Low <<= 1;
  if (bin) {
    m_Low += m_Range;
  }
  m_bitsLeft--;
  if (m_bitsLeft < 12) {
    writeOut();
  }
}

void BinEncoderBase::encodeBinsEP(unsigned bins, unsigned numBins) {
  BinCounter::addEP(numBins);
  if (m_Range == 256) {
    encodeAlignedBinsEP(bins, numBins);
    return;
  }
  while (numBins > 8) {
    numBins -= 8;
    unsigned pattern = bins >> numBins;
    m_Low <<= 8;
    m_Low += m_Range * pattern;
    bins -= pattern << numBins;
    m_bitsLeft -= 8;
    if (m_bitsLeft < 12) {
      writeOut();
    }
  }
  m_Low <<= numBins;
  m_Low += m_Range * bins;
  m_bitsLeft -= numBins;
  if (m_bitsLeft < 12) {
    writeOut();
  }
}

void BinEncoderBase::encodeRemAbsEP(unsigned bins, unsigned goRicePar,
                                    unsigned cutoff,
                                    int maxLog2TrDynamicRange) {
  const unsigned threshold = cutoff << goRicePar;
  if (bins < threshold) {
    const unsigned bitMask = (1 << goRicePar) - 1;
    const unsigned length = (bins >> goRicePar) + 1;
    encodeBinsEP((1 << length) - 2, length);
    encodeBinsEP(bins & bitMask, goRicePar);
  } else {
    const unsigned maxPrefixLength = 32 - cutoff - maxLog2TrDynamicRange;
    unsigned prefixLength = 0;
    unsigned codeValue = (bins >> goRicePar) - cutoff;
    unsigned suffixLength;
    if (codeValue >= ((1 << maxPrefixLength) - 1)) {
      prefixLength = maxPrefixLength;
      suffixLength = maxLog2TrDynamicRange;
    } else {
      while (codeValue > ((2u << prefixLength) - 2u)) {
        prefixLength++;
      }
      suffixLength = prefixLength + goRicePar + 1; //+1 for the separator bit
    }
    const unsigned totalPrefixLength = prefixLength + cutoff;
    const unsigned bitMask = (1 << goRicePar) - 1;
    const unsigned prefix = (1 << totalPrefixLength) - 1;
    const unsigned suffix =
        ((codeValue - ((1 << prefixLength) - 1)) << goRicePar) |
        (bins & bitMask);
    encodeBinsEP(prefix, totalPrefixLength); // prefix
    encodeBinsEP(suffix, suffixLength); // separator, suffix, and rParam bits
  }
}

void BinEncoderBase::encodeBinTrm(unsigned bin) {
  BinCounter::addTrm();
  m_Range -= 2;
  if (bin) {
    m_Low += m_Range;
    m_Low <<= 7;
    m_Range = 2 << 7;
    m_bitsLeft -= 7;
  } else if (m_Range >= 256) {
    return;
  } else {
    m_Low <<= 1;
    m_Range <<= 1;
    m_bitsLeft--;
  }
  if (m_bitsLeft < 12) {
    writeOut();
  }
}

void BinEncoderBase::align() { m_Range = 256; }

unsigned BinEncoderBase::getNumWrittenBits() {
  return (m_Bitstream->getNumberOfWrittenBits() + 8 * m_numBufferedBytes + 23 -
          m_bitsLeft);
}

uint32_t BinEncoderBase::getNumBins() { return BinCounter::getAll(); }

bool BinEncoderBase::isEncoding() { return true; }

void BinEncoderBase::encodeAlignedBinsEP(unsigned bins, unsigned numBins) {
  unsigned remBins = numBins;
  while (remBins > 0) {
    // The process of encoding an EP bin is the same as that of coding a normal
    // bin where the symbol ranges for 1 and 0 are both half the range:
    //
    //  low = (low + range/2) << 1       (to encode a 1)
    //  low =  low            << 1       (to encode a 0)
    //
    //  i.e.
    //  low = (low + (bin * range/2)) << 1
    //
    //  which is equivalent to:
    //
    //  low = (low << 1) + (bin * range)
    //
    //  this can be generalised for multiple bins, producing the following
    //  expression:
    //
    unsigned binsToCode =
        std::min<unsigned>(remBins, 8); // code bytes if able to take advantage
                                        // of the system's byte-write function
    unsigned binMask = (1 << binsToCode) - 1;
    unsigned newBins = (bins >> (remBins - binsToCode)) & binMask;
    m_Low = (m_Low << binsToCode) + (newBins << 8); // range is known to be 256
    remBins -= binsToCode;
    m_bitsLeft -= binsToCode;
    if (m_bitsLeft < 12) {
      writeOut();
    }
  }
}

void BinEncoderBase::writeOut() {
  unsigned leadByte = m_Low >> (24 - m_bitsLeft);
  m_bitsLeft += 8;
  m_Low &= 0xffffffffu >> m_bitsLeft;
  if (leadByte == 0xff) {
    m_numBufferedBytes++;
  } else {
    if (m_numBufferedBytes > 0) {
      unsigned carry = leadByte >> 8;
      unsigned byte = m_bufferedByte + carry;
      m_bufferedByte = leadByte & 0xff;
      m_Bitstream->write(byte, 8);
      byte = (0xff + carry) & 0xff;
      while (m_numBufferedBytes > 1) {
        m_Bitstream->write(byte, 8);
        m_numBufferedBytes--;
      }
    } else {
      m_numBufferedBytes = 1;
      m_bufferedByte = leadByte;
    }
  }
}

template <class BinProbModel>
TBinEncoder<BinProbModel>::TBinEncoder()
    : BinEncoderBase(static_cast<const BinProbModel *>(nullptr)),
      m_Ctx(static_cast<CtxStore<BinProbModel> &>(*this)) {}

template <class BinProbModel>
void TBinEncoder<BinProbModel>::encodeBin(unsigned bin, unsigned ctxId) {
  BinCounter::addCtx(ctxId);
  BinProbModel &rcProbModel = m_Ctx[ctxId];
  uint32_t LPS = rcProbModel.getLPS(m_Range);

  m_Range -= LPS;
  if (bin != rcProbModel.mps()) {
    int numBits = rcProbModel.getRenormBitsLPS(LPS);
    m_bitsLeft -= numBits;
    m_Low += m_Range;
    m_Low = m_Low << numBits;
    m_Range = LPS << numBits;
    if (m_bitsLeft < 12) {
      writeOut();
    }
  } else {
    if (m_Range < 256) {
      int numBits = rcProbModel.getRenormBitsRange(m_Range);
      m_bitsLeft -= numBits;
      m_Low <<= numBits;
      m_Range <<= numBits;
      if (m_bitsLeft < 12) {
        writeOut();
      }
    }
  }
  rcProbModel.update(bin);
  BinEncoderBase::m_BinStore.addBin(bin, ctxId);
}

template <class BinProbModel>
void TBinEncoder<BinProbModel>::setBinStorage(bool b) {
  m_BinStore.setUse(b);
}

template <class BinProbModel>
const BinStore *TBinEncoder<BinProbModel>::getBinStore() const {
  return &m_BinStore;
}

template <class BinProbModel>
BinEncIf *TBinEncoder<BinProbModel>::getTestBinEncoder() const {
  BinEncIf *testBinEncoder = 0;
  if (m_BinStore.inUse()) {
    testBinEncoder = new TBinEncoder<BinProbModel>();
  }
  return testBinEncoder;
}

template <class BinProbModel>
BitEstimatorBase::BitEstimatorBase(const BinProbModel *dummy)
    : BinEncIf(dummy) {
  m_EstFracBits = 0;
}

void BitEstimatorBase::init(OutputBitstream *bitstream) {}

void BitEstimatorBase::uninit() {}

void BitEstimatorBase::start() { m_EstFracBits = 0; }

void BitEstimatorBase::finish() {}

void BitEstimatorBase::restart() {
  m_EstFracBits = (m_EstFracBits >> SCALE_BITS) << SCALE_BITS;
}

void BitEstimatorBase::reset(int qp, int initId) {
  Ctx::init(qp, initId);
  m_EstFracBits = 0;
}

void BitEstimatorBase::resetBits() { m_EstFracBits = 0; }

uint64_t BitEstimatorBase::getEstFracBits() const { return m_EstFracBits; }

unsigned BitEstimatorBase::getNumBins(unsigned) const {
  THROW("not supported for BitEstimator");
  return 0;
}

void BitEstimatorBase::encodeBinEP(unsigned) {
  m_EstFracBits += BinProbModelBase::estFracBitsEP();
}

void BitEstimatorBase::encodeBinsEP(unsigned, unsigned numBins) {
  m_EstFracBits += BinProbModelBase::estFracBitsEP(numBins);
}

uint32_t BitEstimatorBase::getNumBins() {
  THROW("Not supported");
  return 0;
}

unsigned BitEstimatorBase::getNumWrittenBits() { return (uint32_t)0; }

bool BitEstimatorBase::isEncoding() { return false; }

void BitEstimatorBase::encodeRemAbsEP(unsigned bins, unsigned goRicePar,
                                      unsigned cutoff,
                                      int maxLog2TrDynamicRange) {
  const unsigned threshold = cutoff << goRicePar;
  if (bins < threshold) {
    m_EstFracBits +=
        BinProbModelBase::estFracBitsEP((bins >> goRicePar) + 1 + goRicePar);
  } else {
    const unsigned maxPrefixLength = 32 - cutoff - maxLog2TrDynamicRange;
    unsigned prefixLength = 0;
    unsigned codeValue = (bins >> goRicePar) - cutoff;
    unsigned suffixLength;
    if (codeValue >= ((1 << maxPrefixLength) - 1)) {
      prefixLength = maxPrefixLength;
      suffixLength = maxLog2TrDynamicRange;
    } else {
      while (codeValue > ((2u << prefixLength) - 2u)) {
        prefixLength++;
      }
      suffixLength = prefixLength + goRicePar + 1; //+1 for the separator bit
    }
    m_EstFracBits +=
        BinProbModelBase::estFracBitsEP(cutoff + prefixLength + suffixLength);
  }
}

void BitEstimatorBase::align() {
  static const uint64_t add = BinProbModelBase::estFracBitsEP() - 1;
  static const uint64_t mask = ~add;
  m_EstFracBits += add;
  m_EstFracBits &= mask;
}

template <class BinProbModel>
TBitEstimator<BinProbModel>::TBitEstimator()
    : BitEstimatorBase(static_cast<const BinProbModel *>(nullptr)),
      m_Ctx(static_cast<CtxStore<BinProbModel> &>(*this)) {}

template <class BinProbModel>
void TBitEstimator<BinProbModel>::encodeBin(unsigned bin, unsigned ctxId) {
  m_Ctx[ctxId].estFracBitsUpdate(bin, m_EstFracBits);
}

template <class BinProbModel>
void TBitEstimator<BinProbModel>::encodeBinTrm(unsigned bin) {
  m_EstFracBits += BinProbModel::estFracBitsTrm(bin);
}

template <class BinProbModel>
void TBitEstimator<BinProbModel>::setBinStorage(bool b) {}

template <class BinProbModel>
const BinStore *TBitEstimator<BinProbModel>::getBinStore() const {
  return 0;
}

template <class BinProbModel>
BinEncIf *TBitEstimator<BinProbModel>::getTestBinEncoder() const {
  return 0;
}

template class TBinEncoder<BinProbModel_Std>;

template class TBitEstimator<BinProbModel_Std>;
