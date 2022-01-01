#ifndef ENTROPY_CODEC_BIT_STREAM
#define ENTROPY_CODEC_BIT_STREAM

#include <stdint.h>
#include <stdio.h>
#include <vector>

#include "common_def.hpp"

namespace EntropyCoding {

/**
 * Model of a writable bitstream that accumulates bits to produce a
 * bytestream.
 */
class OutputBitstream {
  /**
   * FIFO for storage of bytes.  Use:
   *  - fifo.push_back(x) to append words
   *  - fifo.clear() to empty the FIFO
   *  - &fifo.front() to get a pointer to the data array.
   *    NB, this pointer is only valid until the next push_back()/clear()
   */
  ::std::vector<uint8_t> &m_fifo;

  uint32_t m_num_held_bits; /// number of bits not flushed to bytestream.
  uint8_t m_held_bits;      /// the bits held and not flushed to bytestream.
                            /// this value is always msb-aligned, bigendian.
public:
  // create / destroy
  OutputBitstream(::std::vector<uint8_t> &fifo);
  ~OutputBitstream();

  // interface for encoding
  /**
   * append uiNumberOfBits least significant bits of uiBits to
   * the current bitstream
   */
  void write(uint32_t uiBits, uint32_t uiNumberOfBits);

  /** insert one bits until the bitstream is byte-aligned */
  void writeAlignOne();

  /** insert zero bits until the bitstream is byte-aligned */
  void writeAlignZero();

  // utility functions

  /**
   * Return a pointer to the start of the byte-stream buffer.
   * Pointer is valid until the next write/flush/reset call.
   * NB, data is arranged such that subsequent bytes in the
   * bytestream are stored in ascending addresses.
   */
  uint8_t *getByteStream() const;

  /**
   * Return the number of valid bytes available from  getByteStream()
   */
  uint32_t getByteStreamLength();

  /**
   * Reset all internal state.
   */
  void clear();

  /**
   * returns the number of bits that need to be written to
   * achieve byte alignment.
   */
  int getNumBitsUntilByteAligned() const;

  /**
   * Return the number of bits that have been written since the last clear()
   */
  uint32_t getNumberOfWrittenBits() const;

  void insertAt(const OutputBitstream &src, uint32_t pos);

  /**
   * Return a reference to the internal fifo
   */
  ::std::vector<uint8_t> &getFIFO();

  uint8_t getHeldBits();

  // OutputBitstream& operator= (const OutputBitstream& src);
  /** Return a reference to the internal fifo */
  const ::std::vector<uint8_t> &getFIFO() const;

  void addSubstream(OutputBitstream *pcSubstream);
  void writeByteAlignment();

  //! returns the number of start code emulations contained in the current
  //! buffer
  int countStartCodeEmulations();
};

/**
 * Model of an input bitstream that extracts bits from a predefined
 * bytestream.
 */
class InputBitstream {
public:
  ::std::vector<uint8_t> m_fifo; /// FIFO for storage of complete bytes
  ::std::vector<uint32_t> m_emulationPreventionByteLocation;

  uint32_t m_fifo_idx; /// Read index into m_fifo

  uint32_t m_num_held_bits;
  uint8_t m_held_bits;
  uint32_t m_numBitsRead;

public:
  /**
   * Create a new bitstream reader object that reads from buf.
   */
  InputBitstream();
  virtual ~InputBitstream() = default;
  InputBitstream(const InputBitstream &src);
  InputBitstream(const ::std::vector<uint8_t> &fifo,
                 const ::std::vector<uint32_t> &emulationPreventionByteLocation,
                 const uint32_t fifo_idx, const uint32_t num_held_bits,
                 const uint8_t held_bits, const uint32_t numBitsRead)
      : m_fifo(fifo),
        m_emulationPreventionByteLocation(emulationPreventionByteLocation),
        m_fifo_idx(fifo_idx), m_num_held_bits(num_held_bits),
        m_held_bits(held_bits), m_numBitsRead(numBitsRead) {}

  void resetToStart();

  // interface for decoding
  void pseudoRead(uint32_t uiNumberOfBits, uint32_t &ruiBits);
  void read(uint32_t uiNumberOfBits, uint32_t &ruiBits);
  void readByte(uint32_t &ruiBits);

  void peekPreviousByte(uint32_t &byte);

  uint32_t readOutTrailingBits();
  uint8_t getHeldBits();
  OutputBitstream &operator=(const OutputBitstream &src);
  uint32_t getByteLocation();

  // Peek at bits in word-storage. Used in determining if we have completed
  // reading of current bitstream and therefore slice in LCEC.
  uint32_t peekBits(uint32_t uiBits);

  // utility functions
  uint32_t read(uint32_t numberOfBits);
  uint32_t readByte();
  uint32_t getNumBitsUntilByteAligned();
  uint32_t getNumBitsLeft();
  InputBitstream *
  extractSubstream(uint32_t uiNumBits); // Read the nominated number of bits,
                                        // and return as a bitstream.
  uint32_t getNumBitsRead();
  uint32_t readByteAlignment();

  void pushEmulationPreventionByteLocation(uint32_t pos);
  uint32_t numEmulationPreventionBytesRead();
  const ::std::vector<uint32_t> &getEmulationPreventionByteLocation() const;
  uint32_t getEmulationPreventionByteLocation(uint32_t idx);
  void clearEmulationPreventionByteLocation();
  void setEmulationPreventionByteLocation(const ::std::vector<uint32_t> &vec);

  const ::std::vector<uint8_t> &getFifo() const;
  ::std::vector<uint8_t> &getFifo();
};
} // namespace EntropyCoding

#endif // ENTROPY_CODEC_BIT_STREAM