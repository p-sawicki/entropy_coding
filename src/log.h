#ifndef ENTROPY_CODER_LIB_LOG_H
#define ENTROPY_CODER_LIB_LOG_H

#include <fstream>

enum class SyntaxElement {
  end_of_slice_one_bit = 0x00,
  end_of_tile_one_bit = 0x01,
  end_of_subset_one_bit = 0x02,
};

class Logger {
  std::ofstream fs;

public:
  Logger() : fs("log.txt") {}

#ifdef ENABLE_LOGGING
  void LogElement(const SyntaxElement elem) {
    fs << static_cast<int>(elem) << "\n";
  }

  void LogElementValPair(const SyntaxElement elem, const int val) {
    fs << static_cast<int>(elem) << "\t" << val << "\n";
  }
#else
  void LogElement(const SyntaxElement) {}
  void LogElementValPair(const SyntaxElement, const int) {}
#endif // ENABLE_LOGGING

  Logger(const Logger& rhs) = delete;
  Logger& operator=(const Logger& rhs) = delete;
};

Logger binLogger;

#endif // ENTROPY_CODER_LIB_LOG_H