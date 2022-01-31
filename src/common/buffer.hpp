#ifndef COMMON_BUFFER
#define COMMON_BUFFER

#include <cstring>

#include "common_def.hpp"

namespace Common {

template <typename T> struct AreaBuf : public Size {
  T *buf;
  int stride;

  AreaBuf() : Size(), buf(NULL), stride(0) {}
  AreaBuf(T *_buf, const Size &size)
      : Size(size), buf(_buf), stride(size.width) {}
  AreaBuf(T *_buf, const int &_stride, const Size &size)
      : Size(size), buf(_buf), stride(_stride) {}
  AreaBuf(T *_buf, const SizeType &_width, const SizeType &_height)
      : Size(_width, _height), buf(_buf), stride(_width) {}
  AreaBuf(T *_buf, const int &_stride, const SizeType &_width,
          const SizeType &_height)
      : Size(_width, _height), buf(_buf), stride(_stride) {}

  void fill(const T &val);

  T &at(const int &x, const int &y) { return buf[y * stride + x]; }
  const T &at(const int &x, const int &y) const { return buf[y * stride + x]; }
};

typedef AreaBuf<Pel> PelBuf;
typedef AreaBuf<const Pel> CPelBuf;

typedef AreaBuf<TCoeff> CoeffBuf;
typedef AreaBuf<const TCoeff> CCoeffBuf;

typedef AreaBuf<TCoeff> PLTescapeBuf;
typedef AreaBuf<const TCoeff> CPLTescapeBuf;

typedef AreaBuf<bool> PLTtypeBuf;
typedef AreaBuf<const bool> CPLTtypeBuf;

struct UnitArea;

#define SIZE_AWARE_PER_EL_OP(OP, INC)                                          \
  if ((width & 7) == 0) {                                                      \
    for (int y = 0; y < height; y++) {                                         \
      for (int x = 0; x < width; x += 8) {                                     \
        OP(x + 0);                                                             \
        OP(x + 1);                                                             \
        OP(x + 2);                                                             \
        OP(x + 3);                                                             \
        OP(x + 4);                                                             \
        OP(x + 5);                                                             \
        OP(x + 6);                                                             \
        OP(x + 7);                                                             \
      }                                                                        \
                                                                               \
      INC;                                                                     \
    }                                                                          \
  } else if ((width & 3) == 0) {                                               \
    for (int y = 0; y < height; y++) {                                         \
      for (int x = 0; x < width; x += 4) {                                     \
        OP(x + 0);                                                             \
        OP(x + 1);                                                             \
        OP(x + 2);                                                             \
        OP(x + 3);                                                             \
      }                                                                        \
                                                                               \
      INC;                                                                     \
    }                                                                          \
  } else if ((width & 1) == 0) {                                               \
    for (int y = 0; y < height; y++) {                                         \
      for (int x = 0; x < width; x += 2) {                                     \
        OP(x + 0);                                                             \
        OP(x + 1);                                                             \
      }                                                                        \
                                                                               \
      INC;                                                                     \
    }                                                                          \
  } else {                                                                     \
    for (int y = 0; y < height; y++) {                                         \
      for (int x = 0; x < width; x++) {                                        \
        OP(x);                                                                 \
      }                                                                        \
                                                                               \
      INC;                                                                     \
    }                                                                          \
  }

template <typename T> void AreaBuf<T>::fill(const T &val) {
  if (sizeof(T) == 1) {
    if (width == stride) {
      memset(buf, reinterpret_cast<const signed char &>(val),
               width * height * sizeof(T));
    } else {
      T *dest = buf;
      size_t line = width * sizeof(T);

      for (unsigned y = 0; y < height; y++) {
        memset(dest, reinterpret_cast<const signed char &>(val), line);

        dest += stride;
      }
    }
  } else if (T(0) == val) {
    if (width == stride) {
      memset(buf, 0, width * height * sizeof(T));
    } else {
      T *dest = buf;
      size_t line = width * sizeof(T);

      for (unsigned y = 0; y < height; y++) {
        memset(dest, 0, line);

        dest += stride;
      }
    }
  } else {
    T *dest = buf;

#define FILL_INC dest += stride
#define FILL_OP(ADDR) dest[ADDR] = val

    SIZE_AWARE_PER_EL_OP(FILL_OP, FILL_INC);

#undef FILL_INC
#undef FILL_OP
  }
}
} // namespace Common

#endif