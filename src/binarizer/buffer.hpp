#ifndef __BUFFER__
#define __BUFFER__

#include "common_def.hpp"

namespace EntropyCoding {

template <typename T> struct AreaBuf : public Size {
  T *buf;

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

template <typename T> void AreaBuf<T>::fill(const T &val) {
  if (sizeof(T) == 1) {
    if (width == stride) {
      ::memset(buf, reinterpret_cast<const signed char &>(val),
               width * height * sizeof(T));
    } else {
      T *dest = buf;
      size_t line = width * sizeof(T);

      for (unsigned y = 0; y < height; y++) {
        ::memset(dest, reinterpret_cast<const signed char &>(val), line);

        dest += stride;
      }
    }
  } else if (T(0) == val) {
    if (width == stride) {
      ::memset(buf, 0, width * height * sizeof(T));
    } else {
      T *dest = buf;
      size_t line = width * sizeof(T);

      for (unsigned y = 0; y < height; y++) {
        ::memset(dest, 0, line);

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
} // namespace EntropyCoding

#endif