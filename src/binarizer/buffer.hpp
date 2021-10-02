#include "common_def.hpp"
#include "motion_info.hpp"

template <typename T> struct AreaBuf : public Size {
  T *buf;
  int stride;
  // the proper type causes awful lot of errors
  // ptrdiff_t stride;

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

  operator AreaBuf<const T>() const {
    return AreaBuf<const T>(buf, stride, width, height);
  }

  void fill(const T &val);
  void memset(const int val);

  void copyFrom(const AreaBuf<const T> &other);
  void roundToOutputBitdepth(const AreaBuf<const T> &src, const ClpRng &clpRng);

  void reconstruct(const AreaBuf<const T> &pred, const AreaBuf<const T> &resi,
                   const ClpRng &clpRng);
  void copyClip(const AreaBuf<const T> &src, const ClpRng &clpRng);

  void subtract(const AreaBuf<const T> &other);
  void extendSingleBorderPel();
  void extendBorderPel(unsigned margin);
  void extendBorderPel(unsigned marginX, unsigned marginY);
  void padBorderPel(unsigned marginX, unsigned marginY, int dir);
  void addWeightedAvg(const AreaBuf<const T> &other1,
                      const AreaBuf<const T> &other2, const ClpRng &clpRng,
                      const int8_t bcwIdx);
  void removeWeightHighFreq(const AreaBuf<T> &other, const bool bClip,
                            const ClpRng &clpRng, const int8_t iBcwWeight);
  void addAvg(const AreaBuf<const T> &other1, const AreaBuf<const T> &other2,
              const ClpRng &clpRng);
  void removeHighFreq(const AreaBuf<T> &other, const bool bClip,
                      const ClpRng &clpRng);
  void updateHistogram(std::vector<int32_t> &hist) const;

  T meanDiff(const AreaBuf<const T> &other) const;
  void subtract(const T val);

  void linearTransform(const int scale, const int shift, const int offset,
                       bool bClip, const ClpRng &clpRng);

  void transposedFrom(const AreaBuf<const T> &other);

  void toLast(const ClpRng &clpRng);

  void rspSignal(std::vector<Pel> &pLUT);
  void scaleSignal(const int scale, const bool dir, const ClpRng &clpRng);
  void applyLumaCTI(std::vector<Pel> &pLUTY);
  void applyChromaCTI(Pel *bufY, int strideY, std::vector<Pel> &pLUTUV,
                      int bitDepth, ChromaFormat chrFormat, bool fwdMap);
  T computeAvg() const;

  T &at(const int &x, const int &y) { return buf[y * stride + x]; }
  const T &at(const int &x, const int &y) const { return buf[y * stride + x]; }

  T &at(const Position &pos) { return buf[pos.y * stride + pos.x]; }
  const T &at(const Position &pos) const { return buf[pos.y * stride + pos.x]; }

  T *bufAt(const int &x, const int &y) { return &at(x, y); }
  const T *bufAt(const int &x, const int &y) const { return &at(x, y); }

  T *bufAt(const Position &pos) { return &at(pos); }
  const T *bufAt(const Position &pos) const { return &at(pos); }

  AreaBuf<T> subBuf(const Position &pos, const Size &size) {
    return AreaBuf<T>(bufAt(pos), stride, size);
  }
  AreaBuf<const T> subBuf(const Position &pos, const Size &size) const {
    return AreaBuf<const T>(bufAt(pos), stride, size);
  }
  AreaBuf<T> subBuf(const int &x, const int &y, const unsigned &_w,
                    const unsigned &_h) {
    return AreaBuf<T>(bufAt(x, y), stride, _w, _h);
  }
  AreaBuf<const T> subBuf(const int &x, const int &y, const unsigned &_w,
                          const unsigned &_h) const {
    return AreaBuf<const T>(bufAt(x, y), stride, _w, _h);
  }
};

typedef AreaBuf<Pel> PelBuf;
typedef AreaBuf<const Pel> CPelBuf;

typedef AreaBuf<TCoeff> CoeffBuf;
typedef AreaBuf<const TCoeff> CCoeffBuf;

typedef AreaBuf<MotionInfo> MotionBuf;
typedef AreaBuf<const MotionInfo> CMotionBuf;

typedef AreaBuf<TCoeff> PLTescapeBuf;
typedef AreaBuf<const TCoeff> CPLTescapeBuf;

typedef AreaBuf<bool> PLTtypeBuf;
typedef AreaBuf<const bool> CPLTtypeBuf;

struct UnitArea;

template <typename T> struct UnitBuf {
  typedef static_vector<AreaBuf<T>, MAX_NUM_COMPONENT> UnitBufBuffers;
  typedef static_vector<AreaBuf<const T>, MAX_NUM_COMPONENT>
      ConstUnitBufBuffers;

  ChromaFormat chromaFormat;
  UnitBufBuffers bufs;

  UnitBuf() : chromaFormat(NUM_CHROMA_FORMAT) {}
  UnitBuf(const ChromaFormat &_chromaFormat, const UnitBufBuffers &_bufs)
      : chromaFormat(_chromaFormat), bufs(_bufs) {}
  UnitBuf(const ChromaFormat &_chromaFormat, UnitBufBuffers &&_bufs)
      : chromaFormat(_chromaFormat), bufs(std::forward<UnitBufBuffers>(_bufs)) {
  }
  UnitBuf(const ChromaFormat &_chromaFormat, const AreaBuf<T> &blkY)
      : chromaFormat(_chromaFormat), bufs{blkY} {}
  UnitBuf(const ChromaFormat &_chromaFormat, AreaBuf<T> &&blkY)
      : chromaFormat(_chromaFormat), bufs{std::forward<AreaBuf<T>>(blkY)} {}
  UnitBuf(const ChromaFormat &_chromaFormat, const AreaBuf<T> &blkY,
          const AreaBuf<T> &blkCb, const AreaBuf<T> &blkCr)
      : chromaFormat(_chromaFormat), bufs{blkY, blkCb, blkCr} {}
  UnitBuf(const ChromaFormat &_chromaFormat, AreaBuf<T> &&blkY,
          AreaBuf<T> &&blkCb, AreaBuf<T> &&blkCr)
      : chromaFormat(_chromaFormat), bufs{std::forward<AreaBuf<T>>(blkY),
                                          std::forward<AreaBuf<T>>(blkCb),
                                          std::forward<AreaBuf<T>>(blkCr)} {}

  operator UnitBuf<const T>() const {
    return UnitBuf<const T>(chromaFormat,
                            ConstUnitBufBuffers(bufs.begin(), bufs.end()));
  }

  AreaBuf<T> &get(const ComponentID comp) { return bufs[comp]; }
  const AreaBuf<T> &get(const ComponentID comp) const { return bufs[comp]; }

  AreaBuf<T> &Y() { return bufs[0]; }
  const AreaBuf<T> &Y() const { return bufs[0]; }
  AreaBuf<T> &Cb() { return bufs[1]; }
  const AreaBuf<T> &Cb() const { return bufs[1]; }
  AreaBuf<T> &Cr() { return bufs[2]; }
  const AreaBuf<T> &Cr() const { return bufs[2]; }

  void fill(const T &val);
  void copyFrom(const UnitBuf<const T> &other, const bool lumaOnly = false,
                const bool chromaOnly = false);
  void roundToOutputBitdepth(const UnitBuf<const T> &src,
                             const ClpRngs &clpRngs);
  void reconstruct(const UnitBuf<const T> &pred, const UnitBuf<const T> &resi,
                   const ClpRngs &clpRngs);
  void copyClip(const UnitBuf<const T> &src, const ClpRngs &clpRngs,
                const bool lumaOnly = false, const bool chromaOnly = false);
  void subtract(const UnitBuf<const T> &other);
  void addWeightedAvg(const UnitBuf<const T> &other1,
                      const UnitBuf<const T> &other2, const ClpRngs &clpRngs,
                      const uint8_t bcwIdx = BCW_DEFAULT,
                      const bool chromaOnly = false,
                      const bool lumaOnly = false);
  void addAvg(const UnitBuf<const T> &other1, const UnitBuf<const T> &other2,
              const ClpRngs &clpRngs, const bool chromaOnly = false,
              const bool lumaOnly = false);
  void extendSingleBorderPel();
  void extendBorderPel(unsigned marginX, unsigned marginY);
  void padBorderPel(unsigned margin, int dir);
  void extendBorderPel(unsigned margin);
  void removeHighFreq(const UnitBuf<T> &other, const bool bClip,
                      const ClpRngs &clpRngs,
                      const int8_t bcwWeight = g_BcwWeights[BCW_DEFAULT]);

  UnitBuf<T> subBuf(const UnitArea &subArea);
  const UnitBuf<const T> subBuf(const UnitArea &subArea) const;
  void colorSpaceConvert(const UnitBuf<T> &other, const bool forward,
                         const ClpRng &clpRng);
};

typedef UnitBuf<Pel> PelUnitBuf;
typedef UnitBuf<const Pel> CPelUnitBuf;

typedef UnitBuf<TCoeff> CoeffUnitBuf;
typedef UnitBuf<const TCoeff> CCoeffUnitBuf;

struct PelStorage : public PelUnitBuf {
  PelStorage();
  ~PelStorage();

  void swap(PelStorage &other);
  void createFromBuf(PelUnitBuf buf);
  void create(const UnitArea &_unit);
  void create(const ChromaFormat &_chromaFormat, const Area &_area,
              const unsigned _maxCUSize = 0, const unsigned _margin = 0,
              const unsigned _alignment = 0,
              const bool _scaleChromaMargin = true);
  void destroy();

  PelBuf getBuf(const CompArea &blk);
  const CPelBuf getBuf(const CompArea &blk) const;

  PelBuf getBuf(const ComponentID CompID);
  const CPelBuf getBuf(const ComponentID CompID) const;

  PelUnitBuf getBuf(const UnitArea &unit);
  const CPelUnitBuf getBuf(const UnitArea &unit) const;
  Pel *getOrigin(const int id) const { return m_origin[id]; }

private:
  Pel *m_origin[MAX_NUM_COMPONENT];
};