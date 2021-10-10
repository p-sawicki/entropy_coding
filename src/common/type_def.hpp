#ifndef ENTROPY_CODEC_TYPE_DEF
#define ENTROPY_CODEC_TYPE_DEF

#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>

namespace EntropyCoding {

enum ComponentID {
  COMPONENT_Y = 0,
  COMPONENT_Cb = 1,
  COMPONENT_Cr = 2,
  MAX_NUM_COMPONENT = 3,
  JOINT_CbCr = MAX_NUM_COMPONENT,
  MAX_NUM_TBLOCKS = MAX_NUM_COMPONENT
};

/// supported slice type
enum SliceType {
  B_SLICE = 0,
  P_SLICE = 1,
  I_SLICE = 2,
  NUMBER_OF_SLICE_TYPES = 3
};

enum ChannelType {
  CHANNEL_TYPE_LUMA = 0,
  CHANNEL_TYPE_CHROMA = 1,
  MAX_NUM_CHANNEL_TYPE = 2
};

enum TreeType {
  TREE_D = 0, // default tree status (for single-tree slice, TREE_D means joint
              // tree; for dual-tree I slice, TREE_D means TREE_L for luma and
              // TREE_C for chroma)
  TREE_L = 1, // separate tree only contains luma (may split)
  TREE_C = 2, // separate tree only contains chroma (not split), to avoid small
              // chroma block
};

enum ModeType {
  MODE_TYPE_ALL = 0,   // all modes can try
  MODE_TYPE_INTER = 1, // can try inter
  MODE_TYPE_INTRA = 2, // can try intra, ibc, palette
};

enum ChromaFormat {
  CHROMA_400 = 0,
  CHROMA_420 = 1,
  CHROMA_422 = 2,
  CHROMA_444 = 3,
  NUM_CHROMA_FORMAT = 4
};

/// supported prediction type
enum PredMode {
  MODE_INTER = 0, ///< inter-prediction mode
  MODE_INTRA = 1, ///< intra-prediction mode
  MODE_IBC = 2,   ///< ibc-prediction mode
  MODE_PLT = 3,   ///< plt-prediction mode
  NUMBER_OF_PREDICTION_MODES = 4,
};

enum SbtIdx {
  SBT_OFF_DCT = 0,
  SBT_VER_HALF = 1,
  SBT_HOR_HALF = 2,
  SBT_VER_QUAD = 3,
  SBT_HOR_QUAD = 4,
  NUMBER_SBT_IDX,
  SBT_OFF_MTS, // note: must be after all SBT modes, only used in fast algorithm
               // to discern the best mode is inter EMT
};

/// reference list index
enum RefPicList {
  REF_PIC_LIST_0 = 0, ///< reference list 0
  REF_PIC_LIST_1 = 1, ///< reference list 1
  NUM_REF_PIC_LIST_01 = 2,
  REF_PIC_LIST_X = 100 ///< special mark
};

enum MergeType {
  MRG_TYPE_DEFAULT_N = 0, // 0
  MRG_TYPE_SUBPU_ATMVP,
  MRG_TYPE_IBC,
  NUM_MRG_TYPE // 5
};

enum ApsType {
  ALF_APS = 0,
  LMCS_APS = 1,
  SCALING_LIST_APS = 2,
};

// TODO: Existing names used for the different NAL unit types can be altered to
// better reflect the names in the spec.
//       However, the names in the spec are not yet stable at this point. Once
//       the names are stable, a cleanup effort can be done without use of
//       macros to alter the names used to indicate the different NAL unit
//       types.
enum NalUnitType {
  NAL_UNIT_CODED_SLICE_TRAIL = 0, // 0
  NAL_UNIT_CODED_SLICE_STSA,      // 1
  NAL_UNIT_CODED_SLICE_RADL,      // 2
  NAL_UNIT_CODED_SLICE_RASL,      // 3

  NAL_UNIT_RESERVED_VCL_4,
  NAL_UNIT_RESERVED_VCL_5,
  NAL_UNIT_RESERVED_VCL_6,

  NAL_UNIT_CODED_SLICE_IDR_W_RADL, // 7
  NAL_UNIT_CODED_SLICE_IDR_N_LP,   // 8
  NAL_UNIT_CODED_SLICE_CRA,        // 9
  NAL_UNIT_CODED_SLICE_GDR,        // 10

  NAL_UNIT_RESERVED_IRAP_VCL_11,
  NAL_UNIT_OPI,                   // 12
  NAL_UNIT_DCI,                   // 13
  NAL_UNIT_VPS,                   // 14
  NAL_UNIT_SPS,                   // 15
  NAL_UNIT_PPS,                   // 16
  NAL_UNIT_PREFIX_APS,            // 17
  NAL_UNIT_SUFFIX_APS,            // 18
  NAL_UNIT_PH,                    // 19
  NAL_UNIT_ACCESS_UNIT_DELIMITER, // 20
  NAL_UNIT_EOS,                   // 21
  NAL_UNIT_EOB,                   // 22
  NAL_UNIT_PREFIX_SEI,            // 23
  NAL_UNIT_SUFFIX_SEI,            // 24
  NAL_UNIT_FD,                    // 25

  NAL_UNIT_RESERVED_NVCL_26,
  NAL_UNIT_RESERVED_NVCL_27,

  NAL_UNIT_UNSPECIFIED_28,
  NAL_UNIT_UNSPECIFIED_29,
  NAL_UNIT_UNSPECIFIED_30,
  NAL_UNIT_UNSPECIFIED_31,
  NAL_UNIT_INVALID
};

/// coefficient scanning type used in ACS
enum CoeffScanType {
  SCAN_DIAG = 0, ///< up-right diagonal scan
  SCAN_TRAV_HOR = 1,
  SCAN_TRAV_VER = 2,
  SCAN_NUMBER_OF_TYPES
};

//////////////////////////////////////////////////////////////////////////
// Encoder modes to try out
//////////////////////////////////////////////////////////////////////////

enum EncModeFeature {
  ENC_FT_FRAC_BITS = 0,
  ENC_FT_DISTORTION,
  ENC_FT_RD_COST,
  ENC_FT_ENC_MODE_TYPE,
  ENC_FT_ENC_MODE_OPTS,
  ENC_FT_ENC_MODE_PART,
  NUM_ENC_FEATURES
};

enum MsgLevel {
  SILENT = 0,
  ERROR = 1,
  WARNING = 2,
  INFO = 3,
  NOTICE = 4,
  VERBOSE = 5,
  DETAILS = 6
};

enum SAOModeMergeTypes {
  SAO_MERGE_LEFT = 0,
  SAO_MERGE_ABOVE,
  NUM_SAO_MERGE_TYPES
};

enum SAOModeNewTypes {
  SAO_TYPE_START_EO = 0,
  SAO_TYPE_EO_0 = SAO_TYPE_START_EO,
  SAO_TYPE_EO_90,
  SAO_TYPE_EO_135,
  SAO_TYPE_EO_45,

  SAO_TYPE_START_BO,
  SAO_TYPE_BO = SAO_TYPE_START_BO,

  NUM_SAO_NEW_TYPES
};
constexpr int NUM_SAO_EO_TYPES_LOG2 = 2;

enum SAOEOClasses {
  SAO_CLASS_EO_FULL_VALLEY = 0,
  SAO_CLASS_EO_HALF_VALLEY = 1,
  SAO_CLASS_EO_PLAIN = 2,
  SAO_CLASS_EO_HALF_PEAK = 3,
  SAO_CLASS_EO_FULL_PEAK = 4,
  NUM_SAO_EO_CLASSES,
};
constexpr int NUM_SAO_BO_CLASSES_LOG2 = 5;
constexpr int NUM_SAO_BO_CLASSES = (1 << NUM_SAO_BO_CLASSES_LOG2);

enum ISPType {
  NOT_INTRA_SUBPARTITIONS = 0,
  HOR_INTRA_SUBPARTITIONS = 1,
  VER_INTRA_SUBPARTITIONS = 2,
  NUM_INTRA_SUBPARTITIONS_MODES = 3,
  INTRA_SUBPARTITIONS_RESERVED = 4
};

/// motion vector predictor direction used in AMVP
enum MvpDir {
  MD_LEFT = 0,    ///< MVP of left block
  MD_ABOVE,       ///< MVP of above block
  MD_ABOVE_RIGHT, ///< MVP of above right block
  MD_BELOW_LEFT,  ///< MVP of below left block
  MD_ABOVE_LEFT   ///< MVP of above left block
};

enum CoeffScanGroupType {
  SCAN_UNGROUPED = 0,
  SCAN_GROUPED_4x4 = 1,
  SCAN_NUMBER_OF_GROUP_TYPES = 2
};

enum PLTRunMode { PLT_RUN_INDEX = 0, PLT_RUN_COPY = 1, NUM_PLT_RUN = 2 };

enum SbtPos { SBT_POS0 = 0, SBT_POS1 = 1, NUMBER_SBT_POS };

enum ImvMode { IMV_OFF = 0, IMV_FPEL, IMV_4PEL, IMV_HPEL, NUM_IMV_MODES };

enum MTSIdx {
  MTS_DCT2_DCT2 = 0,
  MTS_SKIP = 1,
  MTS_DST7_DST7 = 2,
  MTS_DCT8_DST7 = 3,
  MTS_DST7_DCT8 = 4,
  MTS_DCT8_DCT8 = 5
};

template <typename T> class dynamic_cache {
public:
  std::vector<T *> m_cache;
  ~dynamic_cache() { deleteEntries(); }

  void deleteEntries() {
    for (auto &p : m_cache) {
      delete p;
      p = nullptr;
    }

    m_cache.clear();
  }

  T *get() {
    T *ret;

    if (!m_cache.empty()) {
      ret = m_cache.back();
      m_cache.pop_back();
    } else {
      ret = new T;
    }

    return ret;
  }

  void cache(T *el) { m_cache.push_back(el); }

  void cache(std::vector<T *> &vel) {
    m_cache.insert(m_cache.end(), vel.begin(), vel.end());
    vel.clear();
  }
};

typedef dynamic_cache<struct CodingUnit> CUCache;
typedef dynamic_cache<struct PredictionUnit> PUCache;
typedef dynamic_cache<struct TransformUnit> TUCache;

struct XUCache {
  CUCache cuCache;
  PUCache puCache;
  TUCache tuCache;
};

#define CH_L CHANNEL_TYPE_LUMA
#define CH_C CHANNEL_TYPE_CHROMA

class Exception : public ::std::exception {
public:
  Exception(const ::std::string &_s) : m_str(_s) {}
  Exception(const Exception &_e) : ::std::exception(_e), m_str(_e.m_str) {}
  virtual ~Exception() noexcept {};
  virtual const char *what() const noexcept { return m_str.c_str(); }
  Exception &operator=(const Exception &_e) {
    ::std::exception::operator=(_e);
    m_str = _e.m_str;
    return *this;
  }
  template <typename T> Exception &operator<<(T t) {
    std::ostringstream oss;
    oss << t;
    m_str += oss.str();
    return *this;
  }

private:
  ::std::string m_str;
};

// if a check fails with THROW or CHECK, please check if ported correctly from
// assert in revision 1196)
#define THROW(x)                                                               \
  throw(Exception("\nERROR: In function \"")                                   \
        << __FUNCTION__ << "\" in " << __FILE__ << ":" << __LINE__ << ": "     \
        << x)
#define CHECK(c, x)                                                            \
  if (c) {                                                                     \
    THROW(x);                                                                  \
  }
#define EXIT(x) throw(Exception("\n") << x << "\n")
#define CHECK_NULLPTR(_ptr)                                                    \
  CHECK(!(_ptr), "Accessing an empty pointer pointer!")

#if !NDEBUG // for non MSVC compiler, define _DEBUG if in debug mode to have
            // same behavior between MSVC and others in debug
#ifndef _DEBUG
#define _DEBUG 1
#endif
#endif

#if defined(_DEBUG)
#define CHECKD(c, x)                                                           \
  if (c) {                                                                     \
    THROW(x);                                                                  \
  }
#else
#define CHECKD(c, x)
#endif // _DEBUG

static inline ChannelType toChannelType(const ComponentID id) {
  return (id == COMPONENT_Y) ? CHANNEL_TYPE_LUMA : CHANNEL_TYPE_CHROMA;
}
static inline bool isLuma(const ComponentID id) { return (id == COMPONENT_Y); }
static inline bool isLuma(const ChannelType id) {
  return (id == CHANNEL_TYPE_LUMA);
}
static inline bool isChroma(const ComponentID id) {
  return (id != COMPONENT_Y);
}
static inline bool isChroma(const ChannelType id) {
  return (id != CHANNEL_TYPE_LUMA);
}
static inline uint32_t getChannelTypeScaleX(const ChannelType id,
                                            const ChromaFormat fmt) {
  return (isLuma(id) || (fmt == CHROMA_444)) ? 0 : 1;
}
static inline uint32_t getChannelTypeScaleY(const ChannelType id,
                                            const ChromaFormat fmt) {
  return (isLuma(id) || (fmt != CHROMA_420)) ? 0 : 1;
}
static inline uint32_t getComponentScaleX(const ComponentID id,
                                          const ChromaFormat fmt) {
  return getChannelTypeScaleX(toChannelType(id), fmt);
}
static inline uint32_t getComponentScaleY(const ComponentID id,
                                          const ChromaFormat fmt) {
  return getChannelTypeScaleY(toChannelType(id), fmt);
}
static inline uint32_t getNumberValidComponents(const ChromaFormat fmt) {
  return (fmt == CHROMA_400) ? 1 : MAX_NUM_COMPONENT;
}
static inline uint32_t getNumberValidChannels(const ChromaFormat fmt) {
  return (fmt == CHROMA_400) ? 1 : MAX_NUM_CHANNEL_TYPE;
}
static inline bool isChromaEnabled(const ChromaFormat fmt) {
  return !(fmt == CHROMA_400);
}
static inline ComponentID getFirstComponentOfChannel(const ChannelType id) {
  return (isLuma(id) ? COMPONENT_Y : COMPONENT_Cb);
}

namespace Level {
enum Tier { MAIN = 0, HIGH = 1, NUMBER_OF_TIERS = 2 };

enum Name {
  // code = (major_level * 16 + minor_level * 3)
  NONE = 0,
  LEVEL1 = 16,
  LEVEL2 = 32,
  LEVEL2_1 = 35,
  LEVEL3 = 48,
  LEVEL3_1 = 51,
  LEVEL4 = 64,
  LEVEL4_1 = 67,
  LEVEL5 = 80,
  LEVEL5_1 = 83,
  LEVEL5_2 = 86,
  LEVEL6 = 96,
  LEVEL6_1 = 99,
  LEVEL6_2 = 102,
  LEVEL6_3 = 105,
  LEVEL15_5 = 255,
};
} // namespace Level

namespace Profile {
enum Name {
  NONE = 0,
#if JVET_W2005_RANGE_EXTENSION_PROFILES
  INTRA = 8,
#endif
  STILL_PICTURE = 64,
  MAIN_10 = 1,
  MAIN_10_STILL_PICTURE = MAIN_10 | STILL_PICTURE,
  MULTILAYER_MAIN_10 = 17,
  MULTILAYER_MAIN_10_STILL_PICTURE = MULTILAYER_MAIN_10 | STILL_PICTURE,
  MAIN_10_444 = 33,
  MAIN_10_444_STILL_PICTURE = MAIN_10_444 | STILL_PICTURE,
  MULTILAYER_MAIN_10_444 = 49,
  MULTILAYER_MAIN_10_444_STILL_PICTURE = MULTILAYER_MAIN_10_444 | STILL_PICTURE,
#if JVET_W2005_RANGE_EXTENSION_PROFILES
  MAIN_12 = 2,
  MAIN_12_444 = 34,
  MAIN_16_444 = 36,
  MAIN_12_INTRA = MAIN_12 | INTRA,
  MAIN_12_444_INTRA = MAIN_12_444 | INTRA,
  MAIN_16_444_INTRA = MAIN_16_444 | INTRA,
  MAIN_12_STILL_PICTURE = MAIN_12 | STILL_PICTURE,
  MAIN_12_444_STILL_PICTURE = MAIN_12_444 | STILL_PICTURE,
  MAIN_16_444_STILL_PICTURE = MAIN_16_444 | STILL_PICTURE,
#endif
};
} // namespace Profile

struct BitDepths {
  ::std::array<int, MAX_NUM_CHANNEL_TYPE> recon; ///< the bit depth as indicated in the SPS
};

#if RExt__HIGH_BIT_DEPTH_SUPPORT
typedef int Pel;              ///< pixel type
typedef int64_t TCoeff;       ///< transform coefficient
typedef int TMatrixCoeff;     ///< transform matrix coefficient
typedef int16_t TFilterCoeff; ///< filter coefficient
typedef int64_t
    Intermediate_Int; ///< used as intermediate value in calculations
typedef uint64_t
    Intermediate_UInt; ///< used as intermediate value in calculations
#else
typedef int16_t Pel;          ///< pixel type
typedef int TCoeff;           ///< transform coefficient
typedef int16_t TMatrixCoeff; ///< transform matrix coefficient
typedef int16_t TFilterCoeff; ///< filter coefficient
typedef int Intermediate_Int; ///< used as intermediate value in calculations
typedef uint32_t
    Intermediate_UInt; ///< used as intermediate value in calculations
#endif

typedef uint64_t SplitSeries; ///< used to encoded the splits that caused a
                              ///< particular CU size
typedef uint64_t
    ModeTypeSeries; ///< used to encoded the ModeType at different split depth

typedef uint64_t Distortion; ///< distortion measurement

template <typename T, size_t N> class static_vector {
  T _arr[N];
  size_t _size;

public:
  typedef T value_type;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef T &reference;
  typedef T const &const_reference;
  typedef T *pointer;
  typedef T const *const_pointer;
  typedef T *iterator;
  typedef T const *const_iterator;

  static const size_type max_num_elements = N;

  static_vector() : _size(0) {}
  static_vector(size_t N_) : _size(N_) {}
  static_vector(size_t N_, const T &_val) : _size(0) { resize(N_, _val); }
  template <typename It> static_vector(It _it1, It _it2) : _size(0) {
    while (_it1 < _it2)
      _arr[_size++] = *_it1++;
  }
  static_vector(std::initializer_list<T> _il) : _size(0) {
    typename std::initializer_list<T>::iterator _src1 = _il.begin();
    typename std::initializer_list<T>::iterator _src2 = _il.end();

    while (_src1 < _src2)
      _arr[_size++] = *_src1++;

    CHECKD(_size > N, "capacity exceeded");
  }
  static_vector &operator=(std::initializer_list<T> _il) {
    _size = 0;

    typename std::initializer_list<T>::iterator _src1 = _il.begin();
    typename std::initializer_list<T>::iterator _src2 = _il.end();

    while (_src1 < _src2)
      _arr[_size++] = *_src1++;

    CHECKD(_size > N, "capacity exceeded");
  }

  void resize(size_t N_) {
    CHECKD(N_ > N, "capacity exceeded");
    while (_size < N_)
      _arr[_size++] = T();
    _size = N_;
  }
  void resize(size_t N_, const T &_val) {
    CHECKD(N_ > N, "capacity exceeded");
    while (_size < N_)
      _arr[_size++] = _val;
    _size = N_;
  }
  void reserve(size_t N_) { CHECKD(N_ > N, "capacity exceeded"); }
  void push_back(const T &_val) {
    CHECKD(_size >= N, "capacity exceeded");
    _arr[_size++] = _val;
  }
  void push_back(T &&val) {
    CHECKD(_size >= N, "capacity exceeded");
    _arr[_size++] = std::forward<T>(val);
  }
  void pop_back() {
    CHECKD(_size == 0, "calling pop_back on an empty vector");
    _size--;
  }
  void pop_front() {
    CHECKD(_size == 0, "calling pop_front on an empty vector");
    _size--;
    for (int i = 0; i < _size; i++)
      _arr[i] = _arr[i + 1];
  }
  void clear() { _size = 0; }
  reference at(size_t _i) {
    CHECKD(_i >= _size, "Trying to access an out-of-bound-element");
    return _arr[_i];
  }
  const_reference at(size_t _i) const {
    CHECKD(_i >= _size, "Trying to access an out-of-bound-element");
    return _arr[_i];
  }
  reference operator[](size_t _i) {
    CHECKD(_i >= _size, "Trying to access an out-of-bound-element");
    return _arr[_i];
  }
  const_reference operator[](size_t _i) const {
    CHECKD(_i >= _size, "Trying to access an out-of-bound-element");
    return _arr[_i];
  }
  reference front() {
    CHECKD(_size == 0, "Trying to access the first element of an empty vector");
    return _arr[0];
  }
  const_reference front() const {
    CHECKD(_size == 0, "Trying to access the first element of an empty vector");
    return _arr[0];
  }
  reference back() {
    CHECKD(_size == 0, "Trying to access the last element of an empty vector");
    return _arr[_size - 1];
  }
  const_reference back() const {
    CHECKD(_size == 0, "Trying to access the last element of an empty vector");
    return _arr[_size - 1];
  }
  pointer data() { return _arr; }
  const_pointer data() const { return _arr; }
  iterator begin() { return _arr; }
  const_iterator begin() const { return _arr; }
  const_iterator cbegin() const { return _arr; }
  iterator end() { return _arr + _size; }
  const_iterator end() const { return _arr + _size; };
  const_iterator cend() const { return _arr + _size; };
  size_type size() const { return _size; };
  size_type byte_size() const { return _size * sizeof(T); }
  bool empty() const { return _size == 0; }

  size_type capacity() const { return N; }
  size_type max_size() const { return N; }
  size_type byte_capacity() const { return sizeof(_arr); }

  iterator insert(const_iterator _pos, const T &_val) {
    CHECKD(_size >= N, "capacity exceeded");
    for (difference_type i = _size - 1; i >= _pos - _arr; i--)
      _arr[i + 1] = _arr[i];
    *const_cast<iterator>(_pos) = _val;
    _size++;
    return const_cast<iterator>(_pos);
  }

  iterator insert(const_iterator _pos, T &&_val) {
    CHECKD(_size >= N, "capacity exceeded");
    for (difference_type i = _size - 1; i >= _pos - _arr; i--)
      _arr[i + 1] = _arr[i];
    *const_cast<iterator>(_pos) = std::forward<T>(_val);
    _size++;
    return const_cast<iterator>(_pos);
  }
  template <class InputIt>
  iterator insert(const_iterator _pos, InputIt first, InputIt last) {
    const difference_type numEl = last - first;
    CHECKD(_size + numEl >= N, "capacity exceeded");
    for (difference_type i = _size - 1; i >= _pos - _arr; i--)
      _arr[i + numEl] = _arr[i];
    iterator it = const_cast<iterator>(_pos);
    _size += numEl;
    while (first != last)
      *it++ = *first++;
    return const_cast<iterator>(_pos);
  }

  iterator insert(const_iterator _pos, size_t numEl,
                  const T &val) { // const difference_type numEl = last - first;
    CHECKD(_size + numEl >= N, "capacity exceeded");
    for (difference_type i = _size - 1; i >= _pos - _arr; i--)
      _arr[i + numEl] = _arr[i];
    iterator it = const_cast<iterator>(_pos);
    _size += numEl;
    for (int k = 0; k < numEl; k++)
      *it++ = val;
    return const_cast<iterator>(_pos);
  }

  void erase(const_iterator _pos) {
    iterator it = const_cast<iterator>(_pos) - 1;
    iterator last = end() - 1;
    while (++it != last)
      *it = *(it + 1);
    _size--;
  }
};

enum SAOMode // mode
{ SAO_MODE_OFF = 0,
  SAO_MODE_NEW,
  SAO_MODE_MERGE,
  NUM_SAO_MODES };

constexpr int MAX_NUM_SAO_CLASSES =
    32; //(NUM_SAO_EO_GROUPS >
        // NUM_SAO_BO_GROUPS)?NUM_SAO_EO_GROUPS:NUM_SAO_BO_GROUPS

struct SAOOffset {
  SAOMode modeIdc; // NEW, MERGE, OFF
  int typeIdc; // union of SAOModeMergeTypes and SAOModeNewTypes, depending on
               // modeIdc.
  int typeAuxInfo; // BO: starting band index
  ::std::array<int, MAX_NUM_SAO_CLASSES> offset;
};

struct SAOBlkParam {
  SAOOffset &operator[](int compIdx) { return offsetParam[compIdx]; }
  const SAOOffset &operator[](int compIdx) const {
    return offsetParam[compIdx];
  }

  ::std::array<SAOOffset, MAX_NUM_COMPONENT> offsetParam;
};

class ChromaCbfs {
public:
  ChromaCbfs() : Cb(true), Cr(true) {}
  ChromaCbfs(bool _cbf) : Cb(_cbf), Cr(_cbf) {}

public:
  bool sigChroma(ChromaFormat chromaFormat) const {
    if (chromaFormat == CHROMA_400) {
      return false;
    }
    return (Cb || Cr);
  }
  bool &cbf(ComponentID compID) {
    bool *cbfs[MAX_NUM_TBLOCKS] = {nullptr, &Cb, &Cr};

    return *cbfs[compID];
  }

public:
  bool Cb;
  bool Cr;
};
} // namespace EntropyCoding

#endif // ENTROPY_CODEC_TYPE_DEF