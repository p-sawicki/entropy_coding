#ifndef COMMON_COMMON_DEF
#define COMMON_COMMON_DEF

#include <limits>

#include "type_def.hpp"

namespace Common {

static const int AFFINE_ME_LIST_SIZE = 4;
static const int AFFINE_ME_LIST_SIZE_LD = 3;
static const double AFFINE_ME_LIST_MVP_TH = 1.0;

// ====================================================================================================================
// Common constants
// ====================================================================================================================
static const uint64_t MAX_UINT64 = 0xFFFFFFFFFFFFFFFFU;
static const uint32_t MAX_UINT =
    0xFFFFFFFFU; ///< max. value of unsigned 32-bit integer
static const int MAX_INT = 2147483647; ///< max. value of signed 32-bit integer
static const uint8_t MAX_UCHAR = 255;
static const uint8_t MAX_SCHAR = 127;
static const double MAX_DOUBLE = 1.7e+308; ///< max. value of double-type value

// ====================================================================================================================
// Coding tool configuration
// ====================================================================================================================
// Most of these should not be changed - they resolve the meaning of otherwise
// magic numbers.

static const int MAX_GOP = 64; ///< max. value of hierarchical GOP size
static const int MAX_NUM_REF_PICS =
    29; ///< max. number of pictures used for reference
static const int MAX_NUM_REF =
    16; ///< max. number of entries in picture reference list
static const int MAX_QP = 63;
static const int NOT_VALID = -1;

static const int AMVP_MAX_NUM_CANDS =
    2; ///< AMVP: advanced motion vector prediction - max number of final
       ///< candidates
static const int AMVP_MAX_NUM_CANDS_MEM =
    3; ///< AMVP: advanced motion vector prediction - max number of candidates
static const int AMVP_DECIMATION_FACTOR = 2;
static const int MRG_MAX_NUM_CANDS = 6;        ///< MERGE
static const int AFFINE_MRG_MAX_NUM_CANDS = 5; ///< AFFINE MERGE
static const int IBC_MRG_MAX_NUM_CANDS = 6;    ///< IBC MERGE

static const int MAX_TLAYER =
    7; ///< Explicit temporal layer QP offset - max number of temporal layer

static const int ADAPT_SR_SCALE =
    1; ///< division factor for adaptive search range

static const int MIN_TB_LOG2_SIZEY = 2;
static const int MAX_TB_LOG2_SIZEY = 6;

static const int MIN_TB_SIZEY = 1 << MIN_TB_LOG2_SIZEY;
static const int MAX_TB_SIZEY = 1 << MAX_TB_LOG2_SIZEY;

static const int MAX_NESTING_NUM_LAYER = 64;

static const int MAX_VPS_LAYERS = 64;
static const int MAX_VPS_SUBLAYERS = 7;
static const int MAX_NUM_OLSS = 256;
static const int MAX_VPS_OLS_MODE_IDC = 2;

static const int MIP_MAX_WIDTH = MAX_TB_SIZEY;
static const int MIP_MAX_HEIGHT = MAX_TB_SIZEY;

static const int MAX_NUM_ALF_ALTERNATIVES_CHROMA = 8;
static const int MAX_NUM_ALF_CLASSES = 25;
static const int MAX_NUM_ALF_LUMA_COEFF = 13;
static const int MAX_NUM_ALF_CHROMA_COEFF = 7;
static const int MAX_ALF_FILTER_LENGTH = 7;
static const int MAX_ALF_PADDING_SIZE = 4;
#define MAX_NUM_CC_ALF_FILTERS 4
static constexpr int MAX_NUM_CC_ALF_CHROMA_COEFF = 8;
static constexpr int CCALF_DYNAMIC_RANGE = 6;
static constexpr int CCALF_BITS_PER_COEFF_LEVEL = 3;

static const int ALF_FIXED_FILTER_NUM = 64;
static const int ALF_CTB_MAX_NUM_APS = 8;
static const int NUM_FIXED_FILTER_SETS = 16;

static const int MAX_BDOF_APPLICATION_REGION = 16;

static const int MAX_CPB_CNT = 32; ///< Upper bound of (cpb_cnt_minus1 + 1)
static const int MAX_NUM_LAYER_IDS = 64;
static const int COEF_REMAIN_BIN_REDUCTION =
    5; ///< indicates the level at which the VLC transitions from Golomb-Rice to
       ///< TU+EG(k)
static const int CU_DQP_TU_CMAX = 5; ///< max number bins for truncated unary
static const int CU_DQP_EG_k = 0;    ///< expgolomb order

static const int SBH_THRESHOLD =
    4; ///< value of the fixed SBH controlling threshold

static const int MAX_NUM_VPS = 16;
static const int MAX_NUM_SPS = 16;
static const int MAX_NUM_PPS = 64;
static const int MAX_NUM_APS = 32;     // Currently APS ID has 5 bits
static const int NUM_APS_TYPE_LEN = 3; // Currently APS Type has 3 bits
static const int MAX_NUM_APS_TYPE =
    8; // Currently APS Type has 3 bits so the max type is 8

static constexpr int MAX_TILE_COLS = 30; // Maximum number of tile columns
static constexpr int MAX_TILES = 990;    // Maximum number of tiles
static constexpr int MAX_SLICES = 1000;  // Maximum number of slices per picture

static const int MLS_GRP_NUM =
    1024; ///< Max number of coefficient groups, max(16, 256)

static const int MLS_CG_SIZE = 4; ///< Coefficient group size of 4x4; =
                                  ///< MLS_CG_LOG2_WIDTH + MLS_CG_LOG2_HEIGHT

static const int RVM_VCEGAM10_M = 4;

static const int MAX_REF_LINE_IDX = 3;  // highest refLine offset in the list
static const int MRL_NUM_REF_LINES = 3; // number of candidates in the array
static const int MULTI_REF_LINE_IDX[4] = {0, 1, 2, 0};

static const int PRED_REG_MIN_WIDTH =
    4; // Minimum prediction region width for ISP subblocks

static const int NUM_LUMA_MODE =
    67; ///< Planar + DC + 65 directional mode (4*16 + 1)
static const int NUM_LMC_MODE = 1 + 2; ///< LMC + MDLM_T + MDLM_L
static const int NUM_INTRA_MODE = (NUM_LUMA_MODE + NUM_LMC_MODE);

static const int NUM_EXT_LUMA_MODE = 28;

static const int NUM_DIR = (((NUM_LUMA_MODE - 3) >> 2) + 1);
static const int PLANAR_IDX = 0; ///< index for intra PLANAR mode
static const int DC_IDX = 1;     ///< index for intra DC     mode
static const int HOR_IDX =
    (1 * (NUM_DIR - 1) + 2); ///< index for intra HORIZONTAL mode
static const int DIA_IDX =
    (2 * (NUM_DIR - 1) + 2); ///< index for intra DIAGONAL   mode
static const int VER_IDX =
    (3 * (NUM_DIR - 1) + 2); ///< index for intra VERTICAL   mode
static const int VDIA_IDX =
    (4 * (NUM_DIR - 1) + 2); ///< index for intra VDIAGONAL  mode
static const int BDPCM_IDX =
    (5 * (NUM_DIR - 1) + 2);             ///< index for intra VDIAGONAL  mode
static const int NOMODE_IDX = MAX_UCHAR; ///< indicating uninitialized elements

static const int NUM_CHROMA_MODE =
    (5 + NUM_LMC_MODE); ///< total number of chroma modes
static const int LM_CHROMA_IDX =
    NUM_LUMA_MODE; ///< chroma mode index for derived from LM mode
static const int MDLM_L_IDX = LM_CHROMA_IDX + 1; ///< MDLM_L
static const int MDLM_T_IDX = LM_CHROMA_IDX + 2; ///< MDLM_T
static const int DM_CHROMA_IDX =
    NUM_INTRA_MODE; ///< chroma mode index for derived from luma intra mode

static const uint32_t NUM_TRAFO_MODES_MTS =
    6; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64, 128
static const uint32_t MTS_INTRA_MAX_CU_SIZE =
    32; ///< Max Intra CU size applying EMT, supported values: 8, 16, 32, 64,
        ///< 128
static const uint32_t MTS_INTER_MAX_CU_SIZE =
    32; ///< Max Inter CU size applying EMT, supported values: 8, 16, 32, 64,
        ///< 128
static const int NUM_MOST_PROBABLE_MODES = 6;
static const int LM_SYMBOL_NUM = (1 + NUM_LMC_MODE);

static const int MAX_NUM_MIP_MODE = 32; ///< maximum number of MIP pred. modes
static const int FAST_UDI_MAX_RDMODE_NUM =
    (NUM_LUMA_MODE + MAX_NUM_MIP_MODE); ///< maximum number of RD comparison in
                                        ///< fast-UDI estimation loop

static const int MAX_LFNST_COEF_NUM = 16;

static const int LFNST_LAST_SIG_LUMA = 1;
static const int LFNST_LAST_SIG_CHROMA = 1;

static const int NUM_LFNST_NUM_PER_SET = 3;

static const int CABAC_INIT_PRESENT_FLAG = 1;

static const int MV_FRACTIONAL_BITS_INTERNAL = 4;
static const int MV_FRACTIONAL_BITS_SIGNAL = 2;
static const int MV_FRACTIONAL_BITS_DIFF =
    MV_FRACTIONAL_BITS_INTERNAL - MV_FRACTIONAL_BITS_SIGNAL;
static const int LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS_SIGNAL =
    1 << MV_FRACTIONAL_BITS_SIGNAL;
static const int LUMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS =
    1 << MV_FRACTIONAL_BITS_INTERNAL;
static const int CHROMA_INTERPOLATION_FILTER_SUB_SAMPLE_POSITIONS =
    1 << (MV_FRACTIONAL_BITS_INTERNAL + 1);

static const int MAX_NUM_SUB_PICS = (1 << 16);
static const int MAX_NUM_LONG_TERM_REF_PICS = 33;
static const int NUM_LONG_TERM_REF_PIC_SPS = 0;

static const int MAX_QP_OFFSET_LIST_SIZE =
    6; ///< Maximum size of QP offset list is 6 entries
static const int MAX_NUM_CQP_MAPPING_TABLES =
    3; ///< Maximum number of chroma QP mapping tables (Cb, Cr and joint Cb-Cr)
static const int MIN_QP_VALUE_FOR_16_BIT =
    -48; ////< Minimum value for QP (-6*(bitdepth - 8) ) for bit depth 16 ;
         /// actual minimum QP value is bit depth dependent
static const int MAX_NUM_QP_VALUES =
    MAX_QP + 1 - MIN_QP_VALUE_FOR_16_BIT; ////< Maximum number of QP values
                                          /// possible - bit depth dependent

// Cost mode support
static const int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP =
    0; ///< QP to use for lossless coding.
static const int LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME =
    4; ///< QP' to use for mixed_lossy_lossless coding.
static const int RExt__GOLOMB_RICE_ADAPTATION_STATISTICS_SETS =
    MAX_NUM_COMPONENT;

static const int RExt__PREDICTION_WEIGHTING_ANALYSIS_DC_PRECISION =
    0; ///< Additional fixed bit precision used during encoder-side weighting
       ///< prediction analysis. Currently only used when
       ///< high_precision_prediction_weighting_flag is set, for backwards
       ///< compatibility reasons.

static const int MAX_TIMECODE_SEI_SETS = 3; ///< Maximum number of time sets

static const int MAX_CU_DEPTH = 7; ///< log2(CTUSize)
static const int MAX_CU_SIZE = 1 << MAX_CU_DEPTH;
static const int MIN_CU_LOG2 = 2;
static const int MIN_PU_SIZE = 4;
static const int MAX_NUM_PARTS_IN_CTU =
    ((MAX_CU_SIZE * MAX_CU_SIZE) >> (MIN_CU_LOG2 << 1));
static const int MAX_NUM_TUS =
    16; ///< Maximum number of TUs within one CU. When max TB size is 32x32, up
        ///< to 16 TUs within one CU (128x128) is supported
static const int MAX_LOG2_DIFF_CU_TR_SIZE = 3;
static const int MAX_CU_TILING_PARTITIONS = 1
                                            << (MAX_LOG2_DIFF_CU_TR_SIZE << 1);

static const int JVET_C0024_ZERO_OUT_TH = 32;

static const int MAX_NUM_PART_IDXS_IN_CTU_WIDTH =
    MAX_CU_SIZE / MIN_PU_SIZE; ///< maximum number of partition indices across
                               ///< the width of a CTU (or height of a CTU)
static const int SCALING_LIST_REM_NUM = 6;

static const int QUANT_SHIFT = 14; ///< Q(4) = 2^14
static const int IQUANT_SHIFT = 6;

static constexpr int SCALE_BITS = 15; // Precision for fractional bit estimates
static constexpr double FRAC_BITS_SCALE = 1.0 / (1 << SCALE_BITS);

static constexpr int SCALING_LIST_PRED_MODES = 2;
static const int SCALING_LIST_NUM =
    MAX_NUM_COMPONENT *
    SCALING_LIST_PRED_MODES; ///< list number for quantization matrix

static const int SCALING_LIST_START_VALUE = 8; ///< start value for dpcm mode
static const int MAX_MATRIX_COEF_NUM =
    64; ///< max coefficient number for quantization matrix
static const int MAX_MATRIX_SIZE_NUM =
    8; ///< max size number for quantization matrix
static const int SCALING_LIST_BITS = 8; ///< bit depth of scaling list entries
static const int LOG2_SCALING_LIST_NEUTRAL_VALUE =
    4; ///< log2 of the value that, when used in a scaling list, has no effect
       ///< on quantisation
static const int SCALING_LIST_DC = 16; ///< default DC value

static const int LAST_SIGNIFICANT_GROUPS = 14;

static const int AFFINE_MIN_BLOCK_SIZE = 4; ///< Minimum affine MC block size

static const int MMVD_REFINE_STEP = 8; ///< max number of distance step
static const int MMVD_MAX_REFINE_NUM =
    (MMVD_REFINE_STEP * 4); ///< max number of candidate from a base candidate
static const int MMVD_BASE_MV_NUM = 2; ///< max number of base candidate
static const int MMVD_ADD_NUM =
    (MMVD_MAX_REFINE_NUM *
     MMVD_BASE_MV_NUM); ///< total number of mmvd candidate
static const int MMVD_MRG_MAX_RD_NUM = MRG_MAX_NUM_CANDS;
static const int MMVD_MRG_MAX_RD_BUF_NUM =
    (MMVD_MRG_MAX_RD_NUM + 1); ///< increase buffer size by 1

static const int MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_LUMA = 28;
static const int MAX_TU_LEVEL_CTX_CODED_BIN_CONSTRAINT_CHROMA = 28;

static const int BIO_EXTEND_SIZE = 1;
static const int BIO_TEMP_BUFFER_SIZE =
    (MAX_CU_SIZE + 2 * BIO_EXTEND_SIZE) * (MAX_CU_SIZE + 2 * BIO_EXTEND_SIZE);

static const int PROF_BORDER_EXT_W = 1;
static const int PROF_BORDER_EXT_H = 1;
static const int BCW_NUM = 5; ///< the number of weight options
static const int BCW_DEFAULT = ((uint8_t)(
    BCW_NUM >> 1)); ///< Default weighting index representing for w=0.5
static const int BCW_SIZE_CONSTRAINT =
    256; ///< disabling Bcw if cu size is smaller than 256
static const int MAX_NUM_HMVP_CANDS =
    (MRG_MAX_NUM_CANDS - 1); ///< maximum number of HMVP candidates to be stored
                             ///< and used in merge list
static const int MAX_NUM_HMVP_AVMPCANDS =
    4; ///< maximum number of HMVP candidates to be used in AMVP list

static const int ALF_VB_POS_ABOVE_CTUROW_LUMA = 4;
static const int ALF_VB_POS_ABOVE_CTUROW_CHMA = 2;

#if W0038_DB_OPT
static const int MAX_ENCODER_DEBLOCKING_QUALITY_LAYERS = 8;
#endif

#if SHARP_LUMA_DELTA_QP
static const uint32_t LUMA_LEVEL_TO_DQP_LUT_MAXSIZE =
    1024; ///< max LUT size for QP offset based on luma

#endif
static const int DMVR_SUBCU_WIDTH = 16;
static const int DMVR_SUBCU_HEIGHT = 16;
static const int DMVR_SUBCU_WIDTH_LOG2 = 4;
static const int DMVR_SUBCU_HEIGHT_LOG2 = 4;
static const int MAX_NUM_SUBCU_DMVR =
    ((MAX_CU_SIZE * MAX_CU_SIZE) >>
     (DMVR_SUBCU_WIDTH_LOG2 + DMVR_SUBCU_HEIGHT_LOG2));
static const int DMVR_NUM_ITERATION = 2;

// QTBT high level parameters
// for I slice luma CTB configuration para.
static const int MAX_BT_DEPTH = 4; ///<  <=7
                                   // for P/B slice CTU config. para.
static const int MAX_BT_DEPTH_INTER =
    4; ///< <=7
       // for I slice chroma CTB configuration para. (in luma samples)
static const int MAX_BT_DEPTH_C = 0; ///< <=7
static const int MIN_DUALTREE_CHROMA_WIDTH = 4;
static const int MIN_DUALTREE_CHROMA_SIZE = 16;

static const int SKIP_DEPTH = 3;
static const int PICTURE_DISTANCE_TH = 1;
static const int FAST_SKIP_DEPTH = 2;

static const double PBINTRA_RATIO = 1.1;
static const int NUM_MRG_SATD_CAND = 4;
static const double MRG_FAST_RATIO = 1.25;
static const int NUM_AFF_MRG_SATD_CAND = 2;

static const double AMAXBT_TH32 = 15.0;
static const double AMAXBT_TH64 = 30.0;

// need to know for static memory allocation
static const int MAX_DELTA_QP = 7; ///< maximum supported delta QP value
static const int MAX_TESTED_QPs =
    (1 + 1 + (MAX_DELTA_QP << 1)); ///< dqp=0 +- max_delta_qp + lossless mode

static const int COM16_C806_TRANS_PREC = 0;

static const int NTAPS_LUMA = 8;   ///< Number of taps for luma
static const int NTAPS_CHROMA = 4; ///< Number of taps for chroma
#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
static const int MAX_LADF_INTERVALS =
    5; /// max number of luma adaptive deblocking filter qp offset intervals
#endif

static const int NTAPS_BILINEAR = 2; ///< Number of taps for bilinear filter

static const int ATMVP_SUB_BLOCK_SIZE = 3; ///< sub-block size for ATMVP
static const int GEO_MAX_NUM_UNI_CANDS = 6;
static const int GEO_MAX_NUM_CANDS =
    GEO_MAX_NUM_UNI_CANDS * (GEO_MAX_NUM_UNI_CANDS - 1);
static const int GEO_MIN_CU_LOG2 = 3;
static const int GEO_MAX_CU_LOG2 = 6;
static const int GEO_MIN_CU_SIZE = 1 << GEO_MIN_CU_LOG2;
static const int GEO_MAX_CU_SIZE = 1 << GEO_MAX_CU_LOG2;
static const int GEO_NUM_CU_SIZE = (GEO_MAX_CU_LOG2 - GEO_MIN_CU_LOG2) + 1;
static const int GEO_NUM_PARTITION_MODE = 64;
static const int GEO_NUM_ANGLES = 32;
static const int GEO_NUM_DISTANCES = 4;
static const int GEO_NUM_PRESTORED_MASK = 6;
static const int GEO_WEIGHT_MASK_SIZE =
    3 * (GEO_MAX_CU_SIZE >> 3) * 2 + GEO_MAX_CU_SIZE;
static const int GEO_MV_MASK_SIZE = GEO_WEIGHT_MASK_SIZE >> 2;
static const int GEO_MAX_TRY_WEIGHTED_SAD = 60;
static const int GEO_MAX_TRY_WEIGHTED_SATD = 8;

static const int SBT_MAX_SIZE = 64; ///< maximum CU size for using SBT
static const int SBT_NUM_SL =
    10; ///< maximum number of historical PU decision saved for a CU
static const int SBT_NUM_RDO = 2; ///< maximum number of SBT mode tried for a PU

static const int NUM_INTER_CU_INFO_SAVE =
    8; ///< maximum number of inter cu information saved for fast algorithm
static const int LDT_MODE_TYPE_INHERIT =
    0; ///< No need to signal mode_constraint_flag, and the modeType of the
       ///< region is inherited from its parent node
static const int LDT_MODE_TYPE_INFER =
    1; ///< No need to signal mode_constraint_flag, and the modeType of the
       ///< region is inferred as MODE_TYPE_INTRA
static const int LDT_MODE_TYPE_SIGNAL =
    2; ///< Need to signal mode_constraint_flag, and the modeType of the region
       ///< is determined by the flag

static const int IBC_MAX_CAND_SIZE = 16; // max block size for ibc search
static const int IBC_NUM_CANDIDATES =
    64; ///< Maximum number of candidates to store/test
static const int CHROMA_REFINEMENT_CANDIDATES =
    8; /// 8 candidates BV to choose from
static const int IBC_FAST_METHOD_NOINTRA_IBCCBF0 = 0x01;
static const int IBC_FAST_METHOD_BUFFERBV = 0X02;
static const int IBC_FAST_METHOD_ADAPTIVE_SEARCHRANGE = 0X04;
static constexpr int MV_EXPONENT_BITCOUNT = 4;
static constexpr int MV_MANTISSA_BITCOUNT = 6;
static constexpr int MV_MANTISSA_UPPER_LIMIT =
    ((1 << (MV_MANTISSA_BITCOUNT - 1)) - 1);
static constexpr int MV_MANTISSA_LIMIT = (1 << (MV_MANTISSA_BITCOUNT - 1));
static constexpr int MV_EXPONENT_MASK = ((1 << MV_EXPONENT_BITCOUNT) - 1);

static constexpr int MV_BITS = 18;
static constexpr int MV_MAX = (1 << (MV_BITS - 1)) - 1;
static constexpr int MV_MIN = -(1 << (MV_BITS - 1));

static const int MVD_MAX = (1 << 17) - 1;
static const int MVD_MIN = -(1 << 17);

static const int PIC_ANALYZE_CW_BINS = 32;
static const int PIC_CODE_CW_BINS = 16;
static const int LMCS_SEG_NUM = 32;
static const int FP_PREC = 11;
static const int CSCALE_FP_PREC = 11;
static const int LOG2_PALETTE_CG_SIZE = 4;
static const int RUN_IDX_THRE = 4;
static const int MAX_CU_BLKSIZE_PLT = 64;
static const int NUM_TRELLIS_STATE = 3;
static const double ENC_CHROMA_WEIGHTING = 0.8;
static const int MAXPLTPREDSIZE = 63;
static const int MAXPLTSIZE = 31;
static const int MAXPLTPREDSIZE_DUALTREE = 31;
static const int MAXPLTSIZE_DUALTREE = 15;
static const double PLT_CHROMA_WEIGHTING = 0.8;
static const int PLT_ENCBITDEPTH = 8;
static const int PLT_FAST_RATIO = 100;
#if RExt__DECODER_DEBUG_TOOL_MAX_FRAME_STATS
static const int EPBIN_WEIGHT_FACTOR = 4;
#endif
static const int ENC_PPS_ID_RPR = 3;
static const int SCALE_RATIO_BITS = 14;
static const int MAX_SCALING_RATIO = 2; // max downsampling ratio for RPR
static const ::std::pair<int, int> SCALE_1X = ::std::pair<int, int>(
    1 << SCALE_RATIO_BITS, 1 << SCALE_RATIO_BITS); // scale ratio 1x
static const int DELTA_QP_ACT[4] = {-5, 1, 3, 1};
static const int MAX_TSRC_RICE = 8; ///< Maximum supported TSRC Rice parameter
static const int MIN_TSRC_RICE = 1; ///< Minimum supported TSRC Rice parameter
static const int MAX_CTI_LUT_SIZE =
    64; ///< Maximum colour transform LUT size for CTI SEI

static const SplitSeries SPLIT_BITS = 5;
static const SplitSeries SPLIT_DMULT = 5;
static const SplitSeries SPLIT_MASK = 31; ///< = (1 << SPLIT_BITS) - 1

template <typename T>
inline T Clip3(const T minVal, const T maxVal, const T a) {
  return std::min<T>(std::max<T>(minVal, a), maxVal);
} ///< general min/max clip

static inline int floorLog2(uint32_t x) {
  if (x == 0) {
    // note: ceilLog2() expects -1 as return value
    return -1;
  }
#ifdef __GNUC__
  return 31 - __builtin_clz(x);
#else
#ifdef _MSC_VER
  unsigned long r = 0;
  _BitScanReverse(&r, x);
  return r;
#else
  int result = 0;
  if (x & 0xffff0000) {
    x >>= 16;
    result += 16;
  }
  if (x & 0xff00) {
    x >>= 8;
    result += 8;
  }
  if (x & 0xf0) {
    x >>= 4;
    result += 4;
  }
  if (x & 0xc) {
    x >>= 2;
    result += 2;
  }
  if (x & 0x2) {
    x >>= 1;
    result += 1;
  }
  return result;
#endif
#endif
}

typedef enum {
  AFFINEMODEL_4PARAM,
  AFFINEMODEL_6PARAM,
  AFFINE_MODEL_NUM
} EAffineModel;

struct ClpRng {
  int min{0};
  int max{0};
  int bd{0};
  int n{0};
};

struct ClpRngs {
  ClpRng comp[MAX_NUM_COMPONENT]; ///< the bit depth as indicated in the SPS
  bool used;
  bool chroma;
};

typedef int PosType;
typedef uint32_t SizeType;
struct Position {
  PosType x;
  PosType y;

  Position() : x(0), y(0) {}
  Position(const PosType _x, const PosType _y) : x(_x), y(_y) {}

  bool operator!=(const Position &other) const {
    return x != other.x || y != other.y;
  }
  bool operator==(const Position &other) const {
    return x == other.x && y == other.y;
  }

  Position offset(const Position pos) const {
    return Position(x + pos.x, y + pos.y);
  }
  Position offset(const PosType _x, const PosType _y) const {
    return Position(x + _x, y + _y);
  }
  void repositionTo(const Position newPos) {
    x = newPos.x;
    y = newPos.y;
  }
  void relativeTo(const Position origin) {
    x -= origin.x;
    y -= origin.y;
  }

  Position operator-(const Position &other) const {
    return {x - other.x, y - other.y};
  }
};

struct Size {
  SizeType width;
  SizeType height;

  Size() : width(0), height(0) {}
  Size(const SizeType _width, const SizeType _height)
      : width(_width), height(_height) {}

  bool operator!=(const Size &other) const {
    return (width != other.width) || (height != other.height);
  }
  bool operator==(const Size &other) const {
    return (width == other.width) && (height == other.height);
  }
  uint32_t area() const { return (uint32_t)width * (uint32_t)height; }
#if REUSE_CU_RESULTS_WITH_MULTIPLE_TUS
  void resizeTo(const Size newSize) {
    width = newSize.width;
    height = newSize.height;
  }
#endif
};

struct Area : public Position, public Size {
  Area() : Position(), Size() {}
  Area(const Position &_pos, const Size &_size) : Position(_pos), Size(_size) {}
  Area(const PosType _x, const PosType _y, const SizeType _w, const SizeType _h)
      : Position(_x, _y), Size(_w, _h) {}

  Position &pos() { return *this; }
  const Position &pos() const { return *this; }
  Size &size() { return *this; }
  const Size &size() const { return *this; }

  const Position &topLeft() const { return *this; }
  Position topRight() const { return {(PosType)(x + width - 1), y}; }
  Position bottomLeft() const { return {x, (PosType)(y + height - 1)}; }
  Position bottomRight() const {
    return {(PosType)(x + width - 1), (PosType)(y + height - 1)};
  }
  Position center() const {
    return {(PosType)(x + width / 2), (PosType)(y + height / 2)};
  }

  bool contains(const Position &_pos) const {
    return (_pos.x >= x) && (_pos.x < (x + width)) && (_pos.y >= y) &&
           (_pos.y < (y + height));
  }
  bool contains(const Area &_area) const {
    return contains(_area.pos()) && contains(_area.bottomRight());
  }

#if GDR_ENABLED
  bool overlap(const Area &_area) const {
    Area thisArea = Area(pos(), size());

    if (contains(_area))
      return false;

    if (_area.contains(thisArea))
      return false;

    bool topLeft = contains(_area.topLeft());
    bool topRight = contains(_area.topRight());
    bool botLeft = contains(_area.bottomLeft());
    bool botRight = contains(_area.bottomRight());

    int sum = (topLeft ? 1 : 0) + (topRight ? 1 : 0) + (botLeft ? 1 : 0) +
              (botRight ? 1 : 0);

    if (0 < sum && sum < 4) {
      return true;
    }

    return false;
  }
#endif

  bool operator!=(const Area &other) const {
    return (Size::operator!=(other)) || (Position::operator!=(other));
  }
  bool operator==(const Area &other) const {
    return (Size::operator==(other)) && (Position::operator==(other));
  }
};

inline Area clipArea(const Area &_area, const Area &boundingBox) {
  Area area = _area;

  if (area.x + area.width > boundingBox.x + boundingBox.width) {
    area.width = boundingBox.x + boundingBox.width - area.x;
  }

  if (area.y + area.height > boundingBox.y + boundingBox.height) {
    area.height = boundingBox.y + boundingBox.height - area.y;
  }

  return area;
}

struct UnitScale {
  UnitScale() : posx(0), posy(0), area(posx + posy) {}
  UnitScale(int sx, int sy) : posx(sx), posy(sy), area(posx + posy) {}
  int posx;
  int posy;
  int area;

  template <typename T> T scaleHor(const T &in) const { return in >> posx; }
  template <typename T> T scaleVer(const T &in) const { return in >> posy; }
  template <typename T> T scaleArea(const T &in) const { return in >> area; }

  Position scale(const Position &pos) const {
    return {pos.x >> posx, pos.y >> posy};
  }
  Size scale(const Size &size) const {
    return {size.width >> posx, size.height >> posy};
  }
  Area scale(const Area &_area) const {
    return Area(scale(_area.pos()), scale(_area.size()));
  }
};

class SizeIndexInfo {
public:
  SizeIndexInfo() {}
  virtual ~SizeIndexInfo() {}
  SizeType numAllWidths() { return (SizeType)m_idxToSizeTab.size(); }
  SizeType numAllHeights() { return (SizeType)m_idxToSizeTab.size(); }
  SizeType numWidths() { return (SizeType)m_numBlkSizes; }
  SizeType numHeights() { return (SizeType)m_numBlkSizes; }
  SizeType sizeFrom(SizeType idx) { return m_idxToSizeTab[idx]; }
  SizeType idxFrom(SizeType size) {
    CHECKD(m_sizeToIdxTab[size] == ::std::numeric_limits<SizeType>::max(),
           "Index of given size does NOT EXIST!");
    return m_sizeToIdxTab[size];
  }
  bool isCuSize(SizeType size) { return m_isCuSize[size]; }
  virtual void init(SizeType maxSize) {}

protected:
  void xInit() {
    m_isCuSize.resize(m_sizeToIdxTab.size(), false);

    ::std::vector<SizeType> grpSizes;

    for (int i = 0, n = 0; i < m_sizeToIdxTab.size(); i++) {
      if (m_sizeToIdxTab[i] != ::std::numeric_limits<SizeType>::max()) {
        m_sizeToIdxTab[i] = n;
        m_idxToSizeTab.push_back(i);
        n++;
      }

      if (m_sizeToIdxTab[i] != ::std::numeric_limits<SizeType>::max() &&
          m_sizeToIdxTab[i >> 1] != ::std::numeric_limits<SizeType>::max() &&
          i >= 4) {
        m_isCuSize[i] = true;
      }

      // collect group sizes (for coefficient group coding)
      SizeType grpSize = i >> ((i & 3) != 0 ? 1 : 2);
      if (m_sizeToIdxTab[i] != ::std::numeric_limits<SizeType>::max() &&
          m_sizeToIdxTab[grpSize] == ::std::numeric_limits<SizeType>::max()) {
        grpSizes.push_back(grpSize);
      }
    }

    m_numBlkSizes = (SizeType)m_idxToSizeTab.size();

    for (SizeType grpSize : grpSizes) {
      if (grpSize > 0 &&
          m_sizeToIdxTab[grpSize] == ::std::numeric_limits<SizeType>::max()) {
        m_sizeToIdxTab[grpSize] = (SizeType)m_idxToSizeTab.size();
        m_idxToSizeTab.push_back(grpSize);
      }
    }
  };

  ::std::vector<bool> m_isCuSize;
  int m_numBlkSizes; // as opposed to number all sizes, which also contains
                     // grouped sizes
  ::std::vector<SizeType> m_sizeToIdxTab;
  ::std::vector<SizeType> m_idxToSizeTab;
};

inline size_t rsAddr(const Position &pos, const uint32_t stride,
                     const UnitScale &unitScale) {
  return (size_t)(stride >> unitScale.posx) *
             (size_t)(pos.y >> unitScale.posy) +
         (size_t)(pos.x >> unitScale.posx);
}

inline size_t rsAddr(const Position &pos, const Position &origin,
                     const uint32_t stride, const UnitScale &unitScale) {
  return (stride >> unitScale.posx) * ((pos.y - origin.y) >> unitScale.posy) +
         ((pos.x - origin.x) >> unitScale.posx);
}

inline size_t rsAddr(const Position &pos, const uint32_t stride) {
  return stride * (size_t)pos.y + (size_t)pos.x;
}

inline size_t rsAddr(const Position &pos, const Position &origin,
                     const uint32_t stride) {
  return stride * (pos.y - origin.y) + (pos.x - origin.x);
}

template <typename T1, typename T2, size_t N>
inline T2 *copy_array(const T1 *src, ::std::array<T2, N> &dest) {
  return ::std::copy(src, src + dest.size(), dest.begin());
}

template <typename T1, typename T2, size_t N>
inline T1 *copy_array(const ::std::array<T2, N> &src, T1 *dest) {
  return ::std::copy(src.begin(), src.end(), dest);
}

class SizeIndexInfoLog2 : public SizeIndexInfo {
public:
  SizeIndexInfoLog2() {}
  ~SizeIndexInfoLog2(){};

  void init(SizeType maxSize) {
    for (int i = 0, n = 0; i <= maxSize; i++) {
      SizeType val = std::numeric_limits<SizeType>::max();
      if (i == (1 << n)) {
        n++;
        val = i;
      }
      m_sizeToIdxTab.push_back(val);
    }
    SizeIndexInfo::xInit();
  }
};
} // namespace Common

#endif // COMMON_COMMON_DEF