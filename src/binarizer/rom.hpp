#include "common_def.hpp"

// ====================================================================================================================
// Scanning order & context mapping table
// ====================================================================================================================
extern int g_riceT[4];
extern int g_riceShift[5];
extern const uint32_t g_groupIdx[MAX_TB_SIZEY];
extern const uint32_t g_minInGroup[LAST_SIGNIFICANT_GROUPS];
extern const uint32_t g_goRiceParsCoeff[32];
inline uint32_t g_goRicePosCoeff0(int st, uint32_t ricePar) {
  return (st < 2 ? 1 : 2) << ricePar;
}

extern SizeIndexInfo *gp_sizeIdxInfo;

extern ScanElement
    *g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES]
                [MAX_CU_SIZE / 2 + 1][MAX_CU_SIZE / 2 + 1];
extern ScanElement g_coefTopLeftDiagScan8x8[MAX_CU_SIZE / 2 + 1][64];

extern const int8_t g_BcwLog2WeightBase;
extern const int8_t g_BcwWeightBase;
extern const int8_t g_BcwWeights[BCW_NUM];
extern const int8_t g_BcwSearchOrder[BCW_NUM];
extern int8_t g_BcwCodingOrder[BCW_NUM];
extern int8_t g_BcwParsingOrder[BCW_NUM];

constexpr uint8_t g_tbMax[257] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8};

// flexible conversion from relative to absolute index
struct ScanElement {
  uint32_t idx;
  uint16_t x;
  uint16_t y;
};

extern       uint32_t   g_log2SbbSize[MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1][2];
extern uint8_t g_paletteRunTopLut[5];
extern uint8_t g_paletteRunLeftLut[5];