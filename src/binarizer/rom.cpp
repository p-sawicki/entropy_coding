#include "rom.hpp"

// ====================================================================================================================
// Scanning order & context model mapping
// ====================================================================================================================
int g_riceT[4] = {32, 128, 512, 2048};
int g_riceShift[5] = {0, 2, 4, 6, 8};

SizeIndexInfo *gp_sizeIdxInfo = NULL;

// scanning order table
ScanElement *g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES]
                        [MAX_CU_SIZE / 2 + 1][MAX_CU_SIZE / 2 + 1];
ScanElement g_coefTopLeftDiagScan8x8[MAX_CU_SIZE / 2 + 1][64];

const uint32_t g_minInGroup[LAST_SIGNIFICANT_GROUPS] = {
    0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96};

const uint32_t g_groupIdx[MAX_TB_SIZEY] = {
    0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,
    8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11};

const uint32_t g_goRiceParsCoeff[32] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                                        1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                                        2, 2, 2, 2, 2, 2, 3, 3, 3, 3};

const int8_t g_BcwLog2WeightBase = 3;
const int8_t g_BcwWeightBase = (1 << g_BcwLog2WeightBase);
const int8_t g_BcwWeights[BCW_NUM] = {-2, 3, 4, 5, 10};
const int8_t g_BcwSearchOrder[BCW_NUM] = {BCW_DEFAULT, BCW_DEFAULT - 2,
                                          BCW_DEFAULT + 2, BCW_DEFAULT - 1,
                                          BCW_DEFAULT + 1};
int8_t g_BcwCodingOrder[BCW_NUM];
int8_t g_BcwParsingOrder[BCW_NUM];