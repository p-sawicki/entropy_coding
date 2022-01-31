#include "rom.hpp"

using namespace Common;

// ====================================================================================================================
// Scanning order & context model mapping
// ====================================================================================================================
int Common::g_riceT[4] = {32, 128, 512, 2048};
int Common::g_riceShift[5] = {0, 2, 4, 6, 8};

SizeIndexInfo *Common::gp_sizeIdxInfo = NULL;

// scanning order table
ScanElement *Common::g_scanOrder[SCAN_NUMBER_OF_GROUP_TYPES][SCAN_NUMBER_OF_TYPES]
                        [MAX_CU_SIZE / 2 + 1][MAX_CU_SIZE / 2 + 1];
ScanElement Common::g_coefTopLeftDiagScan8x8[MAX_CU_SIZE / 2 + 1][64];

const uint32_t Common::g_minInGroup[LAST_SIGNIFICANT_GROUPS] = {
    0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96};

const uint32_t Common::g_groupIdx[MAX_TB_SIZEY] = {
    0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,
    8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11};

const uint32_t Common::g_goRiceParsCoeff[32] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
                                        1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                                        2, 2, 2, 2, 2, 2, 3, 3, 3, 3};

const int8_t Common::g_BcwLog2WeightBase = 3;
const int8_t Common::g_BcwWeightBase = (1 << g_BcwLog2WeightBase);
const int8_t Common::g_BcwWeights[BCW_NUM] = {-2, 3, 4, 5, 10};
const int8_t Common::g_BcwSearchOrder[BCW_NUM] = {BCW_DEFAULT, BCW_DEFAULT - 2,
                                          BCW_DEFAULT + 2, BCW_DEFAULT - 1,
                                          BCW_DEFAULT + 1};
int8_t Common::g_BcwCodingOrder[BCW_NUM];
int8_t Common::g_BcwParsingOrder[BCW_NUM];

uint32_t Common::g_log2SbbSize[MAX_CU_DEPTH + 1][MAX_CU_DEPTH + 1][2] =
    //===== luma/chroma =====
    {{{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 4}, {0, 4}, {0, 4}},
     {{1, 0}, {1, 1}, {1, 1}, {1, 3}, {1, 3}, {1, 3}, {1, 3}, {1, 3}},
     {{2, 0}, {1, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}},
     {{3, 0}, {3, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}},
     {{4, 0}, {3, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}},
     {{4, 0}, {3, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}},
     {{4, 0}, {3, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}},
     {{4, 0}, {3, 1}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}, {2, 2}}};

uint8_t Common::g_paletteRunTopLut[5] = {0, 1, 1, 2, 2};
uint8_t Common::g_paletteRunLeftLut[5] = {0, 1, 2, 3, 4};

class ScanGenerator {
private:
  uint32_t m_line, m_column;
  const uint32_t m_blockWidth, m_blockHeight;
  const uint32_t m_stride;
  const CoeffScanType m_scanType;

public:
  ScanGenerator(uint32_t blockWidth, uint32_t blockHeight, uint32_t stride,
                CoeffScanType scanType)
      : m_line(0), m_column(0), m_blockWidth(blockWidth),
        m_blockHeight(blockHeight), m_stride(stride), m_scanType(scanType) {}

  uint32_t GetCurrentX() const { return m_column; }
  uint32_t GetCurrentY() const { return m_line; }

  uint32_t GetNextIndex(uint32_t blockOffsetX, uint32_t blockOffsetY) {
    const uint32_t rtn =
        ((m_line + blockOffsetY) * m_stride) + m_column + blockOffsetX;

    // advance line and column to the next position
    switch (m_scanType) {
      //------------------------------------------------

    case SCAN_DIAG:

      if ((m_column == m_blockWidth - 1) ||
          (m_line == 0)) // if we reach the end of a rank, go diagonally down to
                         // the next one
      {
        m_line += m_column + 1;
        m_column = 0;

        if (m_line >=
            m_blockHeight) // if that takes us outside the block, adjust so that
                           // we are back on the bottom row
        {
          m_column += m_line - (m_blockHeight - 1);
          m_line = m_blockHeight - 1;
        }
      } else {
        m_column++;
        m_line--;
      }
      break;

    case SCAN_TRAV_HOR:
      if (m_line % 2 == 0) {
        if (m_column == (m_blockWidth - 1)) {
          m_line++;
          m_column = m_blockWidth - 1;
        } else {
          m_column++;
        }
      } else {
        if (m_column == 0) {
          m_line++;
          m_column = 0;
        } else {
          m_column--;
        }
      }
      break;

    case SCAN_TRAV_VER:
      if (m_column % 2 == 0) {
        if (m_line == (m_blockHeight - 1)) {
          m_column++;
          m_line = m_blockHeight - 1;
        } else {
          m_line++;
        }
      } else {
        if (m_line == 0) {
          m_column++;
          m_line = 0;
        } else {
          m_line--;
        }
      }
      break;
      //------------------------------------------------

    default:

      THROW("ERROR: Unknown scan type \""
            << m_scanType << "\"in ScanGenerator::GetNextIndex");
      break;
    }

    return rtn;
  }
};

void Common::initROM() {
  gp_sizeIdxInfo = new SizeIndexInfoLog2();
  gp_sizeIdxInfo->init(MAX_CU_SIZE);

  SizeIndexInfoLog2 sizeInfo;
  sizeInfo.init(MAX_CU_SIZE);

  // initialize scan orders
  for (uint32_t blockHeightIdx = 0; blockHeightIdx < sizeInfo.numAllHeights();
       blockHeightIdx++) {
    for (uint32_t blockWidthIdx = 0; blockWidthIdx < sizeInfo.numAllWidths();
         blockWidthIdx++) {
      const uint32_t blockWidth = sizeInfo.sizeFrom(blockWidthIdx);
      const uint32_t blockHeight = sizeInfo.sizeFrom(blockHeightIdx);
      const uint32_t totalValues = blockWidth * blockHeight;

      //--------------------------------------------------------------------------------------------------

      // non-grouped scan orders

      for (uint32_t scanTypeIndex = 0; scanTypeIndex < SCAN_NUMBER_OF_TYPES;
           scanTypeIndex++) {
        const CoeffScanType scanType = CoeffScanType(scanTypeIndex);
        ScanElement *scan = nullptr;

        if (blockWidthIdx < sizeInfo.numWidths() &&
            blockHeightIdx < sizeInfo.numHeights()) {
          scan = new ScanElement[totalValues];
        }

        g_scanOrder[SCAN_UNGROUPED][scanType][blockWidthIdx][blockHeightIdx] =
            scan;

        if (scan == nullptr) {
          continue;
        }

        ScanGenerator fullBlockScan(blockWidth, blockHeight, blockWidth,
                                    scanType);

        for (uint32_t scanPosition = 0; scanPosition < totalValues;
             scanPosition++) {
          const int rasterPos = fullBlockScan.GetNextIndex(0, 0);
          const int posY = rasterPos / blockWidth;
          const int posX = rasterPos - (posY * blockWidth);

          scan[scanPosition].idx = rasterPos;
          scan[scanPosition].x = posX;
          scan[scanPosition].y = posY;
        }
      }

      //--------------------------------------------------------------------------------------------------

      // grouped scan orders
      const uint32_t *log2Sbb =
          g_log2SbbSize[floorLog2(blockWidth)][floorLog2(blockHeight)];
      const uint32_t log2CGWidth = log2Sbb[0];
      const uint32_t log2CGHeight = log2Sbb[1];

      const uint32_t groupWidth = 1 << log2CGWidth;
      const uint32_t groupHeight = 1 << log2CGHeight;
      const uint32_t widthInGroups =
          std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, blockWidth) >> log2CGWidth;
      const uint32_t heightInGroups =
          std::min<unsigned>(JVET_C0024_ZERO_OUT_TH, blockHeight) >>
          log2CGHeight;

      const uint32_t groupSize = groupWidth * groupHeight;
      const uint32_t totalGroups = widthInGroups * heightInGroups;

      for (uint32_t scanTypeIndex = 0; scanTypeIndex < SCAN_NUMBER_OF_TYPES;
           scanTypeIndex++) {
        const CoeffScanType scanType = CoeffScanType(scanTypeIndex);

        ScanElement *scan = new ScanElement[totalValues];

        g_scanOrder[SCAN_GROUPED_4x4][scanType][blockWidthIdx][blockHeightIdx] =
            scan;

        if (blockWidth > JVET_C0024_ZERO_OUT_TH ||
            blockHeight > JVET_C0024_ZERO_OUT_TH) {
          for (uint32_t i = 0; i < totalValues; i++) {
            scan[i].idx = totalValues - 1;
            scan[i].x = blockWidth - 1;
            scan[i].y = blockHeight - 1;
          }
        }

        ScanGenerator fullBlockScan(widthInGroups, heightInGroups, groupWidth,
                                    scanType);

        for (uint32_t groupIndex = 0; groupIndex < totalGroups; groupIndex++) {
          const uint32_t groupPositionY = fullBlockScan.GetCurrentY();
          const uint32_t groupPositionX = fullBlockScan.GetCurrentX();
          const uint32_t groupOffsetX = groupPositionX * groupWidth;
          const uint32_t groupOffsetY = groupPositionY * groupHeight;
          const uint32_t groupOffsetScan = groupIndex * groupSize;

          ScanGenerator groupScan(groupWidth, groupHeight, blockWidth,
                                  scanType);

          for (uint32_t scanPosition = 0; scanPosition < groupSize;
               scanPosition++) {
            const int rasterPos =
                groupScan.GetNextIndex(groupOffsetX, groupOffsetY);
            const int posY = rasterPos / blockWidth;
            const int posX = rasterPos - (posY * blockWidth);

            scan[groupOffsetScan + scanPosition].idx = rasterPos;
            scan[groupOffsetScan + scanPosition].x = posX;
            scan[groupOffsetScan + scanPosition].y = posY;
          }

          fullBlockScan.GetNextIndex(0, 0);
        }
      }

      //--------------------------------------------------------------------------------------------------
    }
  }

  // initialize CoefTopLeftDiagScan8x8 for LFNST
  for (uint32_t blockWidthIdx = 0; blockWidthIdx < sizeInfo.numAllWidths();
       blockWidthIdx++) {
    const uint32_t blockWidth = sizeInfo.sizeFrom(blockWidthIdx);

    const static uint8_t g_auiXYDiagScan8x8[64][2] = {
        {0, 0}, {0, 1}, {1, 0}, {0, 2}, {1, 1}, {2, 0}, {0, 3}, {1, 2},
        {2, 1}, {3, 0}, {1, 3}, {2, 2}, {3, 1}, {2, 3}, {3, 2}, {3, 3},
        {0, 4}, {0, 5}, {1, 4}, {0, 6}, {1, 5}, {2, 4}, {0, 7}, {1, 6},
        {2, 5}, {3, 4}, {1, 7}, {2, 6}, {3, 5}, {2, 7}, {3, 6}, {3, 7},
        {4, 0}, {4, 1}, {5, 0}, {4, 2}, {5, 1}, {6, 0}, {4, 3}, {5, 2},
        {6, 1}, {7, 0}, {5, 3}, {6, 2}, {7, 1}, {6, 3}, {7, 2}, {7, 3},
        {4, 4}, {4, 5}, {5, 4}, {4, 6}, {5, 5}, {6, 4}, {4, 7}, {5, 6},
        {6, 5}, {7, 4}, {5, 7}, {6, 6}, {7, 5}, {6, 7}, {7, 6}, {7, 7}};
    for (int i = 0; i < 64; i++) {
      g_coefTopLeftDiagScan8x8[blockWidthIdx][i].idx =
          g_auiXYDiagScan8x8[i][0] + g_auiXYDiagScan8x8[i][1] * blockWidth;
      g_coefTopLeftDiagScan8x8[blockWidthIdx][i].x = g_auiXYDiagScan8x8[i][0];
      g_coefTopLeftDiagScan8x8[blockWidthIdx][i].y = g_auiXYDiagScan8x8[i][1];
    }
  }
}

void Common::destroyROM() {
  unsigned numWidths = gp_sizeIdxInfo->numAllWidths();
  unsigned numHeights = gp_sizeIdxInfo->numAllHeights();

  for (uint32_t groupTypeIndex = 0; groupTypeIndex < SCAN_NUMBER_OF_GROUP_TYPES;
       groupTypeIndex++) {
    for (uint32_t scanOrderIndex = 0; scanOrderIndex < SCAN_NUMBER_OF_TYPES;
         scanOrderIndex++) {
      for (uint32_t blockWidthIdx = 0; blockWidthIdx <= numWidths;
           blockWidthIdx++) {
        for (uint32_t blockHeightIdx = 0; blockHeightIdx <= numHeights;
             blockHeightIdx++) {
          delete[] g_scanOrder[groupTypeIndex][scanOrderIndex][blockWidthIdx]
                              [blockHeightIdx];
          g_scanOrder[groupTypeIndex][scanOrderIndex][blockWidthIdx]
                     [blockHeightIdx] = nullptr;
        }
      }
    }
  }

  delete gp_sizeIdxInfo;
  gp_sizeIdxInfo = nullptr;
}
