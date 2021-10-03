#include "coding_structure.hpp"

uint32_t getCtuAddr        (const Position& pos, const PreCalcValues &pcv);


// CS tools
namespace CS {
UnitArea getArea(const CodingStructure &cs, const UnitArea &area,
                 const ChannelType chType);
bool isDualITree(const CodingStructure &cs);
} // namespace CS

// CU tools
namespace CU {
bool isIntra(const CodingUnit &cu);
bool isInter(const CodingUnit &cu);
bool isIBC(const CodingUnit &cu);
bool isPLT(const CodingUnit &cu);

bool isSameCtu(const CodingUnit &cu, const CodingUnit &cu2);
bool isSameSliceAndTile(const CodingUnit &cu, const CodingUnit &cu2);
bool isLastSubCUOfCtu(const CodingUnit &cu);
uint32_t getCtuAddr(const CodingUnit &cu);
int predictQP(const CodingUnit &cu, const int prevQP);

uint32_t getNumPUs(const CodingUnit &cu);
PartSplit getSplitAtDepth(const CodingUnit &cu, const unsigned depth);
ModeType getModeTypeAtDepth(const CodingUnit &cu, const unsigned depth);

bool isBcwIdxCoded(const CodingUnit &cu);
uint8_t getValidBcwIdx(const CodingUnit &cu);
void setBcwIdx(CodingUnit &cu, uint8_t uh);
bool bdpcmAllowed(const CodingUnit &cu, const ComponentID compID);
bool isMTSAllowed(const CodingUnit &cu, const ComponentID compID);

bool divideTuInRows(const CodingUnit &cu);
PartSplit getISPType(const CodingUnit &cu, const ComponentID compID);
bool isISPFirst(const CodingUnit &cu, const CompArea &tuArea,
                const ComponentID compID);
bool canUseISP(const CodingUnit &cu, const ComponentID compID);
bool canUseISP(const int width, const int height,
               const int maxTrSize = MAX_TB_SIZEY);
bool canUseLfnstWithISP(const CompArea &cuArea, const ISPType ispSplitType);
bool canUseLfnstWithISP(const CodingUnit &cu, const ChannelType chType);
uint32_t getISPSplitDim(const int width, const int height,
                        const PartSplit ispType);

PUTraverser traversePUs(CodingUnit &cu);
TUTraverser traverseTUs(CodingUnit &cu);
cPUTraverser traversePUs(const CodingUnit &cu);
cTUTraverser traverseTUs(const CodingUnit &cu);

bool hasSubCUNonZeroMVd(const CodingUnit &cu);
bool hasSubCUNonZeroAffineMVd(const CodingUnit &cu);

uint8_t getSbtIdx(const uint8_t sbtInfo);
uint8_t getSbtPos(const uint8_t sbtInfo);
uint8_t targetSbtAllowed(uint8_t idx, uint8_t sbtAllowed);
} // namespace CU

namespace PU {
int getLMSymbolList(const PredictionUnit &pu, int *modeList);
int getIntraMPMs(const PredictionUnit &pu, unsigned *mpm,
                 const ChannelType &channelType = CHANNEL_TYPE_LUMA);
bool isMIP(const PredictionUnit &pu,
           const ChannelType &chType = CHANNEL_TYPE_LUMA);
bool isDMChromaMIP(const PredictionUnit &pu);
uint32_t getIntraDirLuma(const PredictionUnit &pu);
void getIntraChromaCandModes(const PredictionUnit &pu,
                             unsigned modeList[NUM_CHROMA_MODE]);
const PredictionUnit &getCoLocatedLumaPU(const PredictionUnit &pu);
bool isBipredRestriction(const PredictionUnit &pu);

bool isLMCMode(unsigned mode);
} // namespace PU

// TU tools
namespace TU {
bool getCbf(const TransformUnit &tu, const ComponentID &compID);
bool getCbfAtDepth(const TransformUnit &tu, const ComponentID &compID,
                   const unsigned &depth);
void setCbfAtDepth(TransformUnit &tu, const ComponentID &compID,
                   const unsigned &depth, const bool &cbf);
bool isTSAllowed(const TransformUnit &tu, const ComponentID compID);

TransformUnit *getPrevTU(const TransformUnit &tu, const ComponentID compID);
bool getPrevTuCbfAtDepth(const TransformUnit &tu, const ComponentID compID,
                         const int trDepth);
} // namespace TU

bool allowLfnstWithMip(const Size &block) {
  if (block.width >= 16 && block.height >= 16) {
    return true;
  }
  return false;
}

int getNumModesMip(const Size &block);
int getMipSizeId(const Size &block);