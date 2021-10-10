#ifndef ENTROPY_CODEC_UNIT_PARTITIONER
#define ENTROPY_CODEC_UNIT_PARTITIONER

#include <vector>

#include "unit.hpp"

namespace EntropyCoding {

static_assert(
    MAX_CU_TILING_PARTITIONS >= 4,
    "Minimum required number of partitions for the Partitioning type is 4!");
typedef ::std::vector<UnitArea> Partitioning;

enum PartSplit {
  CTU_LEVEL = 0,
  CU_QUAD_SPLIT,

  CU_HORZ_SPLIT,
  CU_VERT_SPLIT,
  CU_TRIH_SPLIT,
  CU_TRIV_SPLIT,
  TU_MAX_TR_SPLIT,
  TU_NO_ISP,
  TU_1D_HORZ_SPLIT,
  TU_1D_VERT_SPLIT,
  SBT_VER_HALF_POS0_SPLIT,
  SBT_VER_HALF_POS1_SPLIT,
  SBT_HOR_HALF_POS0_SPLIT,
  SBT_HOR_HALF_POS1_SPLIT,
  SBT_VER_QUAD_POS0_SPLIT,
  SBT_VER_QUAD_POS1_SPLIT,
  SBT_HOR_QUAD_POS0_SPLIT,
  SBT_HOR_QUAD_POS1_SPLIT,
  NUM_PART_SPLIT,
  CU_MT_SPLIT =
      1000, ///< dummy element to indicate the MT (multi-type-tree) split
  CU_BT_SPLIT = 1001,  ///< dummy element to indicate the BT split
  CU_DONT_SPLIT = 2000 ///< dummy element to indicate no splitting
};

struct PartLevel {
  PartSplit split;
  Partitioning parts;
  unsigned idx;
  bool checkdIfImplicit;
  bool isImplicit;
  PartSplit implicitSplit;
  PartSplit firstSubPartSplit;
  bool canQtSplit;
  bool qgEnable;
  bool qgChromaEnable;
  int modeType;

  PartLevel();
  PartLevel(const PartSplit _split, const Partitioning &_parts);
  PartLevel(const PartSplit _split, Partitioning &&_parts);
  PartLevel(const PartSplit _split, const Partitioning &_parts,
            const unsigned _idx, const bool _checkdIfImplicit,
            const bool _isImplicit, const PartSplit _implicitSplit,
            const PartSplit _firstSubPartSplit, const bool _canQtSplit,
            const bool _qgEnable, const bool _qgChromaEnable,
            const int _modeType)
      : split(_split), parts(_parts), idx(_idx),
        checkdIfImplicit(_checkdIfImplicit), isImplicit(_isImplicit),
        implicitSplit(_implicitSplit), firstSubPartSplit(_firstSubPartSplit),
        canQtSplit(_canQtSplit), qgEnable(_qgEnable),
        qgChromaEnable(_qgChromaEnable), modeType(_modeType) {}
};

// set depending on max QT / BT possibilities
typedef static_vector<PartLevel, 2 * MAX_CU_DEPTH + 1> PartitioningStack;

class Partitioner {
protected:
  PartitioningStack m_partStack;
#if _DEBUG
  UnitArea m_currArea;
#endif

public:
  unsigned currDepth;
  unsigned currQtDepth;
  unsigned currTrDepth;
  unsigned currBtDepth;
  unsigned currMtDepth;
  unsigned currSubdiv;
  Position currQgPos;
  Position currQgChromaPos;

  unsigned currImplicitBtDepth;
  ChannelType chType;
  TreeType treeType;
  ModeType modeType;

  Partitioner() {}
  Partitioner(const PartitioningStack &partStack, const unsigned _currDepth,
              const unsigned _currQtDepth, const unsigned _currTrDepth,
              const unsigned _currBtDepth, const unsigned _currMtDepth,
              const unsigned _currSubdiv, const Position &_currQgPos,
              const Position &_currQgChromaPos,
              const unsigned _currImplicitBtDepth, const ChannelType _chType,
              const TreeType _treeType, const ModeType _modeType)
      : m_partStack(partStack), currDepth(_currDepth),
        currQtDepth(_currQtDepth), currTrDepth(_currTrDepth),
        currBtDepth(_currBtDepth), currMtDepth(_currMtDepth),
        currSubdiv(_currSubdiv), currQgPos(_currQgPos),
        currQgChromaPos(_currQgChromaPos),
        currImplicitBtDepth(_currImplicitBtDepth), chType(_chType),
        treeType(_treeType), modeType(_modeType) {}

  virtual ~Partitioner() {}

  const PartLevel &currPartLevel() const { return m_partStack.back(); }
  const UnitArea &currArea() const {
    return currPartLevel().parts[currPartIdx()];
  }
  const unsigned currPartIdx() const { return currPartLevel().idx; }
  const PartitioningStack &getPartStack() const { return m_partStack; }
  const bool currQgEnable() const { return currPartLevel().qgEnable; }
  const bool currQgChromaEnable() const {
    return currPartLevel().qgChromaEnable;
  }

  SplitSeries getSplitSeries() const;
  ModeTypeSeries getModeTypeSeries() const;

  virtual void initCtu(const UnitArea &ctuArea, const ChannelType _chType,
                       const Slice &slice) = 0;
  virtual void splitCurrArea(const PartSplit split,
                             const CodingStructure &cs) = 0;
  virtual void exitCurrSplit() = 0;
  virtual bool nextPart(const CodingStructure &cs, bool autoPop = false) = 0;

  virtual void setCUData(CodingUnit &cu);

public:
  virtual void canSplit(const CodingStructure &cs, bool &canNo, bool &canQt,
                        bool &canBh, bool &canBv, bool &canTh, bool &canTv) = 0;
  virtual bool canSplit(const PartSplit split, const CodingStructure &cs) = 0;
  virtual bool isSplitImplicit(const PartSplit split,
                               const CodingStructure &cs) = 0;
  virtual PartSplit getImplicitSplit(const CodingStructure &cs) = 0;
  bool isSepTree(const CodingStructure &cs);
  bool isConsInter() { return modeType == MODE_TYPE_INTER; }
  bool isConsIntra() { return modeType == MODE_TYPE_INTRA; }
};

class QTBTPartitioner : public Partitioner {
public:
  QTBTPartitioner() : Partitioner() {}
  
  QTBTPartitioner(const PartitioningStack &partStack, const unsigned _currDepth,
                  const unsigned _currQtDepth, const unsigned _currTrDepth,
                  const unsigned _currBtDepth, const unsigned _currMtDepth,
                  const unsigned _currSubdiv, const Position &_currQgPos,
                  const Position &_currQgChromaPos,
                  const unsigned _currImplicitBtDepth,
                  const ChannelType _chType, const TreeType _treeType,
                  const ModeType _modeType)
      : Partitioner(partStack, _currDepth, _currQtDepth, _currTrDepth,
                    _currBtDepth, _currMtDepth, _currSubdiv, _currQgPos,
                    _currQgChromaPos, _currImplicitBtDepth, _chType, _treeType,
                    _modeType) {}

  void initCtu(const UnitArea &ctuArea, const ChannelType _chType,
               const Slice &slice);
  void splitCurrArea(const PartSplit split, const CodingStructure &cs);
  void exitCurrSplit();
  bool nextPart(const CodingStructure &cs, bool autoPop = false);

  void canSplit(const CodingStructure &cs, bool &canNo, bool &canQt,
                bool &canBh, bool &canBv, bool &canTh, bool &canTv);
  bool canSplit(const PartSplit split, const CodingStructure &cs);
  bool isSplitImplicit(const PartSplit split, const CodingStructure &cs);
  PartSplit getImplicitSplit(const CodingStructure &cs);
};

class TUIntraSubPartitioner : public Partitioner {
public:
  TUIntraSubPartitioner(Partitioner &_initialState) {
    // we copy the input partitioner data
    m_partStack.push_back(PartLevel(TU_NO_ISP, {_initialState.currArea()}));

    currDepth = _initialState.currDepth;
    currQtDepth = _initialState.currQtDepth;
    currTrDepth = _initialState.currTrDepth;
    currBtDepth = _initialState.currBtDepth;
    currMtDepth = _initialState.currMtDepth;
    chType = _initialState.chType;
#if _DEBUG
    m_currArea = _initialState.currArea();
#endif
    treeType = _initialState.treeType;
    modeType = _initialState.modeType;
  }
  TUIntraSubPartitioner(const PartitioningStack &partStack,
                        const unsigned _currDepth, const unsigned _currQtDepth,
                        const unsigned _currTrDepth,
                        const unsigned _currBtDepth,
                        const unsigned _currMtDepth, const unsigned _currSubdiv,
                        const Position &_currQgPos,
                        const Position &_currQgChromaPos,
                        const unsigned _currImplicitBtDepth,
                        const ChannelType _chType, const TreeType _treeType,
                        const ModeType _modeType)
      : Partitioner(partStack, _currDepth, _currQtDepth, _currTrDepth,
                    _currBtDepth, _currMtDepth, _currSubdiv, _currQgPos,
                    _currQgChromaPos, _currImplicitBtDepth, _chType, _treeType,
                    _modeType) {}

  void initCtu(const UnitArea &ctuArea, const ChannelType chType,
               const Slice &slice){}; // not needed
  void splitCurrArea(const PartSplit split, const CodingStructure &cs);
  void exitCurrSplit();
  bool nextPart(const CodingStructure &cs, bool autoPop = false);
  void canSplit(const CodingStructure &cs, bool &canNo, bool &canQt,
                bool &canBh, bool &canBv, bool &canTh, bool &canTv){};
  bool canSplit(const PartSplit split, const CodingStructure &cs);
  bool isSplitImplicit(const PartSplit split, const CodingStructure &cs) {
    return false;
  }; // not needed
  PartSplit getImplicitSplit(const CodingStructure &cs) {
    return CU_DONT_SPLIT;
  }; // not needed
};

//////////////////////////////////////////////////////////////////////////
// Partitioner namespace - contains methods calculating the actual splits
//////////////////////////////////////////////////////////////////////////

namespace PartitionerImpl {
Partitioning getCUSubPartitions(const UnitArea &cuArea,
                                const CodingStructure &cs,
                                const PartSplit splitType = CU_QUAD_SPLIT);
Partitioning getMaxTuTiling(const UnitArea &curArea, const CodingStructure &cs);
void getTUIntraSubPartitions(Partitioning &sub, const UnitArea &tuArea,
                             const CodingStructure &cs,
                             const PartSplit splitType);
Partitioning getSbtTuTiling(const UnitArea &curArea, const CodingStructure &cs,
                            const PartSplit splitType);
}; // namespace PartitionerImpl
} // namespace EntropyCoding

#endif
