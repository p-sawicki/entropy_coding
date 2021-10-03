/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2021, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     CodingStructure.h
 *  \brief    A class managing the coding information for a specific image part
 */

#ifndef __CODINGSTRUCTURE__
#define __CODINGSTRUCTURE__

#include <vector>

#include "buffer.hpp"
#include "slice.hpp"
#include "unit.hpp"
#include "unit_partitioner.hpp"

struct Picture;

enum PictureType {
  PIC_RECONSTRUCTION = 0,
  PIC_ORIGINAL,
  PIC_TRUE_ORIGINAL,
  PIC_FILTERED_ORIGINAL,
  PIC_PREDICTION,
  PIC_RESIDUAL,
  PIC_ORG_RESI,
  PIC_RECON_WRAP,
  PIC_ORIGINAL_INPUT,
  PIC_TRUE_ORIGINAL_INPUT,
  PIC_FILTERED_ORIGINAL_INPUT,
  NUM_PIC_TYPES
};

// ---------------------------------------------------------------------------
// coding structure
// ---------------------------------------------------------------------------

class CodingStructure {
public:
  UnitArea area;

  Picture *picture;
  CodingStructure *parent;
  Slice *slice;

  UnitScale unitScale[MAX_NUM_COMPONENT];
  int chromaQpAdj;
  const SPS *sps;
  const PPS *pps;
  PicHeader *picHeader;
  const PreCalcValues *pcv;

  const CodingUnit *getCU(const Position &pos, const ChannelType _chType) const;
  const PredictionUnit *getPU(const Position &pos,
                              const ChannelType _chType) const;
  const TransformUnit *getTU(const Position &pos, const ChannelType _chType,
                             const int subTuIdx = -1) const;

  CodingUnit *getCU(const Position &pos, const ChannelType _chType);
  CodingUnit *getLumaCU(const Position &pos);
  PredictionUnit *getPU(const Position &pos, const ChannelType _chType);
  TransformUnit *getTU(const Position &pos, const ChannelType _chType,
                       const int subTuIdx = -1);

  const CodingUnit *getCURestricted(const Position &pos, const Position curPos,
                                    const unsigned curSliceIdx,
                                    const unsigned curTileIdx,
                                    const ChannelType _chType) const;
  const CodingUnit *getCURestricted(const Position &pos,
                                    const CodingUnit &curCu,
                                    const ChannelType _chType) const;
  const PredictionUnit *getPURestricted(const Position &pos,
                                        const PredictionUnit &curPu,
                                        const ChannelType _chType) const;

  CodingUnit &addCU(const UnitArea &unit, const ChannelType _chType);
  PredictionUnit &addPU(const UnitArea &unit, const ChannelType _chType);
  TransformUnit &addTU(const UnitArea &unit, const ChannelType _chType);
  void addEmptyTUs(Partitioner &partitioner);

  TreeType treeType; // because partitioner can not go deep to tu and cu coding
                     // (e.g., addCU()), need another variable for indicating
                     // treeType
  ModeType modeType;

  const int signalModeCons(const PartSplit split, Partitioner &partitioner,
                           const ModeType modeTypeParent) const;
  std::vector<CodingUnit *> cus;
  std::vector<PredictionUnit *> pus;
  std::vector<TransformUnit *> tus;

  PLTBuf prevPLT;
  void reorderPrevPLT(PLTBuf &prevPLT, uint8_t curPLTSize[MAX_NUM_CHANNEL_TYPE],
                      Pel curPLT[MAX_NUM_COMPONENT][MAXPLTSIZE],
                      bool reuseflag[MAX_NUM_CHANNEL_TYPE][MAXPLTPREDSIZE],
                      uint32_t compBegin, uint32_t numComp, bool jointPLT);

 private:
   // needed for TU encoding
  bool m_isTuEnc;

  unsigned *m_cuIdx[MAX_NUM_CHANNEL_TYPE];
  unsigned *m_puIdx[MAX_NUM_CHANNEL_TYPE];
  unsigned *m_tuIdx[MAX_NUM_CHANNEL_TYPE];

  unsigned m_numCUs;
  unsigned m_numPUs;
  unsigned m_numTUs;

  CUCache &m_cuCache;
  PUCache &m_puCache;
  TUCache &m_tuCache;

  TCoeff *m_coeffs[MAX_NUM_COMPONENT];
  Pel *m_pcmbuf[MAX_NUM_COMPONENT];
  bool *m_runType[MAX_NUM_CHANNEL_TYPE];
  int m_offsets[MAX_NUM_COMPONENT];
};

static inline uint32_t getNumberValidTBlocks(const PreCalcValues &pcv) {
  return (pcv.chrFormat == CHROMA_400)
             ? 1
             : (pcv.multiBlock422 ? MAX_NUM_TBLOCKS : MAX_NUM_COMPONENT);
}

#endif
