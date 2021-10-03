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

/** \file     Picture.h
 *  \brief    Description of a coded picture
 */

#ifndef __PICTURE__
#define __PICTURE__

#include <deque>

#include "buffer.hpp"
#include "coding_structure.hpp"
#include "hash.hpp"
#include "mcts.hpp"
#include "sei_color_transform.hpp"
#include "slice.hpp"
#include "unit.hpp"

class SEI;
class AQpLayer;

typedef std::list<SEI *> SEIMessages;

#define M_BUFS(JID, PID) m_bufs[PID]

struct Picture : public UnitArea {
  CodingStructure *cs;

  SAOBlkParam *getSAO(int id = 0) { return &m_sao[id][0]; };

  std::vector<SAOBlkParam> m_sao[2];

  std::vector<uint8_t> m_alfCtuEnableFlag[MAX_NUM_COMPONENT];
  uint8_t *getAlfCtuEnableFlag(int compIdx) {
    return m_alfCtuEnableFlag[compIdx].data();
  }

  std::vector<short> m_alfCtbFilterIndex;
  short *getAlfCtbFilterIndex() { return m_alfCtbFilterIndex.data(); }

  std::vector<uint8_t> m_alfCtuAlternative[MAX_NUM_COMPONENT];
  uint8_t *getAlfCtuAlternativeData(int compIdx) {
    return m_alfCtuAlternative[compIdx].data();
  }
};

#endif
