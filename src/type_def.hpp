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

#ifndef TYPE_DEF_HPP
#define TYPE_DEF_HPP

#include <sstream>
#include <stdexcept>


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

class Exception : public std::exception {
public:
  Exception(const std::string &_s) : m_str(_s) {}
  Exception(const Exception &_e) : std::exception(_e), m_str(_e.m_str) {}
  virtual ~Exception() noexcept {};
  virtual const char *what() const noexcept { return m_str.c_str(); }
  Exception &operator=(const Exception &_e) {
    std::exception::operator=(_e);
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
  std::string m_str;
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

#endif // TYPE_DEF_HPP