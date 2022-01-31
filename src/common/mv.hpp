#ifndef COMMON_MV
#define COMMON_MV

#include "common_def.hpp"

namespace Common {

// ====================================================================================================================
// Class definition
// ====================================================================================================================

enum MvPrecision {
  MV_PRECISION_4PEL = 0,    // 4-pel
  MV_PRECISION_INT = 2,     // 1-pel, shift 2 bits from 4-pel
  MV_PRECISION_HALF = 3,    // 1/2-pel
  MV_PRECISION_QUARTER = 4, // 1/4-pel (the precision of regular MV difference
                            // signaling), shift 4 bits from 4-pel
  MV_PRECISION_SIXTEENTH =
      6, // 1/16-pel (the precision of internal MV), shift 6 bits from 4-pel
  MV_PRECISION_INTERNAL = 2 + MV_FRACTIONAL_BITS_INTERNAL,
};

/// basic motion vector class
class Mv {
private:
  static const MvPrecision m_amvrPrecision[4];
  static const MvPrecision m_amvrPrecAffine[3];
  static const MvPrecision m_amvrPrecIbc[3];

public:
  int hor; ///< horizontal component of motion vector
  int ver; ///< vertical component of motion vector

  //   //
  //   ------------------------------------------------------------------------------------------------------------------
  //   // constructors
  //   //
  //   ------------------------------------------------------------------------------------------------------------------

  Mv() : hor(0), ver(0) {}
  Mv(int iHor, int iVer) : hor(iHor), ver(iVer) {}

  //   //
  //   ------------------------------------------------------------------------------------------------------------------
  //   // set
  //   //
  //   ------------------------------------------------------------------------------------------------------------------

  void set(int iHor, int iVer) {
    hor = iHor;
    ver = iVer;
  }
  void setZero() { hor = ver = 0; }

  //   //
  //   ------------------------------------------------------------------------------------------------------------------
  //   // get
  //   //
  //   ------------------------------------------------------------------------------------------------------------------

  int getHor() const { return hor; }
  int getVer() const { return ver; }

  const Mv &operator<<=(const int i) {
    hor <<= i;
    ver <<= i;
    return *this;
  }

  void changePrecision(const MvPrecision &src, const MvPrecision &dst) {
    const int shift = (int)dst - (int)src;
    if (shift >= 0) {
      *this <<= shift;
    } else {
      const int rightShift = -shift;
      const int nOffset = 1 << (rightShift - 1);
      hor = hor >= 0 ? (hor + nOffset - 1) >> rightShift
                     : (hor + nOffset) >> rightShift;
      ver = ver >= 0 ? (ver + nOffset - 1) >> rightShift
                     : (ver + nOffset) >> rightShift;
    }
  }

  // translational MV
  void changeTransPrecInternal2Amvr(const int amvr) {
    changePrecision(MV_PRECISION_INTERNAL, m_amvrPrecision[amvr]);
  }

  // affine MV
  void changeAffinePrecInternal2Amvr(const int amvr) {
    changePrecision(MV_PRECISION_INTERNAL, m_amvrPrecAffine[amvr]);
  }

  // IBC block vector
  void changeIbcPrecInternal2Amvr(const int amvr) {
    changePrecision(MV_PRECISION_INTERNAL, m_amvrPrecIbc[amvr]);
  }
}; // END CLASS DEFINITION MV
} // namespace Common

#endif // COMMON_MV
