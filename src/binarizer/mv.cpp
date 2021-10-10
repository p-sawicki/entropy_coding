#include "mv.hpp"

const EntropyCoding::MvPrecision EntropyCoding::Mv::m_amvrPrecision[4] = {
    MV_PRECISION_QUARTER, MV_PRECISION_INT, MV_PRECISION_4PEL,
    MV_PRECISION_HALF}; // for cu.imv=0, 1, 2 and 3
const EntropyCoding::MvPrecision EntropyCoding::Mv::m_amvrPrecAffine[3] = {
    MV_PRECISION_QUARTER, MV_PRECISION_SIXTEENTH,
    MV_PRECISION_INT}; // for cu.imv=0, 1 and 2
const EntropyCoding::MvPrecision EntropyCoding::Mv::m_amvrPrecIbc[3] = {
    MV_PRECISION_INT, MV_PRECISION_INT,
    MV_PRECISION_4PEL}; // for cu.imv=0, 1 and 2
