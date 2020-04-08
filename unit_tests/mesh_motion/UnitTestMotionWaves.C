#include <gtest/gtest.h>
#include <limits>

#include "mesh_motion/MotionRotation.h"
#include "mesh_motion/MotionScaling.h"
#include "mesh_motion/MotionTranslation.h"

#include "UnitTestRealm.h"

namespace {

  const double testTol = 1e-14;

  std::vector<double> transform(
    const sierra::nalu::MotionBase::TransMatType& transMat,
    const double* xyz )
  {
    std::vector<double> transCoord(3,0.0);

    // perform matrix multiplication between transformation matrix
    // and original coordinates to obtain transformed coordinates
    for (int d = 0; d < sierra::nalu::MotionBase::threeDVecSize; d++) {
      transCoord[d] = transMat[d][0]*xyz[0]
                     +transMat[d][1]*xyz[1]
                     +transMat[d][2]*xyz[2]
                     +transMat[d][3];
    }

    return transCoord;
  }

}

TEST(meshMotion, airy_wave)
{
  // create a yaml node describing rotation
  const std::string Airy_Wave_Info =
    "omega: 3.0              \n"
    "centroid: [0.3,0.5,0.0] \n"
    ;
}
