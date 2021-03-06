// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//



#ifndef PstabErrorIndicatorEdgeAlgorithm_h
#define PstabErrorIndicatorEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class PstabErrorIndicatorEdgeAlgorithm : public Algorithm
{
public:

  PstabErrorIndicatorEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    VectorFieldType *Gpdx,
    const bool simpleGradApproach = false);
  ~PstabErrorIndicatorEdgeAlgorithm();

  void execute();

  // extract fields; nodal
  ScalarFieldType *pressure_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  VectorFieldType *edgeAreaVec_;
  GenericFieldType *pstabEI_;

  const double simpleGradApproachScale_;
};

} // namespace nalu
} // namespace Sierra

#endif
