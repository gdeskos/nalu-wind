// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//



#ifndef PostProcessingInfo_h
#define PostProcessingInfo_h

#include <NaluParsedTypes.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class PostProcessingData;

class PostProcessingInfo
{

public:
  
  PostProcessingInfo();
  ~PostProcessingInfo();
  
  void load(const YAML::Node & y_node);
  
  std::vector<PostProcessingData *> ppDataVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
