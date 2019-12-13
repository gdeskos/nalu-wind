#include "mesh_motion/MotionWaves.h"

#include <NaluParsing.h>
#include "utils/ComputeVectorDivergence.h"

// stk_mesh/base/fem
#include <stk_mesh/base/FieldBLAS.hpp>

#include <cmath>

namespace sierra{
namespace nalu{

MotionWaves::MotionWaves(
	stk::mesh::MetaData& meta,
	const YAML::Node& node)
	: MotionBase()
{
  load(node);

  // declare divergence of mesh velocity for this motion
  ScalarFieldType *divV = &(meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_mesh_velocity"));
  stk::mesh::put_field_on_mesh(*divV, meta.universal_part(), nullptr);
  stk::mesh::field_fill(0.0, *divV);
}

void MotionWaves::load(const YAML::Node& node)
{
  // Get type of input prescribed wave
	
	get_if_present(node,"wave_model", waveModel_, waveModel_);

  if (waveModel_=="Linear_precibed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "wave", frequency_, frequency_);	 
	}
	else if (waveModel_ == "HOS"){

	}
	else {
    throw std::runtime_error("invalid wave_type model specified ");
	} 	     
  double eps = std::numeric_limits<double>::epsilon();
	
  // Time parameters
  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_-eps;
  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_+eps;

}


} // nalu
} // sierra

