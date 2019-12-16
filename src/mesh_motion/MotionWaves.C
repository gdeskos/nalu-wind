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
	
	get_if_present(node,"wave_motion_model", waveModel_, waveModel_);

  if (waveModel_=="Linear_prescribed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_); 
	}
	else if (waveModel_ == "HOS"){

	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	} 	     
  double eps = std::numeric_limits<double>::epsilon();
	
  // Time parameters
  get_if_present(node, "start_time", startTime_, startTime_);
  startTime_ = startTime_-eps;
  get_if_present(node, "end_time", endTime_, endTime_);
  endTime_ = endTime_+eps;
	get_if_present(node, "sea_level_z",sealevelz_,sealevelz_); 
  // get sea-level based on if it was defined

}
void MotionWaves::build_transformation(
  const double time,
  const double* xyz)
{
  if(time < (startTime_)) return;

  double motionTime = (time < endTime_)? time : endTime_;

	ThreeDVecType curr_disp={};
	curr_disp[0]=0.;
	curr_disp[1]=0.;

	if(waveModel_== "Linear_prescribed"){
	curr_disp[2]=sealevelz_+amplitude_*std::cos(2*M_PI/wavelength_*xyz[0]-2*M_PI/waveperiod_*time)*std::exp(-0.1*xyz[2]/amplitude_);
  }
	else if (waveModel_ =="HOS"){
	
	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	}
	translation_mat(curr_disp);
}

void MotionWaves::translation_mat(const ThreeDVecType& curr_disp)
{
  reset_mat(transMat_);

  // Build matrix for translating object
  transMat_[0][3] = curr_disp[0];
  transMat_[1][3] = curr_disp[1];
  transMat_[2][3] = curr_disp[2];
}

MotionBase::ThreeDVecType MotionWaves::compute_velocity(
  const double time,
  const TransMatType&  /* compTrans */,
  const double* mxyz,
  const double* /* cxyz */ )
{
  ThreeDVecType vel = {};

  if( (time < startTime_) || (time > endTime_) ) return vel;

  double WaveVelocity = amplitude_*2*M_PI/waveperiod_*std::sin(2.*M_PI/wavelength_*mxyz[0]-2*M_PI/wavelength_*time);

	vel[0] = 0.;
	vel[1] = 0.;
  vel[2] = WaveVelocity ;

  return vel;
}

void MotionWaves::post_compute_geometry(
  stk::mesh::BulkData& bulk,
  stk::mesh::PartVector& partVec,
  stk::mesh::PartVector& partVecBc,
  bool& computedMeshVelDiv)
{
  if(computedMeshVelDiv) return;

  // compute divergence of mesh velocity
  VectorFieldType* meshVelocity = bulk.mesh_meta_data().get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "mesh_velocity");

  ScalarFieldType* meshDivVelocity = bulk.mesh_meta_data().get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "div_mesh_velocity");

  compute_vector_divergence(bulk, partVec, partVecBc, meshVelocity, meshDivVelocity, true);
  computedMeshVelDiv = true;
}


} // nalu
} // sierra

