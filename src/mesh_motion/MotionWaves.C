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

  if (waveModel_=="Sinusoidal_full_domain"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_); 
	// Compute the wave number
  }
  else if (waveModel_=="Linear_prescribed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_); 
	}
  else if (waveModel_=="StokesSecondOrder_prescribed"){
  get_if_present(node, "amplitude", amplitude_, amplitude_);
  get_if_present(node, "waveperiod", waveperiod_, waveperiod_);	
	get_if_present(node, "wavelength",wavelength_,wavelength_);  
  get_if_present(node, "waterdepth_", waterdepth_, waterdepth_);	
  }
  else if (waveModel_=="StokesThirdOrder_prescribed"){
  
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
  // Compute wave-related parameters
  wavenumber_=2.*M_PI/wavelength_;
  wavefrequency_=2.*M_PI/waveperiod_;
  //Define dispersion
	dispersion_= std::sqrt(std::tanh(wavenumber_*waterdepth_));

}
void MotionWaves::build_transformation(
  const double time,
  const double* xyz)
{
  if(time < (startTime_)) return;

  double motionTime = (time < endTime_)? time : endTime_;

  double phase=wavenumber_*xyz[0]-wavefrequency_*motionTime;
	ThreeDVecType curr_disp={};
	curr_disp[0]=0.;
	curr_disp[1]=0.;
	
	
	if(waveModel_== "Sinusoidal_full_domain"){
	curr_disp[2]=sealevelz_+amplitude_*std::cos(phase);
	}
	else if(waveModel_== "Linear_prescribed"){
	curr_disp[2]=sealevelz_+amplitude_*std::cos(phase)*std::exp(-0.1*xyz[2]/amplitude_);
	}
  else if (waveModel_ == "StokesSecondOrder_prescribed"){
  curr_disp[2]=sealevelz_+amplitude_*(std::cos(phase)+wavenumber_*amplitude_*(3-dispersion_*dispersion_)/(4.*dispersion_*dispersion_*dispersion_)*std::cos(2.*phase))*std::exp(-0.1*xyz[2]/amplitude_);
  }
	else if (waveModel_ == "StokesThirdOrder_prescribed"){
	//curr_disp[2]=sealevelz_*((1.0-1.0/16.0*(wavenumber_*amplitude_)*(wavenumber_*amplitude_))*std::cos(phase)+1.0/2.0*(wavenumber_*amplitude_)*std::cos(2.*phase)+3./8.*(wavenumber_*amplitude_)*(wavenumber_*amplitude_)*std::cos(3*phase));
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
  
  double motionTime = (time < endTime_)? time : endTime_;

  double VerticalWaveVelocity=0;
  double HorizontalWaveVelocity=0;
  double phase=wavenumber_*mxyz[0]-wavefrequency_*motionTime;
	
  if(waveModel_== "Sinusoidal_full_domain"){
  VerticalWaveVelocity = amplitude_*wavefrequency_*std::sin(phase);  
  HorizontalWaveVelocity = 0.;
  }
  else if(waveModel_== "Linear_prescribed"){
  VerticalWaveVelocity = amplitude_*wavefrequency_*std::sin(phase);  
  HorizontalWaveVelocity = amplitude_*wavefrequency_*std::cos(phase);
  }
  else if (waveModel_ == "StokesSecondOrder_prescribed"){
  VerticalWaveVelocity=amplitude_*wavefrequency_/std::tanh(wavenumber_*waterdepth_)*std::sin(phase)
													+3./4.*wavefrequency_*wavenumber_*std::sinh(2*wavenumber_*waterdepth_)/std::pow(std::sinh(wavenumber_*waterdepth_),4)*std::sin(2.*phase);	
  HorizontalWaveVelocity= amplitude_*wavefrequency_/std::tanh(wavenumber_*waterdepth_)*std::cos(phase)
													+3./4.*wavefrequency_*wavenumber_*std::cosh(2*wavenumber_*waterdepth_)/std::pow(std::sinh(wavenumber_*waterdepth_),4)*std::cos(2.*phase);	
	}
	else if (waveModel_ == "StokesThirdOrder_prescribed"){
  VerticalWaveVelocity=0.;
  HorizontalWaveVelocity=0.;	
  }
	else if (waveModel_ =="HOS"){
  VerticalWaveVelocity=0.;
  HorizontalWaveVelocity=0.;		
	}
	else {
    throw std::runtime_error("invalid wave_motion model specified ");
	}

  double eps = std::numeric_limits<double>::epsilon();
  
 
  if(mxyz[2] < sealevelz_ + eps){
  vel[0] = HorizontalWaveVelocity;
	vel[2] = VerticalWaveVelocity ;
  }else{
  vel[0] = 0.;
  vel[2] = 0.;
  }
	vel[1] = 0.;
  
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

