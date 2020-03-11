#ifndef MOTIONWAVES_H
#define MOTIONWAVES_H

#include <iostream>
#include <cmath>
#include "MotionBase.h"

namespace sierra{
namespace nalu{

class MotionWaves : public MotionBase
{
public:
  MotionWaves(
    stk::mesh::MetaData&,
    const YAML::Node&);

  virtual ~MotionWaves()
  {
  }

virtual void build_transformation(const double, const double*);

  /** Function to compute motion-specific velocity
   *
   * @param[in] time           Current time
   * @param[in] compTrans      Transformation matrix
   *                           for points other than xyz
   * @param[in] mxyz           Model coordinates
   * @param[in] mxyz           Transformed coordinates
   */
  virtual ThreeDVecType compute_velocity(
    const double time,
    const TransMatType& compTrans,
    const double* mxyz,
    const double* cxyz );

  /** perform post compute geometry work for this motion
   *
   * @param[in] computedMeshVelDiv flag to denote if divergence of
   *                               mesh velocity already computed
   */
	void post_compute_geometry(
		stk::mesh::BulkData&,
    stk::mesh::PartVector&,
    stk::mesh::PartVector&,
    bool& computedMeshVelDiv );

private:
  MotionWaves() = delete;
  MotionWaves(const MotionWaves&) = delete;

  void load(const YAML::Node&);

  void translation_mat(const ThreeDVecType&);

  std::string waveModel_{"Sinusoidal_full_domain"};   
  double amplitude_{0.1};
  double waveperiod_{1.0};
  double wavefrequency_{2.*M_PI}; // Angular frequency omega=2*pi/tau (tau being the period)
  double wavelength_{1.0};
  double wavenumber_{2.*M_PI}; // Angular wavenumber k=2*pi/lambda (lambda being the wavenumber)
  double sealevelz_{0.0};
  // Deformation damping function
  double meshdampinglength_{1000};
  int meshdampingcoeff_{3};
  double dispersion_{1.0};
  double waterdepth_{50.};

};

} // nalu
} // sierra

#endif /* MOTIONWAVES_H */
