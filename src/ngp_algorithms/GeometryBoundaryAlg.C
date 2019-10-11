/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "ngp_algorithms/GeometryBoundaryAlg.h"
#include "BuildTemplates.h"
#include "master_element/MasterElement.h"
#include "master_element/MasterElementFactory.h"
#include "ngp_algorithms/ViewHelper.h"
#include "ngp_utils/NgpLoopUtils.h"
#include "ngp_utils/NgpFieldOps.h"
#include "Realm.h"
#include "ScratchViews.h"
#include "SolutionOptions.h"
#include "utils/StkHelpers.h"


namespace sierra {
namespace nalu {

template <typename AlgTraits>
GeometryBoundaryAlg<AlgTraits>::GeometryBoundaryAlg(
  Realm& realm, stk::mesh::Part* part)
  : Algorithm(realm, part),
    dataNeeded_(realm_.meta_data()),
    exposedAreaVec_(get_field_ordinal(
      realm.meta_data(), "exposed_area_vector", realm.meta_data().side_rank())),
    meSCS_(MasterElementRepo::get_surface_master_element<AlgTraits>())
{
  dataNeeded_.add_cvfem_surface_me(meSCS_);
  const auto coordID = get_field_ordinal(
    realm_.meta_data(), realm_.solutionOptions_->get_coordinates_name());
  dataNeeded_.add_coordinates_field(coordID, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataNeeded_.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
}

template<typename AlgTraits>
void GeometryBoundaryAlg<AlgTraits>::execute()
{
  using ElemSimdDataType = sierra::nalu::nalu_ngp::ElemSimdData<ngp::Mesh>;

  const auto& meshInfo = realm_.mesh_info();
  const auto& meta = meshInfo.meta();
  const auto ngpMesh = meshInfo.ngp_mesh();
  const auto& fieldMgr = meshInfo.ngp_field_manager();
  auto exposedAreaVec = fieldMgr.template get_field<double>(exposedAreaVec_);
  const auto areaVecOps = nalu_ngp::simd_elem_field_updater(ngpMesh, exposedAreaVec);

  const stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(partVec_);

  sierra::nalu::nalu_ngp::run_elem_algorithm(
    meshInfo, meta.side_rank(), dataNeeded_, sel,
    KOKKOS_LAMBDA(ElemSimdDataType & edata) {
      auto& scrViews = edata.simdScrView;
      const auto& meViews = scrViews.get_me_views(sierra::nalu::CURRENT_COORDINATES);
      const auto& v_area = meViews.scs_areav;

      for (int ip = 0; ip < AlgTraits::numFaceIp_; ++ip)
        for (int d=0; d < AlgTraits::nDim_; ++d)
          areaVecOps(edata, ip * AlgTraits::nDim_ + d) = v_area(ip, d);
    });

  exposedAreaVec.modify_on_device();
}

INSTANTIATE_KERNEL_FACE(GeometryBoundaryAlg)

}  // nalu
}  // sierra