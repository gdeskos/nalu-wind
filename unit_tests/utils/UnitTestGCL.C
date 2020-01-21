// Copyright 2017 National Technology & Engineering Solutions of Sandia, LLC
// (NTESS), National Renewable Energy Laboratory, University of Texas Austin,
// Northwest Research Associates. Under the terms of Contract DE-NA0003525
// with NTESS, the U.S. Government retains certain rights in this software.
//
// This software is released under the BSD 3-clause license. See LICENSE file
// for more details.
//

#include "utils/UnitTestGCL.h"

namespace {

  namespace hex8_golds_x_rot {
    namespace ngp_mesh_velocity {
      static constexpr double swept_vol[12] = {
        0.0, -0.0625, 0.0, -0.0625, 0.0, 0.0625,
        0.0, 0.0625, 0.0625, 0.0625, -0.0625, -0.0625
      };

      static constexpr double face_vel_mag[12] = {
        0.0, -0.125, 0.0, -0.125, 0.0, 0.125,
        0.0, 0.125, 0.125, 0.125, -0.125, -0.125
      };
    }
  }

  namespace hex8_golds_y_rot {
    namespace ngp_mesh_velocity {
      static constexpr double swept_vol[12] = {
        0.0625, 0.0, -0.0625, 0.0, -0.0625, 0.0,
        0.0625, 0.0, -0.0625, 0.0625, 0.0625, -0.0625
      };

      static constexpr double face_vel_mag[12] = {
        0.125, 0.0, -0.125, 0.0, -0.125, 0.0,
        0.125, 0.0, -0.125, 0.125, 0.125, -0.125
      };
    }
  }
}

TEST_F(GCLTest, rigid_rotation_elem)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: 1.5707963267948966                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}

TEST_F(GCLTest, rigid_rotation_edge)
{
  if (bulk_.parallel_size() > 1) return;

  realm_.realmUsesEdges_ = true;
  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: 1.5707963267948966                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}

TEST_F(GCLTest, rigid_scaling_elem)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: scaling                                                 \n"
    "        rate: [1.0, 1.0, 1.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_relative_error();
}

TEST_F(GCLTest, rigid_scaling_edge)
{
  if (bulk_.parallel_size() > 1) return;

  realm_.realmUsesEdges_ = true;
  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: scaling                                                 \n"
    "        rate: [1.0, 1.0, 1.0]                                         \n"
    "        centroid: [0.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_relative_error();
}


TEST_F(GCLTest, rigid_translation)
{
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "3x3x3|offset:0,65,0";
  const bool secondOrder = false;
  const double deltaT = 0.003; // approx 0.25 deg motion for given omega
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: translation                                             \n"
    "        velocity: [1.0, 0.0, 0.0]                                     \n";

  fill_mesh_and_init_fields(meshDims);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();
  compute_div_mesh_vel();
  compute_dvoldt();
  compute_absolute_error();
}

TEST_F(GCLTest, NGP_mesh_velocity_x_rot)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "1x1x1";
  const bool secondOrder = false;
  const double deltaT = 0.5; 
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: -3.141592653589793                                     \n"
    "        axis: [1.0, 0.0, 0.0]                                         \n"
    "        centroid: [0.5, 0.5, 0.5]                                     \n";
  
  fill_mesh_and_init_fields(meshDims, false);
  init_time_integrator(secondOrder, deltaT, 1);
  register_algorithms(mesh_motion);
  init_states();

  const double tol = 1.0e-15;
  namespace gold_values = ::hex8_golds_x_rot::ngp_mesh_velocity;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);
    int counter=0;
    for (const auto* b: bkts) {
      const double *sv = stk::mesh::field_data(*sweptVol_, *b, 0);
      const double *fvm = stk::mesh::field_data(*faceVelMag_, *b, 0);
      for (int i=0; i < 12; i++) {
        EXPECT_NEAR(gold_values::swept_vol[i], sv[i], tol);
        counter++;
        EXPECT_NEAR(gold_values::face_vel_mag[i], fvm[i], tol);
        counter++;
        }
    }
    EXPECT_EQ(counter,24);
  }
}


TEST_F(GCLTest, NGP_mesh_velocity_y_rot)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "1x1x1";
  const bool secondOrder = false;
  const double deltaT = 0.5; 
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: -3.141592653589793                                     \n"
    "        axis: [0.0,1.0,0.0]                                           \n"
    "        centroid: [0.5, 0.5, 0.5]                                     \n";
  
  fill_mesh_and_init_fields(meshDims, false);
  init_time_integrator(secondOrder, deltaT, 1);
  register_algorithms(mesh_motion);
  init_states();

  const double tol = 1.0e-15;
  namespace gold_values = ::hex8_golds_y_rot::ngp_mesh_velocity;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);
    int counter=0;
    for (const auto* b: bkts) {
      const double *sv = stk::mesh::field_data(*sweptVol_, *b, 0);
      const double *fvm = stk::mesh::field_data(*faceVelMag_, *b, 0);
      for (int i=0; i < 12; i++) {
        EXPECT_NEAR(gold_values::swept_vol[i], sv[i], tol);
        counter++;
        EXPECT_NEAR(gold_values::face_vel_mag[i], fvm[i], tol);
        counter++;
        }
    }
    EXPECT_EQ(counter,24);
  }
}

TEST_F(GCLTest, NGP_mesh_velocity_y_rot_scs_center)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  const std::string meshDims = "1x1x1";
  const bool secondOrder = false;
  const double deltaT = 0.25; 
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: -3.141592653589793                                     \n"
    "        axis: [0.0,1.0,0.0]                                           \n"
    "        centroid: [0.5, 0.5, 0.75]                                    \n";
  
  fill_mesh_and_init_fields(meshDims, false);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();

  const double tol = 1.0e-15;
  namespace gold_values = ::hex8_golds_y_rot::ngp_mesh_velocity;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::ELEM_RANK, sel);
    int counter=0;
    for (const auto* b: bkts) {
      const double *sv = stk::mesh::field_data(*sweptVol_, *b, 0);
      const double *fvm = stk::mesh::field_data(*faceVelMag_, *b, 0);
      for (int i=0; i < 12; i++) {
        // check only for scs through the center of which rotation axis passes
        // in addition to all scs perpendicular to rotation axis - total 6 scs
        if((i==1) || (i==3) || (i==4) || (i==5) || (i==6) || (i==7))
        {
          EXPECT_NEAR(0.0, sv[i], tol);
          counter++;
          EXPECT_NEAR(0.0, fvm[i], tol);
          counter++;
        }
      }
    }
    EXPECT_EQ(counter,12);
  }
}

TEST_F(GCLTest, NGP_mesh_velocity_linear_waves)
{
  // Analytical values computed by hand (well actually a python script)
  double VNm1_analytical=1.;
	double VN_analytical=1.0327665042944956;
  double VNp1_analytical=1.0655330085889911;
  double dVdt=0.3276650429449557; //based on BDF1
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;
		const double tol = 1.0e-14;
  const std::string meshDims = "1x1x1";
  const bool secondOrder = false;
  const double deltaT = 0.25; 
  const std::string mesh_motion =
    "mesh_motion:                                                          \n"
    "  - name: interior                                                    \n"
    "    frame: non_inertial                                               \n"
    "    mesh_parts: [ block_1 ]                                           \n"
    "    motion:                                                           \n"
    "      - type: rotation                                                \n"
    "        omega: -3.141592653589793                                     \n"
    "        axis: [0.0,1.0,0.0]                                           \n"
    "        centroid: [0.5, 0.5, 0.75]                                    \n";
 /* 
  fill_mesh_and_init_fields(meshDims, false);
  init_time_integrator(secondOrder, deltaT);
  register_algorithms(mesh_motion);
  init_states();

  MeshVelocityKernelHex8Mesh::fill_mesh_and_init_fields("generated:2x2x2",false,false);

  unit_test_utils::HelperObjects helperObjs(
    bulk_, stk::topology::HEX_8, 1, partVec_[0]);

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.currentTime_ = 0.2;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.;
  helperObjs.realm.timeIntegrator_ = &timeIntegrator;
  helperObjs.realm.realmUsesEdges_ = false;
  helperObjs.realm.solutionOptions_->meshMotion_ = true;
  
  sierra::nalu::GeometryAlgDriver geomAlgDriver(helperObjs.realm);
  geomAlgDriver.register_elem_algorithm<sierra::nalu::GeometryInteriorAlg>(
    sierra::nalu::INTERIOR, partVec_[0], "geometry");

  sierra::nalu::MeshVelocityAlg<sierra::nalu::AlgTraitsHex8> mvAlg(helperObjs.realm, partVec_[0]);
  
  //First set the mesh displacement corresponding to translation of the z-axis
  // create a yaml node describing rotation
  VectorFieldType *meshDispNp1 = &(meshDisp_->field_of_state(stk::mesh::StateNP1));
  VectorFieldType *meshDispN = &(meshDisp_->field_of_state(stk::mesh::StateN));
  VectorFieldType *meshDispNm1 = &(meshDisp_->field_of_state(stk::mesh::StateNM1));

  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
    for (const auto node: *b) {
      double* mcoord = stk::mesh::field_data(*coordinates_, node);
      double* ccoord = stk::mesh::field_data(*cCoords_,node);
      for(int j=0; j<3; j++)
	ccoord[j] = mcoord[j];
      }
    }
  }

  
  double full_dnvNm1 = 0;
  double full_dnvN = 0;
  double full_dnvNp1 = 0;
  geomAlgDriver.execute();
  bool foundNode = false; 
  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
    for (const auto node: *b) {
      double* mcoord = stk::mesh::field_data(*coordinates_, node);
      double* dnv = stk::mesh::field_data(*dnvField_,node);
      if (bulk_.num_elements(node) == 8) {
        full_dnvNm1 = *dnv;
        foundNode = true;
        break;
      }
    }
  }
  }
  EXPECT_NEAR(VNm1_analytical, full_dnvNm1, tol);
  
  //Compute intermediate volume at N
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
      for (const auto node: *b) {
          double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
          double* dispN = stk::mesh::field_data(*meshDispN, node);
          double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
          double* ccoord = stk::mesh::field_data(*cCoords_, node);
          double* mcoord = stk::mesh::field_data(*coordinates_, node);
          // Coordinates at n-1
          dispNm1[0] = 0.0;
          dispNm1[1] = 0.0;
          // Coordinates at n
          dispN[0] = 0.0;
          dispN[1] = 0.0;
          // Coordinates at n+1
          dispNp1[0] = 0.0;
          dispNp1[1] = 0.0;
         
					if(mcoord[2]>0+tol){
          dispNm1[2]= 0.0;//0.2*std::cos(2.*M_PI*mcoord[0]/3.-2.*M_PI/0.3*0.);
          dispNp1[2]= 0.1*std::cos(2.*M_PI*mcoord[0]/8.);
          }else{
					dispNm1[2]=0.;
					dispNp1[2]=0.;
					}
					ccoord[0]=dispNp1[0]+mcoord[0];
					ccoord[1]=dispNp1[1]+mcoord[1];
					ccoord[2]=dispNp1[2]+mcoord[2];
			}
  }
  geomAlgDriver.execute();
  foundNode = false; 
  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
    for (const auto node: *b) {
      double* mcoord = stk::mesh::field_data(*coordinates_, node);
      double* dnv = stk::mesh::field_data(*dnvField_,node);
      if (bulk_.num_elements(node) == 8) {
        full_dnvN = *dnv;
        foundNode = true;
        break;
      }
    }
  }
  }
  EXPECT_NEAR(VN_analytical, full_dnvN, tol);
  //Compute final volume at Np1
  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
      for (const auto node: *b) {
          double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
          double* dispN = stk::mesh::field_data(*meshDispN, node);
          double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
          double* ccoord = stk::mesh::field_data(*cCoords_, node);
          double* mcoord = stk::mesh::field_data(*coordinates_, node);
          // Coordinates at n-1
          dispNm1[0] = 0.0;
          dispNm1[1] = 0.0;
					dispNm1[2]=0.;
          // Coordinates at n
          dispN[0] = 0.0;
          dispN[1] = 0.0;
          dispN[2] = 0.0;
          // Coordinates at n+1
          dispNp1[0] = 0.0;
          dispNp1[1] = 0.0;
         
					if(mcoord[2]>0+tol){
          dispNp1[2]= 0.2*std::cos(2.*M_PI*mcoord[0]/8.);
          }else{
					dispNp1[2]=0.;
					}
					ccoord[0]=dispNp1[0]+mcoord[0];
					ccoord[1]=dispNp1[1]+mcoord[1];
					ccoord[2]=dispNp1[2]+mcoord[2];
			}
  }
	}
  geomAlgDriver.execute();
  foundNode = false; 
  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
    for (const auto node: *b) {
      double* mcoord = stk::mesh::field_data(*coordinates_, node);
      double* dnv = stk::mesh::field_data(*dnvField_,node);
      if (bulk_.num_elements(node) == 8) {
        full_dnvNp1 = *dnv;
        foundNode = true;
        break;
      }
    }
  }
  }
  EXPECT_NEAR(VNp1_analytical, full_dnvNp1, tol);
  //Do the actual mesh motion calculation
  {
  stk::mesh::Selector sel = meta_.universal_part();
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
  for (const auto* b: bkts) {
      for (const auto node: *b) {
          double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
          double* dispN = stk::mesh::field_data(*meshDispN, node);
          double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
          double* ccoord = stk::mesh::field_data(*cCoords_, node);
          double* mcoord = stk::mesh::field_data(*coordinates_, node);
          // Coordinates at n-1
          dispNm1[0] = 0.0;
          dispNm1[1] = 0.0;
          // Coordinates at n
          dispN[0] = 0.0;
          dispN[1] = 0.0;
          // Coordinates at n+1
          dispNp1[0] = 0.0;
          dispNp1[1] = 0.0;
         
					if(mcoord[2]>0+tol){
          dispNm1[2]= 0.0;
          dispN[2]= 0.1*std::cos(2.*M_PI*mcoord[0]/8.);
          dispNp1[2]= 0.2*std::cos(2.*M_PI*mcoord[0]/8.);
          }else{
					dispNm1[2]=0.;
					dispN[2]=0.;
					dispNp1[2]=0.;
					}
					ccoord[0]=dispNp1[0]+mcoord[0];
					ccoord[1]=dispNp1[1]+mcoord[1];
					ccoord[2]=dispNp1[2]+mcoord[2];
			}
  }
  }
  geomAlgDriver.execute();

  mvAlg.execute();
  stk::mesh::PartVector bndyPartVec;
  sierra::nalu::compute_scalar_divergence(bulk_, partVec_, bndyPartVec, faceVelMag_, divMeshVelField_);

  double full_dnv_mdv = 1e15;
  foundNode = false;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    for (const auto* b: bkts) {
      for (const auto node: *b) {
        const double* meshVelDiv = stk::mesh::field_data(*divMeshVelField_,node);
        if (bulk_.num_elements(node) == 8) {
          full_dnv_mdv = *meshVelDiv;
          foundNode = true;
          break;
        }
        if (foundNode)
          break;
      }
      if (foundNode)
        break;
    }
  }
  EXPECT_NEAR(dVdt, full_dnv_mdv, tol);
}

TEST_F(MeshVelocityKernelHex8Mesh, NGP_mesh_vel_div_scaling)
{
  // Only execute for 1 processor runs
  if (bulk_.parallel_size() > 1) return;

  MeshVelocityKernelHex8Mesh::fill_mesh_and_init_fields("generated:2x2x2",true,false);

  unit_test_utils::HelperObjects helperObjs(
    bulk_, stk::topology::HEX_8, 1, partVec_[0]);

  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.currentTime_ = 0.2;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.0;
  helperObjs.realm.timeIntegrator_ = &timeIntegrator;
  // Force computation of edge area vector
  helperObjs.realm.realmUsesEdges_ = false;
  helperObjs.realm.solutionOptions_->meshMotion_ = true;
  
  sierra::nalu::GeometryAlgDriver geomAlgDriver(helperObjs.realm);
  geomAlgDriver.register_elem_algorithm<sierra::nalu::GeometryInteriorAlg>(
    sierra::nalu::INTERIOR, partVec_[0], "geometry");
  geomAlgDriver.register_elem_algorithm<
    sierra::nalu::MeshVelocityAlg>(sierra::nalu::INTERIOR, partVec_[0], "mesh_vel");

  bool foundNode = false;
  stk::mesh::Entity full_dnv_node;
  {
    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    for (const auto* b: bkts) {
      for (const auto node: *b) {
        if (bulk_.num_elements(node) == 8) {
          full_dnv_node = node;
          foundNode = true;
          break;
        }
      }
      if (foundNode)
        break;
    }
  }
  
  //First set the mesh displacement corresponding to rotation about x-axis
  // create a yaml node describing rotation
  const std::string scaleInfo =
      "rate: [1.0,1.0,1.0]     \n"
      "centroid: [0.0,0.0,0.0] \n"
      ;
  YAML::Node scaleNode = YAML::Load(scaleInfo);
  sierra::nalu::MotionScaling scaleClass(meta_, scaleNode);
  VectorFieldType *meshDispNp1 = &(meshDisp_->field_of_state(stk::mesh::StateNP1));
  VectorFieldType *meshDispN = &(meshDisp_->field_of_state(stk::mesh::StateN));
  VectorFieldType *meshDispNm1 = &(meshDisp_->field_of_state(stk::mesh::StateNM1));

  //Get dual nodal volume at time N
  {
      stk::mesh::Selector sel = meta_.universal_part();
      const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
      for (const auto* b: bkts) {
          for (const auto node: *b) {
              double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
              double* ccoord = stk::mesh::field_data(*cCoords_, node);
              double* mcoord = stk::mesh::field_data(*coordinates_, node);
              
              scaleClass.build_transformation(0.1,mcoord);
              std::vector<double> rot_xyz = transform(scaleClass.get_trans_mat(), mcoord);
              dispNp1[0] = rot_xyz[0] - mcoord[0];
              dispNp1[1] = rot_xyz[1] - mcoord[1];
              dispNp1[2] = rot_xyz[2] - mcoord[2];
              
              ccoord[0] = rot_xyz[0];
              ccoord[1] = rot_xyz[1];
              ccoord[2] = rot_xyz[2];
          }
      }
  }
  geomAlgDriver.execute();
  double dnv_n = *(stk::mesh::field_data(*dnvField_,full_dnv_node));

  //Now do the full thing at time N+1
  {
      stk::mesh::Selector sel = meta_.universal_part();
      const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
      for (const auto* b: bkts) {
          for (const auto node: *b) {
              double* dispNp1 = stk::mesh::field_data(*meshDispNp1, node);
              double* dispN = stk::mesh::field_data(*meshDispN, node);
              double* dispNm1 = stk::mesh::field_data(*meshDispNm1, node);
              double* ccoord = stk::mesh::field_data(*cCoords_, node);
              double* mcoord = stk::mesh::field_data(*coordinates_, node);
              dispNm1[0] = 0.0;
              dispNm1[1] = 0.0;
              dispNm1[2] = 0.0;

              scaleClass.build_transformation(0.1,mcoord);
              std::vector<double> rot_xyz = transform(scaleClass.get_trans_mat(), mcoord);
              dispN[0] = rot_xyz[0] - mcoord[0];
              dispN[1] = rot_xyz[1] - mcoord[1];
              dispN[2] = rot_xyz[2] - mcoord[2];
              
              scaleClass.build_transformation(0.2,mcoord);
              rot_xyz = transform(scaleClass.get_trans_mat(), mcoord);
              dispNp1[0] = rot_xyz[0] - mcoord[0];
              dispNp1[1] = rot_xyz[1] - mcoord[1];
              dispNp1[2] = rot_xyz[2] - mcoord[2];
              
              ccoord[0] = rot_xyz[0];
              ccoord[1] = rot_xyz[1];
              ccoord[2] = rot_xyz[2];
          }
      }
  }
  geomAlgDriver.execute();
  double dnv_np1 = *(stk::mesh::field_data(*dnvField_,full_dnv_node));
  double dvdt = (dnv_np1 - dnv_n)/0.1;
  
  stk::mesh::PartVector bndyPartVec;
  sierra::nalu::compute_scalar_divergence(bulk_, partVec_, bndyPartVec, faceVelMag_, divMeshVelField_);
  double full_dnv_mdv = *(stk::mesh::field_data(*divMeshVelField_,full_dnv_node));
  
  const double tol = 1.0e-14;
  EXPECT_NEAR(dvdt, full_dnv_mdv, tol);
*/
}
