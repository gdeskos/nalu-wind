target_sources(${utest_ex_name} PRIVATE
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestContinuityAdvElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestContinuityMassElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestContinuityInflowElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestContinuityOpenElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestEnthalpyTGradBCElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestFaceBasic.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestFaceElemBasic.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestKernelUtils.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumActuatorSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumAdvDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumBuoyancyBoussinesqSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumBuoyancySrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumCoriolisSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumSSTTAMSForcingSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumSSTTAMSDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumMassElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumNSOElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumNSOKeElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestMomentumUpwAdvDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSSTSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSSTTAMSSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarAdvDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarMassElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarFluxBCElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarUpwAdvDiffElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestScalarOpenElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestSteadyThermal3dContactSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestTurbKineticEnergyKsgsSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestTurbKineticEnergyKsgsDesignOrderSrcElem.C
   ${CMAKE_CURRENT_SOURCE_DIR}/UnitTestWallDistElem.C
)
