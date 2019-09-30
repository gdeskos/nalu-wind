/*------------------------------------------------------------------------*/
/*  Copyright 2019 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <array>

#include "gtest/gtest.h"
#include "UnitTestUtils.h"
#include "ngp_utils/NgpTypes.h"

namespace sierra{
namespace nalu{
namespace nalu_ngp {

struct BaseClass
{
  KOKKOS_FUNCTION virtual ~BaseClass() {} 
  KOKKOS_FUNCTION virtual void Warning(int,
    std::array<DoubleType,3> &,
    std::array<DoubleType,3> &){}
  KOKKOS_FUNCTION virtual void NoWarning(
    std::array<DoubleType,3> &,
    std::array<DoubleType,3> &){}
};

class NgpCompileTest : public ::testing::Test {};

void show_cuda_compiler_warning()
{
  const std::string debuggingName("show_cuda_compiler_warning");
  BaseClass* base = kokkos_malloc_on_device<BaseClass>(debuggingName);
  Kokkos::parallel_for(debuggingName, 1, KOKKOS_LAMBDA (const int) {
    new (base) BaseClass();
  });

  using ShmemType = typename NGPMeshTraits<ngp::Mesh>::ShmemType;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const size_t&) {
    std::array<DoubleType,3> D0;
#define TRIGGER_WARNING
#ifdef TRIGGER_WARNING
    base->Warning(0, D0, D0);
#endif
    base->NoWarning(D0, D0);
  });
}

TEST_F(NgpCompileTest, Show_Compiler_Warning)
{
  show_cuda_compiler_warning();
}

}}}
