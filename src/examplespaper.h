#pragma once

#include "VoxelVolume.h"
#include "distancefield.h"
#include "poremorphology.h"
#include <Eigen/Dense>
#include <iostream>
#include <utility>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using namespace std;
using namespace Eigen;

// VoxelVolume<int> two_circles()
//{
//   VoxelVolume<int> image;
//   image.s = Vector3l(9,13,1);
//   image.set_spacing_and_voxelValues_from_s();
//   image.voxelValues.clear();
//   image.voxelValues.resize(image.s.cast<size_t>().prod(),0);

//  Vector2f m0(4,4), m1(4,7); float r0 = 3.5*3.5, r1 = 2.5*2.5;

//  ofstream myCircles("output/myCircles0");
//  myCircles << m0.transpose() << " " << sqrt(r0) << endl
//            << m1.transpose() << " " << sqrt(r1) << endl;

//#pragma omp parallel for
//  for(size_t pxID=0; pxID<image.s.cast<size_t>().prod(); ++pxID)
//  {
//    int sGrid = 16;
//    size_t hits=0;
//    Vector3f px0 =
//    (image.vxID_to_vx(pxID)).cast<float>()+Vector3f(0.5/(2*sGrid),
//                                                                   0.5/(2*sGrid),
//                                                                   0);
//    for(int pxX = -sGrid; pxX<sGrid; ++pxX)
//      for(int pxY = -sGrid; pxY<sGrid; ++pxY)
//      {
//        Vector2f px(px0(0)+float(pxX)/(2*sGrid), px0(1)+float(pxY)/(2*sGrid));
//        if((px-m0).squaredNorm()>r0 &&
//           (px-m1).squaredNorm()>r1)
//          ++hits;
//      }

//    image[pxID] = (hits*256)/(4*sGrid*sGrid);
//  }

//  image.hasValues = true;
//  return image;
//}

////------------------------------------------------------------------------------

// VoxelVolume<int> two_circles_one_box()
//{
//   VoxelVolume<int> image;
//   image.s = Vector3l(16,24,1);
//   image.set_spacing_and_voxelValues_from_s();
//   image.voxelValues.clear();
//   image.voxelValues.resize(image.s.cast<size_t>().prod(),0);

//  Vector2f m0(7,7), m1(9,17); float r0 = 6.5*6.5, r1 = 5.5*5.5;

//  ofstream myCircles("output/myCircles0");
//  myCircles << m0.transpose() << " " << sqrt(r0) << endl
//            << m1.transpose() << " " << sqrt(r1) << endl;

//#pragma omp parallel for
//  for(size_t pxID=0; pxID<image.s.cast<size_t>().prod(); ++pxID)
//  {
//    int sGrid = 16;
//    size_t hits=0;
//    Vector3f px0 =
//    (image.vxID_to_vx(pxID)).cast<float>()+Vector3f(0.5/(2*sGrid),
//                                                                   0.5/(2*sGrid),
//                                                                   0);
//    for(int pxX = -sGrid; pxX<sGrid; ++pxX)
//      for(int pxY = -sGrid; pxY<sGrid; ++pxY)
//      {
//        Vector2f px(px0(0)+float(pxX)/(2*sGrid), px0(1)+float(pxY)/(2*sGrid));
//        if((px-m0).squaredNorm()>r0 &&
//           (px-m1).squaredNorm()>r1)
//          ++hits;
//      }

//    image[pxID] = (hits*256)/(4*sGrid*sGrid);
//  }

//  image.hasValues = true;
//  return image;
//}

////------------------------------------------------------------------------------

// void create_example()
//{
//   DistanceField distanceField;
//   VoxelVolume<int> exampleVolume = two_circles();
//   exampleVolume.export_stack_for_gp(0,"output/greyValues");
//   distanceField.create_distance_field<int>(exampleVolume,127.5);
//   distanceField.export_stack_for_gp(0,"output/distanceField");

//  ofstream centralCicles("output/myCircles");
//  for(size_t j=0; j<distanceField.s(1); ++j)
//  {
//    Vector3l vx(4,j,0);
//    centralCicles << vx(0) << " " << vx(1) << " " <<
//    distanceField[distanceField.vx_to_vxID(vx)] << endl;
//  }

//}

////------------------------------------------------------------------------------

// void create_example2()
//{
//   DistanceField distanceField;
//   VoxelVolume<int> exampleVolume = two_circles_one_box();
//   exampleVolume.export_stack_for_gp(0,"output/greyValues");
//   distanceField.create_distance_field<int>(exampleVolume,127.5);
//   distanceField.export_pgm_stacks("output/distance_field/");
//   distanceField.export_stack_for_gp(0,"output/distanceField");

//  ofstream centralCicles("output/myCircles");
//  for(size_t j=0; j<distanceField.s(1); ++j)
//  {
//    Vector3l vx(5,j,0);
//    centralCicles << vx(0) << " " << vx(1) << " " <<
//    distanceField[distanceField.vx_to_vxID(vx)] << endl;
//  }

//  PoreMorphology poreMorphology(distanceField);
//  poreMorphology.create_pore_morphology();
//  poreMorphology.export_ppm_stacks("output/morphology/before_throat_reduction/");
//  poreMorphology.reduce_throat_volume();
//  poreMorphology.export_statistics("output/statistics/");
//  poreMorphology.export_ppm_stacks("output/morphology/");

//}
//------------------------------------------------------------------------------
} // namespace fred
