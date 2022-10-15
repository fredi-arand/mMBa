#ifndef VOXELVOLUME_H
#define VOXELVOLUME_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <math.h>
#include <map>
#include <type_traits>


// ************************************************************************************************



template <typename T>
struct VoxelVolume
{

    VoxelVolume() : s(Eigen::Vector3i(0,0,0)) {}
    template<typename T2>
    VoxelVolume(VoxelVolume<T2> voxelVolume){
      if(std::is_same<T,T2>::value){
        *this=voxelVolume; return;
      }
      s = voxelVolume.s; this->voxelValues.reserve(s.prod());
      for( auto const & element : voxelVolume.voxelValues ){
        voxelValues.push_back(T(element));
      }
    }

    Eigen::Vector3i to_vx(size_t vxID) const
    { return Eigen::Vector3i( vxID%s(0), (vxID/s(0))%s(1), vxID/(s(0)*s(1)) ); }

    size_t to_vxID(Eigen::Vector3i const &vx) const;


    T operator[](size_t vxID) const {return voxelValues[vxID];}

    T & operator[](size_t vxID) {return voxelValues[vxID];}

    T operator()(Eigen::Vector3i x) const {return this->operator [](to_vxID(x));}
    T operator()(long x0, long x1, long x2) const
    { return this->operator ()(Eigen::Vector3i(x0,x1,x2)); }

    T & operator()(Eigen::Vector3i x)  {return this->operator [](to_vxID(x));}
    T & operator()(long x0, long x1, long x2)
    { return this->operator ()(Eigen::Vector3i(x0,x1,x2)); }

    // MEMBERS
    std::vector<T> voxelValues;
    Eigen::Vector3i s;
};


//*************************************************************************************************

template <typename T>
size_t VoxelVolume<T>::to_vxID (Eigen::Vector3i const &vx) const
{

  Eigen::Vector3i spacing(1,s(0),s(0)*s(1));

  if((vx.array()>=0).all() &&
      (vx.array()<s.array()).all() )
    return spacing.dot(vx);

  Eigen::Vector3i vxNew = vx;
  for(int dim=0; dim<3; ++dim)
  {
    if(vxNew(dim)<0) {vxNew(dim)=0;}
    if(vxNew(dim)>=s(dim)) {vxNew(dim)=s(dim)-1;}
  }

  return spacing.dot(vxNew);
}

//*************************************************************************************************



#endif // VOXELVOLUME_H
