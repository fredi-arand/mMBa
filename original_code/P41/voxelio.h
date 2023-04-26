#ifndef VOXELIO_H
#define VOXELIO_H

#include "voxelvolume.h"
#include <Eigen/Dense>
#include <string>
#include <fstream>

namespace voxelIO {

  //***********************************************************************************************

  bool check_nr_of_arguments(int argc, int argCDesired){
    if(argc!=argCDesired)
    {
      std::cout << "\nWrong number of arguments\n";
      return false;}
    return true;
  }

  //***********************************************************************************************

  template <typename T>
  VoxelVolume<T> to_voxelVolume(std::string filename, Eigen::Vector3i const & s)
  {

    std::ifstream myfile(filename);

    if(!myfile.good())
    {std::cout << "\nWARNING: Can't import volume, no file found!\n"; return VoxelVolume<T>();}

    if((s.array()<1).any())
    {std::cout << "\nWARNING: Can't import volume, bad dimensions!\n"; return VoxelVolume<T>();}

    std::cout << "\nImporting Raw Volume from \"" << filename << "\" ...\n";

    VoxelVolume<T> voxelVolume;
    voxelVolume.s = s;
    std::vector<T> & voxelValues = voxelVolume.voxelValues;

    voxelValues.resize(s.cast<size_t>().prod(),T(0));

    // myfile.read only handles files up to 2 GB (the daily arschgeburt)
    size_t totalByteSize = sizeof(T) * s.cast<size_t>().prod();
    const size_t oneGB  = 30 << size_t(1);

    char* writePosition =reinterpret_cast<char *>(&voxelValues[0]);
    while( totalByteSize != 0 )
    {
      size_t bytesToRead = std::max(totalByteSize, oneGB);
      myfile.read(writePosition, bytesToRead);
      totalByteSize -= bytesToRead;
      writePosition += bytesToRead;
    }

    if(myfile.eof() || myfile.peek() != EOF)
    {std::cout << "\nWARNING: Couldnt't import volume, error while reading in file\n";
     return VoxelVolume<T>();}

    return voxelVolume;

  }

  //***********************************************************************************************

  template <typename T>
  void to_file(VoxelVolume<T> voxelVolume, std::string fileName)
  {

    std::ofstream myfile(fileName);

    std::cout << "\nExporting Raw Volume to \"" << fileName << "\" ...\n";

    auto const & voxelValues = voxelVolume.voxelValues;
    auto const & s = voxelVolume.s;

    myfile.write((char*)&voxelValues[0], sizeof(T)*s.prod());

  }

  //***********************************************************************************************

}

#endif // VOXELIO_H
