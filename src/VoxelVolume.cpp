#include "VoxelVolume.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <string>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::resize(Vector3l const &s, T const &value) {
  if ((s.array() <= 0).any())
    throw runtime_error("s must be positive");
  if (!data.empty())
    throw runtime_error("resize must be called on empty VoxelVolume");
  this->s = s;
  spacing << 1, s(0), s(0) * s(1);
  data.resize(s.cast<size_t>().prod(), value);
}
//------------------------------------------------------------------------------
template <typename T>
size_t VoxelVolume<T>::vx_to_vxID(Vector3l const &vx) const {
  if ((vx.array() >= 0).all() && (vx.array() < s.array()).all())
    return spacing.dot(vx);

  Vector3l vxNew = vx;
  for (int dim = 0; dim < 3; ++dim) {
    if (vxNew(dim) < 0) {
      vxNew(dim) = 0;
    }
    if (vxNew(dim) >= s(dim)) {
      vxNew(dim) = s(dim) - 1;
    }
  }

  return spacing.dot(vxNew);
}
//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::import_raw_volume(Vector3l const &s,
                                       const char *filename) {
  ifstream myFile(filename);

  if (!myFile.good()) {
    cout << "\nWARNING: Can't import volume, no file found!\n";
    return;
  }

  resize(s);

  cout << "\nImporting Raw Volume from \"" << filename << "\" ...\n";

  if constexpr (!is_same_v<T, bool>)
    myFile.read((char *)&data[0], sizeof(T) * data.size());

  if (myFile.eof() || myFile.peek() != EOF) {
    cout << "\nWARNING: Couldnt't import volume, error while reading in file\n";
    return;
  }
}
//------------------------------------------------------------------------------
template <>
void VoxelVolume<MorphologyValue>::export_pgm_stacks(const char *) const;
//------------------------------------------------------------------------------
template <>
void VoxelVolume<Eigen::Vector3<uint8_t>>::export_pgm_stacks(
    const char *) const;
//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::export_pgm_stacks(const char *foldername) const {
  if (data.empty()) {
    cout << "\nWARNING: Can't export pgm stacks, empty Voxel Volume!\n";
    return;
  }

  cout << "\nExporting as pgm stacks to " << foldername << " ...\n";

  auto minMaxElement = minmax_element(data.begin(), data.end());
  T minElement = *minMaxElement.first;
  T maxElement = *minMaxElement.second;

  if (minElement == maxElement) {
    cout << "\nWARNING: Can't export image, minValue = maxValue!\n";
    return;
  }

  for (int k = 0; k < s(2); ++k) {
    vector<uint8_t> currImage(s(1) * s(0), 0);

#pragma omp parallel for
    for (int j = 0; j < s(1); ++j) {
      auto pxIt = currImage.begin() + spacing(1) * (s(1) - 1 - j);

      for (auto vxIt = data.begin() + spacing.dot(Vector3l(0, j, k));
           vxIt != data.begin() + spacing.dot(Vector3l(0, j + 1, k));
           ++vxIt, ++pxIt)
        *pxIt = ((*vxIt - minElement) * 255) / (maxElement - minElement);
    }

    char numberBuffer[64];
    snprintf(numberBuffer, 64, "%06i", k);
    ofstream myFile(string(foldername) + "stack" + numberBuffer + ".pgm");
    myFile << "P5\n" << s(0) << " " << s(1) << endl << 255 << endl;
    myFile.write((const char *)currImage.data(), currImage.size());
  }
}
//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::export_raw(const char *filename) const {

  ofstream myFile(filename);

  cout << "\nExporting Raw Volume to \"" << filename << "\" ...\n";

  if constexpr (!is_same_v<T, bool>)
    myFile.write((char *)&data[0], sizeof(T) * s.cast<size_t>().prod());
}
//------------------------------------------------------------------------------
template <>
void VoxelVolume<MorphologyValue>::export_stack_for_gp(long,
                                                       const char *) const;
//------------------------------------------------------------------------------
template <>
void VoxelVolume<Eigen::Vector3<uint8_t>>::export_stack_for_gp(
    long, const char *) const;
//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::export_stack_for_gp(long stackID,
                                         const char *filename) const {
  std::ofstream myFile(filename);
  for (long n = stackID * spacing(2); n < (stackID + 1) * spacing(2); ++n) {
    myFile << n % s(0) << " " << (n / s(0)) % s(1) << " " << float(data[n])
           << std::endl;
    if (n % spacing(1) == spacing(1) - 1)
      myFile << std::endl;
  }
}
//------------------------------------------------------------------------------
template struct VoxelVolume<uint8_t>;
template struct VoxelVolume<bool>;
template struct VoxelVolume<float>;
template struct VoxelVolume<int>;
template struct VoxelVolume<uint32_t>;
template struct VoxelVolume<MorphologyValue>;
template struct VoxelVolume<Eigen::Vector3<uint8_t>>;
//------------------------------------------------------------------------------
} // namespace fred
//------------------------------------------------------------------------------
