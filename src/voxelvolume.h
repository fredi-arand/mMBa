#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <vector>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using namespace std;
using namespace Eigen;
using Vector3l = Vector3<long>;
//------------------------------------------------------------------------------
template <typename T> struct VoxelVolume {

  VoxelVolume()
      : s(Vector3l(0, 0, 0)), spacing(Vector3l(1, 1, 1)), hasValues(false) {}

  // standard functions to translate between ID and Coord
  Vector3l vxID_to_vx(size_t vxID) const {
    return Vector3l(vxID % s(0), (vxID / s(0)) % s(1), vxID / (s(0) * s(1)));
  }

  // use this to always get vxID even for illegal values of vx
  size_t vx_to_vxID(Vector3l const &vx) const;

  void set_spacing_and_voxelValues_from_s() {
    // quickly set spacing and voxelValues for given s
    spacing << 1, s(0), s(0) * s(1);
    voxelValues.resize(s.cast<size_t>().prod());
  }

  // switch x and y
  void switch_xy();

  // import volume from raw file
  void import_raw_volume(Vector3l const &s, string filename);

  // binarize volume
  T binarize_volume(VoxelVolume<bool> &isFict, VoxelVolume<bool> &isPhys) const;
  void binarize_volume(VoxelVolume<bool> &isFict, VoxelVolume<bool> &isPhys,
                       T isoValue) const;

  // create from multiplication of two other volumes
  void create_from_two_volumes(VoxelVolume<T> const &volA,
                               VoxelVolume<T> const &volB);

  // export volume as pgm stacks
  void export_pgm_stacks(string foldername) const;
  void export_raw(string filename) const;

  // export for gnuplot
  void export_stack_for_gp(long stackID, string filename) const {
    ofstream myFile(filename);
    for (long n = stackID * spacing(2); n < (stackID + 1) * spacing(2); ++n) {
      myFile << n % s(0) << " " << (n / s(0)) % s(1) << " "
             << float(voxelValues[n]) << endl;
      if (n % spacing(1) == spacing(1) - 1)
        myFile << endl;
    }
  }

  T operator[](size_t vxID) const { return voxelValues[vxID]; }

  T operator()(Vector3l x) const { return this->operator[](vx_to_vxID(x)); }
  T operator()(long x0, long x1, long x2) const {
    return this->operator()(Vector3l(x0, x1, x2));
  }

  size_t size() { return voxelValues.size(); }

  // voxelValues: value of voxel (i,j,k) is stored at voxelValues(spacing(2)*k +
  // spacing(1)*j + i) s:           size of voxel grid spacing:     spacing in
  // 0/1/2-direction
  vector<T> voxelValues;
  Vector3l s, spacing;
  bool hasValues;
};

//------------------------------------------------------------------------------

template <typename T> void VoxelVolume<T>::switch_xy() {
  VoxelVolume<T> switchedVolume(*this);
  switchedVolume.voxelValues.clear();
  switchedVolume.voxelValues.resize(s.cast<size_t>().prod());

  Vector3l xi;
  size_t i = 0;
  for (xi(2) = 0; xi(2) < s(2); ++xi(2))
    for (xi(0) = 0; xi(0) < s(0); ++xi(0))
      for (xi(1) = 0; xi(1) < s(1); ++xi(1), ++i)
        switchedVolume.voxelValues[i] = voxelValues[xi.dot(spacing)];

  voxelValues = switchedVolume.voxelValues;
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
T VoxelVolume<T>::binarize_volume(VoxelVolume<bool> &isFict,
                                  VoxelVolume<bool> &isPhys) const {
  auto minMaxValue = minmax_element(voxelValues.begin(), voxelValues.end());
  T isoValue = (*(minMaxValue.first) + *(minMaxValue.second)) / 2;
  binarize_volume(isFict, isPhys, isoValue);
  return isoValue;
}

//------------------------------------------------------------------------------

template <typename T>
void VoxelVolume<T>::binarize_volume(VoxelVolume<bool> &isFict,
                                     VoxelVolume<bool> &isPhys,
                                     T isoValue) const {
  if (!hasValues) {
    cout << "\nWARNING: Cannot binarize volume, empty Voxel Volume!\n";
  }

  cout << "\nCreating Binary Volumes ...\n";

  isFict.s = s.array() - 1;
  isFict.spacing = Vector3l(1, isFict.s(0), isFict.s(0) * isFict.s(1));
  isFict.voxelValues.resize(isFict.s.cast<size_t>().prod());

  isPhys.s = isFict.s;
  isPhys.spacing = isFict.spacing;
  isPhys.voxelValues.resize(isFict.voxelValues.size());

  // make two runs to use parallel architecture
  for (size_t runID = 0; runID < 2; ++runID) {

    // can't use parallel processing on bool vector for some cases
    if (spacing(2) < 8)
      omp_set_num_threads(1);

#pragma omp parallel
    {

#pragma omp for
      for (long k = runID; k < s(2) - 1; k += 2) {
        size_t vxID = static_cast<size_t>(isFict.spacing(2)) * k;

        for (int j = 0; j < s(1) - 1; ++j) {
          vector<size_t> checkCube(8, 0);
          for (int n = 0; n < 8; ++n)
            checkCube[n] =
                spacing.dot(Vector3l(n % 2, j + (n / 2) % 2, k + n / 4));

          for (int i = 0; i < s(0) - 1; ++i) {
            size_t currSum = 0;
            for (int n = 0; n < 8; ++n)
              currSum += this->voxelValues[checkCube[n]] >= isoValue;

            if (currSum != 0)
              isPhys.voxelValues[vxID] = true;

            if (currSum != 8)
              isFict.voxelValues[vxID] = true;

            for (int n = 0; n < 8; ++n)
              ++checkCube[n];

            ++vxID;
          }
        }
      }
    }

    if (spacing(2) < 8)
      omp_set_num_threads(0);
  }
}

//------------------------------------------------------------------------------

template <typename T>
void VoxelVolume<T>::create_from_two_volumes(const VoxelVolume<T> &volA,
                                             const VoxelVolume<T> &volB) {
  if ((volA.s.array() != volB.s.array()).any() || !volA.hasValues ||
      volB.hasValues) {
    cout << "\nWARNING: can't combine Volumes\n";
    return;
  }

  cout << "\nCombining Volumes ...\n";
  s = volA.s;
  spacing = volA.spacing;
  voxelValues.resize(s.cast<size_t>().prod(), 0);
  transform(volA.voxelValues.begin(), volA.voxelValues.end(),
            volB.voxelValues.begin(), voxelValues.begin(), multiplies<T>());

  hasValues = true;
}

//------------------------------------------------------------------------------
template <typename T>
void VoxelVolume<T>::import_raw_volume(Vector3l const &s, string filename) {
  ifstream myFile(filename);

  if (!myFile.good()) {
    cout << "\nWARNING: Can't import volume, no file found!\n";
    return;
  }

  if ((s.array() < 1).any()) {
    cout << "\nWARNING: Can't import volume, bad dimensions!\n";
    return;
  }

  cout << "\nImporting Raw Volume from \"" << filename << "\" ...\n";

  this->s = s;
  spacing = Vector3l(1, s(0), s(0) * s(1));

  voxelValues.resize(s.cast<int32_t>().prod());

  if constexpr (!is_same_v<T, bool>)
    myFile.read((char *)&voxelValues[0], sizeof(T) * voxelValues.size());

  if (myFile.eof() || myFile.peek() != EOF) {
    cout << "\nWARNING: Couldnt't import volume, error while reading in file\n";
    return;
  }

  hasValues = true;
}

//------------------------------------------------------------------------------

template <typename T>
void VoxelVolume<T>::export_pgm_stacks(string foldername) const {
  if (!hasValues) {
    cout << "\nWARNING: Can't export pgm stacks, empty Voxel Volume!\n";
    return;
  }

  cout << "\nExporting as pgm stacks to " << foldername << " ...\n";

  auto minMaxElement = minmax_element(voxelValues.begin(), voxelValues.end());
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

      for (auto vxIt = voxelValues.begin() + spacing.dot(Vector3l(0, j, k));
           vxIt != voxelValues.begin() + spacing.dot(Vector3l(0, j + 1, k));
           ++vxIt, ++pxIt)
        *pxIt = ((*vxIt - minElement) * 255) / (maxElement - minElement);
    }

    char numberBuffer[64];
    sprintf(numberBuffer, "%06i", k);
    ofstream myFile(foldername + "stack" + numberBuffer + ".pgm");
    myFile << "P5\n" << s(0) << " " << s(1) << endl << 255 << endl;
    myFile.write((const char *)currImage.data(), currImage.size());
  }
}

//------------------------------------------------------------------------------
template <typename T> void VoxelVolume<T>::export_raw(string filename) const {

  ofstream myFile(filename);

  cout << "\nExporting Raw Volume to \"" << filename << "\" ...\n";

  if constexpr (!is_same_v<T, bool>)
    myFile.write((char *)&voxelValues[0], sizeof(T) * s.cast<size_t>().prod());
}
//------------------------------------------------------------------------------
template struct VoxelVolume<uint8_t>;
template struct VoxelVolume<bool>;
template struct VoxelVolume<float>;
//------------------------------------------------------------------------------
} // namespace fred
