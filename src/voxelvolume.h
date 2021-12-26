#ifndef VOXELVOLUME_H
#define VOXELVOLUME_H

#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

//------------------------------------------------------------------------------

template <typename T> struct VoxelVolume {

  VoxelVolume()
      : s(Vector3i(0, 0, 0)), spacing(Vector3i(1, 1, 1)), hasValues(false) {}

  // standard functions to translate between ID and Coord
  Vector3i vxID_to_vx(size_t vxID) const {
    return Vector3i(vxID % s(0), (vxID / s(0)) % s(1), vxID / (s(0) * s(1)));
  }

  // use this to always get vxID even for illegal values of vx
  size_t vx_to_vxID(Vector3i const &vx) const;

  // quickly set spacing and voxelValues for given s
  void set_spacing_and_voxelValues_from_s() {
    spacing << 1, s(0), s(0) * s(1);
    voxelValues.resize(s.prod());
  }

  // switch x and y
  void switch_xy();

  // import volume from raw file
  void import_raw_volume(Vector3i const &s, string filename);

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
  void export_stack_for_gp(size_t stackID, string filename) const {
    ofstream myfile(filename);
    for (size_t n = stackID * spacing(2); n < (stackID + 1) * spacing(2); ++n) {
      myfile << n % s(0) << " " << (n / s(0)) % s(1) << " "
             << float(voxelValues[n]) << endl;
      if (n % spacing(1) == spacing(1) - 1)
        myfile << endl;
    }
  }

  T operator[](size_t vxID) const { return voxelValues[vxID]; }
  T &operator[](size_t vxID) { return voxelValues[vxID]; }

  T operator()(Vector3i x) const { return this->operator[](vx_to_vxID(x)); }
  T operator()(long x0, long x1, long x2) const {
    return this->operator()(Vector3i(x0, x1, x2));
  }

  T &operator()(Vector3i x) { return this->operator[](vx_to_vxID(x)); }
  T &operator()(long x0, long x1, long x2) {
    return this->operator()(Vector3i(x0, x1, x2));
  }

  size_t size() { return voxelValues.size(); }

  // voxelValues: value of voxel (i,j,k) is stored at voxelValues(spacing(2)*k +
  // spacing(1)*j + i) s:           size of voxel grid spacing:     spacing in
  // 0/1/2-direction
  vector<T> voxelValues;
  Vector3i s, spacing;
  bool hasValues;
};

//------------------------------------------------------------------------------

template <typename T> void VoxelVolume<T>::switch_xy() {
  VoxelVolume<T> switchedVolume(*this);
  switchedVolume.voxelValues.clear();
  switchedVolume.voxelValues.resize(s.prod());

  Vector3i xi;
  size_t i = 0;
  for (xi(2) = 0; xi(2) < s(2); ++xi(2))
    for (xi(0) = 0; xi(0) < s(0); ++xi(0))
      for (xi(1) = 0; xi(1) < s(1); ++xi(1), ++i)
        switchedVolume[i] = voxelValues[xi.dot(spacing)];

  voxelValues = switchedVolume.voxelValues;
}

//------------------------------------------------------------------------------

template <typename T>
size_t VoxelVolume<T>::vx_to_vxID(Vector3i const &vx) const {
  if ((vx.array() >= 0).all() && (vx.array() < s.array()).all())
    return spacing.dot(vx);

  Vector3i vxNew = vx;
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
  isFict.spacing = Vector3i(1, isFict.s(0), isFict.s(0) * isFict.s(1));
  isFict.voxelValues.resize(isFict.s.prod());

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
      for (size_t k = runID; k < s(2) - 1; k += 2) {
        size_t vxID = k * isFict.spacing(2);

        for (int j = 0; j < s(1) - 1; ++j) {
          vector<size_t> checkCube(8, 0);
          for (int n = 0; n < 8; ++n)
            checkCube[n] =
                spacing.dot(Vector3i(n % 2, j + (n / 2) % 2, k + n / 4));

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
  voxelValues.resize(s.prod(), 0);
  transform(volA.voxelValues.begin(), volA.voxelValues.end(),
            volB.voxelValues.begin(), voxelValues.begin(), multiplies<T>());

  hasValues = true;
}

//------------------------------------------------------------------------------

template <typename T>
void VoxelVolume<T>::import_raw_volume(Vector3i const &s, string filename) {
  ifstream myfile(filename);

  if (!myfile.good()) {
    cout << "\nWARNING: Can't import volume, no file found!\n";
    return;
  }

  if ((s.array() < 1).any()) {
    cout << "\nWARNING: Can't import volume, bad dimensions!\n";
    return;
  }

  cout << "\nImporting Raw Volume from \"" << filename << "\" ...\n";

  this->s = s;
  spacing = Vector3i(1, s(0), s(0) * s(1));

  voxelValues.resize(s.cast<int32_t>().prod());

  myfile.read((char *)&voxelValues[0], sizeof(T) * voxelValues.size());

  if (myfile.eof() || myfile.peek() != EOF) {
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

      for (auto vxIt = voxelValues.begin() + spacing.dot(Vector3i(0, j, k));
           vxIt != voxelValues.begin() + spacing.dot(Vector3i(0, j + 1, k));
           ++vxIt, ++pxIt)
        *pxIt = ((*vxIt - minElement) * 255) / (maxElement - minElement);
    }

    char numberBuffer[64];
    sprintf(numberBuffer, "%06i", k);
    ofstream myfile(foldername + "stack" + numberBuffer + ".pgm");
    myfile << "P5\n" << s(0) << " " << s(1) << endl << 255 << endl;
    myfile.write((const char *)currImage.data(), currImage.size());
  }
}

//------------------------------------------------------------------------------

template <typename T> void VoxelVolume<T>::export_raw(string filename) const {

  ofstream myfile(filename);

  cout << "\nExporting Raw Volume to \"" << filename << "\" ...\n";

  myfile.write((char *)&voxelValues[0], sizeof(T) * s.prod());
}

#endif // VOXELVOLUME_H
