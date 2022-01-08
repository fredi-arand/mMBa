#include "DistanceField.h"
#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
//------------------------------------------------------------------------------
namespace fred {
//------------------------------------------------------------------------------
using namespace std;
using namespace std::chrono;
using namespace Eigen;
//------------------------------------------------------------------------------
template DistanceField
DistanceField::create<float>(Vector3l const &s, const char *filename,
                             std::optional<float> const &isoValue);
template DistanceField
DistanceField::create<uint8_t>(Vector3l const &s, const char *filename,
                               std::optional<float> const &isoValue);
//------------------------------------------------------------------------------
template <typename T>
DistanceField DistanceField::create(Vector3l const &s, const char *filename,
                                    std::optional<float> const &isoValue) {
  VoxelVolume<T> voxelVolume;
  voxelVolume.import_raw_volume(s, filename);
  return create(voxelVolume, isoValue);
}
//------------------------------------------------------------------------------
template <typename T>
DistanceField DistanceField::create(VoxelVolume<T> const &voxelVolume,
                                    std::optional<float> const &optValue) {

  if (voxelVolume.data.empty()) {
    cout << "\nWARNING: Can't create distance field, empty Voxel Volume!\n";
    return DistanceField();
  }

  float isoValue;
  if (!optValue) {
    auto minMax =
        minmax_element(voxelVolume.data.begin(), voxelVolume.data.end());
    isoValue = (float(*(minMax.first)) + float(*(minMax.second))) / 2.0;
  } else
    isoValue = *optValue;

  high_resolution_clock::time_point tStart = high_resolution_clock::now();

  cout << "\nCreating Distance Field (isoValue: " << float(isoValue)
       << ") ...\n";

  Vector3l const &s = voxelVolume.s;

  float posInf = s.squaredNorm(), negInf = -posInf;
  DistanceField distanceField;
  distanceField.resize(s, posInf);

  // run thrice
  for (size_t dimStart = 0; dimStart < 3; ++dimStart) {
    VoxelVolume<float> testDistanceField;
    testDistanceField.resize(s, posInf);
#pragma omp parallel for
    for (size_t n = 0; n < distanceField.size(); ++n)
      if (voxelVolume[n] >= isoValue)
        testDistanceField[n] = 0;

    cout << "\rIteration " << dimStart + 1 << " of 3" << flush;

    for (size_t dim = dimStart; dim < dimStart + 3; ++dim)
    // run multiple lines in parallel
#pragma omp parallel for
      for (int k = 0; k < s((dim + 2) % 3); ++k)
        for (int j = 0; j < s((dim + 1) % 3); ++j) {

          Vector3l vx(0, 0, 0);
          vx((dim + 2) % 3) = k;
          vx((dim + 1) % 3) = j;
          vector<size_t> vxIDs;
          vxIDs.reserve(s(dim % 3));
          for (vx(dim % 3) = 0; vx(dim % 3) < s(dim % 3); ++vx(dim % 3))
            vxIDs.push_back(distanceField.vx_to_vxID(vx));

          // float coordinates of points on grid line, and value at minimum
          vector<float> xMins;
          xMins.reserve(2 * s(dim % 3));
          vector<float> fMins;
          fMins.reserve(2 * s(dim % 3));
          vector<bool> onJunction;
          onJunction.reserve(2 * s(dim % 3));

          for (int i = 0; i < s(dim % 3) - 1; ++i) {
            xMins.push_back(i);
            onJunction.push_back(true);
            fMins.push_back(testDistanceField[vxIDs[i]]);

            if (dim == dimStart &&
                (voxelVolume[vxIDs[i]] >= isoValue) !=
                    (voxelVolume[vxIDs[i + 1]] >= isoValue) &&
                voxelVolume[vxIDs[i]] != isoValue &&
                voxelVolume[vxIDs[i + 1]] != isoValue) {
              float isoIntersection =
                  float(voxelVolume[vxIDs[i]] - isoValue) /
                  float(voxelVolume[vxIDs[i]] - voxelVolume[vxIDs[i + 1]]);
              if (isoIntersection > 0 && isoIntersection < 1) {
                xMins.push_back((float)i + isoIntersection);
                onJunction.push_back(false);
                fMins.push_back(0);
              }
            }
          }
          xMins.push_back(s(dim % 3) - 1);
          onJunction.push_back(true);
          fMins.push_back(testDistanceField[vxIDs[s(dim % 3) - 1]]);

          // find and save lower envelope
          long int currParabola = 0;
          vector<size_t> envelopeParabolas(1, 0);
          envelopeParabolas.reserve(xMins.size());
          vector<float> envelopeRange(2, posInf);
          envelopeRange.reserve(xMins.size() + 1);
          envelopeRange[0] = negInf;

          for (size_t testParabola = 1; testParabola < xMins.size();
               ++testParabola) {
            bool leftIntersection = false;
          weLoveGoTo:
            float intersection;
            if (currParabola == -1)
              intersection = negInf;
            else
              intersection = (fMins[testParabola] -
                              fMins[envelopeParabolas[currParabola]] +
                              xMins[testParabola] * xMins[testParabola] -
                              xMins[envelopeParabolas[currParabola]] *
                                  xMins[envelopeParabolas[currParabola]]) /
                             (2 * xMins[testParabola] -
                              2 * xMins[envelopeParabolas[currParabola]]);
            if (currParabola != -1 &&
                intersection <= envelopeRange[currParabola]) {
              leftIntersection = true;
              --currParabola;
              goto weLoveGoTo;
            }

            ++currParabola;

            if (leftIntersection) {
              envelopeParabolas[currParabola] = testParabola;
              envelopeParabolas.resize(currParabola + 1);
            } else {
              envelopeParabolas.push_back(testParabola);
            }

            envelopeRange[currParabola] = intersection;
            if (leftIntersection) {
              envelopeRange[currParabola + 1] = posInf;
              envelopeRange.resize(currParabola + 2);
            } else {
              envelopeRange.push_back(posInf);
            }
          }

          // write values to testDistanceField
          currParabola = 0;
          size_t vxID = 0;
          for (size_t testParabola = 0; testParabola < xMins.size();
               ++testParabola) {

            while (envelopeRange[currParabola + 1] < xMins[testParabola])
              ++currParabola;

            if (onJunction[testParabola]) {
              testDistanceField[vxIDs[vxID]] =
                  (xMins[testParabola] -
                   xMins[envelopeParabolas[currParabola]]) *
                      (xMins[testParabola] -
                       xMins[envelopeParabolas[currParabola]]) +
                  fMins[envelopeParabolas[currParabola]];
              ++vxID;
            }
          }

        } // end of 1D line iteration
          // end of plane iteration
          // end of volume iteration in dim-direction

        // update minDistanceField
#pragma omp parallel for // NOLINT
    for (size_t n = 0; n < testDistanceField.size(); ++n)
      if (testDistanceField[n] < distanceField[n])
        distanceField[n] = testDistanceField[n];

  } // end of distance field iteration (exact in dimension dimStart)

#pragma omp parallel for
  for (size_t n = 0; n < distanceField.size(); ++n)
    distanceField[n] = sqrt(distanceField[n]);

  cout << endl;

  high_resolution_clock::time_point tEnd = high_resolution_clock::now();
  cout << "Duration: "
       << double(duration_cast<milliseconds>(tEnd - tStart).count()) / 1000.0
       << " s" << endl;

  return distanceField;
}
//------------------------------------------------------------------------------
void DistanceField::calculate_porosity() const {
  size_t counter = 0;
  for (size_t n = 0; n != data.size(); ++n)
    if (data[n] >= 0.5)
      ++counter;

  cout << endl << float(counter) / float(data.size()) << endl;
}
//------------------------------------------------------------------------------
} // namespace fred
