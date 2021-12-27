#ifndef IMPORTHOMBERG_H
#define IMPORTHOMBERG_H

#include "Eigen/Dense"
#include "voxelvolume.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace Eigen;

VoxelVolume<uint8_t> create_coverage_rep(string path) {

  // initalize: voxel resolution, volume covered: [A, A+dAB], subcell points
  long int resolution1D = 512;
  Vector3i resolution(resolution1D, resolution1D, resolution1D);
  Vector3f A(-0.5, -0.5, -0.5);
  Vector3f dAB(resolution.cast<float>());
  size_t subcellPoints = 8;

  // fine volume stores binary volumes for subcell points: in or out
  VoxelVolume<bool> fineVolume;
  fineVolume.s = resolution * subcellPoints;
  //  fineVolume.set_spacing_and_voxelValues_from_s();
  fineVolume.spacing =
      Vector3i(1, size_t(fineVolume.s(0)),
               size_t(fineVolume.s(0)) * size_t(fineVolume.s(1)));
  fineVolume.voxelValues.clear();
  fineVolume.voxelValues.resize(size_t(fineVolume.s(0)) *
                                    size_t(fineVolume.s(1)) *
                                    size_t(fineVolume.s(2)),
                                0);

  cout << fineVolume.voxelValues.size() << endl;

  // import points from path
  vector<Vector3f> sphereCenters;
  vector<float> sphereRadii;
  ifstream myfile(path);
  size_t dummyCounter = 0;
  string currentLine;
  while (getline(myfile, currentLine)) {

    if (currentLine.at(0) == '#' || currentLine.size() == 0)
      continue;

    float x0, x1, x2, r;
    istringstream(currentLine) >> x0 >> x1 >> x2 >> r;

    cout << endl << x0 << " " << x1 << " " << x2 << " " << r;

    ++dummyCounter;
    //    if(dummyCounter<16481)
    //      continue;

    Vector3f M(x0, x1, x2);
    M *= double(resolution1D);
    r *= double(resolution1D);
    //    if((M.array()<0).any() ||
    //       (M.array()>=dAB.array()-1).any())
    //      continue;

    //    M*=20.0;
    //    r*=20.0;

    sphereCenters.push_back(M);
    sphereRadii.push_back(r);
  }

  // use balls to erase 1s
  for (size_t sphereID = 0; sphereID < sphereCenters.size(); ++sphereID) {
    cout << endl
         << "sphere " << sphereID + 1 << " of " << sphereCenters.size() << endl;

    // radius and Coordinate w.r.t. coordinates of fine volume
    float rFine = float(subcellPoints) * sphereRadii[sphereID];
    Vector3f mFine = A + float(subcellPoints) * (sphereCenters[sphereID] - A);

    // bounds for check
    Vector3i xFineMin(ceil(mFine(0) - rFine), ceil(mFine(1) - rFine),
                      ceil(mFine(2) - rFine));
    Vector3i xFineMax(floor(mFine(0) + rFine), floor(mFine(1) + rFine),
                      floor(mFine(2) + rFine));

    // rsquared
    float rFineSquared = rFine * rFine;

    long int KMin = max(xFineMin(2), 0);
    long int KMax = min(xFineMax(2), fineVolume.s(2) - 1);
    if (KMax - KMin > 64) {
#pragma omp parallel for
      for (long int K = KMin; K <= KMax; ++K)
        for (long int J = max(xFineMin(1), 0);
             J <= min(xFineMax(1), fineVolume.s(1) - 1); ++J)
          for (long int I = max(xFineMin(0), 0);
               I <= min(xFineMax(0), fineVolume.s(0) - 1); ++I)
            if ((Vector3f(I, J, K) - mFine).squaredNorm() <= rFineSquared)
              fineVolume
                  .voxelValues[size_t(I) * size_t(fineVolume.spacing(0)) +
                               size_t(J) * size_t(fineVolume.spacing(1)) +
                               size_t(K) * size_t(fineVolume.spacing(2))] =
                  true;
    } else {
      for (long int K = KMin; K <= KMax; ++K)
        for (long int J = max(xFineMin(1), 0);
             J <= min(xFineMax(1), fineVolume.s(1) - 1); ++J)
          for (long int I = max(xFineMin(0), 0);
               I <= min(xFineMax(0), fineVolume.s(0) - 1); ++I)
            if ((Vector3f(I, J, K) - mFine).squaredNorm() <= rFineSquared)
              fineVolume
                  .voxelValues[size_t(I) * size_t(fineVolume.spacing(0)) +
                               size_t(J) * size_t(fineVolume.spacing(1)) +
                               size_t(K) * size_t(fineVolume.spacing(2))] =
                  true;
    }

    cout << "\r" << float(sphereID + 1) * 100.0 / sphereCenters.size()
         << " %     " << flush;
    //    cout << "\r" << sphereCenters[sphereID].transpose() << "      " <<
    //    sphereRadii[sphereID] << "     " <<
    //    float(sphereID+1)*100.0/sphereCenters.size() << " %     " << flush;
  }

  cout << endl;

  // write to file voxelVolume
  VoxelVolume<uint8_t> voxelVolume;
  voxelVolume.s = resolution;
  voxelVolume.set_spacing_and_voxelValues_from_s();
  voxelVolume.voxelValues.clear();
  voxelVolume.voxelValues.resize(voxelVolume.s.cast<size_t>().prod(), 0);

  size_t possibleCounts = subcellPoints * subcellPoints * subcellPoints;
#pragma omp parallel for
  for (size_t vxID = 0; vxID < voxelVolume.voxelValues.size(); ++vxID) {
    size_t counts = 0;

    Vector3i vx0Sub = subcellPoints * voxelVolume.vxID_to_vx(vxID);

    for (size_t K = 0; K < subcellPoints; ++K)
      for (size_t J = 0; J < subcellPoints; ++J)
        for (size_t I = 0; I < subcellPoints; ++I)
          counts += fineVolume.voxelValues
                        [size_t(vx0Sub(0) + I) * size_t(fineVolume.spacing(0)) +
                         size_t(vx0Sub(1) + J) * size_t(fineVolume.spacing(1)) +
                         size_t(vx0Sub(2) + K) * size_t(fineVolume.spacing(2))];

    voxelVolume[vxID] = (counts * 255) / possibleCounts;

    //    if(vxID%1024 == 0)
    //      cout << "\r" << float(vxID+1)*100.0/voxelVolume.voxelValues.size()
    //      << " %     " << flush;
  }
  cout << "\r100 %     " << endl;

  return voxelVolume;
}

//------------------------------------------------------------------------------

VoxelVolume<uint8_t> import_homberg(string path) {
  ifstream myfile(path);

  //  Vector3f minX(1000,1000,1000), maxX(-1000,-1000,-1000);
  //  Vector3f minX(0,0,0);
  //  Vector3f maxX(0.05,0.05,0.05);
  Vector3f minX(0, 0, 0);
  Vector3f maxX(1, 1, 1);

  VoxelVolume<uint8_t> greyVolume;
  //  size_t cubeLength = 750;
  size_t cubeLength = 512;
  greyVolume.s = Vector3i(cubeLength, cubeLength, cubeLength);
  greyVolume.set_spacing_and_voxelValues_from_s();
  greyVolume.voxelValues.clear();
  greyVolume.voxelValues.resize(greyVolume.s.cast<size_t>().prod(), 0);

  cout << endl << cubeLength << endl;

  //  size_t numThreads = 16;
  size_t numThreads = omp_get_max_threads();

  vector<map<size_t, VoxelVolume<uint8_t>>> subVoxelMaps(numThreads);

  cout << endl;

  size_t someCounter = 0;

  while (true) {
    string dummy;
    float x0, x1, x2, r;
    myfile >> dummy >> x0 >> x1 >> x2 >> r;

    if (myfile.eof())
      break;

    cout << someCounter++ << endl;

    Vector3f x(x0, x1, x2);

    // uncomment to find minX and maxX before rerunning program
    //    for(int dim=0; dim<3; ++dim)
    //      if(x(dim)-r<minX(dim))
    //        minX(dim)=x(dim)-r;
    //    for(int dim=0; dim<3; ++dim)
    //      if(x(dim)+r>maxX(dim))
    //        maxX(dim)=x(dim)+r;
    //    continue;

    Vector3f vxCenter =
        (x - minX).array() * float(cubeLength) / (maxX - minX).array();
    Vector3i vxMin, vxMax;
    float rVx = r * float(cubeLength) / (maxX(0) - minX(0));

    bool sphereInside = true;
    if ((vxCenter.array() + rVx < 0).any())
      sphereInside = false;
    if ((vxCenter.array() - rVx >= cubeLength).any())
      sphereInside = false;

    if (!sphereInside)
      continue;

    for (size_t dim = 0; dim < 3; ++dim) {
      vxMin(dim) = floor(vxCenter(dim) - rVx);
      if (vxMin(dim) < 0) {
        vxMin(dim) = 0;
      }
      vxMax(dim) = ceil(vxCenter(dim) + rVx);
      if (vxMax(dim) > static_cast<int>(cubeLength)) {
        vxMax(dim) = cubeLength;
      }
    }

    vector<vector<Vector3i>> checkVxs(numThreads);
    for (int k = vxMin(2); k < vxMax(2); ++k)
      for (int j = vxMin(1); j < vxMax(1); ++j)
        for (int i = vxMin(0); i < vxMax(0); ++i)
          checkVxs[greyVolume.vx_to_vxID(Vector3i(i, j, k)) % numThreads]
              .push_back(Vector3i(i, j, k));

    float rVxSquared = rVx * rVx;

#pragma omp parallel for
    for (size_t threadID = 0; threadID < checkVxs.size(); ++threadID)
      for (size_t n = 0; n < checkVxs[threadID].size(); ++n) {
        Vector3i checkVx = checkVxs[threadID][n];
        size_t checkVxID = greyVolume.vx_to_vxID(checkVx);

        if (greyVolume[checkVxID] == 255)
          continue;

        bool processedBefore = false;
        VoxelVolume<uint8_t> subVoxelVolume;
        subVoxelVolume.s = Vector3i(8, 8, 8);
        subVoxelVolume.set_spacing_and_voxelValues_from_s();
        subVoxelVolume.voxelValues.clear();
        subVoxelVolume.voxelValues.resize(8 * 8 * 8, 0);

        // quickly check all corner voxels to save time
        Vector3f cornerCheck(checkVx(0) + 0.5 / 8, checkVx(1) + 0.5 / 8,
                             checkVx(2) + 0.5 / 8);
        bool isInside = (cornerCheck - vxCenter).squaredNorm() <= rVxSquared;
        bool allSame = true;
        for (int cornerID = 1; cornerID < 8; ++cornerID) {
          cornerCheck =
              Vector3f(checkVx(0) + (0.5 + 8 * (cornerID % 2)) / 8.0,
                       checkVx(1) + (0.5 + 8 * ((cornerID / 2) % 2)) / 8.0,
                       checkVx(2) + (0.5 + 8 * (cornerID / 4)) / 8.0);
          if (((cornerCheck - vxCenter).squaredNorm() <= rVxSquared) !=
              isInside) {
            allSame = false;
            break;
          }
        }

        if (allSame) {
          if (isInside)
            greyVolume[checkVxID] = 255;
          continue;
        }

        auto it = subVoxelMaps[threadID].find(checkVxID);
        if (it != subVoxelMaps[threadID].end()) {
          processedBefore = true;
          subVoxelVolume = it->second;
          //#pragma omp critical
          //          cout << "\ncopied from map " << threadID << endl;
        }

        size_t hits = 0;
        for (int K = 0; K < 8; ++K)
          for (int J = 0; J < 8; ++J)
            for (int I = 0; I < 8; ++I) {
              Vector3f xCheck(checkVx(0) + (0.5 + I) / 8.0,
                              checkVx(1) + (0.5 + J) / 8.0,
                              checkVx(2) + (0.5 + K) / 8.0);

              if ((vxCenter - xCheck).squaredNorm() <= rVxSquared) {
                ++hits;
                subVoxelVolume[subVoxelVolume.vx_to_vxID(Vector3i(I, J, K))] =
                    1;
              }
            }

        if (hits == 8 * 8 * 8 && !processedBefore) {
          greyVolume[checkVxID] = 255;
          continue;
        }

        if (hits == 0 && !processedBefore) {
          //          greyVolume[checkVxID] = 0;
          continue;
        }

        subVoxelMaps[threadID][checkVxID] = subVoxelVolume;
      }
  }

  cout << "-\n";
#pragma omp parallel for
  for (size_t threadID = 0; threadID < numThreads; ++threadID)
    for (auto it = subVoxelMaps[threadID].begin();
         it != subVoxelMaps[threadID].end(); ++it) {
      size_t counter = 0;
      for (size_t vxID = 0; vxID < it->second.voxelValues.size(); ++vxID)
        counter += (it->second)[vxID];

      size_t value = (counter * 256) / (8 * 8 * 8 + 1);
      if (value > greyVolume[it->first])
        greyVolume[it->first] = value;
    }

  return greyVolume;
}

VoxelVolume<uint8_t> import_homberg_negative(string path) {
  ifstream myfile(path);

  Vector3f minX(-0.5, -0.5, -0.5);
  Vector3f maxX(511.5, 511.5, 511.5);

  VoxelVolume<uint8_t> greyVolume;
  size_t cubeLength = 512;
  greyVolume.s = Vector3i(cubeLength, cubeLength, cubeLength);
  greyVolume.set_spacing_and_voxelValues_from_s();
  greyVolume.voxelValues.clear();
  greyVolume.voxelValues.resize(greyVolume.s.cast<size_t>().prod(), 255);

  cout << endl << cubeLength << endl;

  //  size_t numThreads = 16;
  size_t numThreads = omp_get_max_threads();

  vector<map<size_t, VoxelVolume<uint8_t>>> subVoxelMaps(numThreads);

  cout << endl;

  while (true) {
    string dummy;
    float x0, x1, x2, r;
    myfile >> dummy >> x0 >> x1 >> x2 >> r;

    if (myfile.eof())
      break;

    cout << ".\n";

    Vector3f x(x0, x1, x2);

    Vector3f vxCenter =
        (x - minX).array() * float(cubeLength) / (maxX - minX).array();
    Vector3i vxMin, vxMax;
    float rVx = r * float(cubeLength) / (maxX(0) - minX(0));

    for (size_t dim = 0; dim < 3; ++dim) {
      vxMin(dim) = floor(vxCenter(dim) - rVx);
      if (vxMin(dim) < 0) {
        vxMin(dim) = 0;
      }
      vxMax(dim) = ceil(vxCenter(dim) + rVx);
      if (vxMax(dim) > static_cast<int>(cubeLength)) {
        vxMax(dim) = cubeLength;
      }
    }

    vector<vector<Vector3i>> checkVxs(numThreads);
    for (int k = vxMin(2); k < vxMax(2); ++k)
      for (int j = vxMin(1); j < vxMax(1); ++j)
        for (int i = vxMin(0); i < vxMax(0); ++i)
          checkVxs[greyVolume.vx_to_vxID(Vector3i(i, j, k)) % numThreads]
              .push_back(Vector3i(i, j, k));

    float rVxSquared = rVx * rVx;

#pragma omp parallel for
    for (size_t threadID = 0; threadID < checkVxs.size(); ++threadID)
      for (size_t n = 0; n < checkVxs[threadID].size(); ++n) {
        Vector3i checkVx = checkVxs[threadID][n];
        size_t checkVxID = greyVolume.vx_to_vxID(checkVx);

        if (greyVolume[checkVxID] == 0)
          continue;

        bool processedBefore = false;
        VoxelVolume<uint8_t> subVoxelVolume;
        subVoxelVolume.s = Vector3i(8, 8, 8);
        subVoxelVolume.set_spacing_and_voxelValues_from_s();
        subVoxelVolume.voxelValues.clear();
        subVoxelVolume.voxelValues.resize(8 * 8 * 8, 1);

        // quickly check all corner voxels to save time
        Vector3f cornerCheck(checkVx(0) + 0.5 / 8, checkVx(1) + 0.5 / 8,
                             checkVx(2) + 0.5 / 8);
        bool isInside = (cornerCheck - vxCenter).squaredNorm() <= rVxSquared;
        bool allSame = true;
        for (int cornerID = 1; cornerID < 8; ++cornerID) {
          cornerCheck =
              Vector3f(checkVx(0) + (0.5 + 8 * (cornerID % 2)) / 8.0,
                       checkVx(1) + (0.5 + 8 * ((cornerID / 2) % 2)) / 8.0,
                       checkVx(2) + (0.5 + 8 * (cornerID / 4)) / 8.0);
          if (((cornerCheck - vxCenter).squaredNorm() <= rVxSquared) !=
              isInside) {
            allSame = false;
            break;
          }
        }

        if (allSame) {
          if (isInside)
            greyVolume[checkVxID] = 0;
          continue;
        }

        auto it = subVoxelMaps[threadID].find(checkVxID);
        if (it != subVoxelMaps[threadID].end()) {
          processedBefore = true;
          subVoxelVolume = it->second;
          //#pragma omp critical
          //          cout << "\ncopied from map " << threadID << endl;
        }

        size_t hits = 0;
        for (int K = 0; K < 8; ++K)
          for (int J = 0; J < 8; ++J)
            for (int I = 0; I < 8; ++I) {
              Vector3f xCheck(checkVx(0) + (0.5 + I) / 8.0,
                              checkVx(1) + (0.5 + J) / 8.0,
                              checkVx(2) + (0.5 + K) / 8.0);

              if ((vxCenter - xCheck).squaredNorm() <= rVxSquared) {
                ++hits;
                subVoxelVolume[subVoxelVolume.vx_to_vxID(Vector3i(I, J, K))] =
                    0;
              }
            }

        if (hits == 8 * 8 * 8 && !processedBefore) {
          greyVolume[checkVxID] = 0;
          continue;
        }

        if (hits == 0 && !processedBefore) {
          //          greyVolume[checkVxID] = 255;
          continue;
        }

        subVoxelMaps[threadID][checkVxID] = subVoxelVolume;
      }
  }

  cout << "-\n";
#pragma omp parallel for
  for (size_t threadID = 0; threadID < numThreads; ++threadID)
    for (auto it = subVoxelMaps[threadID].begin();
         it != subVoxelMaps[threadID].end(); ++it) {
      size_t counter = 0;
      for (size_t vxID = 0; vxID < it->second.voxelValues.size(); ++vxID)
        counter += (it->second)[vxID];

      size_t value = (counter * 256) / (8 * 8 * 8 + 1);
      if (value < greyVolume[it->first])
        greyVolume[it->first] = value;
    }

  return greyVolume;
}

#endif // IMPORTHOMBERG_H
