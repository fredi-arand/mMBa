#pragma once

#include "DistanceField.h"
#include "PoreMorphology.h"
#include <fstream>
#include <iostream>
#include <set>

namespace fred {
//------------------------------------------------------------------------------
inline void
gnuplot_distance_field_and_maximal_balls(string const &folderName,
                                         DistanceField const &distanceField) {

  ofstream file_gp(folderName + "distance_field.txt");
  set<pair<float, pair<unsigned, unsigned>>> maximalBalls;

  Vector3l const &s = distanceField.s;

  srand(0);
  map<size_t, float> color;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    color[n] = float(rand()) / RAND_MAX;

  for (int i = 0; i < s(0); ++i)
    for (int j = 0; j < s(1); ++j) {
      file_gp << i << " " << j << " " << distanceField(i, j, 0) << " "
              << color.at(distanceField.vx_to_vxID(Vector3l(i, j, 0))) << endl;
      if (distanceField(i, j, 0) == 0.0)
        continue;
      maximalBalls.insert(
          make_pair(distanceField(i, j, 0), make_pair(s(1) - j, s(0) - i)));
    }

  ofstream file2_gp(folderName + "maximal_balls.txt");
  for (auto it = maximalBalls.end(); it-- != maximalBalls.begin();) {
    float d = it->first;
    unsigned j = (it->second).first;
    j = s(1) - j;
    unsigned i = (it->second).second;
    i = s(0) - i;
    file2_gp << i << " " << j << " " << d << " "
             << color.at(distanceField.vx_to_vxID(Vector3l(i, j, 0))) << endl;
  }
}
//------------------------------------------------------------------------------
template <typename T>
void gnuplot_volume_to_images(string const &folderName,
                              VoxelVolume<T> const &voxelVolume) {

  Vector3l const &s = voxelVolume.s;

  for (int k = 0; k < s(2); ++k) {
    ofstream myFile(folderName + "stack_" + to_string(k) + ".txt");
    for (int j = 0; j < s(1); ++j)
      for (int i = 0; i < s(0); ++i) {
        myFile << i << " " << j << " " << float(voxelVolume(i, j, k)) << endl;
      }
  }
}
//------------------------------------------------------------------------------
static unsigned fileCounter = 0;
//------------------------------------------------------------------------------
inline void gnuplot_palette_file(string fileName,
                                 PoreMorphology const &poreMorphology) {

  auto const &morphologyVolume = poreMorphology.morphologyVolume;
  Vector3l const &s = morphologyVolume.s;
  auto const &parentToVoxelIndex = poreMorphology.parentToVoxelIndex;

  srand(0);
  map<size_t, float> color;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    color[n] = float(rand()) / RAND_MAX;

  ofstream image_palette_file(fileName);
  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    MorphologyValue morphologyValue = morphologyVolume[voxelIndex];
    uint32_t flag = morphologyValue.state;
    uint32_t parent = morphologyValue.parentId;
    Vector3l position = morphologyVolume.vxID_to_vx(voxelIndex);

    if (flag == MorphologyValue::BACKGROUND) {
      image_palette_file << position.transpose() << " -1.0\n";
      continue;
    }

    if (flag == MorphologyValue::THROAT) {
      image_palette_file << position.transpose() << " -2.0 \n";
      continue;
    }

    if (parentToVoxelIndex.count(parent) == 0) {
      image_palette_file << position.transpose() << " -3.0 \n";
      continue;
    }

    size_t parentIndex = parentToVoxelIndex.at(parent);
    float x = color[parentIndex];

    image_palette_file << position.transpose() << " " << x << endl;
  }
}
//------------------------------------------------------------------------------
inline void
update_neighbors_box(DistanceField const &distanceField,
                     VoxelVolume<MorphologyValue> &morphologyVolume,
                     size_t const &voxelIndex_i,
                     map<uint32_t, size_t> const &parentToVoxelIndex) {

  auto const &s = morphologyVolume.s;

  double padding = 0.5;

  Vector3l const voxelCoordinate_i = morphologyVolume.vxID_to_vx(voxelIndex_i);

  MorphologyValue const &morphologyValue_i = morphologyVolume[voxelIndex_i];
  uint32_t parent_i = morphologyValue_i.parentId;

  float const &r_i = distanceField[voxelIndex_i] + padding;
  long const roundedR_i = floor(r_i);
  float const &r_i_squared = r_i * r_i;

  for (long K = -roundedR_i; K <= roundedR_i; ++K)
    for (long J = -roundedR_i; J <= roundedR_i; ++J)
      for (long I = -roundedR_i; I <= roundedR_i; ++I) {
        if (K == 0 && J == 0 && I == 0)
          continue;

        Vector3l const voxelCoordinate_j =
            voxelCoordinate_i + Vector3l(I, J, K);
        if ((voxelCoordinate_j.array() < 0).any() ||
            (voxelCoordinate_j.array() >= s.array()).any())
          continue;

        float r_ij_squared =
            (voxelCoordinate_i - voxelCoordinate_j).cast<float>().squaredNorm();
        if (r_ij_squared > r_i_squared)
          continue;

        size_t const voxelIndex_j =
            morphologyVolume.vx_to_vxID(voxelCoordinate_j);

        MorphologyValue &morphologyValue_j =
            morphologyVolume.voxelValues[voxelIndex_j];
        uint32_t flag_j = morphologyValue_j.state;
        uint32_t parent_j = morphologyValue_j.parentId;

        if (flag_j != MorphologyValue::INIT)
          continue;

        float const &r_j = distanceField[voxelIndex_j] + padding;
        if (r_j > r_i) { /*cout << "\nblub\n";*/
          continue;
        }

        float r_ij = sqrt(r_ij_squared);

        // change to slave if possible
        if (parent_j == 0) {
          morphologyValue_j.parentId = parent_i;
          parent_j = parent_i;
        }

        if (parent_j == parent_i) {

          // try to enclose
          if (r_ij + r_j <= r_i + .2 * r_j) {
            morphologyValue_j.state = MorphologyValue::ENCLOSED;
            flag_j = MorphologyValue::ENCLOSED;
          }

          continue;
        }

        // some value other than the current master has been written
        // --> mark as throat
        morphologyValue_j.state = MorphologyValue::THROAT;
        flag_j = MorphologyValue::ENCLOSED;
      }

  // paint

  ofstream activeBall("thesis/steps/" + to_string(fileCounter) + "ball.txt");
  ofstream allBalls("thesis/steps/" + to_string(fileCounter) + "all.txt");

  srand(0);
  map<size_t, float> color;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    color[n] = float(rand()) / RAND_MAX;

  vector<size_t> processingOrder;

  for (size_t index = 0; index < s.cast<size_t>().prod(); ++index)
    if (distanceField[index] > 0)
      processingOrder.push_back(index);

#ifdef ENABLE_GNU_PARALLEL
  __gnu_parallel::sort(processingOrder.begin(), processingOrder.end(),
                       [&](size_t const &i, size_t const &j) {
                         return distanceField[i] == distanceField[j]
                                    ? i < j // for reproducibility
                                    : distanceField[i] > distanceField[j];
                       });
#else
  sort(processingOrder.begin(), processingOrder.end(),
       [&](size_t const &i, size_t const &j) {
         return distanceField[i] == distanceField[j]
                    ? i < j // for reproducibility
                    : distanceField[i] > distanceField[j];
       });
#endif

  for (size_t n = 0; n < processingOrder.size(); ++n) {
    size_t voxelIndex = processingOrder[n];
    MorphologyValue morphologyValue = morphologyVolume[voxelIndex];
    uint32_t flag = morphologyValue.state;
    uint32_t parent = morphologyValue.parentId;
    Vector3l position = morphologyVolume.vxID_to_vx(voxelIndex);

    if (flag == MorphologyValue::INIT) {
      if (parentToVoxelIndex.count(parent) == 0)
        allBalls << position.transpose() << " " << distanceField[voxelIndex]
                 << " " << -2.0 << endl;
      else {
        auto parentIndex = parentToVoxelIndex.at(parent);
        allBalls << position.transpose() << " " << distanceField[voxelIndex]
                 << " " << -1.0 + color.at(parentIndex) << endl;
      }

      continue;
    }

    if (flag == MorphologyValue::THROAT) {
      allBalls << position.transpose() << " " << distanceField[voxelIndex]
               << " " << 2.0 << endl;
      continue;
    }

    auto parentIndex = parentToVoxelIndex.at(parent);
    allBalls << position.transpose() << " " << distanceField[voxelIndex] << " "
             << color.at(parentIndex) << endl;
  }

  activeBall << voxelCoordinate_i.transpose() << " "
             << distanceField[voxelIndex_i] << " "
             << color.at(parentToVoxelIndex.at(parent_i)) << endl;

  vector<Vector3l> colors2;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    colors2.push_back(Vector3l(rand() % 255, rand() % 255, rand() % 255));

  ofstream image_file("thesis/steps/" + to_string(fileCounter) + "img.txt");

  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    MorphologyValue morphologyValue = morphologyVolume[voxelIndex];
    uint32_t flag = morphologyValue.state;
    uint32_t parent = morphologyValue.parentId;
    Vector3l position = morphologyVolume.vxID_to_vx(voxelIndex);

    if (flag == MorphologyValue::BACKGROUND) {
      image_file << position.transpose() << " 255 255 255 0\n";
      continue;
    }

    if (flag == MorphologyValue::THROAT) {
      image_file << position.transpose() << " 127 127 127 255\n";
      continue;
    }

    if (parentToVoxelIndex.count(parent) == 0) {
      image_file << position.transpose() << " 0 0 0 255\n";
      continue;
    }

    image_file << position.transpose() << " "
               << colors2[parentToVoxelIndex.at(parent)].transpose()
               << " 255\n";
  }

  // in fig10.gp, "set palette rgb 33,13,10" is used.
  // gnuplot: "show palette rgbformulae" explains how values from x=0 to 1 are
  // mapped: 33 -> r = |2*x-0.5| 13 -> g = sin(180*x) 10 -> b = cos(90*x)

  ofstream image_palette_file("thesis/steps/" + to_string(fileCounter) +
                              "img_palette.txt");
  for (size_t voxelIndex = 0; voxelIndex < s.cast<size_t>().prod();
       ++voxelIndex) {
    MorphologyValue morphologyValue = morphologyVolume[voxelIndex];
    uint32_t flag = morphologyValue.state;
    uint32_t parent = morphologyValue.parentId;
    Vector3l position = morphologyVolume.vxID_to_vx(voxelIndex);

    if (flag == MorphologyValue::BACKGROUND) {
      image_palette_file << position.transpose() << " -1.0\n";
      continue;
    }

    if (flag == MorphologyValue::THROAT) {
      image_palette_file << position.transpose() << " -2.0 \n";
      continue;
    }

    if (parentToVoxelIndex.count(parent) == 0) {
      image_palette_file << position.transpose() << " -3.0 \n";
      continue;
    }

    size_t parentIndex = parentToVoxelIndex.at(parent);
    float x = color[parentIndex];

    image_palette_file << position.transpose() << " " << x << endl;
  }

  ++fileCounter;
}
//------------------------------------------------------------------------------
inline void mb_step_by_step(DistanceField const &distanceField,
                            float rMinMaster, float rMinBall) {

  uint32_t parentCounter{0};
  map<uint32_t, size_t> parentToVoxelIndex;

  Vector3l const &s = distanceField.s;

  cout << "\nCreating Pore Morphology:\n";

  VoxelVolume<MorphologyValue> morphologyVolume;

  morphologyVolume.s = distanceField.s;
  morphologyVolume.spacing = distanceField.spacing;

  morphologyVolume.voxelValues.clear();
  morphologyVolume.voxelValues.resize(morphologyVolume.s.cast<size_t>().prod(),
                                      {MorphologyValue::BACKGROUND, 0});

  // each voxel in the void space is its own master
  size_t voidVoxels = 0;
  for (size_t n = 0; n < s.cast<size_t>().prod(); ++n)
    if (distanceField[n] > rMinBall) {
      ++voidVoxels;
      morphologyVolume.voxelValues[n] = {MorphologyValue::INIT, 0};
    }

  if (voidVoxels == 0) {
    cout << "\nnothing to do\n";
    return;
  }

  vector<size_t> processingOrder;
  processingOrder.clear();
  processingOrder.reserve(distanceField.voxelValues.size() / 16);

#ifdef ENABLE_GNU_PARALLEL
  float r_max = *(__gnu_parallel::max_element(distanceField.voxelValues.begin(),
                                              distanceField.voxelValues.end()));
#else
  float r_max = *(max_element(distanceField.voxelValues.begin(),
                              distanceField.voxelValues.end()));
#endif

  float r_infimum = r_max;

  while (r_infimum != 0.0) {
    if (r_max <= 2.0 || r_max <= 2.0 * rMinBall)
      r_infimum = rMinBall;
    else
      r_infimum = r_max / 2.0;

    //    r_infimum = 0.0;

    processingOrder.clear();
    for (size_t index = 0; index < s.cast<size_t>().prod(); ++index)
      if (distanceField[index] > r_infimum && distanceField[index] <= r_max &&
          morphologyVolume[index].state == MorphologyValue::INIT)
        processingOrder.push_back(index);

    //    processingOrder.shrink_to_fit();

    cout << "Pores: " << parentToVoxelIndex.size() << endl;
    cout << scientific << r_infimum << " < r <= " << r_max << endl;

#ifdef ENABLE_GNU_PARALLEL
    __gnu_parallel::sort(processingOrder.begin(), processingOrder.end(),
                         [&](size_t const &i, size_t const &j) {
                           return distanceField[i] == distanceField[j]
                                      ? i < j // for reproducibility
                                      : distanceField[i] > distanceField[j];
                         });
#else
    sort(processingOrder.begin(), processingOrder.end(),
         [&](size_t const &i, size_t const &j) {
           return distanceField[i] == distanceField[j]
                      ? i < j // for reproducibility
                      : distanceField[i] > distanceField[j];
         });
#endif

    if (processingOrder.size() == 0) {
      r_max = r_infimum;
      continue;
    }

    // size_t progressCounter = 0;
    // size_t forLoopCounter = 0;
    for (auto const &voxelIndex_i : processingOrder) {

      float const &r_i = distanceField[voxelIndex_i];

      //      while(progressCounter <=
      //      (forLoopCounter*100)/processingOrder.size())
      //      {
      //        cout << scientific << progressCounter << " %,\tr=" << r_i
      //             << ",\tpores: " << parentToVoxelIndex.size() << endl;
      //        ++progressCounter;
      //      }
      //      ++forLoopCounter;

      //    if(roundedR_i<omp_get_num_threads())
      //      omp_set_num_threads(1);

      MorphologyValue &morphologyValue_i =
          morphologyVolume.voxelValues[voxelIndex_i];
      uint32_t flag_i = morphologyValue_i.state;
      uint32_t parent_i = morphologyValue_i.parentId;

      // cases: throat, enclosed
      if (flag_i != MorphologyValue::INIT)
        continue;

      // case: not allowed to be parent
      if (r_i < rMinMaster && flag_i == MorphologyValue::INIT && parent_i == 0)
        continue;

      // case: parent.
      if (parent_i == 0) {
        ++parentCounter;
        parentToVoxelIndex[parentCounter] = voxelIndex_i;
        morphologyValue_i.parentId = parentCounter;
        parent_i = morphologyValue_i.parentId;
      }

      // ball always encloses itself. Morphology is fixed at this point.
      morphologyValue_i.state = MorphologyValue::ENCLOSED;
      flag_i = morphologyValue_i.state;

      // check and update neighborhood
      update_neighbors_box(distanceField, morphologyVolume, voxelIndex_i,
                           parentToVoxelIndex);
      //    update_neighbors_flood(voxelIndex_i);
    }

    r_max = r_infimum;
  }

  //  omp_set_num_threads(0);

  // count changed voxels
  size_t ignoredVoxels = 0;
  for (auto &morphologyValue : morphologyVolume.voxelValues)
    if (morphologyValue.state == MorphologyValue::INIT) {
      morphologyValue.state = MorphologyValue::BACKGROUND;
      ++ignoredVoxels;
    }

  cout << "\nIgnored Void Voxel Fraction: " << float(ignoredVoxels) / voidVoxels
       << endl;
  cout << "Pores: " << parentToVoxelIndex.size() << endl;
}

} // namespace fred
